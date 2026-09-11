"""
Build per-run SpectraCodec payload messages from a directory of supporting files.

Every LCMS run gets a message of the form:

    {"unique_file_id": <run id from your manifest>, "payload": {...}}

- Files in PAYLOAD_DIR are embedded identically in every run ("common" scope),
  EXCEPT files matched by a run's RUN_SPECIFIC_COLUMN value, which are embedded
  only in that run ("run_specific" scope). Set RUN_SPECIFIC_COLUMN = None to
  embed everything in every run.
- Each file entry records file type / mime type / encoding / sha256 so the
  original files can be reconstructed and verified on recovery.

Outputs (written to EXPERIMENT_DIR):
  payloads/<mzml_stem>.json   one message per LCMS run
  file_id_map.csv             mzML filename -> unique id mapping
  payload_stats.txt           size accounting per run

See README.md in this directory for the manifest schema and capacity limits.
"""
import os
import json
import base64
import zlib
import hashlib
from datetime import date

import pandas as pd

# ---------------- configuration ----------------
PAYLOAD_DIR = "/path/to/your/payload_directory"      # files to embed
EXPERIMENT_DIR = "/path/to/your/experiment_directory"  # where the mzML files live
MANIFEST_XLSX = os.path.join(PAYLOAD_DIR, "manifest.xlsx")  # one row per LCMS run
MANIFEST_SHEET = "runs"
ID_COLUMN = "lcms_run_uuid"        # unique id per run (survives file renaming)
FILENAME_COLUMN = "mzml_stem"      # mzML filename without .mzML extension
# Column whose value prefixes a run-specific file in PAYLOAD_DIR
# (e.g. strain name matching "<strain>.proteins.faa"), or None to disable:
RUN_SPECIFIC_COLUMN = None
EXPERIMENT_DESCRIPTION = "Describe your experiment here"
# ------------------------------------------------

# file types treated as plain text (embedded as utf-8 strings, compress well)
TEXT_EXTENSIONS = {".faa", ".md", ".txt", ".csv", ".json", ".fasta", ".fa"}

FILE_TYPE_INFO = {
    ".xlsx": ("Microsoft Excel workbook",
              "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
    ".docx": ("Microsoft Word document",
              "application/vnd.openxmlformats-officedocument.wordprocessingml.document"),
    ".pdf":  ("PDF document", "application/pdf"),
    ".faa":  ("FASTA amino acid sequences (plain text)", "text/x-fasta"),
    ".md":   ("Markdown text document", "text/markdown"),
    ".csv":  ("Comma-separated values (plain text)", "text/csv"),
}

RECOVERY_INSTRUCTIONS = (
    "Each entry in payload['files'] is one original file. To reconstruct: if "
    "content_encoding == 'base64', base64-decode 'content' and write the bytes to "
    "'filename' (binary mode); if content_encoding == 'utf-8', write 'content' as "
    "UTF-8 text to 'filename'. Verify each file by comparing its SHA-256 hex digest "
    "to 'sha256'. Entries with scope == 'common' are identical in every LCMS run of "
    "this experiment; entries with scope == 'run_specific' appear only in their "
    "matching runs, so recovering the complete file set may require decoding "
    "several runs. 'unique_file_id' identifies the run this message is embedded "
    "in, even if the mzML file is renamed."
)


def file_entry(name, scope, extra=None):
    path = os.path.join(PAYLOAD_DIR, name)
    ext = os.path.splitext(name)[1].lower()
    with open(path, "rb") as f:
        raw = f.read()
    ftype, mime = FILE_TYPE_INFO.get(ext, (f"{ext} file", "application/octet-stream"))
    entry = {
        "filename": name,
        "scope": scope,
        "file_type": ftype,
        "mime_type": mime,
        "size_bytes": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest(),
    }
    if extra:
        entry.update(extra)
    if ext in TEXT_EXTENSIONS:
        entry["content_encoding"] = "utf-8"
        entry["content"] = raw.decode("utf-8")
    else:
        entry["content_encoding"] = "base64"
        entry["content"] = base64.b64encode(raw).decode("ascii")
    return entry


def main():
    payload_files = sorted(
        f for f in os.listdir(PAYLOAD_DIR)
        if os.path.isfile(os.path.join(PAYLOAD_DIR, f)) and not f.startswith(".")
    )
    df = pd.read_excel(MANIFEST_XLSX, sheet_name=MANIFEST_SHEET)

    # split payload files into common vs run-specific
    specific_entries = {}
    if RUN_SPECIFIC_COLUMN:
        for key in df[RUN_SPECIFIC_COLUMN].dropna().unique():
            matches = [f for f in payload_files if f.startswith(str(key) + ".")]
            assert len(matches) == 1, f"expected one run-specific file for {key}: {matches}"
            specific_entries[key] = file_entry(
                matches[0], "run_specific", {RUN_SPECIFIC_COLUMN: str(key)})
        specific_names = {e["filename"] for e in specific_entries.values()}
    else:
        specific_names = set()
    common_entries = [file_entry(n, "common")
                      for n in payload_files if n not in specific_names]
    print("common files:", [e["filename"] for e in common_entries])
    print("run-specific files:", sorted(specific_names))

    payload_out = os.path.join(EXPERIMENT_DIR, "payloads")
    os.makedirs(payload_out, exist_ok=True)

    stats, id_rows = [], []
    for _, row in df.iterrows():
        key = row[RUN_SPECIFIC_COLUMN] if RUN_SPECIFIC_COLUMN else None
        files = list(common_entries)
        if pd.notna(key) and key in specific_entries:
            files.append(specific_entries[key])
        message = {
            "unique_file_id": row[ID_COLUMN],
            "payload": {
                "payload_format_version": "1.1",
                "created": str(date.today()),
                "experiment": EXPERIMENT_DESCRIPTION,
                "recovery_instructions": RECOVERY_INSTRUCTIONS,
                "n_files": len(files),
                "files": files,
            },
        }
        text = json.dumps(message)
        out_path = os.path.join(payload_out, row[FILENAME_COLUMN] + ".json")
        with open(out_path, "w") as f:
            f.write(text)

        # size accounting through the SpectraCodec pipeline stages
        comp = zlib.compress(text.encode("utf-8"))
        n_bits = len(base64.b64encode(comp)) * 7
        order = next(o for o in range(15) if 4**o > n_bits + 2)
        stats.append(
            f"{row[FILENAME_COLUMN]}: json={len(text)/1e6:.2f}MB "
            f"compressed={len(comp)/1e6:.2f}MB bits={n_bits/1e6:.1f}M order={order}"
        )
        id_rows.append({
            "mzml_filename": row[FILENAME_COLUMN] + ".mzML",
            "unique_file_id": row[ID_COLUMN],
            "payload_json": os.path.basename(out_path),
        })

    pd.DataFrame(id_rows).to_csv(os.path.join(EXPERIMENT_DIR, "file_id_map.csv"), index=False)
    with open(os.path.join(EXPERIMENT_DIR, "payload_stats.txt"), "w") as f:
        f.write("\n".join(stats) + "\n")

    orders = [int(s.rsplit("order=", 1)[1]) for s in stats]
    print(f"\nwrote {len(id_rows)} payload messages to {payload_out}")
    print(f"hilbert orders required: min={min(orders)} max={max(orders)} (max supported: 14)")
    if max(orders) >= 14:
        print("WARNING: order 14 payloads need >20 GB RAM to encode/decode and add "
              "~0.5-1 GB of peaks per mzML. Consider shrinking the payload "
              "(convert docx to md, drop or shard large files).")


if __name__ == "__main__":
    main()
