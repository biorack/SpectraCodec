"""
Decode every encoded mzML and reconstruct the original payload files.

Deliberately self-contained: reconstruction uses ONLY information embedded in
the mzML messages (filename, content_encoding, sha256) — no reference to the
original payload directory — demonstrating that the files are fully recoverable
from the spectra alone.

Reads:  EXPERIMENT_DIR/encoded/<stem>.mzML
Writes: EXPERIMENT_DIR/decoded/     reconstructed files
        EXPERIMENT_DIR/decode_log.csv
"""
import os
import sys
import json
import time
import base64
import hashlib
import traceback
from multiprocessing import Pool

import pandas as pd

# make spectra_codec importable when run from this directory
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

# ---------------- configuration ----------------
EXPERIMENT_DIR = "/path/to/your/experiment_directory"
N_WORKERS = 2
# ------------------------------------------------

ENC = os.path.join(EXPERIMENT_DIR, "encoded")
DEC = os.path.join(EXPERIMENT_DIR, "decoded")


def reconstruct_entry(entry):
    """Rebuild one embedded file per its own recovery metadata; verify sha256.
    Returns (filename, sha256, 'written'|'exists'|'MISMATCH...')."""
    if entry["content_encoding"] == "base64":
        data = base64.b64decode(entry["content"])
    else:
        data = entry["content"].encode("utf-8")
    digest = hashlib.sha256(data).hexdigest()
    if digest != entry["sha256"]:
        return (entry["filename"], digest, "SHA MISMATCH vs embedded sha256")
    out = os.path.join(DEC, entry["filename"])
    if os.path.exists(out):
        with open(out, "rb") as f:
            existing = hashlib.sha256(f.read()).hexdigest()
        if existing != digest:
            return (entry["filename"], digest, "CONFLICT with previously decoded copy")
        return (entry["filename"], digest, "exists")
    tmp = out + f".tmp{os.getpid()}"
    with open(tmp, "wb") as f:
        f.write(data)
    os.rename(tmp, out)
    return (entry["filename"], digest, "written")


def decode_one(mzml_filename):
    from spectra_codec import SpectraCodec
    stem = mzml_filename[:-len(".mzML")]
    t0 = time.time()
    try:
        message = SpectraCodec().decode_message_from_file(os.path.join(ENC, mzml_filename))
        obj = json.loads(message)
        uuid = obj["unique_file_id"]
        outcomes = [reconstruct_entry(e) for e in obj["payload"]["files"]]
        bad = [o for o in outcomes if o[2] not in ("written", "exists")]
        status = "ok" if not bad else "FAILED: " + "; ".join(f"{o[0]}: {o[2]}" for o in bad)
        return (stem, uuid, status, len(outcomes), time.time() - t0)
    except Exception:
        return (stem, "", "FAILED: " + traceback.format_exc(limit=3).replace("\n", " | "),
                0, time.time() - t0)


def main():
    os.makedirs(DEC, exist_ok=True)
    id_map = pd.read_csv(os.path.join(EXPERIMENT_DIR, "file_id_map.csv"))
    expected = dict(zip(id_map.mzml_filename, id_map.unique_file_id))
    files = sorted(expected)
    print(f"{len(files)} encoded files, {N_WORKERS} workers", flush=True)

    results = []
    t_start = time.time()
    with Pool(N_WORKERS) as pool:
        for i, res in enumerate(pool.imap_unordered(decode_one, files), 1):
            stem, uuid, status, n_files, dt = res
            uuid_ok = expected.get(stem + ".mzML") == uuid
            results.append((stem, uuid, uuid_ok, status, n_files, dt))
            flag = "ok" if status == "ok" and uuid_ok else "PROBLEM"
            print(f"[{i}/{len(files)}] {flag:7.7s} {dt:5.1f}s files={n_files} "
                  f"uuid_match={uuid_ok}  {stem}", flush=True)

    log = pd.DataFrame(results, columns=["mzml_stem", "decoded_uuid", "uuid_match",
                                         "status", "n_embedded_files", "seconds"])
    log.to_csv(os.path.join(EXPERIMENT_DIR, "decode_log.csv"), index=False)

    n_ok = int(((log.status == "ok") & log.uuid_match).sum())
    print(f"\ndone in {(time.time()-t_start)/60:.1f} min: {n_ok}/{len(files)} fully ok",
          flush=True)
    for _, r in log[(log.status != "ok") | (~log.uuid_match)].iterrows():
        print("PROBLEM:", r.mzml_stem, "| uuid_match:", r.uuid_match, "|", r.status, flush=True)
    print("reconstructed files in", DEC, ":", sorted(os.listdir(DEC)), flush=True)


if __name__ == "__main__":
    main()
