"""
Encode every payload message into its mzML file.

Reads:  EXPERIMENT_DIR/<stem>.mzML            (original files, left untouched)
        EXPERIMENT_DIR/payloads/<stem>.json   (from make_payload.py)
Writes: EXPERIMENT_DIR/encoded/<stem>.mzML

encode_message_to_file() internally asserts decode(encoded) == message, so every
output file is round-trip verified. Already-encoded files are skipped, making the
script safe to re-run after an interruption.
"""
import os
import sys
import time
import traceback
from multiprocessing import Pool

import pandas as pd

# make spectra_codec importable when run from this directory
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

# ---------------- configuration ----------------
EXPERIMENT_DIR = "/path/to/your/experiment_directory"
N_WORKERS = 2  # each worker peaks at a few GB of RAM for multi-MB payloads
# ------------------------------------------------

PAYLOADS = os.path.join(EXPERIMENT_DIR, "payloads")
OUT = os.path.join(EXPERIMENT_DIR, "encoded")


def encode_one(stem):
    from spectra_codec import SpectraCodec
    src = os.path.join(EXPERIMENT_DIR, stem + ".mzML")
    dst = os.path.join(OUT, stem + ".mzML")
    tmp = dst + ".part"
    if os.path.exists(dst):
        return (stem, "skipped", 0.0, os.path.getsize(dst))
    t0 = time.time()
    try:
        with open(os.path.join(PAYLOADS, stem + ".json")) as f:
            message = f.read()
        encoder = SpectraCodec()
        encoder.encode_message_to_file(message, src, tmp)  # asserts round-trip
        os.rename(tmp, dst)
        return (stem, "ok", time.time() - t0, os.path.getsize(dst))
    except Exception:
        if os.path.exists(tmp):
            os.remove(tmp)
        return (stem, "FAILED: " + traceback.format_exc(limit=3).replace("\n", " | "),
                time.time() - t0, 0)


def main():
    os.makedirs(OUT, exist_ok=True)
    df = pd.read_csv(os.path.join(EXPERIMENT_DIR, "file_id_map.csv"))
    stems = [f[:-len(".mzML")] for f in df["mzml_filename"]]
    print(f"{len(stems)} runs, {N_WORKERS} workers", flush=True)

    results = []
    t_start = time.time()
    with Pool(N_WORKERS) as pool:
        for i, res in enumerate(pool.imap_unordered(encode_one, stems), 1):
            stem, status, dt, size = res
            results.append(res)
            print(f"[{i}/{len(stems)}] {status:8.8s} {dt:6.1f}s {size/1e6:7.1f}MB  {stem}",
                  flush=True)

    ok = sum(1 for r in results if r[1] == "ok")
    skipped = sum(1 for r in results if r[1] == "skipped")
    failed = [r for r in results if r[1].startswith("FAILED")]
    print(f"\ndone in {(time.time()-t_start)/60:.1f} min: "
          f"{ok} encoded, {skipped} skipped, {len(failed)} failed", flush=True)
    for r in failed:
        print("FAILED:", r[0], "|", r[1], flush=True)

    log = pd.DataFrame(results, columns=["mzml_stem", "status", "seconds", "bytes"])
    log.to_csv(os.path.join(EXPERIMENT_DIR, "encode_log.csv"), index=False)
    print("log written to", os.path.join(EXPERIMENT_DIR, "encode_log.csv"), flush=True)


if __name__ == "__main__":
    main()
