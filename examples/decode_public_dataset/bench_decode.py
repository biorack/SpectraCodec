"""
Benchmark: SpectraCodec message decode vs SQLite read of the same metadata.
Uses the single mzML file already downloaded from the HF dataset.
"""
import time
import os
import sys
import glob
import json
import sqlite3
import tempfile

import numpy as np
import pymzml

# make spectra_codec importable when run from this directory
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from spectra_codec import SpectraCodec

# locate the cached file downloaded by check_hf_message.py
cache = os.path.expanduser("~/.cache/huggingface/hub")
mzmls = glob.glob(os.path.join(cache, "**", "*.mzML"), recursive=True)
path = mzmls[0]
size_mb = os.path.getsize(path) / 1e6
print(f"File: {os.path.basename(path)}")
print(f"Size on disk: {size_mb:.1f} MB\n")

def timeit(fn, n=5):
    ts = []
    out = None
    for _ in range(n):
        t0 = time.perf_counter()
        out = fn()
        ts.append(time.perf_counter() - t0)
    return min(ts), sum(ts)/len(ts), out

# ---- Full SpectraCodec decode (file open + spectrum read + hilbert decode) ----
def full_decode():
    return SpectraCodec().decode_message_from_file(path)

t_min, t_avg, msg = timeit(full_decode)
print(f"[SpectraCodec] full decode_message_from_file: min={t_min*1000:.1f} ms  avg={t_avg*1000:.1f} ms")
print(f"               message length: {len(msg)} chars\n")

# ---- Decompose: (a) read first spectrum via pymzml, (b) hilbert decode only ----
def read_first_spectrum():
    run = pymzml.run.Reader(path)
    for spec in run:
        return np.array(spec.mz), np.array(spec.i)

t_min_read, t_avg_read, (mz, inten) = timeit(read_first_spectrum)
print(f"[breakdown] (a) pymzml open + read first spectrum: min={t_min_read*1000:.1f} ms  avg={t_avg_read*1000:.1f} ms")

def hilbert_only():
    c = SpectraCodec()
    c.handle_decoding_details(mz.copy(), inten.copy(), method='hilbert')
    return c.decoded_message

t_min_h, t_avg_h, _ = timeit(hilbert_only)
print(f"[breakdown] (b) hilbert decode only (coords->text): min={t_min_h*1000:.1f} ms  avg={t_avg_h*1000:.1f} ms\n")

# ---- SQLite: store the same JSON, then read it back ----
parsed = json.loads(msg)
db = os.path.join(tempfile.gettempdir(), "meta_bench.sqlite")
if os.path.exists(db):
    os.remove(db)
con = sqlite3.connect(db)
con.execute("CREATE TABLE meta (id INTEGER PRIMARY KEY, filename TEXT, json TEXT)")
con.execute("INSERT INTO meta (filename, json) VALUES (?, ?)",
            (os.path.basename(path), msg))
con.commit()
con.close()

def sqlite_read_and_parse():
    con = sqlite3.connect(db)
    row = con.execute("SELECT json FROM meta WHERE id=1").fetchone()
    con.close()
    return json.loads(row[0])

def sqlite_read_persistent():
    con = sqlite3.connect(db)
    row = con.execute("SELECT json FROM meta WHERE id=1").fetchone()
    d = json.loads(row[0])
    con.close()
    return d

t_min_sql, t_avg_sql, _ = timeit(sqlite_read_and_parse, n=50)
print(f"[SQLite] connect + SELECT + json.loads: min={t_min_sql*1000:.3f} ms  avg={t_avg_sql*1000:.3f} ms")

# with a warm/open connection (amortized, e.g. reading many rows)
con = sqlite3.connect(db)
def sqlite_warm():
    row = con.execute("SELECT json FROM meta WHERE id=1").fetchone()
    return json.loads(row[0])
t_min_warm, t_avg_warm, _ = timeit(sqlite_warm, n=50)
con.close()
print(f"[SQLite] warm connection SELECT + json.loads: min={t_min_warm*1000:.3f} ms  avg={t_avg_warm*1000:.3f} ms\n")

print("==== SUMMARY (min times) ====")
print(f"SpectraCodec full decode:      {t_min*1000:9.1f} ms")
print(f"  of which file read (pymzml): {t_min_read*1000:9.1f} ms")
print(f"  of which hilbert decode:     {t_min_h*1000:9.1f} ms")
print(f"SQLite cold (connect+query):   {t_min_sql*1000:9.3f} ms")
print(f"SQLite warm (query only):      {t_min_warm*1000:9.3f} ms")
print(f"\nSlowdown vs SQLite cold: {t_min/t_min_sql:8.0f}x")
print(f"Slowdown vs SQLite warm: {t_min/t_min_warm:8.0f}x")
