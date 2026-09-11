"""End-to-end test on a single run before committing to the full batch:
encode one payload into its mzML, decode it back, and verify the embedded
files reconstruct with correct SHA-256. Run this first — it surfaces payload
size problems (Hilbert order / RAM) in one file instead of one hundred."""
import os
import sys
import json
import time
import base64
import hashlib
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from spectra_codec import SpectraCodec

# ---------------- configuration ----------------
EXPERIMENT_DIR = "/path/to/your/experiment_directory"
STEM = "your_mzml_filename_without_extension"
# ------------------------------------------------

src = os.path.join(EXPERIMENT_DIR, STEM + ".mzML")
dst = os.path.join(tempfile.gettempdir(), STEM + ".encoded.mzML")

with open(os.path.join(EXPERIMENT_DIR, "payloads", STEM + ".json")) as f:
    message = f.read()
print(f"message: {len(message)/1e6:.2f} MB  |  source mzML: {os.path.getsize(src)/1e6:.1f} MB")

t0 = time.time()
enc = SpectraCodec()
enc.encode_message_to_file(message, src, dst)  # asserts decode == original internally
print(f"encode+internal-verify: {time.time()-t0:.1f} s  |  "
      f"encoded mzML: {os.path.getsize(dst)/1e6:.1f} MB "
      f"(+{(os.path.getsize(dst)-os.path.getsize(src))/1e6:.1f} MB)")

t0 = time.time()
dec = SpectraCodec()  # fresh instance = independent decode
recovered = dec.decode_message_from_file(dst)
print(f"independent decode: {time.time()-t0:.1f} s  |  match: {recovered == message}")

obj = json.loads(recovered)
print("unique_file_id:", obj["unique_file_id"])
ok = 0
for e in obj["payload"]["files"]:
    data = (base64.b64decode(e["content"]) if e["content_encoding"] == "base64"
            else e["content"].encode("utf-8"))
    assert hashlib.sha256(data).hexdigest() == e["sha256"], f"sha mismatch: {e['filename']}"
    ok += 1
print(f"all {ok} embedded files reconstruct with correct sha256")
os.remove(dst)
print("PASS")
