# Embedding experiment files into LCMS runs

This example embeds an experiment's supporting documents — manifest, methods,
protocols, reference sequences, anything — directly into the mzML files of the
runs they describe. The files travel with the data forever: anyone who receives
an mzML (even renamed) can recover every document and verify it byte-for-byte.

## Real-world example

This exact workflow produced the public dataset
[lbnl-metabolomics/ExoW3-NLDM](https://huggingface.co/datasets/lbnl-metabolomics/ExoW3-NLDM):
54 LCMS runs, each carrying the experiment's manifest, methods, protocols, and
compound identifications — plus, in each biological run, the proteome fasta of
that run's strain (`scope: "run_specific"`). Decode any one file and rebuild
the documents:

```python
import json, base64, hashlib
from huggingface_hub import hf_hub_download
from spectra_codec import SpectraCodec

repo = "lbnl-metabolomics/ExoW3-NLDM"
fname = "20260203_EB_MdR_101544-059_ExoW3_20251007_QE119_HILICZ_USHXG03396_NEG_MS2_025_TxCtrl-NLDM-NA-24hr-NA_1__146.mzML"
path = hf_hub_download(repo, fname, repo_type="dataset")

msg = json.loads(SpectraCodec().decode_message_from_file(path))
for e in msg["payload"]["files"]:
    data = (base64.b64decode(e["content"]) if e["content_encoding"] == "base64"
            else e["content"].encode("utf-8"))
    assert hashlib.sha256(data).hexdigest() == e["sha256"]
    open(e["filename"], "wb").write(data)   # working .xlsx / .md / .pdf files
```

## Workflow

```
1. make_payload.py       payload directory + manifest  ->  payloads/<run>.json
2. test_encode_one.py    validate ONE run end-to-end before the batch
3. encode_all.py         mzML + payload  ->  encoded/<run>.mzML   (round-trip verified)
4. decode_all.py         encoded mzML    ->  decoded/<original files>  (sha256 verified)
```

Edit the configuration block at the top of each script.

## Manifest schema

An Excel workbook with one row per LCMS run. Required columns (names
configurable in `make_payload.py`):

| column | purpose |
|---|---|
| `lcms_run_uuid` | unique id per run — embedded as `unique_file_id`, identifies the run even if someone renames the mzML |
| `mzml_stem` | mzML filename without the `.mzML` extension |
| *(optional)* run-specific column | value that prefixes a payload file to embed only in that run (e.g. a strain name matching `<strain>.proteins.faa`) |

## Message format

```json
{
  "unique_file_id": "<run id from manifest>",
  "payload": {
    "payload_format_version": "1.1",
    "recovery_instructions": "...plain-English reconstruction steps...",
    "files": [
      {
        "filename": "methods.md",
        "scope": "common",              // identical in every run
        "file_type": "Markdown text document",
        "mime_type": "text/markdown",
        "content_encoding": "utf-8",    // or "base64" for binary files
        "size_bytes": 5703,
        "sha256": "…",
        "content": "…"
      }
    ]
  }
}
```

Recovery needs nothing but the mzML and `spectra_codec.py`: decode the message,
follow `recovery_instructions`, verify each rebuilt file against its `sha256`.
Files with `scope: "common"` arrive in every run, giving free N-fold redundancy.

## Capacity limits (measured)

The message is zlib-compressed, base64-encoded (7 bits/char), and laid onto a
Hilbert curve; capacity is `4^order` bits, `max_order = 14`:

| compressed message | Hilbert order | peaks added | mzML growth | encode RAM | practical? |
|---|---|---|---|---|---|
| ~0.1 MB | 10 | ~0.5M | a few MB | trivial | yes |
| ~1.4 MB | 12 | ~6M | ~15 MB | ~2 GB | yes |
| ~2.6 MB | 13 | ~12M | ~40 MB | ~4 GB | yes (~20 s/file) |
| ~16 MB | 14 (max) | ~73M | ~0.5–1 GB | >20 GB | avoid |

Tips for staying small:
- Convert Word/PDF documents to Markdown (often 100x smaller).
- Plain-text files (fasta, csv, md) compress well; already-compressed binaries
  (docx, xlsx, pdf) do not.
- Put large run-associated files (e.g. per-strain proteomes) in only their own
  runs via the run-specific column instead of every file.
