# SpectraCodec

Encodes arbitrary text (typically JSON metadata) as (m/z, intensity) peaks in
the first spectrum of an mzML file, using a Hilbert space-filling curve. The
message travels inside the mass-spec data itself: it survives copying, renaming,
and re-upload, and is recovered from the spectra alone.

## Core API

Everything lives in `spectra_codec.py` (single-file library, no package install):

```python
from spectra_codec import SpectraCodec

codec = SpectraCodec()
codec.encode_message_to_file(message, "input.mzML", "output.mzML")
message = codec.decode_message_from_file("output.mzML")
```

`encode_message_to_file` internally asserts `decode(output) == message`, so a
successful encode is already round-trip verified. Use a fresh `SpectraCodec()`
instance per file — encoding mutates instance state (`self.order`, curve arrays).

## Decoding: recovering embedded files (Word docs, spreadsheets, PDFs, ...)

The most common task: someone hands you an mzML and you want the documents
embedded in it. The decoded message is JSON; each entry in `payload["files"]`
carries everything needed to rebuild the original file and prove it's intact:

```python
import json, base64, hashlib
from spectra_codec import SpectraCodec

msg = json.loads(SpectraCodec().decode_message_from_file("run.mzML"))
print(msg["unique_file_id"])          # identifies the run even if renamed

for entry in msg["payload"]["files"]:
    # entry["file_type"] / entry["mime_type"] say what it is
    # (e.g. "Microsoft Word document", "Microsoft Excel workbook")
    if entry["content_encoding"] == "base64":      # binary: docx, xlsx, pdf
        data = base64.b64decode(entry["content"])
    else:                                          # utf-8 text: md, csv, fasta
        data = entry["content"].encode("utf-8")
    assert hashlib.sha256(data).hexdigest() == entry["sha256"], entry["filename"]
    with open(entry["filename"], "wb") as f:       # writes a working .docx/.xlsx/...
        f.write(data)
```

Key facts:
- `content_encoding` tells you how to rebuild: `base64` → decode to bytes and
  write binary; `utf-8` → the string IS the file content. A rebuilt `.docx` or
  `.xlsx` opens directly in Word/Excel — it is byte-identical to the original.
- Always verify `sha256` after rebuilding; it is the integrity proof.
- `scope: "common"` files are identical in every run of an experiment;
  `scope: "run_specific"` files exist only in their run — recovering the full
  set may require decoding several runs.
- Batch version with logging and cross-run consistency checks:
  `examples/embed_experiment_files/decode_all.py`.
- Decode is read-only and fast (~3–10 s for multi-MB payloads, dominated by
  the Hilbert-coordinate merge; small metadata-only messages take ~50 ms).

## Signing: tamper evidence (two Ed25519 signatures per message)

Every encode signs automatically when a key is available ($SPECTRACODEC_SIGNING_KEY
or `~/secrets/spectracodec_private_key`); it adds a top-level `provenance` block
right after `unique_file_id` with two signatures:

- `payload_signature` — over the canonical JSON of the message minus
  `provenance`: proves the metadata/embedded documents are authentic.
- `spectra_signature` — over `unique_file_id` + the spectral digest: proves
  the acquired spectra are untampered.

The spectral digest (`scd-1`) is a chained SHA-256 over every spectrum
**except the first** (the carrier), hashing position, native ID, and the
decoded m/z + intensity arrays as little-endian float64 — so it is identical
before and after encoding. That exclusion is what makes sign-then-embed
possible; do not change it without versioning a new spec (`scd-2`, ...).

Verify (works for any third party with the published verification key):

```bash
python spectra_codec.py verify run.mzML --key spectracodec_verification_key.pem
```

The repo's published verification key is `spectracodec_verification_key.pem`
(fingerprint in README.md).

or `spectra_codec.verify_signed_file(path, verification_key_path=...)` →
report dict with `valid` + per-check booleans. Omitting the key verifies
against the key embedded in the file — internal consistency only, NOT origin
(report says `embedded_untrusted`). Unsigned/legacy files report
`signed: False` gracefully. `python spectra_codec.py keygen <priv> <pub>`
mints a keypair. Sha256 hashes and the mzML checksum are unkeyed and can be
recomputed by an attacker; only the signatures prove origin. The private key
must never enter this (public) repo.

## How the encoding works (pipeline order matters)

message → UTF-8 → `zlib.compress` → base64 → 7 bits per base64 char (the
65-symbol alphabet incl. `=` needs 7 bits, not 6) → bits laid along a Hilbert
curve → coordinates of 1-bits become peaks:
`mz = x*0.0001 + 5`, `intensity = y*0.001 + 100`.
Decoding reads peaks in the m/z 5–12 / intensity 100–1000 window of the first
spectrum and reverses the pipeline.

## Hard-won constraints

- **Capacity**: a Hilbert curve of order N holds 4^N bits; `max_order=14`
  (≈268M bits ≈ 28 MB of compressed payload) is the ceiling. Practical limits
  are much lower: order 13 (~2.6 MB compressed) ≈ 20 s and ~4 GB RAM per file;
  order 14 needs >20 GB RAM and adds ~0.5–1 GB of peaks per mzML. Keep payloads
  at order ≤13.
- **lxml huge_tree**: the three `etree.parse()` calls in `spectra_codec.py` use
  `XMLParser(huge_tree=True)`. Without it, encoded peak arrays >10 MB make lxml
  refuse the file ("Text node too long"). Do not remove.
- Payloads are zlib-compressed as a whole: plain text (fasta/md/csv) shrinks a
  lot; already-compressed binaries (docx/xlsx/pdf) do not. Converting docs to
  Markdown before embedding often shrinks payloads ~100x.
- The commented-out class at the bottom of `spectra_codec.py` is a legacy
  encoding (different parameters); files encoded with it will not decode with
  the current defaults.

## Environment

```bash
mamba env create -f environment.yml   # creates spectra_codec_env
```
Core runtime deps are just: numpy, pandas, lxml, pymzml, pyteomics; examples
add huggingface_hub (and openpyxl for xlsx manifests).

## Layout

- `spectra_codec.py` — the library (this is the whole product)
- `examples/decode_public_dataset/` — download one file from the public HF
  dataset (`bpbowen/20210915_JGI-AK_MK_506588_SoilWaterRep_...`) and decode its
  embedded message; plus a decode-speed benchmark
- `examples/embed_experiment_files/` — full workflow embedding an experiment's
  documents (manifest, methods, sequences) into every run, with per-run unique
  IDs and sha256-verified recovery; see its README for the payload schema.
  Real output of this workflow: the public dataset
  `lbnl-metabolomics/ExoW3-NLDM` (54 runs, docs + per-strain proteome fasta
  embedded in each)
- `examples/manuscript/` — notebooks and figure code for the SpectraCodec paper
  (`dataset_utils.get_peak_height` additionally needs github.com/biorack/metatlas)

## Conventions

- Messages are JSON with a top-level `unique_file_id` (identifies the run even
  if the mzML is renamed) and a `payload` object; embedded files carry
  `content_encoding` (utf-8|base64) + `sha256` for verified reconstruction.
- Testing = round-trip: encode asserts internally; independent verification is
  a fresh instance decoding the written file and comparing.
