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
  IDs and sha256-verified recovery; see its README for the payload schema
- `examples/manuscript/` — notebooks and figure code for the SpectraCodec paper
  (`dataset_utils.get_peak_height` additionally needs github.com/biorack/metatlas)

## Conventions

- Messages are JSON with a top-level `unique_file_id` (identifies the run even
  if the mzML is renamed) and a `payload` object; embedded files carry
  `content_encoding` (utf-8|base64) + `sha256` for verified reconstruction.
- Testing = round-trip: encode asserts internally; independent verification is
  a fresh instance decoding the written file and comparing.
