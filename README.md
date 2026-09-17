  <img src="logo.png" alt="SpectraCodec Logo" width="100%"/>
# A Hilbert curve-based method for encoding metadata as data 

If you have questions about your rights to use, distribute this software, or use for commercial purposes, please contact Berkeley Lab's Intellectual Property Office at IPO@lbl.gov.


## Step 1: Define a message.  Typically this will be JSON text, but here is something more fun.

```python
# # Test the conversion
message = """
🌐 SpectraCodec: Hilbert curve metadata encoding 🌐
العربية: سبيكترا كودك: ترميز البيانات بمنحنى هيلبرت
हिंदी: स्पेक्ट्राकोडेक: हिल्बर्ट वक्र डेटा एन्कोडिंग
ქართული: სპექტრაკოდეკი: ჰილბერტის მრუდი
ไทย: สเปกตราโคเดก: การเข้ารหัสข้อมูลด้วยเส้นโค้งฮิลเบิร์ต
日本語: スペクトラコーデック：ヒルベルト曲線符号化
አማርኛ: ስፔክትራኮዴክ፡ ሂልበርት ኩርቭ ኢንኮዲንግ
தமிழ்: ஸ்பெக்ட்ராகோடெக்: ஹில்பர்ட் வளைவு குறியீடு
Монгол: СпектраКодек: Хилбертийн муруй
Ελληνικά: ΣπέκτραΚόντεκ: Κωδικοποίηση καμπύλης Χίλμπερτ
עברית: ספקטראקודק: קידוד עקומת הילברט
🔬💻🌊🎨📊✨
"""
```

## Step 2: Encode the message directly into an mzML file as (m/z, intensity) peaks

```python
from spectra_codec import SpectraCodec
encoder = SpectraCodec()
original_filename = 'your_file.mzML'
output_filename = 'output.mzML'
encoder.encode_message_to_file(message, original_filename, output_filename)
```

## Step 3: Extract the SpectraCodec message (metadata) from the mzML file

```python
message = encoder.decode_message_from_file(output_filename)#, parser='pymzml', method='hilbert')
print("Decoded message:", message)
```

## Try it on real data

Two public datasets carry SpectraCodec messages in every LCMS run:

- [lbnl-metabolomics/20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680](https://huggingface.co/datasets/lbnl-metabolomics/20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680)
  — each run's full metadata (and the manuscript it belongs to)
- [lbnl-metabolomics/ExoW3-NLDM](https://huggingface.co/datasets/lbnl-metabolomics/ExoW3-NLDM)
  — each run embeds the experiment's manifest, methods, protocols, compound
  identifications, and its strain's proteome fasta; the original documents
  (Excel, PDF, Markdown, fasta) are recoverable byte-for-byte from the spectra

Decode one:

```bash
python examples/decode_public_dataset/check_hf_message.py
```

## Signed messages: proving a run is authentic

Hashes alone can't stop tampering — anyone can edit spectra and recompute
every checksum. SpectraCodec therefore embeds two Ed25519 signatures in a
`provenance` block of each signed message: one over the message payload (the
metadata and embedded documents), one over a digest of all acquired spectra
(everything except the carrier spectrum that holds the message). An attacker
who alters either the data or the metadata cannot re-create the signatures
without the signer's private key.

Verify a run against a published verification key:

```bash
python spectra_codec.py verify run.mzML --key spectracodec_verification_key.pem
```

`spectracodec_verification_key.pem` in this repository is the LBNL
SpectraCodec verification key. Its fingerprint (SHA-256 of the raw 32-byte
public key) is:

```
SHA256:lbrclU9WYbrXjl/AWl8VOWFUEm5EQgGybgdSKZzTwfA
```

A valid signature from this key means the run was signed by the holder of the
corresponding private key and has not been modified since. Verifying without
`--key` uses the key embedded in the file itself, which only proves internal
consistency — not origin. To sign your own runs, mint a keypair with
`python spectra_codec.py keygen`, keep the private key out of any repository,
and publish the verification key and its fingerprint somewhere you control.

## Installation

No package install needed — `spectra_codec.py` is a single-file library.
For a ready-made environment:

```bash
mamba env create -f environment.yml   # creates spectra_codec_env
```

Core dependencies: numpy, pandas, lxml, pymzml, pyteomics, psims.
(`openpyxl` is needed to open recovered `.xlsx` files in the examples.)

## Examples

- `examples/decode_public_dataset/` — download a file from the HF dataset and
  decode its embedded message; plus a decode-speed benchmark.
- `examples/embed_experiment_files/` — embed an experiment's supporting
  documents (manifest, methods, protocols, sequences) into every run of an
  experiment, with a unique ID per run and sha256-verified recovery of the
  original files. Start with its README.
- `examples/manuscript/` — notebooks and figure code from the SpectraCodec
  manuscript.

