"""
Download a single mzML file from the manuscript's HuggingFace dataset and
check whether its SpectraCodec-embedded message decodes cleanly with the
current version of spectra_codec.py.
"""
import os
import sys
import json
import traceback

# make spectra_codec importable when run from this directory
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from huggingface_hub import list_repo_files, hf_hub_download
from spectra_codec import SpectraCodec

REPO_ID = "lbnl-metabolomics/20210915_JGI-AK_MK_506588_SoilWaterRep_final_QE-HF_C18_USDAY63680"


def main():
    print(f"Listing files in dataset: {REPO_ID}")
    files = list_repo_files(REPO_ID, repo_type="dataset")
    mzml_files = [f for f in files if f.lower().endswith(".mzml")]
    print(f"Total files: {len(files)}  |  .mzML files: {len(mzml_files)}")
    if not mzml_files:
        print("No .mzML files found. First 20 files:")
        for f in files[:20]:
            print("  -", f)
        sys.exit(1)

    target = sorted(mzml_files)[0]
    print(f"\nDownloading a single file: {target}")
    local_path = hf_hub_download(REPO_ID, target, repo_type="dataset")
    print(f"Downloaded to: {local_path}")

    print("\n--- Decoding with current SpectraCodec defaults (hilbert) ---")
    try:
        codec = SpectraCodec()
        msg = codec.decode_message_from_file(local_path)
        print("Decode succeeded. Message length:", len(msg) if msg else 0)
        print("\n===== DECODED MESSAGE (repr) =====")
        print(repr(msg))
        print("\n===== DECODED MESSAGE (raw) =====")
        print(msg)

        # Integrity check: is it valid JSON (expected metadata format)?
        print("\n===== INTEGRITY CHECK =====")
        try:
            parsed = json.loads(msg)
            print("Valid JSON. Top-level keys:")
            if isinstance(parsed, dict):
                for k in parsed:
                    print("  -", k)
            else:
                print("  (JSON is not an object; type =", type(parsed).__name__, ")")
        except json.JSONDecodeError as e:
            print("NOT valid JSON:", e)
            print("(May still be an intentional non-JSON message.)")
    except Exception:
        print("Decode FAILED with current defaults:")
        traceback.print_exc()
        print("\nThis is consistent with an old-version encoding. "
              "Next step would be trying legacy parameters.")


if __name__ == "__main__":
    main()
