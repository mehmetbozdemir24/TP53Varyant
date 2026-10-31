"""Tum hatti sirayla calistirir.

    python run_pipeline.py              # build + dogrulama + ML
    python run_pipeline.py --download   # once veri setlerini yeniden indir
    python run_pipeline.py --esm2       # ESM-2 LLR hesapla (GPU) ve ML'i onunla yeniden egit
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def run(*args):
    print(f"\n>>> {' '.join(args)}", flush=True)
    subprocess.run([sys.executable, *args], cwd=ROOT, check=True)


if __name__ == "__main__":
    if "--download" in sys.argv:
        run("scripts/download.py")
    run("-m", "tp53.build")
    run("scripts/validate_vcep.py")
    if "--esm2" in sys.argv:
        run("scripts/esm2_llr.py")
    run("-m", "tp53.ml")
    print("\nTamam. Arayuz icin:  streamlit run app.py")
