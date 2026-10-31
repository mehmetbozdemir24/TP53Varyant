"""LLM modulu testi: ClinVar gonderimleri + (istege bagli) Ollama ile kanit cikarimi.
    python scripts/test_llm.py 12374          # sadece gonderimleri getir
    python scripts/test_llm.py 12374 --llm    # + LLM cikarimi
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tp53 import llm

vid = sys.argv[1] if len(sys.argv) > 1 else "12374"
subs = llm.clinvar_submissions(vid)
print(f"{len(subs)} gonderim, {sum(1 for s in subs if s['comment'])} metinli")
for s in subs[:4]:
    print("-", s["scv"], s["submitter"], s["classification"], "|", s["comment"][:160].replace("\n", " "))
if "--llm" in sys.argv:
    t = time.time()
    ev = llm.extract_evidence(f"ClinVar {vid}", subs)
    print(f"\nLLM ({ev.get('_model')}, {time.time() - t:.0f} sn):")
    print(json.dumps({k: v for k, v in ev.items() if k not in ("evidence_quotes",)}, ensure_ascii=False, indent=1)[:3000])
    print("\nOneri:", llm.to_clinical_suggestion(ev))
