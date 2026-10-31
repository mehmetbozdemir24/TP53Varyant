"""ESM-2 maskeli marjinal log-olabilirlik orani (LLR) - TP53 (P04637) icin tum 393 x 19 missense.

LLR(pos, mut) = log p(mut | maskelenmis dizi) - log p(wt | maskelenmis dizi)   (negatif = zararli)
Ayrica her pozisyondaki tahmin entropisi (dusuk entropi = korunmus pozisyon).

Calistirma (GPU onerilir):
    python scripts/esm2_llr.py                         # varsayilan facebook/esm2_t33_650M_UR50D
    python scripts/esm2_llr.py --model facebook/esm2_t36_3B_UR50D
Cikti: data/processed/esm2_llr.csv  (pos, wt, mut, esm2_llr, esm2_entropy)
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import torch
from transformers import AutoTokenizer, EsmForMaskedLM

ROOT = Path(__file__).resolve().parents[1]
AAS = "ACDEFGHIKLMNPQRSTVWY"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--batch", type=int, default=8)
    a = ap.parse_args()
    seq = json.load(open(ROOT / "data" / "raw" / "uniprot_P04637.json", encoding="utf-8"))["sequence"]["value"]
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"model={a.model} cihaz={dev} uzunluk={len(seq)}", flush=True)
    tok = AutoTokenizer.from_pretrained(a.model)
    model = EsmForMaskedLM.from_pretrained(a.model, torch_dtype=torch.float16 if dev == "cuda" else torch.float32).to(dev).eval()
    ids = tok(seq, return_tensors="pt")["input_ids"][0]           # <cls> seq <eos>
    aa_ids = [tok.convert_tokens_to_ids(x) for x in AAS]
    rows = []
    with torch.no_grad():
        for start in range(0, len(seq), a.batch):
            pos = list(range(start, min(start + a.batch, len(seq))))
            batch = ids.repeat(len(pos), 1)
            for k, p in enumerate(pos):
                batch[k, p + 1] = tok.mask_token_id
            logits = model(batch.to(dev)).logits.float()
            for k, p in enumerate(pos):
                lp = torch.log_softmax(logits[k, p + 1, aa_ids], dim=-1).cpu()
                ent = float(-(lp.exp() * lp).sum())
                wt = seq[p]
                lwt = float(lp[AAS.index(wt)])
                for j, m in enumerate(AAS):
                    if m != wt:
                        rows.append((p + 1, wt, m, float(lp[j]) - lwt, ent))
            print(f"  {pos[-1] + 1}/{len(seq)}", end="\r", flush=True)
    out = pd.DataFrame(rows, columns=["pos", "wt", "mut", "esm2_llr", "esm2_entropy"])
    p = ROOT / "data" / "processed" / "esm2_llr.csv"
    p.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(p, index=False)
    r = out[(out.pos == 175) & (out.mut == "H")]
    print(f"\nkaydedildi: {p} ({len(out)} satir). R175H LLR = {float(r.esm2_llr.iloc[0]):.3f}")


if __name__ == "__main__":
    sys.exit(main())
