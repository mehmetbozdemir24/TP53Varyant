"""VCEP motorunu ClinGen eRepo uzman siniflamalariyla karsilastirir (otomatik kriter bazinda da)."""
import sys
from collections import Counter
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd
from tp53.acmg import TP53Classifier, Clinical, parse_code
from tp53.common import PROC

df = pd.read_parquet(PROC / "variants.parquet")
clf = TP53Classifier(df)
AUTO = ["PVS1", "PS1", "PS3", "BS3", "PM1", "PM2", "PM5", "PP3", "BP4", "BP7", "BA1", "BS1"]
rows, code_cmp = [], Counter()
for c, r in clf.erepo.iterrows():
    res = clf.vcep(c, Clinical(use_expert=True))
    ours = {parse_code(l)[0]: l for l in res["met"]}
    exp = {parse_code(l)[0]: l for l in (r["erepo_met"] or "").split(";") if l}
    for k in AUTO:
        a, b = ours.get(k), exp.get(k)
        if a == b:
            code_cmp[(k, "ayni")] += 1
        else:
            code_cmp[(k, "farkli")] += 1
            rows.append(dict(c=c, p1=r["p1"], code=k, bizim=a, uzman=b, ver=r["erepo_version"]))
    rows.append(dict(c=c, p1=r["p1"], code="SINIF", bizim=res["classification"], uzman=r["erepo_class"],
                     bizim_puan=res["points"], uzman_puan=r["erepo_points"], ver=r["erepo_version"]))
out = pd.DataFrame(rows)
cls = out[out.code == "SINIF"]
agree = (cls.bizim == cls.uzman).mean()
print(f"Sinif uyumu: {agree:.1%} ({(cls.bizim == cls.uzman).sum()}/{len(cls)})")
print(pd.crosstab(cls.uzman, cls.bizim))
print("\nKriter bazinda uyum (otomatik kriterler):")
for k in AUTO:
    s, d = code_cmp[(k, 'ayni')], code_cmp[(k, 'farkli')]
    print(f"  {k:5s} ayni={s:3d} farkli={d:3d}  ({s/(s+d):.0%})")
out.to_csv(PROC / "validation_vcep_vs_erepo.csv", index=False)
print("\nFarkli siniflar:")
print(cls[cls.bizim != cls.uzman][["c", "p1", "uzman", "bizim", "uzman_puan", "bizim_puan", "ver"]].to_string(index=False))
