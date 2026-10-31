"""Otomatik VCEP ve Franklin-benzeri modu, ClinVar'da >= 2 yildizli (eRepo disi) kayitlarla karsilastirir."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd
from tp53.acmg import TP53Classifier, Clinical
from tp53.common import PROC

clf = TP53Classifier(pd.read_parquet(PROC / "variants.parquet"))
d = clf.df
stars = pd.to_numeric(d.clinvar_stars, errors="coerce").fillna(0)
g = d.clinvar_germline.fillna("")
sel = d[(stars >= 2) & d.erepo_class.isna() & g.str.contains("athogenic|enign") & ~g.str.contains("onflict")]
grp = lambda x: "P/LP" if "athogenic" in x and "enign" not in x else ("B/LB" if "enign" in x and "athogenic" not in x else "VUS")
rows = []
for c, r in sel.iterrows():
    truth = grp(r.clinvar_germline)
    if truth == "VUS":
        continue
    rows.append(dict(c=c, p1=r.p1, consequence=r.consequence, clinvar=truth,
                     vcep_auto=grp(clf.vcep(c, Clinical(use_expert=False))["classification"]),
                     franklin_mode=grp(clf.franklin(c)["classification"])))
x = pd.DataFrame(rows)
print(f"{len(x)} varyant (ClinVar >= 2 yildiz, eRepo disi)")
for col in ("vcep_auto", "franklin_mode"):
    print(f"\n{col}:\n", pd.crosstab(x.clinvar, x[col]))
    wrong = (((x.clinvar == "P/LP") & (x[col] == "B/LB")) | ((x.clinvar == "B/LB") & (x[col] == "P/LP"))).sum()
    print(f"ayni yon: {(x.clinvar == x[col]).mean():.1%}  ters yon hata: {wrong}")
x.to_csv(PROC / "validation_clinvar2star.csv", index=False)
