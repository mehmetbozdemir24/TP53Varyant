"""Makine ogrenmesi modulleri.

Model A - "TP53 Aggregated Prediction": Franklin'in toplu tahmin skorunun TP53'e ozel karsiligi.
          Girdi: yalnizca in-silico skorlar (+ istege bagli ESM-2 LLR). Etiket: ClinGen eRepo uzman
          siniflamasi + ClinVar (>= 1 yildiz, celiskisiz) P/LP vs B/LB missense varyantlar.
Model B - "Fonksiyon / hipomorf": Kato 2003 transaktivasyon sinifini (Functional / Partially functional /
          Non-functional) tahmin eder. 'Partially functional' olasiligi hipomorf aday skoru olarak kullanilir.
          B1: sekans + yapi + in-silico ozellikleri;  B2: B1 + diger fonksiyonel testler (Giacomelli, Funk, Kotler).

Dogrulama: kodona gore gruplu 5 katli capraz dogrulama (ayni kodondaki varyantlar ayni katta).
Calistirma: python -m tp53.ml
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (average_precision_score, balanced_accuracy_score, confusion_matrix, f1_score,
                             roc_auc_score)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .common import CHARGE, DOMAINS, HYDRO, MODELS, PROC, VOLUME

INSILICO = ["revel", "alphamissense", "bayesdel", "varity_r", "varity_er", "sift", "sift4g", "polyphen2_hdiv",
            "polyphen2_hvar", "primateai", "fathmm", "fathmm_xf", "mutationtaster", "mutationassessor", "provean",
            "dann", "metalr", "metasvm", "metarnn", "clinpred", "vest4", "mvp", "mpc", "esm1b", "eve", "gmvp",
            "deogen2", "cadd_phred", "gerp", "phylop100", "phastcons100", "lrt", "agvgd_num"]
# PS3/BS3 icin yayimlanmis esik veya gosterim icin kullanilan araclar (karsilastirma tablosu)
SINGLE_TOOLS = ["revel", "alphamissense", "bayesdel", "varity_r", "cadd_phred", "esm1b", "primateai", "eve",
                "polyphen2_hvar", "sift"]
LOWER_IS_DAMAGING = {"sift", "sift4g", "fathmm", "provean", "esm1b", "lrt", "esm2_llr"}


def load():
    df = pd.read_parquet(PROC / "variants.parquet")
    m = df[df["consequence"] == "missense"].copy()
    for c in INSILICO[:-1] + ["grantham", "codon", "hotspot_aa_count", "nci_tumor_n", "tcga_icgc_genie_count"]:
        if c in m:
            m[c] = pd.to_numeric(m[c], errors="coerce")
    m["agvgd_num"] = pd.to_numeric(m["agvgd"].astype(str).str.extract(r"(\d+)")[0], errors="coerce")
    esm = PROC / "esm2_llr.csv"
    if esm.exists():
        e = pd.read_csv(esm)
        e["p1"] = e["wt"] + e["pos"].astype(str) + e["mut"]
        m = m.merge(e[["p1", "esm2_llr", "esm2_entropy"]], on="p1", how="left")
    return df, m


def structural(m):
    X = pd.DataFrame(index=m.index)
    X["codon"] = m["codon"]
    for name, a, b in DOMAINS:
        X[f"dom_{name[:6]}"] = ((m["codon"] >= a) & (m["codon"] <= b)).astype(int)
    for k, tab in (("hyd", HYDRO), ("vol", VOLUME), ("chg", CHARGE)):
        w = m["aa_ref"].map(lambda a: tab.get(a, 0.0)); u = m["aa_alt"].map(lambda a: tab.get(a, 0.0))
        X[f"{k}_wt"], X[f"{k}_mut"], X[f"{k}_delta"] = w, u, u - w
    X["grantham"] = m["grantham"]
    X["to_pro"] = (m["aa_alt"] == "P").astype(int)
    X["from_gly"] = (m["aa_ref"] == "G").astype(int)
    X["to_gly"] = (m["aa_alt"] == "G").astype(int)
    X["cys_change"] = ((m["aa_ref"] == "C") | (m["aa_alt"] == "C")).astype(int)
    # cinko baglayan ve DNA temas kalintilari (UniProt/Cho 1994): C176, H179, C238, C242, K120, S241, R248, R273, A276, C277, R280
    X["zinc_site"] = m["codon"].isin([176, 179, 238, 242]).astype(int)
    X["dna_contact"] = m["codon"].isin([120, 241, 248, 273, 276, 277, 280, 283]).astype(int)
    return X


def features(m, kind):
    cols = [c for c in INSILICO if c in m]
    if "esm2_llr" in m:
        cols += ["esm2_llr", "esm2_entropy"]
    X = m[cols].copy()
    if kind in ("B1", "B2"):
        X = X.join(structural(m))
    if kind == "B2":
        for a in ("giacomelli", "funk", "kotler"):
            v = m[a].astype(str)
            X[f"{a}_lof"] = np.where(v.str.upper().str.startswith("LOF"), 1.0, np.where(v.str.lower().str.startswith("nolof"), 0.0, np.nan))
        v = m["kawaguchi"].astype(str).str.lower()
        X["kawaguchi_abn"] = np.where(v.str.contains("abn|mono|dim"), 1.0, np.where(v.str.contains("norm|tetra"), 0.0, np.nan))
    return X


# ------------------------------------------------------------------ etiketler
def path_labels(m):
    """1 = P/LP, 0 = B/LB. Oncelik: ClinGen eRepo > ClinVar (>= 1 yildiz, celiskisiz)."""
    y = pd.Series(np.nan, index=m.index)
    src = pd.Series(None, index=m.index, dtype=object)
    cv = m["clinvar_germline"].fillna("")
    stars = pd.to_numeric(m["clinvar_stars"], errors="coerce").fillna(0)
    ok = stars >= 1
    y[ok & cv.str.contains("athogenic") & ~cv.str.contains("onflict|enign")] = 1
    y[ok & cv.str.contains("enign") & ~cv.str.contains("onflict|athogenic")] = 0
    src[y.notna()] = "ClinVar"
    er = m["erepo_class"].fillna("")
    y[er.str.contains("athogenic")] = 1
    y[er.str.contains("enign")] = 0
    src[er.str.contains("athogenic|enign")] = "ClinGen eRepo"
    return y, src


def kato_labels(m):
    mp = {"Functional": 0, "Partially functional": 1, "Non-functional": 2}
    return m["kato"].map(mp)


def _cv(model_fn, X, y, groups, n_splits=5, seed=42):
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    oof = np.zeros((len(y), len(np.unique(y)))) if len(np.unique(y)) > 2 else np.zeros(len(y))
    for tr, te in cv.split(X, y, groups):
        mdl = model_fn()
        mdl.fit(X.iloc[tr], y.iloc[tr])
        p = mdl.predict_proba(X.iloc[te])
        oof[te] = p if oof.ndim == 2 else p[:, 1]
    return oof


def lr():
    return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), LogisticRegression(C=0.5, max_iter=3000))


def hgb(classes=2):
    return HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, l2_regularization=1.0,
                                          min_samples_leaf=10, random_state=0)


# ------------------------------------------------------------------ Model A
def train_pathogenicity(m):
    y, src = path_labels(m)
    lab = y.notna()
    X = features(m, "A")
    keep = X.columns[X[lab].notna().mean() >= 0.5]       # cok eksik araclari at
    X = X[keep]
    Xl, yl, g = X[lab], y[lab].astype(int), m.loc[lab, "codon"]
    res = {"n": int(lab.sum()), "n_pos": int(yl.sum()), "n_neg": int((1 - yl).sum()),
           "label_sources": src[lab].value_counts().to_dict(), "features": list(keep), "models": {}, "single_tools": {}}
    best, best_auc = None, -1
    for name, fn in (("logistic", lr), ("gboost", hgb)):
        aucs = []
        for seed in (1, 2, 3):
            oof = _cv(fn, Xl, yl, g, seed=seed)
            aucs.append((roc_auc_score(yl, oof), average_precision_score(yl, oof)))
        a = np.mean(aucs, axis=0)
        res["models"][name] = {"roc_auc": round(float(a[0]), 4), "pr_auc": round(float(a[1]), 4)}
        if a[0] > best_auc:
            best, best_auc = name, a[0]
    for t in SINGLE_TOOLS + (["esm2_llr"] if "esm2_llr" in m else []):
        if t not in m:
            continue
        s = pd.to_numeric(m.loc[lab, t], errors="coerce")
        ok = s.notna()
        if ok.sum() < 30:
            continue
        sc = -s[ok] if t in LOWER_IS_DAMAGING else s[ok]
        res["single_tools"][t] = {"roc_auc": round(float(roc_auc_score(yl[ok], sc)), 4), "coverage": round(float(ok.mean()), 3)}
    fn = lr if best == "logistic" else hgb
    final = CalibratedClassifierCV(fn(), method="isotonic", cv=5)
    final.fit(Xl, yl)
    score = final.predict_proba(X)[:, 1]
    res["best"] = best
    return final, keep, score, res


# ------------------------------------------------------------------ Model B
def train_function(m, kind="B1"):
    y = kato_labels(m)
    lab = y.notna()
    X = features(m, kind)
    keep = X.columns[X[lab].notna().mean() >= 0.3]
    X = X[keep]
    Xl, yl, g = X[lab], y[lab].astype(int), m.loc[lab, "codon"]
    oof = _cv(lambda: hgb(3), Xl, yl, g)
    pred = oof.argmax(1)
    res = {"kind": kind, "n": int(lab.sum()), "class_counts": yl.value_counts().sort_index().to_dict(),
           "macro_f1": round(float(f1_score(yl, pred, average="macro")), 4),
           "balanced_acc": round(float(balanced_accuracy_score(yl, pred)), 4),
           "auc_partial_vs_rest": round(float(roc_auc_score((yl == 1).astype(int), oof[:, 1])), 4),
           "auc_nonfunc_vs_rest": round(float(roc_auc_score((yl == 2).astype(int), oof[:, 2])), 4),
           "confusion": confusion_matrix(yl, pred).tolist(), "features": list(keep)}
    final = hgb(3).fit(Xl, yl)
    proba = final.predict_proba(X)
    # egitimdeki varyantlar icin sizinti olmamasi adina CV disi (oof) tahmini kullan
    proba[np.where(lab)[0]] = oof
    return final, keep, proba, res


def run():
    import joblib
    df, m = load()
    out = m[["c", "p1"]].copy()
    mdlA, colsA, sA, rA = train_pathogenicity(m)
    out["ml_path_score"] = sA
    y, _ = path_labels(m)
    # egitim varyantlari icin capraz dogrulama disi skor
    lab = y.notna()
    XA = features(m, "A")[colsA]
    oofA = _cv(lr if rA["best"] == "logistic" else hgb, XA[lab], y[lab].astype(int), m.loc[lab, "codon"])
    out.loc[lab, "ml_path_score"] = oofA
    out["ml_path_pct"] = out["ml_path_score"].rank(pct=True)
    metrics = {"A_pathogenicity": rA}
    for kind in ("B1", "B2"):
        mdl, cols, P, r = train_function(m, kind)
        out[f"ml_{kind}_functional"], out[f"ml_{kind}_partial"], out[f"ml_{kind}_nonfunctional"] = P[:, 0], P[:, 1], P[:, 2]
        metrics[f"{kind}_function"] = r
        joblib.dump({"model": mdl, "cols": list(cols)}, MODELS / f"function_{kind}.joblib")
    joblib.dump({"model": mdlA, "cols": list(colsA)}, MODELS / "pathogenicity_A.joblib")
    out.to_parquet(PROC / "ml_predictions.parquet", index=False)
    (MODELS / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk not in ("features", "confusion")} for k, v in metrics.items()},
                     indent=1))
    return metrics


if __name__ == "__main__":
    run()
