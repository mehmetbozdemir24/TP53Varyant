"""Ham veri setlerinden birlesik TP53 varyant tablosunu kurar.

Cikti: data/processed/variants.parquet   (her satir bir varyant, anahtar: NM_000546.6 c. notasyonu)
       data/processed/gene.json          (CDS dizisi, protein, domainler, hotspotlar, ozet sayimlar)
       data/processed/erepo.parquet      (ClinGen uzman siniflamalari, kod listeleri)
       data/processed/nci_function.parquet, nci_tumor_summary.parquet, civic.parquet

Calistirma: python -m tp53.build
"""
from __future__ import annotations

import json
import math
import re
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from scipy.stats import chi2

from .common import (AA1, AA3, COMP, ENST, EXONS, HOTSPOT_CODONS, PROC, RAW, consequence, domain_of,
                     exon_of_cpos, grantham, is_nan, norm_c, p3_to_p1, parse_cpos, MISSENSE_RE)

FOUNDER_POPS = {"asj", "fin", "ami", "mid", "remaining"}   # VCEP: frekans kurallarinda dislanir
NUM = lambda x: (None if x is None or x == "" or (isinstance(x, float) and math.isnan(x)) else float(x))


def log(*a):
    print("[build]", *a, flush=True)


# =================================================================== kaynak okuyucular
def read_vcep():
    s2 = pd.read_excel(RAW / "vcep_TableS2.xlsx", header=1).iloc[1:]
    s2.columns = ["c", "p3", "agvgd", "bayesdel_vcep", "s2_code", "s2_spliceai", "s2_splice_note"][:len(s2.columns)]
    s2 = s2[s2["c"].astype(str).str.startswith("c.")].copy()
    s2["c"] = s2["c"].map(norm_c)
    s2["agvgd"] = s2["agvgd"].astype(str).str.replace("Class ", "", regex=False)
    s3 = pd.read_excel(RAW / "vcep_TableS3.xlsx", header=1).iloc[1:]
    s3.columns = ["p1", "s3_seq_change", "kato", "funk", "giacomelli", "kotler", "kawaguchi", "s3_code"]
    s3 = s3[s3["p1"].astype(str).str.match(r"^[A-Z]\d+[A-Z]$")].copy()
    s1 = pd.read_excel(RAW / "vcep_TableS1.xlsx", header=1)
    s1 = s1.rename(columns={s1.columns[2]: "c", s1.columns[15]: "s1_pvs1", s1.columns[16]: "s1_rationale",
                            s1.columns[11]: "s1_pred1", s1.columns[13]: "s1_rna"})
    s1 = s1[s1["c"].astype(str).str.startswith("c.")][["c", "s1_pvs1", "s1_rationale", "s1_pred1", "s1_rna"]]
    s1["c"] = s1["c"].map(norm_c)
    # Tablo her donor/akseptor blogu (WT satiri + 6 varyant) icin gucu birlesik hucrede veriyor:
    # blok = bir 'WT donor/acceptor' satirindan sonrakine kadar.
    raw = pd.read_excel(RAW / "vcep_TableS1.xlsx", header=1)
    raw["block"] = raw.iloc[:, 2].astype(str).str.startswith("WT").cumsum()
    for i, col in ((15, "s1_pvs1"), (16, "s1_rationale"), (13, "s1_rna")):
        raw[col] = raw.groupby("block")[raw.columns[i]].transform(lambda s: s.ffill().bfill())
    raw["c"] = raw.iloc[:, 2].map(norm_c)
    s1 = s1.drop(columns=["s1_pvs1", "s1_rationale", "s1_rna"]).merge(
        raw.dropna(subset=["c"])[["c", "s1_pvs1", "s1_rationale", "s1_rna"]], on="c", how="left")
    log(f"VCEP S1={len(s1)} S2={len(s2)} S3={len(s3)}")
    return s1, s2, s3


def _agg(v, how="max"):
    """dbNSFP transkript-bazli listelerini tek degere indirger (yuksek=zararli ise max, dusuk=zararli ise min)."""
    if v is None:
        return None
    if isinstance(v, (list, tuple)):
        vals = [float(x) for x in v if x not in (None, ".", "") and not (isinstance(x, float) and math.isnan(x))]
        if not vals:
            return None
        return max(vals) if how == "max" else min(vals)
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _mt_oriented(m):
    """MutationTaster 'score' tahminin dogruluk olasiligidir; D/A (hastalik) ise oldugu gibi, N/P (polimorfizm) ise
    1-skor alinarak 'zararlilik olasiligina' cevrilir. Transkriptler arasi en yuksek deger dondurulur."""
    if not isinstance(m, dict) or m.get("score") is None:
        return None
    sc, pr = m["score"], m.get("pred")
    sc = sc if isinstance(sc, list) else [sc]
    pr = pr if isinstance(pr, list) else [pr] * len(sc)
    vals = []
    for s, p in zip(sc, pr):
        try:
            s = float(s)
        except (TypeError, ValueError):
            continue
        vals.append(s if p in ("D", "A") else 1 - s if p in ("N", "P") else None)
    vals = [v for v in vals if v is not None]
    return max(vals) if vals else None


def _g(d, *path):
    for p in path:
        if not isinstance(d, dict):
            return None
        d = d.get(p)
    return d


def read_dbnsfp():
    """MyVariant.info: dbNSFP in-silico skorlari + CADD + dbscSNV + gnomAD v2/ExAC (hg19 _id)."""
    mv = json.load(open(RAW / "myvariant_TP53_dbnsfp.json", encoding="utf-8"))
    rows = []
    for h in mv:
        d = h.get("dbnsfp") or {}
        m = re.match(r"chr17:g\.(\d+)([ACGT])>([ACGT])", h["_id"])
        if not m:
            continue
        hg38 = _g(d, "hg38", "start")
        cadd = h.get("cadd")
        if isinstance(cadd, list):
            cadd = cadd[0]
        dbsc = h.get("dbscsnv") or {}
        if isinstance(dbsc, list):
            dbsc = dbsc[0]
        ge = h.get("gnomad_exome") or {}
        gg = h.get("gnomad_genome") or {}
        aa = d.get("aa") or {}
        if isinstance(aa, list):
            aa = aa[0] if aa else {}
        rows.append(dict(
            hg19_pos=int(m.group(1)), g_ref=m.group(2), g_alt=m.group(3), hg38_pos=hg38,
            rsid=_g(h, "dbsnp", "rsid") or d.get("rsid"),
            aa_ref=aa.get("ref"), aa_alt=aa.get("alt"),
            revel=_agg(_g(d, "revel", "score")),
            alphamissense=_agg(_g(d, "alphamissense", "score")),
            bayesdel=_agg(_g(d, "bayesdel", "no_af", "score")),
            varity_r=_agg(_g(d, "varity", "r", "score")),
            varity_er=_agg(_g(d, "varity", "er", "score")),
            sift=_agg(_g(d, "sift", "score"), "min"),
            sift4g=_agg(_g(d, "sift4g", "score"), "min"),
            polyphen2_hdiv=_agg(_g(d, "polyphen2", "hdiv", "score")),
            polyphen2_hvar=_agg(_g(d, "polyphen2", "hvar", "score")),
            primateai=_agg(_g(d, "primateai", "score")),
            fathmm=_agg(_g(d, "fathmm", "score"), "min"),
            fathmm_xf=_agg(_g(d, "fathmm-xf", "coding_score")),
            mutationtaster=_mt_oriented(_g(d, "mutationtaster")),
            mutationassessor=_agg(_g(d, "mutationassessor", "score")),
            provean=_agg(_g(d, "provean", "score"), "min"),
            dann=_agg(_g(d, "dann", "score")),
            metalr=_agg(_g(d, "metalr", "score")),
            metasvm=_agg(_g(d, "metasvm", "score")),
            metarnn=_agg(_g(d, "metarnn", "score")),
            clinpred=_agg(_g(d, "clinpred", "score")),
            vest4=_agg(_g(d, "vest4", "score")),
            mvp=_agg(_g(d, "mvp", "score")),
            mpc=_agg(_g(d, "mpc", "score")),
            esm1b=_agg(_g(d, "esm1b", "score"), "min"),
            eve=_agg(_g(d, "eve", "score")),
            gmvp=_agg(_g(d, "gmvp", "score")),
            deogen2=_agg(_g(d, "deogen2", "score")),
            mutformer=_agg(_g(d, "mutformer", "score")),
            phactboost=_agg(_g(d, "phactboost", "score")),
            mutpred=_agg(_g(d, "mutpred", "score")),
            lrt=_agg(_g(d, "lrt", "score"), "min"),
            gerp=_agg(_g(d, "gerp++", "rs")),
            phylop100=_agg(_g(d, "phylop", "100way_vertebrate", "score")),
            phastcons100=_agg(_g(d, "phastcons", "100way_vertebrate", "score")),
            genocanyon=_agg(_g(d, "genocanyon", "score")),
            fitcons=_agg(_g(d, "fitcons", "integrated", "score")),
            mutpred_mech="; ".join(sorted({x.get("mechanism", "") for x in (_g(d, "mutpred", "pred") or [])
                                          if isinstance(x, dict) and x.get("p_val", 1) < 0.05}))[:300] or None,
            cadd_phred=_agg(_g(cadd, "phred")) if isinstance(cadd, dict) else None,
            dbscsnv_ada=_agg(dbsc.get("ada_score")) if isinstance(dbsc, dict) else None,
            dbscsnv_rf=_agg(dbsc.get("rf_score")) if isinstance(dbsc, dict) else None,
            gnomad2_exome_af=_agg(_g(ge, "af", "af")), gnomad2_genome_af=_agg(_g(gg, "af", "af")),
            exac_af=_agg(_g(d, "exac", "adj_af")),
        ))
    df = pd.DataFrame(rows)
    log(f"dbNSFP {len(df)} SNV")
    return df


def faf95(ac, an):
    """gnomAD filtreleme alel frekansi (Whiffin 2017): Poisson %95 tek yonlu alt guven siniri / AN."""
    if not an or not ac:
        return 0.0
    return float(chi2.ppf(0.05, 2 * ac) / 2.0 / an)


def read_gnomad():
    g = json.load(open(RAW / "gnomad_v4_TP53.json", encoding="utf-8"))
    rows = []
    for v in g:
        pops = defaultdict(lambda: [0, 0])
        ac = an = 0
        filt = []
        for src in ("exome", "genome"):
            s = v.get(src)
            if not s:
                continue
            ac += s.get("ac") or 0
            an += s.get("an") or 0
            filt += [f"{src}:{f}" for f in (s.get("filters") or [])]
            for p in s.get("populations") or []:
                pid = p["id"]
                if "_" in pid or pid in ("XX", "XY"):
                    continue
                pops[pid][0] += p.get("ac") or 0
                pops[pid][1] += p.get("an") or 0
        best = (0.0, None)
        pop_tab = {}
        for pid, (pac, pan) in pops.items():
            pop_tab[pid] = {"ac": pac, "an": pan, "af": (pac / pan if pan else 0.0)}
            if pid in FOUNDER_POPS or pan < 2000 or pac < 2:
                continue
            f = faf95(pac, pan)
            if f > best[0]:
                best = (f, pid)
        # PM2 icin: kurucu olmayan ve birden fazla alel gorulen gruplarin en yuksek AF'si
        multi = [d["af"] for pid, d in pop_tab.items() if pid not in FOUNDER_POPS and d["ac"] > 1]
        rows.append(dict(
            gnomad_id=v["variant_id"], hg38_pos=v["pos"], g_ref=v["ref"], g_alt=v["alt"],
            c_gnomad=norm_c(v.get("hgvsc")) if v.get("transcript_id") == ENST else None,
            p_gnomad=v.get("hgvsp") if v.get("transcript_id") == ENST else None,
            gnomad_csq=v.get("consequence"),
            gnomad_ac=ac, gnomad_an=an, gnomad_af=(ac / an if an else 0.0),
            gnomad_hom=sum((v.get(s) or {}).get("ac_hom") or 0 for s in ("exome", "genome")),
            gnomad_faf95_popmax=best[0], gnomad_faf95_pop=best[1],
            gnomad_max_pop_af_multi=max(multi) if multi else 0.0,
            gnomad_filters=";".join(filt) or None, gnomad_flags=";".join(v.get("flags") or []) or None,
            gnomad_pops=json.dumps(pop_tab),
        ))
    df = pd.DataFrame(rows)
    log(f"gnomAD v4 {len(df)} varyant")
    return df


def read_nci():
    """NCI TP53 Database r21: somatik tumor, germline aile ve fonksiyonel test kayitlari."""
    t = pd.read_csv(RAW / "nci_TumorVariantDownload_r21.csv", low_memory=False, encoding="latin1")
    t["c"] = t["c_description"].map(norm_c)
    t = t[t["c"].notna()]
    sp = ["SpliceAI_DS_AG", "SpliceAI_DS_AL", "SpliceAI_DS_DG", "SpliceAI_DS_DL"]
    for c in sp:
        t[c] = pd.to_numeric(t[c], errors="coerce")
    t["spliceai_max"] = t[sp].max(axis=1)
    first = ["MUT_ID", "hg19_Chr17_coordinates", "hg38_Chr17_coordinates", "ExonIntron", "ProtDescription",
             "Effect", "Type", "CpG_site", "WT_nucleotide", "Mutant_nucleotide", "AGVGDClass", "SIFTClass",
             "Polyphen2", "TransactivationClass", "DNE_LOFclass", "DNEclass", "StructureFunctionClass", "Hotspot",
             "Structural_motif", "TCGA_ICGC_GENIE_count", "COSMIClink", "CLINVARlink", "spliceai_max",
             "g_description_GRCh38"]
    agg = t.groupby("c").agg(**{k: (k, "first") for k in first}, nci_tumor_n=("c", "size"))
    agg = agg.rename(columns={"MUT_ID": "nci_mut_id", "ExonIntron": "exon_intron", "ProtDescription": "p_nci",
                              "Effect": "nci_effect", "AGVGDClass": "agvgd_nci", "TransactivationClass": "kato_nci",
                              "DNE_LOFclass": "dne_lof_class", "DNEclass": "dne_class",
                              "StructureFunctionClass": "structure_function_class", "Hotspot": "nci_hotspot",
                              "Structural_motif": "structural_motif", "TCGA_ICGC_GENIE_count": "tcga_icgc_genie_count",
                              "COSMIClink": "cosmic_id", "CLINVARlink": "clinvar_id_nci", "Type": "mut_type",
                              "CpG_site": "cpg", "SIFTClass": "sift_class_nci", "Polyphen2": "polyphen_class_nci"})
    topo = (t.groupby(["c", "Short_topo"]).size().rename("n").reset_index()
              .sort_values(["c", "n"], ascending=[True, False]))
    topo.to_parquet(PROC / "nci_tumor_topography.parquet", index=False)

    gl = pd.read_csv(RAW / "nci_GermlineDownload_r21.csv", low_memory=False, encoding="latin1")
    gl["c"] = gl["c_description"].map(norm_c)
    gl = gl[gl["c"].notna()]
    gagg = gl.groupby("c").agg(nci_germ_families=("Family_ID", "nunique"),
                               nci_germ_individuals=("Individual_ID", "nunique"),
                               nci_germ_classes=("Class", lambda s: "; ".join(f"{k} ({v})" for k, v in Counter(s.dropna()).most_common())))
    # Li-Fraumeni fenotipini karsilayan aileler (LFS, LFL, Chompret); FH / noFH / Other sayilmaz
    lfs_fam = gl[gl["Class"].isin(["LFS", "LFL", "TP53_Chompret"])].groupby("c")["Family_ID"].nunique()
    gagg["nci_germ_lfs_families"] = lfs_fam.reindex(gagg.index).fillna(0)
    keep = ["c", "Family_ID", "Class", "Country", "Individual_ID", "Sex", "Germline_carrier", "Topography",
            "Short_topo", "Morphology", "Age_at_diagnosis", "Ref_ID", "PubMedLink"]
    gl[[k for k in keep if k in gl.columns]].to_parquet(PROC / "nci_germline.parquet", index=False)

    fn = pd.read_csv(RAW / "nci_FunctionDownload_r21.csv", low_memory=False, encoding="latin1")
    fn = fn.rename(columns={"AAchange": "p1"})
    fn.to_parquet(PROC / "nci_function.parquet", index=False)
    fcount = fn.groupby("p1").size().rename("nci_functional_assays")
    log(f"NCI tumor={len(t)} satir/{len(agg)} varyant, germline={len(gl)} satir/{len(gagg)} varyant, fonksiyon={len(fn)}")
    return agg.join(gagg, how="outer"), fcount


CV_RANK = {"practice guideline": 4, "reviewed by expert panel": 3,
           "criteria provided, multiple submitters, no conflicts": 2,
           "criteria provided, conflicting classifications": 1, "criteria provided, single submitter": 1}


def read_clinvar():
    docs = json.load(open(RAW / "clinvar_TP53_esummary.json", encoding="utf-8"))
    rows = []
    for d in docs:
        title = d.get("title", "")
        m = re.match(r"^NM_000546\.\d+\(TP53\):(c\.[^ ]+)(?: \((p\.[^)]+)\))?", title)
        if not m:
            continue
        gc = d.get("germline_classification") or {}
        oc = d.get("oncogenicity_classification") or {}
        ci = d.get("clinical_impact_classification") or {}
        traits = [t.get("trait_name") for t in gc.get("trait_set") or []]
        vs = (d.get("variation_set") or [{}])[0]
        spdi = vs.get("canonical_spdi") or ""
        rows.append(dict(
            c=norm_c(m.group(1)), p_clinvar=m.group(2), clinvar_vcv=d.get("accession"), clinvar_id=d.get("uid"),
            clinvar_germline=gc.get("description") or None, clinvar_review=gc.get("review_status") or None,
            clinvar_stars=CV_RANK.get(gc.get("review_status") or "", 0),
            clinvar_last_eval=(gc.get("last_evaluated") or "")[:10] or None,
            clinvar_traits="; ".join(t for t in traits if t)[:500] or None,
            clinvar_n_scv=len((d.get("supporting_submissions") or {}).get("scv") or []),
            clinvar_oncogenicity=oc.get("description") or None,
            clinvar_somatic_tier=ci.get("description") or None,
            clinvar_spdi=spdi or None,
            clinvar_mc=", ".join(d.get("molecular_consequence_list") or []) or None,
        ))
    df = pd.DataFrame(rows).sort_values("clinvar_stars", ascending=False).drop_duplicates("c")
    log(f"ClinVar {len(df)} TP53 NM_000546 varyant")
    return df.set_index("c")


def code_points(code):
    """VCEP v2.4 Bayes puanlari: 'PS3_Moderate' -> 2, 'BS1' -> -4, 'BA1' -> -8 (tek basina benign)."""
    code = code.strip()
    base, _, mod = code.partition("_")
    mod = mod.lower()
    if base == "BA1":
        return -8
    if base.startswith("PVS"):
        dflt = 8
    elif base.startswith("PS"):
        dflt = 4
    elif base.startswith("PM"):
        dflt = 2
    elif base.startswith("PP"):
        dflt = 1
    elif base.startswith("BS"):
        dflt = -4
    elif base.startswith("BP"):
        dflt = -1
    else:
        return 0
    sign = 1 if dflt > 0 else -1
    lev = {"verystrong": 8, "very strong": 8, "very_strong": 8, "strong": 4, "moderate": 2, "supporting": 1}
    for k, v in lev.items():
        if mod.replace("_", " ").startswith(k.replace("_", " ")):
            return sign * v
    return dflt


def read_erepo():
    e = json.load(open(RAW / "erepo_TP53.json", encoding="utf-8"))["variantInterpretations"]
    rows = []
    for v in e:
        c = next((norm_c(h) for h in v["hgvs"] if h.startswith("NM_000546.")), None)
        gl = v["guidelines"][0]
        codes = [x for a in gl.get("agents", []) for x in a.get("evidenceCodes", [])]
        met = [x["label"] for x in codes if x.get("status") == "Met"]
        rows.append(dict(c=c, erepo_class=gl["outcome"]["label"], erepo_version=gl.get("version"),
                         erepo_date=v.get("publishedDate"), erepo_uuid=v.get("uuid"), erepo_caid=v.get("caid"),
                         erepo_met=";".join(met), erepo_points=sum(code_points(x) for x in met),
                         erepo_not_met=";".join(x["label"] for x in codes if x.get("status") != "Met"),
                         erepo_url=f"https://erepo.clinicalgenome.org/evrepo/ui/classification/{v.get('uuid')}"))
    df = pd.DataFrame(rows).dropna(subset=["c"]).drop_duplicates("c")
    df.to_parquet(PROC / "erepo.parquet", index=False)
    log(f"eRepo {len(df)} uzman siniflamasi")
    return df.set_index("c")


def read_uniprot():
    u = json.load(open(RAW / "uniprot_P04637.json", encoding="utf-8"))
    seq = u["sequence"]["value"]
    feats = []
    for f in u["features"]:
        if f["type"] in ("Region", "Motif", "DNA binding", "Binding site", "Site", "Modified residue", "Helix",
                         "Beta strand", "Turn", "Cross-link", "Mutagenesis"):
            feats.append(dict(type=f["type"], start=f["location"]["start"]["value"], end=f["location"]["end"]["value"],
                              description=f.get("description", "") or ", ".join(
                                  (x.get("name") or "") for x in [f.get("ligand") or {}])))
    # UniProt'un kendi (kurate) aciklamasi: "in LFS; germline mutation ..." gibi
    curated = {}
    for f in u["features"]:
        if f["type"] == "Natural variant" and f.get("featureId") and f.get("alternativeSequence"):
            alt = (f["alternativeSequence"].get("alternativeSequences") or [""])[0]
            p1 = f"{f['alternativeSequence'].get('originalSequence', '')}{f['location']['start']['value']}{alt}"
            curated[p1] = f.get("description") or ""
    nat = defaultdict(list)
    vv = json.load(open(RAW / "uniprot_variation_P04637.json", encoding="utf-8"))
    for f in vv["features"]:
        if f.get("begin") != f.get("end") or not f.get("alternativeSequence"):
            continue
        p1 = f"{f.get('wildType', '')}{f['begin']}{f['alternativeSequence']}"
        # Proteins API, ClinVar/Ensembl'den derlenmis klinik anlam ekler; kanit kriteri olmayan (0 yildiz) kayitlar atlanir
        sig = sorted({s.get("type") for s in f.get("clinicalSignificances") or []
                      if s.get("type") and s.get("reviewStatus") != "no assertion criteria provided"})
        if f.get("ftId", "").startswith("VAR_"):
            nat[p1].append(dict(ftId=f["ftId"], sig="; ".join(sig), desc=curated.get(p1, "")))
    up = pd.DataFrame([dict(p1=k, uniprot_var=";".join(x["ftId"] for x in v),
                            uniprot_sig="; ".join(sorted({x["sig"] for x in v if x["sig"]})) or None,
                            uniprot_desc=(v[0]["desc"] or None)) for k, v in nat.items()]).set_index("p1")
    log(f"UniProt dizi={len(seq)} aa, ozellik={len(feats)}, VAR_ kaydi={len(up)}")
    return seq, feats, up


def read_hotspots():
    h = json.load(open(RAW / "cancerhotspots_single_TP53.json", encoding="utf-8"))
    res, aa = {}, {}
    for r in h:
        pos = r["aminoAcidPosition"]["start"]
        res[pos] = dict(residue=r["residue"], tumorCount=r.get("tumorCount"), qValue=r.get("qValue"),
                        composition=r.get("tumorTypeComposition") or {}, variants=r.get("variantAminoAcid") or {})
        for a, n in (r.get("variantAminoAcid") or {}).items():
            aa[f"{r['residue']}{a}"] = n
    log(f"cancerhotspots {len(res)} TP53 kalintisi")
    return res, aa


def read_civic():
    out = []
    for k in ("ClinicalEvidenceSummaries", "VariantSummaries"):
        p = RAW / f"civic_{k}_TP53.tsv"
        if p.exists():
            out.append(pd.read_csv(p, sep="\t", low_memory=False))
    if out:
        out[0].to_parquet(PROC / "civic_evidence.parquet", index=False)
        out[1].to_parquet(PROC / "civic_variants.parquet", index=False)
        log(f"CIViC kanit={len(out[0])} varyant={len(out[1])}")


# =================================================================== birlestirme
def _ckey(c):
    """'c.524G>A' -> ('524', 'G', 'A') ; 'c.375+1G>T' -> ('375+1','G','T'); SNV degilse None."""
    m = re.match(r"^c\.([\d+\-*]+)([ACGT])>([ACGT])$", c if isinstance(c, str) else "")
    return (m.group(1), m.group(2), m.group(3)) if m else None


def build_position_map(gn, nci, cv):
    """c. konumu <-> hg38 konumu eslemesi (SNV ciftlerinden) ve kodlayan dizi (c. zincirinde)."""
    c2g, cref = {}, {}
    pairs = []
    pairs += [(c, p) for c, p in zip(gn["c_gnomad"], gn["hg38_pos"]) if c]
    pairs += [(c, NUM(p)) for c, p in zip(nci.index, nci.get("hg38_Chr17_coordinates", pd.Series(dtype=float)))]
    for c, spdi in zip(cv.index, cv["clinvar_spdi"]):
        m = re.match(r"NC_000017\.\d+:(\d+):[ACGT]:[ACGT]$", spdi if isinstance(spdi, str) else "")
        if m:
            pairs.append((c, int(m.group(1)) + 1))
    for c, p in pairs:
        k = _ckey(c)
        if not k or p is None or (isinstance(p, float) and math.isnan(p)):
            continue
        c2g.setdefault(k[0], int(p))
        cref.setdefault(k[0], k[1])
    # kodlayan konumlarda eksikleri ayni ekson icinden doldur (eksi iplik: c +1 => g -1)
    for ex, a, b in EXONS:
        known = [(i, c2g[str(i)]) for i in range(a, b + 1) if str(i) in c2g]
        if not known:
            continue
        i0, g0 = known[0]
        for i in range(a, b + 1):
            c2g.setdefault(str(i), g0 - (i - i0))
    g2c = {v: k for k, v in c2g.items()}
    return c2g, g2c, cref


def mutate_ptc(c, cds):
    """Cerceve kaymasi icin yeni durma kodonunun konumu (kodon). CDS icinde durma yoksa 'beyond'."""
    m = re.match(r"^c\.(\d+)(?:_(\d+))?(del|dup|ins|delins)([ACGTN]*)$", c or "")
    if not m or not cds:
        return None
    a = int(m.group(1)); b = int(m.group(2) or a); op = m.group(3); seq = m.group(4)
    s = cds
    if op == "del":
        new = s[:a - 1] + s[b:]
    elif op == "dup":
        new = s[:b] + s[a - 1:b] + s[b:]
    elif op == "ins":
        if not seq:
            return None
        new = s[:a] + seq + s[a:]
    else:
        if not seq:
            return None
        new = s[:a - 1] + seq + s[b:]
    if (len(new) - len(s)) % 3 == 0:
        return "inframe"
    start = (a - 1) // 3
    for k in range(start, len(new) // 3):
        if new[3 * k:3 * k + 3] in ("TAA", "TAG", "TGA"):
            return k + 1
    return "beyond"


def build():
    from .common import HG19_HG38_OFFSET, translate
    s1, s2, s3 = read_vcep()
    db = read_dbnsfp()
    gn = read_gnomad()
    nci, fcount = read_nci()
    cv = read_clinvar()
    er = read_erepo()
    seq, feats, up = read_uniprot()
    hs_res, hs_aa = read_hotspots()
    read_civic()

    c2g, g2c, cref = build_position_map(gn, nci, cv)
    g_ref_at = dict(zip(db["hg38_pos"].dropna().astype(int), db.loc[db["hg38_pos"].notna(), "g_ref"]))
    for i in range(1, 1183):
        if str(i) not in cref and c2g.get(str(i)) in g_ref_at:
            cref[str(i)] = COMP[g_ref_at[c2g[str(i)]]]
    unknown = [i for i in range(1, 1183) if str(i) not in cref]
    # kalanlar 4 katli dejenere kodonlarin 3. bazi (missense yok); protein degismez, yer tutucu 'C'
    cds = "".join(cref.get(str(i), "C") for i in range(1, 1183))
    if unknown:
        log(f"CDS'de kaynaklarda gorulmeyen {len(unknown)} baz (4 katli dejenere 3. pozisyon): {unknown}")
    prot = translate(cds)
    ok = prot.rstrip("*") == seq
    log(f"CDS {len(cds)} nt, eksik baz={cds.count('N')}, ceviri UniProt ile {'AYNI' if ok else 'FARKLI'}")
    if not ok:
        diff = [i + 1 for i, (x, y) in enumerate(zip(prot, seq)) if x != y][:10]
        log(f"  farkli kodonlar: {diff}")

    # dbNSFP satirlarina c. ata
    db["c"] = [
        (f"c.{g2c[p]}{COMP[r]}>{COMP[a]}" if p in g2c else None)
        for p, r, a in zip(db["hg38_pos"].fillna(-1).astype(int), db["g_ref"], db["g_alt"])
    ]
    unmapped = db["c"].isna().sum()
    db = db.dropna(subset=["c"]).drop_duplicates("c").set_index("c")
    log(f"dbNSFP c. eslendi: {len(db)} (eslenemeyen {unmapped})")

    gn_c = gn.dropna(subset=["c_gnomad"]).drop_duplicates("c_gnomad").set_index("c_gnomad")
    gn_c = gn_c.drop(columns=["hg38_pos", "g_ref", "g_alt"])
    s2i = s2.drop_duplicates("c").set_index("c")
    s1i = s1.drop_duplicates("c").set_index("c")

    allc = set(db.index) | set(gn_c.index) | set(nci.index) | set(cv.index) | set(er.index) | set(s2i.index)
    allc = {c for c in allc if c and c.startswith("c.")}
    df = pd.DataFrame(index=sorted(allc))
    df.index.name = "c"
    df = (df.join(db.drop(columns=["aa_ref", "aa_alt"]), how="left")
            .join(gn_c, how="left").join(nci, how="left").join(cv, how="left").join(er, how="left")
            .join(s2i, how="left").join(s1i, how="left"))

    # ---- protein degisimi (tek harf) : S2 > gnomAD > ClinVar > NCI > ceviri
    def pick_p(row):
        for col in ("p3", "p_gnomad", "p_clinvar"):
            v = row.get(col)
            if isinstance(v, str) and v.startswith("p."):
                return p3_to_p1(v)
        v = row.get("p_nci")
        if isinstance(v, str) and v.startswith("p.") and v != "p.?":
            return p3_to_p1(v)
        return None
    df["p1"] = [pick_p(r) for r in df.to_dict("records")]
    for c in df.index[df["p1"].isna()]:
        k = _ckey(c)
        if k and k[0].isdigit() and 1 <= int(k[0]) <= 1182:
            pos = int(k[0]); codon = (pos - 1) // 3 + 1
            wt = cds[3 * (codon - 1):3 * codon]
            mt = list(wt); mt[(pos - 1) % 3] = k[2]; mt = "".join(mt)
            a, b = translate(wt), translate(mt)
            df.at[c, "p1"] = f"{a}{codon}{'=' if a == b else b}"
    df["p1"] = df["p1"].str.replace(r"^([A-Z])(\d+)\1$", r"\1\2=", regex=True)

    # ---- temel anotasyonlar
    df["consequence"] = [consequence(c, p) for c, p in zip(df.index, df["p1"])]
    cod = df["p1"].str.extract(r"^[A-Z\*](\d+)")[0]
    df["codon"] = pd.to_numeric(cod, errors="coerce")
    mm = df["p1"].str.extract(r"^([A-Z\*])\d+([A-Z\*])$")
    df["aa_ref"], df["aa_alt"] = mm[0], mm[1]
    df["domain"] = df["codon"].map(lambda x: domain_of(int(x)) if not is_nan(x) else None)
    df["grantham"] = [grantham(a, b) if isinstance(a, str) and isinstance(b, str) else None
                      for a, b in zip(df["aa_ref"], df["aa_alt"])]
    cpos = [parse_cpos(c) for c in df.index]
    df["c_pos"] = [x[0] if x else None for x in cpos]
    df["c_offset"] = [x[1] if x else None for x in cpos]
    df["exon"] = [exon_of_cpos(x[0]) if x and x[1] == 0 else None for x in cpos]
    df["hotspot_codon"] = df["codon"].isin(list(HOTSPOT_CODONS)) & (df["consequence"] == "missense")
    df["hotspot_residue_tumors"] = df["codon"].map(lambda x: (hs_res.get(int(x)) or {}).get("tumorCount") if not is_nan(x) else None)
    df["hotspot_aa_count"] = [hs_aa.get(p[:-1] + p[-1]) if isinstance(p, str) and MISSENSE_RE.match(p) else None
                              for p in df["p1"]]
    df["ptc_codon"] = [mutate_ptc(c, cds) if cq in ("frameshift", "indel_other") else None
                       for c, cq in zip(df.index, df["consequence"])]
    for c in df.index[df["consequence"] == "indel_other"]:
        v = df.at[c, "ptc_codon"]
        if v == "inframe":
            df.at[c, "consequence"] = "inframe_indel"
        elif v is not None and not is_nan(v):
            df.at[c, "consequence"] = "frameshift"
    df.loc[df["consequence"] == "nonsense", "ptc_codon"] = df.loc[df["consequence"] == "nonsense", "codon"]
    df["ptc_codon"] = df["ptc_codon"].astype(object)
    # p. notasyonu olmayan cerceve kaymalari icin p. uret (orn. P72fs*51)
    for c in df.index[(df["consequence"] == "frameshift") & df["p1"].isna()]:
        m = re.match(r"^c\.(\d+)", c)
        if not m:
            continue
        codon = (int(m.group(1)) - 1) // 3 + 1
        if 1 <= codon <= len(seq):
            ptc = df.at[c, "ptc_codon"]
            tail = f"*{int(ptc) - codon + 1}" if isinstance(ptc, (int, np.integer)) else ""
            df.at[c, "p1"] = f"{seq[codon - 1]}{codon}fs{tail}"
            df.at[c, "codon"] = codon
    df["domain"] = df["codon"].map(lambda x: domain_of(int(x)) if not is_nan(x) else None)
    df = df.join(s3.set_index("p1"), on="p1", how="left")
    # S3'te ayni protein degisimi farkli nukleotid degisimleri icin ayri satir olabilir (Kotler/Funk farki):
    # c. ile uyusan satiri tut, yoksa ilkini.
    def _nt(c):
        m = re.match(r"^c\.\d+([ACGT]+)>([ACGT]+)$", c) or re.match(r"^c\.\d+_\d+delins([ACGT]+)$", c)
        if not m:
            return None
        return f"{COMP.get(m.group(1), '')}>{m.group(2)}" if len(m.groups()) == 2 else None
    dup = df.index.duplicated(keep=False)
    if dup.any():
        tmp = df[dup].copy()
        tmp["_ok"] = [(_nt(c) is not None and isinstance(sc, str) and sc == f"{_s1}>{_s2}")
                      for c, sc, _s1, _s2 in zip(tmp.index, tmp["s3_seq_change"],
                                                 [re.sub(r'^c\.\d+([ACGT]).*', r'\1', c) for c in tmp.index],
                                                 [re.sub(r'^.*>([ACGT])$', r'\1', c) for c in tmp.index])]
        keep = tmp.sort_values("_ok", ascending=False).groupby(level=0).head(1).drop(columns="_ok")
        df = pd.concat([df[~dup], keep]).sort_index()
    df = df.join(up, on="p1", how="left")
    df = df.join(fcount, on="p1", how="left")

    # ---- koordinatlar
    def g38(c):
        k = _ckey(c)
        return c2g.get(k[0]) if k else None
    df["hg38_pos"] = df["hg38_pos"].fillna(pd.Series({c: g38(c) for c in df.index}))
    df["hg19_pos"] = df["hg19_pos"].fillna(df["hg38_pos"] - HG19_HG38_OFFSET)
    snv = [(_ckey(c)) for c in df.index]
    df["g_ref"] = df["g_ref"].fillna(pd.Series({c: COMP[k[1]] for c, k in zip(df.index, snv) if k}))
    df["g_alt"] = df["g_alt"].fillna(pd.Series({c: COMP[k[2]] for c, k in zip(df.index, snv) if k}))
    df["franklin_id"] = [f"chr17-{int(p)}-{r}-{a}" if not is_nan(p) and isinstance(r, str) else None
                         for p, r, a in zip(df["hg19_pos"], df["g_ref"], df["g_alt"])]
    df["hg38_id"] = [f"17-{int(p)}-{r}-{a}" if not is_nan(p) and isinstance(r, str) else None
                     for p, r, a in zip(df["hg38_pos"], df["g_ref"], df["g_alt"])]
    df["in_gnomad"] = df["gnomad_an"].notna()
    df["spliceai"] = df[["spliceai_max", "s2_spliceai"]].apply(pd.to_numeric, errors="coerce").max(axis=1)
    df["agvgd"] = df["agvgd"].fillna(df["agvgd_nci"])
    df["bayesdel_vcep"] = pd.to_numeric(df["bayesdel_vcep"], errors="coerce").fillna(df["bayesdel"])

    df = df.reset_index()
    for col in df.columns:   # parquet icin karisik tipli kolonlari metne cevir
        if df[col].dtype == object:
            df[col] = df[col].map(lambda v: None if is_nan(v) else (v if isinstance(v, str) else str(v)))
    df.to_parquet(PROC / "variants.parquet", index=False)

    gene = dict(cds=cds, protein=seq, uniprot_features=feats,
                hotspots={str(k): v for k, v in hs_res.items()},
                domains=[dict(name=n, start=a, end=b) for n, a, b in __import__("tp53.common", fromlist=["DOMAINS"]).DOMAINS],
                counts=dict(variants=len(df), missense=int((df.consequence == "missense").sum()),
                            erepo=int(df.erepo_class.notna().sum()), clinvar=int(df.clinvar_vcv.notna().sum())))
    (PROC / "gene.json").write_text(json.dumps(gene), encoding="utf-8")
    log(f"variants.parquet: {len(df)} varyant, {df.shape[1]} kolon")
    log("sonuc tipleri:", dict(Counter(df["consequence"])))
    return df


if __name__ == "__main__":
    build()
