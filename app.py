"""TP53-Franklin: Franklin (Genoox) benzeri, TP53'e odakli varyant yorumlama arayuzu.

Calistirma:  streamlit run app.py
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import altair as alt
import pandas as pd
import streamlit as st

from tp53.acmg import (CATEGORY_ORDER, LFS_MODERATE, LFS_STRONG, Clinical, TP53Classifier)
from tp53.common import DOMAINS, PROC, RAW, MODELS, p1_to_p3
from tp53.help_texts import HELP, TOOL_INFO

st.set_page_config(page_title="TP53-Franklin", page_icon="🧬", layout="wide")

NAVY, RED, GREEN, GREY, TEAL = "#1b2448", "#d93a3a", "#2a9d8f", "#8a8f9c", "#1f6f8b"
CLASS_COLOR = {"Pathogenic": RED, "Likely Pathogenic": "#ef7f45", "Uncertain Significance": GREY,
               "Likely Benign": "#57b8a6", "Benign": GREEN}
CLASS_TR = {"Pathogenic": "Patojenik", "Likely Pathogenic": "Olası patojenik", "Uncertain Significance": "VUS",
            "Likely Benign": "Olası benign", "Benign": "Benign"}

st.markdown(f"""
<style>
.block-container {{padding-top: 1.2rem; max-width: 1500px;}}
.vhead {{background:{NAVY}; color:white; padding:14px 22px; border-radius:10px; margin-bottom:10px}}
.vhead h1 {{color:white; font-size:1.9rem; margin:0; letter-spacing:.5px}}
.vhead .sub {{color:#c9cfe3; font-size:.92rem; margin-top:4px}}
.vhead a {{color:#9fd3ff; text-decoration:underline; margin-right:10px}}
.card {{border:1px solid #e3e6ee; border-radius:10px; padding:12px 16px; margin-bottom:10px; background:#fbfcfe}}
.ev {{display:flex; gap:14px; align-items:flex-start; border:1px solid #e3e6ee; border-radius:10px; padding:10px 14px;
      margin:6px 0; background:#f4f6fa}}
.badge {{min-width:104px; text-align:center; border-radius:7px; padding:6px 4px; color:white; font-weight:700}}
.badge small {{display:block; font-weight:500; font-size:.72rem; opacity:.95}}
.unmet {{color:#5b6070; font-size:.86rem; padding:2px 0 8px 4px}}
.cls {{font-size:1.6rem; font-weight:800}}
.pill {{display:inline-block; border-radius:12px; padding:2px 10px; margin:2px 4px 2px 0; font-size:.82rem;
        background:#eef1f7; border:1px solid #d9deea}}
.muted {{color:#6b7080; font-size:.85rem}}
.help {{border-left:4px solid #4a7bd0; background:#f1f6ff; border-radius:6px; padding:10px 14px; margin:4px 0 12px 0;
        font-size:.9rem; color:#26324d}}
.help b {{color:#1b2448}}
</style>""", unsafe_allow_html=True)


# ====================================================================== veri
@st.cache_resource(show_spinner="Veri yukleniyor...")
def load():
    df = pd.read_parquet(PROC / "variants.parquet")
    ml = PROC / "ml_predictions.parquet"
    if ml.exists():
        df = df.merge(pd.read_parquet(ml).drop(columns=["p1"]), on="c", how="left")
    gene = json.loads((PROC / "gene.json").read_text(encoding="utf-8"))
    extra = {}
    for k in ("nci_tumor_topography", "nci_germline", "nci_function", "civic_evidence", "erepo"):
        p = PROC / f"{k}.parquet"
        extra[k] = pd.read_parquet(p) if p.exists() else pd.DataFrame()
    metrics = json.loads((MODELS / "metrics.json").read_text()) if (MODELS / "metrics.json").exists() else {}
    manifest = json.loads((RAW / "MANIFEST.json").read_text(encoding="utf-8")) if (RAW / "MANIFEST.json").exists() else {}
    val = PROC / "validation_vcep_vs_erepo.csv"
    extra["validation"] = pd.read_csv(val) if val.exists() else pd.DataFrame()
    return df, TP53Classifier(df), gene, extra, metrics, manifest


df, clf, gene, X, METRICS, MANIFEST = load()


def f(x, nd=3):
    try:
        v = float(x)
        if math.isnan(v):
            return "—"
        if nd == 0:
            return f"{v:,.0f}".replace(",", ".")
        if abs(v) < 0.001 and v != 0:
            return f"{v:.{nd}g}"
        return f"{v:.{nd}f}".rstrip("0").rstrip(".")
    except (TypeError, ValueError):
        return "—" if x is None else str(x)


def s(x):
    return None if x is None or (isinstance(x, float) and math.isnan(x)) or str(x) in ("nan", "None", "") else str(x)


# ====================================================================== kenar cubugu
with st.sidebar:
    st.markdown("### 🔎 TP53 varyant ara")
    q = st.text_input("Varyant", value=st.session_state.get("q", "R175H"),
                      help="Ornekler: R175H, p.Arg175His, c.524G>A, TP53:c.524G>A, chr17-7578406-C-T (hg19), "
                           "17-7675088-C-T (hg38), rs28934578")
    cols = st.columns(3)
    for i, ex in enumerate(["R175H", "R248Q", "Y220C", "R337H", "P72R", "c.375+1G>A", "R213*", "c.215del", "G245S"]):
        if cols[i % 3].button(ex, width="stretch"):
            st.session_state["q"] = ex
            st.rerun()
    st.divider()
    st.caption(f"Veri: {gene['counts']['variants']:,} varyant · {gene['counts']['missense']:,} missense · "
               f"{gene['counts']['erepo']} uzman sınıflaması · {gene['counts']['clinvar']:,} ClinVar kaydı")
    st.caption("Kurallar: ClinGen TP53 VCEP v2.4 (CSpec GN009). Araştırma amaçlıdır; tanı/tedavi tavsiyesi değildir.")
    SHOW_HELP = st.toggle("Açıklamaları göster", value=True, help="Her bölümün başındaki 'Bu bölüm ne anlatıyor?' kutuları")

hits = clf.find(q)
if not len(hits):
    st.error(f"'{q}' bulunamadı. TP53 (NM_000546.6) c., p., hg19/hg38 koordinat veya rsID girin.")
    st.stop()
if len(hits) > 1:
    choice = st.sidebar.selectbox("Birden fazla eşleşme", hits.index.tolist(),
                                  format_func=lambda c: f"{c} ({hits.loc[c, 'p1']})")
else:
    choice = hits.index[0]
r = clf.df.loc[choice]
C = r["c"]
P1 = s(r.get("p1"))
P3 = p1_to_p3(P1) if P1 else None

# ====================================================================== baslik
links = []
if s(r.get("rsid")):
    links.append(f'<a href="https://www.ncbi.nlm.nih.gov/snp/{r["rsid"]}" target="_blank">{r["rsid"]}</a>')
if s(r.get("clinvar_id")):
    links.append(f'<a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/{r["clinvar_id"]}/" target="_blank">ClinVar</a>')
if s(r.get("hg38_id")):
    links.append(f'<a href="https://gnomad.broadinstitute.org/variant/{r["hg38_id"]}?dataset=gnomad_r4" target="_blank">gnomAD</a>')
if s(r.get("franklin_id")):
    links.append(f'<a href="https://franklin.genoox.com/clinical-db/variant/snp/{r["franklin_id"]}" target="_blank">Franklin</a>')
if s(r.get("erepo_url")):
    links.append(f'<a href="{r["erepo_url"]}" target="_blank">ClinGen eRepo</a>')
if P1:
    links.append(f'<a href="https://tp53.cancer.gov/results_somatic_mutation_list?protein_description=p.{P1}" target="_blank">NCI TP53</a>')
coord = " | ".join(x for x in [
    f"hg19 {s(r.get('franklin_id'))}" if s(r.get("franklin_id")) else None,
    f"hg38 {s(r.get('hg38_id'))}" if s(r.get("hg38_id")) else None,
    P3, "NM_000546.6", f"ekson {int(float(r['exon']))}" if s(r.get("exon")) else None, s(r.get("domain"))] if x)
st.markdown(f"""<div class="vhead"><h1>TP53:{C}</h1>
<div class="sub">{coord}</div><div class="sub">{' '.join(links)}</div></div>""", unsafe_allow_html=True)

# ====================================================================== klinik form durumu
FIELDS = dict(dn_strong_conf=0, dn_strong_assumed=0, dn_mod_conf=0, dn_mod_assumed=0, chompret=0, classic=0, her2=0,
              meioses=0, families=1, vaf525=0, vaf2535=0, bs2=0, bs4=False, bp2=False, use_expert=True)
pending = st.session_state.pop("pending_clinical", None)
for k, v in FIELDS.items():
    key = f"{C}|{k}"
    if pending and k in pending:
        st.session_state[key] = pending[k]
    st.session_state.setdefault(key, v)


def clinical_from_state():
    g = lambda k: st.session_state[f"{C}|{k}"]
    dn = ([(LFS_STRONG, True)] * g("dn_strong_conf") + [(LFS_STRONG, False)] * g("dn_strong_assumed") +
          [(LFS_MODERATE, True)] * g("dn_mod_conf") + [(LFS_MODERATE, False)] * g("dn_mod_assumed"))
    return Clinical(denovo=dn, chompret=g("chompret"), classic_lfs=g("classic"), her2_lt40=g("her2"), meioses=g("meioses"),
                    families=g("families"), vaf_5_25=g("vaf525"), vaf_25_35=g("vaf2535"), bs2_females60=g("bs2"),
                    nonsegregation=g("bs4"), trans_or_cis=g("bp2"), use_expert=g("use_expert"))


tabs = st.tabs(["ACMG Sınıflama", "Varyant Değerlendirme", "Fonksiyon & ML", "Yayınlar", "Gen Değerlendirme",
                "İlişkili Hastalıklar", "Somatik Klinik Kanıt", "LLM Asistanı", "Veri Kaynakları"])


def class_strip(cls, title="Önerilen sınıflandırma", points=None):
    order = ["Benign", "Likely Benign", "Uncertain Significance", "Likely Pathogenic", "Pathogenic"]
    idx = order.index(cls)
    dots = "".join(
        f'<div style="flex:1;text-align:center"><div style="height:8px;background:{CLASS_COLOR[o]};'
        f'border-radius:4px;margin:0 3px"></div><div style="font-size:.78rem;margin-top:4px;'
        f'{"font-weight:800;color:" + CLASS_COLOR[o] if i == idx else "color:#777"}">{CLASS_TR[o]}</div></div>'
        for i, o in enumerate(order))
    pts = f'<span class="muted"> · {points:+d} puan</span>' if points is not None else ""
    st.markdown(f'<div class="card"><div class="muted">{title}</div><div class="cls" style="color:{CLASS_COLOR[cls]}">'
                f'{cls}{pts}</div><div style="display:flex;margin-top:8px">{dots}</div></div>', unsafe_allow_html=True)


def explain(key, title="Bu bölüm ne anlatıyor?"):
    """Konuyu bilmeyen kullanici icin bolum basi aciklama kutusu (kenar cubugundan kapatilabilir)."""
    if not SHOW_HELP:
        return
    txt = HELP[key].replace("{v}", P1 or C)
    with st.expander(f"ℹ️ {title}", expanded=True):
        st.markdown(txt)


def evidence_cards(res):
    evs = pd.DataFrame(res["evidence"])
    for cat in CATEGORY_ORDER:
        sub = evs[evs.category == cat]
        if sub.empty:
            continue
        met = sub[sub.met & (sub.points != 0)]
        noinf = sub[sub.met & (sub.points == 0)]
        unmet = sub[~sub.met]
        st.markdown(f"##### {cat}")
        if met.empty and noinf.empty:
            st.markdown('<div class="unmet">Karşılanan kriter yok</div>', unsafe_allow_html=True)
        for _, e in pd.concat([met, noinf]).iterrows():
            col = RED if e.code.startswith("P") else GREEN
            if e.points == 0:
                col = GREY
            kind = ("Patojenik" if e.code.startswith("P") else "Benign") + f" {e.strength or ''}"
            src = f'<div class="muted">Kaynak: {e.source}{"" if e.auto else " · klinik/kullanıcı kanıtı"}</div>'
            st.markdown(f'<div class="ev"><div class="badge" style="background:{col}">{e.label}<small>{e.strength or ""}'
                        f' · {e.points:+d}</small></div><div><b>{kind}:</b> {e.description}<br>{e.reason}{src}</div></div>',
                        unsafe_allow_html=True)
        if not unmet.empty:
            items = "; ".join(f"<b>{u.code}</b>" + ("" if u.applicable else " (TP53'te uygulanmaz)") +
                              f": <span title='{u.reason}'>{u.reason[:140]}</span>" for _, u in unmet.iterrows())
            st.markdown(f'<div class="unmet">KARŞILANMAYAN: {items}</div>', unsafe_allow_html=True)


# ====================================================================== 1) ACMG
with tabs[0]:
    explain("acmg_intro")
    with st.expander("➕ Klinik kanıt ekle (PS2, PS4, PP1, PP4, BS2, BS4, BP2)", expanded=False):
        if SHOW_HELP:
            st.info(HELP["acmg_form"])
        st.caption("Puan tabloları: Genome Medicine 2025 (TP53 VCEP v2) Tablo 2-3. LLM Asistanı sekmesinden ClinVar "
                   "metinlerinden öneri alıp buraya aktarabilirsiniz.")
        k = lambda n: f"{C}|{n}"
        st.checkbox("Varyant ClinGen eRepo'da ise uzman panelin klinik kodlarını kullan (girilen değer yoksa)", key=k("use_expert"))
        a, b, c3, d = st.columns(4)
        with a:
            st.markdown("**PS2 – de novo probandlar**")
            st.number_input("Güçlü ilişkili kanser, anne+baba doğrulanmış (4 p)", 0, 20, key=k("dn_strong_conf"))
            st.number_input("Güçlü ilişkili, varsayılan (2 p)", 0, 20, key=k("dn_strong_assumed"))
            st.number_input("Orta ilişkili, doğrulanmış (2 p)", 0, 20, key=k("dn_mod_conf"))
            st.number_input("Orta ilişkili, varsayılan (1 p)", 0, 20, key=k("dn_mod_assumed"))
            st.caption("Güçlü: meme ca <31y, koroid pleksus ca, adrenokortikal <18y, rabdo/osteosarkom <46y. "
                       "Orta: meme 31-49y, beyin tm <46y, ACC 18-49y, akciğer <46y, diğer sarkom <60y, hipodiploid ALL, SHH-MB.")
        with b:
            st.markdown("**PS4 – probandlar**")
            st.number_input("Revised Chompret (0.5 p)", 0, 50, key=k("chompret"))
            st.number_input("Classic LFS (1 p)", 0, 50, key=k("classic"))
            st.number_input("HER2+ meme ca <40y (0.5 p)", 0, 50, key=k("her2"))
        with c3:
            st.markdown("**PP1 / PP4**")
            st.number_input("Kosegregasyon meiozu", 0, 30, key=k("meioses"))
            st.number_input("Aile sayısı", 1, 20, key=k("families"))
            st.number_input("VAF %5-25 gözlem (PP4)", 0, 20, key=k("vaf525"))
            st.number_input("VAF %25-35 gözlem (PP4)", 0, 20, key=k("vaf2535"))
        with d:
            st.markdown("**Benign klinik**")
            st.number_input("BS2: ≥60y kansersiz kadın taşıyıcı", 0, 50, key=k("bs2"))
            st.checkbox("BS4: segregasyon yok", key=k("bs4"))
            st.checkbox("BP2: trans/cis P/LP varyant", key=k("bp2"))
    clinical = clinical_from_state()
    res = clf.vcep(C, clinical)
    fr = clf.franklin(C)

    left, right = st.columns([2.2, 1])
    with left:
        class_strip(res["classification"], "Önerilen sınıflandırma — ClinGen TP53 VCEP v2.4 (Bayes puan)", res["points"])
    with right:
        st.markdown(f'<div class="card"><div class="muted">Yalnızca otomatik (varyant) kanıtlar</div>'
                    f'<b style="color:{CLASS_COLOR[res["auto_classification"]]}">{res["auto_classification"]}</b> '
                    f'({res["auto_points"]:+d})<br><div class="muted" style="margin-top:6px">Franklin-benzeri genel ACMG</div>'
                    f'<b style="color:{CLASS_COLOR[fr["classification"]]}">{fr["classification"]}</b> ({fr["points"]:+d})<br>'
                    f'<span class="muted">{", ".join(fr["met"])}</span></div>', unsafe_allow_html=True)
        if s(r.get("erepo_class")):
            st.markdown(f'<div class="card"><div class="muted">ClinGen uzman paneli (eRepo, v{s(r.get("erepo_version")) or "?"}, '
                        f'{s(r.get("erepo_date"))})</div><b style="color:{CLASS_COLOR.get(r["erepo_class"], GREY)}">'
                        f'{r["erepo_class"]}</b> ({f(r.get("erepo_points"), 0)} puan)<br>'
                        f'<span class="muted">{s(r.get("erepo_met")).replace(";", ", ")}</span></div>', unsafe_allow_html=True)
    st.markdown("**Karşılanan kriterler:** " + " ".join(f'<span class="pill">{m}</span>' for m in res["met"]),
                unsafe_allow_html=True)
    for n in res["notes"]:
        st.info(n)
    explain("acmg_cards", "Kriter kartları nasıl okunur?")
    mode = st.radio("Kanıt görünümü", ["ClinGen TP53 VCEP v2.4", "Franklin-benzeri genel ACMG"], horizontal=True)
    evidence_cards(res if mode.startswith("ClinGen") else fr)

# ====================================================================== 2) Varyant degerlendirme
TOOLS = [  # (kolon, ad, yon, esik-zararli, esik-benign, aciklama)
    ("ml_path_score", "TP53 Aggregated (bu araç, ML)", 1, 0.5, 0.5, "Uzman+ClinVar etiketli TP53 meta-modeli"),
    ("revel", "REVEL", 1, 0.644, 0.29, "ClinGen SVI 2022: >=0.644 PP3, >=0.773 Mod, >=0.932 Strong"),
    ("alphamissense", "AlphaMissense", 1, 0.564, 0.34, "Yazar esikleri"),
    ("bayesdel_vcep", "BayesDel (noAF)", 1, 0.16, -0.008, "TP53 VCEP esikleri"),
    ("varity_r", "VARITY_R", 1, 0.5, 0.5, ""),
    ("esm1b", "ESM1b", -1, -7.5, -7.5, "Brandes 2023"),
    ("eve", "EVE", 1, 0.5, 0.5, ""),
    ("primateai", "PrimateAI", 1, 0.803, 0.803, "Yazar esigi"),
    ("cadd_phred", "CADD", 1, 25.3, 22.7, "ClinGen SVI 2022"),
    ("polyphen2_hvar", "PolyPhen-2 HVAR", 1, 0.909, 0.447, ""),
    ("sift", "SIFT", -1, 0.05, 0.05, "<=0.05 zararli"),
    ("fathmm", "FATHMM", -1, -1.5, -1.5, "<=-1.5 zararli"),
    ("mutationassessor", "MutationAssessor", 1, 1.935, 0.8, ">3.5 yuksek"),
    ("mutationtaster", "MutationTaster", 1, 0.5, 0.5, ""),
    ("metalr", "MetaLR", 1, 0.5, 0.5, ""),
    ("metasvm", "MetaSVM", 1, 0.0, 0.0, ""),
    ("dann", "DANN", 1, 0.96, 0.9, ""),
    ("provean", "PROVEAN", -1, -2.5, -2.5, ""),
    ("vest4", "VEST4", 1, 0.764, 0.449, "ClinGen SVI 2022"),
    ("mutpred", "MutPred2", 1, 0.737, 0.391, "ClinGen SVI 2022"),
    ("clinpred", "ClinPred", 1, 0.5, 0.5, ""),
    ("gerp", "GERP++ (korunmuşluk)", 1, 4.4, 3.0, ">4.4 yuksek kisit"),
    ("phylop100", "phyloP100 (korunmuşluk)", 1, 2.0, 0.0, ""),
    ("genocanyon", "GenoCanyon", 1, 0.5, 0.5, ""),
    ("fitcons", "fitCons", 1, 0.5, 0.5, ""),
    ("spliceai", "SpliceAI (max Δ)", 1, 0.2, 0.1, "VCEP: >=0.2 PP3, <=0.1 BP4"),
    ("dbscsnv_ada", "dbscSNV Ada", 1, 0.6, 0.15, ""),
    ("dbscsnv_rf", "dbscSNV RF", 1, 0.6, 0.15, ""),
]


@st.cache_data
def missense_distributions():
    m = df[df.consequence == "missense"]
    return {t[0]: pd.to_numeric(m[t[0]], errors="coerce").dropna().values for t in TOOLS if t[0] in m}


def norm_trait(t: str) -> str:
    """ClinVar durum adlarini birlestirir: 'Li-Fraumeni syndrome 1', 'Li-fraumeni syndrome 1 (lfs1)' -> 'Li-Fraumeni syndrome'."""
    import re as _re
    t = _re.sub(r"\s*\([^)]*\)\s*$", "", t.strip())          # sondaki kisaltma: (LFS1)
    t = _re.sub(r",?\s+\d+$", "", t)                           # sondaki alt tip numarasi: ' 1', ', 1'
    return t[:1].upper() + t[1:] if t else t


def traits_of(x) -> list:
    out = []
    for t in (s(x) or "").split(";"):
        t = t.strip()
        if t and not re.search(r"not provided|not specified", t, re.I):
            n = norm_trait(t)
            if n.lower() not in [o.lower() for o in out]:
                out.append(n)
    return out


def verdict(v, d, hi, lo):
    if v is None:
        return "—"
    if d == 1:
        return "Zararlı" if v >= hi else "Benign" if v <= lo else "Belirsiz"
    return "Zararlı" if v <= hi else "Benign" if v > lo else "Belirsiz"


with tabs[1]:
    explain("va_metrics")
    a, b, c3, d, e = st.columns(5)
    a.metric("Sonuç tipi", s(r.get("consequence")) or "—")
    b.metric("ClinVar (germline)", s(r.get("clinvar_germline")) or "kayıt yok")
    c3.metric("gnomAD v4 AC/AN", f"{f(r.get('gnomad_ac'), 0)}/{f(r.get('gnomad_an'), 0)}")
    d.metric("NCI somatik tümör", f(r.get("nci_tumor_n"), 0))
    e.metric("TP53 ML skoru", f(r.get("ml_path_score")))

    st.markdown("#### Bölge görüntüleyici (ClinVar ≥1★ ve ClinGen uzman sınıflamaları)")
    explain("va_region", "Bu grafik nasıl okunur?")
    flt = st.radio("Filtre", ["Missense", "LOF (nonsense/frameshift/splice)", "Tümü"], horizontal=True, key="rv")
    reg = df[(pd.to_numeric(df.clinvar_stars, errors="coerce").fillna(0) >= 1) | df.erepo_class.notna()].copy()
    if flt.startswith("Missense"):
        reg = reg[reg.consequence == "missense"]
    elif flt.startswith("LOF"):
        reg = reg[reg.consequence.isin(["nonsense", "frameshift", "splice_canonical", "start_lost"])]
    sig = reg.erepo_class.fillna(reg.clinvar_germline).fillna("")
    reg["grup"] = pd.Series("VUS/çelişkili", index=reg.index)
    reg.loc[sig.str.contains("athogenic") & ~sig.str.contains("onflict|enign"), "grup"] = "P/LP"
    reg.loc[sig.str.contains("enign") & ~sig.str.contains("onflict|athogenic"), "grup"] = "B/LB"
    reg["codon"] = pd.to_numeric(reg.codon, errors="coerce")
    agg = reg.dropna(subset=["codon"]).groupby(["codon", "grup"]).size().rename("n").reset_index()
    cmap = alt.Scale(domain=["P/LP", "VUS/çelişkili", "B/LB"], range=[RED, "#b9bdc8", GREEN])
    ymax = float(agg.groupby("codon")["n"].sum().max()) if len(agg) else 1.0
    dom = pd.DataFrame([dict(ad=n, start=a_, end=b_, y0=-0.16 * ymax, y1=-0.05 * ymax) for n, a_, b_ in DOMAINS])
    xs = alt.Scale(domain=[1, 393], nice=False, clamp=True)
    bars = alt.Chart(agg).mark_bar(width=2).encode(
        x=alt.X("codon:Q", scale=xs, title="Kodon (alttaki şerit: VCEP domainleri)"),
        y=alt.Y("n:Q", title="Varyant sayısı", stack=True),
        color=alt.Color("grup:N", scale=cmap, title="Sınıf"), tooltip=["codon", "grup", "n"])
    track = alt.Chart(dom).mark_rect(opacity=.85).encode(
        x=alt.X("start:Q", scale=xs), x2="end:Q", y="y0:Q", y2="y1:Q",
        color=alt.Color("ad:N", title="Domain", scale=alt.Scale(scheme="tableau10")), tooltip=["ad", "start", "end"])
    here = pd.DataFrame({"codon": [float(r["codon"])], "etiket": [f"▼ {P1 or C}"], "y": [ymax * 1.12]}) \
        if s(r.get("codon")) else pd.DataFrame({"codon": [], "etiket": [], "y": []})
    rule = alt.Chart(here).mark_rule(color="#7b2cbf", strokeDash=[5, 3], size=2.5).encode(x=alt.X("codon:Q", scale=xs))
    label = alt.Chart(here).mark_text(color="#7b2cbf", fontWeight="bold", fontSize=13, dy=-4).encode(
        x=alt.X("codon:Q", scale=xs), y=alt.Y("y:Q"), text="etiket:N")
    st.altair_chart(alt.layer(bars, track, rule, label).resolve_scale(color="independent").properties(height=280),
                    width="stretch")
    st.caption("Mor kesik çizgi ve etiket: aradığın varyantın proteindeki yeri.")
    up_pos = [x for x in gene["uniprot_features"] if x["type"] in ("Binding site", "Site") and s(r.get("codon"))
              and abs(float(x["start"]) - float(r["codon"])) <= 3]
    if up_pos:
        st.caption("Yakın UniProt fonksiyonel bölgeleri: " + "; ".join(f"{x['type']} {x['start']} {x['description']}" for x in up_pos))

    st.markdown("#### Klinik kanıt")
    explain("va_clinical", "Bu kartlar ne anlatıyor?")
    c1, c2, c4 = st.columns(3)
    with c1:
        st.markdown(f"""<div class="card"><b>ClinGen eRepo (uzman paneli)</b><br>
        {s(r.get('erepo_class')) or 'Kayıt yok'}<br><span class="muted">{(s(r.get('erepo_met')) or '').replace(';', ', ')}</span>
        {'<br><a href="' + r['erepo_url'] + '" target="_blank">kaydı aç</a>' if s(r.get('erepo_url')) else ''}</div>""",
                    unsafe_allow_html=True)
    with c2:
        st.markdown(f"""<div class="card"><b>ClinVar</b> {s(r.get('clinvar_vcv')) or ''}<br>
        Germline: {s(r.get('clinvar_germline')) or '—'} <span class="muted">({s(r.get('clinvar_review')) or '—'},
        {f(r.get('clinvar_n_scv'), 0)} gönderim, son: {s(r.get('clinvar_last_eval')) or '—'})</span><br>
        Onkojenite: {s(r.get('clinvar_oncogenicity')) or '—'} · Somatik: {s(r.get('clinvar_somatic_tier')) or '—'}<br>
        <span class="muted">{(s(r.get('clinvar_traits')) or '')[:300]}</span></div>""", unsafe_allow_html=True)
    with c4:
        st.markdown(f"""<div class="card"><b>UniProt P04637</b><br>{s(r.get('uniprot_var')) or 'VAR_ kaydı yok'}<br>
        <span class="muted">UniProt notu: {(s(r.get('uniprot_desc')) or '—')[:250]}</span><br>
        <span class="muted">Derlenen klinik anlam (ClinVar/Ensembl, kriterli kayıtlar): {s(r.get('uniprot_sig')) or '—'}</span><br>
        <b>NCI germline</b>: {f(r.get('nci_germ_families'), 0)} aile, {f(r.get('nci_germ_individuals'), 0)} birey<br>
        <span class="muted">Birey kayıt sınıfları: {s(r.get('nci_germ_classes')) or '—'}</span></div>""",
                    unsafe_allow_html=True)

    st.markdown("#### Tahminler")
    explain("va_pred", "Tahminler grafiği nasıl okunur?")
    dist = missense_distributions()
    rows = []
    for col, name, dirn, hi, lo, note in TOOLS:
        v = pd.to_numeric(r.get(col), errors="coerce") if col in r else None
        v = None if v is None or pd.isna(v) else float(v)
        pct = None
        if v is not None and col in dist and len(dist[col]):
            arr = dist[col]
            # orta sira (midrank): esit degerler yuzdeligi sisirmesin
            lo_, le_ = ((arr < v).mean(), (arr <= v).mean()) if dirn == 1 else ((arr > v).mean(), (arr >= v).mean())
            pct = float((lo_ + le_) / 2 * 100)
        rows.append(dict(Araç=name, Skor=v, Yorum=verdict(v, dirn, hi, lo), **{"TP53 missense yüzdelik": pct},
                         **{"Ne ölçer?": TOOL_INFO.get(col, "")},
                         Eşik=f"{'>=' if dirn == 1 else '<='}{hi}", Not=note))
    pt = pd.DataFrame(rows)
    if s(r.get("agvgd")):
        st.caption(f"Align-GVGD sınıfı: {r['agvgd']} · VCEP Tablo S2 kodu: {s(r.get('s2_code')) or '—'} · "
                   f"MutPred2 mekanizma: {s(r.get('mutpred_mech')) or '—'}")
    chart_df = pt.dropna(subset=["TP53 missense yüzdelik"])
    if len(chart_df):
        ch = alt.Chart(chart_df).mark_bar().encode(
            y=alt.Y("Araç:N", sort=None, title=None), x=alt.X("TP53 missense yüzdelik:Q", scale=alt.Scale(domain=[0, 100]),
                                                             title="Sıralama: bu varyant TP53 varyantlarının yüzde kaçından daha zararlı bulundu (0-100)"),
            color=alt.Color("Yorum:N", scale=alt.Scale(domain=["Zararlı", "Belirsiz", "Benign"], range=[RED, "#c9cdd6", GREEN])),
            tooltip=["Araç", "Ne ölçer?", alt.Tooltip("Skor:Q", format=".3f"), "Yorum",
                     alt.Tooltip("TP53 missense yüzdelik:Q", format=".0f"), "Eşik"]).properties(height=26 * len(chart_df))
        ch = ch.configure_axisY(labelLimit=260, labelOverlap=False)
        st.altair_chart(ch, width="stretch")
        vc = pt[pt.Yorum.isin(["Zararlı", "Belirsiz", "Benign"])].Yorum.value_counts()
        st.markdown(f"**Özet:** skoru olan {int(vc.sum())} aracın **{vc.get('Zararlı', 0)}** tanesi zararlı, "
                    f"**{vc.get('Benign', 0)}** tanesi zararsız, **{vc.get('Belirsiz', 0)}** tanesi kararsız diyor. "
                    f"ACMG'de kullanılan tek tahmin: VCEP Tablo S2 kodu **{s(r.get('s2_code')) or 'yok'}**.")
    st.dataframe(pt, hide_index=True, width="stretch",
                 column_config={"Skor": st.column_config.NumberColumn(format="%.4f"),
                                "TP53 missense yüzdelik": st.column_config.ProgressColumn(min_value=0, max_value=100, format="%.0f")})

    st.markdown("#### Popülasyon frekansları")
    explain("va_pop", "Popülasyon frekansı ne demek?")
    pops = json.loads(s(r.get("gnomad_pops")) or "{}")
    names = dict(afr="Afrika/Afro-Amerikan", amr="Latin Amerika", asj="Aşkenazi Yahudi*", eas="Doğu Asya", fin="Fin*",
                 mid="Orta Doğu*", nfe="Avrupa (Fin olmayan)", sas="Güney Asya", ami="Amish*", remaining="Diğer*")
    ptab = pd.DataFrame([dict(Popülasyon=names.get(k, k), AC=v["ac"], AN=v["an"], AF=v["af"]) for k, v in pops.items()])
    c1, c2 = st.columns([1.3, 1])
    with c1:
        if len(ptab):
            st.dataframe(ptab.sort_values("AF", ascending=False), hide_index=True, width="stretch",
                         column_config={"AF": st.column_config.NumberColumn(format="%.2e")})
            st.caption("* VCEP frekans kurallarında dışlanan kurucu/karışık gruplar. gnomAD v4.1 exome+genome toplamı.")
        else:
            st.info("gnomAD v4.1'de bu varyant gözlenmemiş.")
    with c2:
        st.markdown(f"""<div class="card"><b>gnomAD v4.1 toplam AF:</b> {f(r.get('gnomad_af'))}<br>
        <b>FAF95 popmax:</b> {f(r.get('gnomad_faf95_popmax'))} ({s(r.get('gnomad_faf95_pop')) or '—'})<br>
        <b>Homozigot:</b> {f(r.get('gnomad_hom'), 0)} · <b>Filtre:</b> {s(r.get('gnomad_filters')) or 'PASS'}<br>
        <b>gnomAD v2 exome / genome:</b> {f(r.get('gnomad2_exome_af'))} / {f(r.get('gnomad2_genome_af'))}<br>
        <b>ExAC:</b> {f(r.get('exac_af'))}<hr style="margin:6px 0">
        <span class="muted">VCEP eşikleri: PM2_Supporting AF &lt; 0.00003 (grup &lt; 0.00004) · BS1 FAF 0.0003–0.001 ·
        BA1 FAF ≥ 0.001 (grup AN ≥ 2000, AC ≥ 2)</span></div>""", unsafe_allow_html=True)

# ====================================================================== 3) Fonksiyon & ML
with tabs[2]:
    st.markdown("#### Fonksiyonel test sonuçları (VCEP Tablo S3)")
    explain("fn_tests", "Fonksiyonel test nedir?")
    fa = pd.DataFrame([
        dict(Test="Kato 2003 – maya transaktivasyonu (8 promotör)", Sonuç=s(r.get("kato")) or s(r.get("kato_nci")) or "—",
             Ölçüt="≤%20 fonksiyonsuz, %20-75 kısmen, >%75 fonksiyonel"),
        dict(Test="Funk 2025 – CRISPR saturasyon mutajenezi", Sonuç=s(r.get("funk")) or "—", Ölçüt="RFS ≥ 0 LOF"),
        dict(Test="Giacomelli 2018 – A549 büyüme baskılama (etoposid)", Sonuç=s(r.get("giacomelli")) or "—", Ölçüt="Z ≤ -0.21 LOF"),
        dict(Test="Kotler 2018 – H1299 göreli uygunluk", Sonuç=s(r.get("kotler")) or "—", Ölçüt="RFS ≥ -1.0 LOF"),
        dict(Test="Kawaguchi 2005 – oligomerizasyon (maya)", Sonuç=s(r.get("kawaguchi")) or "—", Ölçüt="monomer/dimer anormal"),
        dict(Test="NCI: dominant negatif / LOF sınıfı", Sonuç=s(r.get("dne_lof_class")) or "—", Ölçüt="NCI TP53 r21"),
        dict(Test="NCI: yapı-fonksiyon sınıfı", Sonuç=s(r.get("structure_function_class")) or "—", Ölçüt="NCI TP53 r21"),
    ])
    st.dataframe(fa, hide_index=True, width="stretch")
    st.caption(f"VCEP ön-kodu: {s(r.get('s3_code')) or '—'}")
    nf = X["nci_function"]
    if P1 and len(nf):
        sub = nf[nf["p1"] == P1]
        if len(sub):
            with st.expander(f"NCI TP53 fonksiyonel deney kayıtları ({len(sub)})"):
                keep = [c for c in ["Loss_of_Function", "Dominant_Negative_Activity", "Gain_of_Function", "Temperature_Sensitivity",
                                    "Cell_assay", "Cell_lines", "Assay_design", "Method", "Authors", "Year", "Title", "PubMed"] if c in sub]
                st.dataframe(sub[keep], hide_index=True, width="stretch")

    st.markdown("#### Makine öğrenmesi")
    if SHOW_HELP:
        h1, h2 = st.columns(2)
        h1.info(HELP["ml_a"])
        h2.info(HELP["ml_b"])
    if s(r.get("ml_path_score")):
        c1, c2 = st.columns(2)
        with c1:
            st.markdown(f"""<div class="card"><b>TP53 Aggregated Prediction (Model A)</b><br>
            <span class="cls" style="color:{RED if float(r['ml_path_score']) >= .5 else GREEN}">{float(r['ml_path_score']):.3f}</span>
            <span class="muted"> · TP53 missense yüzdeliği %{100 * float(r['ml_path_pct']):.0f}</span><br>
            <span class="muted">Franklin'in toplu tahmin skorunun TP53'e özel karşılığı: in-silico skorlar → kalibre lojistik/GB model.
            Eğitimdeki varyantlar için çapraz doğrulama dışı (sızıntısız) skor gösterilir.</span></div>""", unsafe_allow_html=True)
        with c2:
            pb = pd.DataFrame([dict(Model=m_, Fonksiyonel=float(r.get(f"ml_{m_}_functional") or 0),
                                    **{"Kısmen fonksiyonel (hipomorf)": float(r.get(f"ml_{m_}_partial") or 0)},
                                    Fonksiyonsuz=float(r.get(f"ml_{m_}_nonfunctional") or 0)) for m_ in ("B1", "B2")])
            st.markdown("**Model B: Kato sınıfı olasılıkları**")
            order = ["Fonksiyonel", "Kısmen fonksiyonel (hipomorf)", "Fonksiyonsuz"]
            pl = pb.melt(id_vars="Model", var_name="sinif", value_name="olasilik")
            pl["sira"] = pl["sinif"].map({k_: i_ for i_, k_ in enumerate(order)})
            st.altair_chart(alt.Chart(pl).mark_bar().encode(
                y=alt.Y("Model:N", title=None),
                x=alt.X("olasilik:Q", stack="zero", scale=alt.Scale(domain=[0, 1]), title="Olasılık (toplam 1)"),
                color=alt.Color("sinif:N", scale=alt.Scale(domain=order, range=[GREEN, "#e8a33d", RED]), legend=None),
                order=alt.Order("sira:Q"),
                tooltip=["Model", alt.Tooltip("sinif:N", title="Sınıf"), alt.Tooltip("olasilik:Q", title="Olasılık", format=".2f")])
                .properties(height=alt.Step(30)), width="stretch")
            st.caption("Yeşil: fonksiyonel · turuncu: kısmen fonksiyonel (hipomorf) · kırmızı: fonksiyonsuz")
            st.dataframe(pb.round(2), hide_index=True, width="stretch")
            st.caption("B1: sekans + yapı + in-silico. B2: + Giacomelli/Funk/Kotler/Kawaguchi. 'Kısmen fonksiyonel' olasılığı "
                       "hipomorf aday skorudur.")
    else:
        st.info("ML tahmini yalnızca missense varyantlar için üretilir.")
    if METRICS:
        with st.expander("Model performansı (kodon gruplu 5 katlı çapraz doğrulama)", expanded=False):
            if SHOW_HELP:
                st.info(HELP["ml_perf"])
            A = METRICS["A_pathogenicity"]
            st.markdown(f"**Model A** – {A['n']} etiketli missense ({A['n_pos']} P/LP, {A['n_neg']} B/LB; "
                        f"{A['label_sources']}). En iyi: **{A['best']}**")
            perf = [dict(Yöntem=f"TP53 ML ({k})", ROC_AUC=v["roc_auc"], Kapsam=1.0) for k, v in A["models"].items()]
            perf += [dict(Yöntem=k, ROC_AUC=v["roc_auc"], Kapsam=v["coverage"]) for k, v in A["single_tools"].items()]
            st.dataframe(pd.DataFrame(perf).sort_values("ROC_AUC", ascending=False), hide_index=True)
            st.caption("Tek araç AUC'leri yalnızca skoru olan alt kümede hesaplanmıştır. Not: REVEL/VARITY gibi araçlar "
                       "ClinVar etiketleriyle eğitildiği için döngüsellik (circularity) riski vardır.")
            for k in ("B1_function", "B2_function"):
                if k in METRICS:
                    B = METRICS[k]
                    st.markdown(f"**Model {B['kind']}** – n={B['n']}, macro-F1 {B['macro_f1']}, dengeli doğruluk {B['balanced_acc']}, "
                                f"AUC(kısmen vs diğer) {B['auc_partial_vs_rest']}, AUC(fonksiyonsuz vs diğer) {B['auc_nonfunc_vs_rest']}")
                    st.dataframe(pd.DataFrame(B["confusion"], index=["gerçek F", "gerçek Kısmen", "gerçek Fsuz"],
                                              columns=["tahmin F", "tahmin Kısmen", "tahmin Fsuz"]))
            st.caption("ESM-2 maskeli marjinal LLR hesaplanmışsa (scripts/esm2_llr.py) özellik olarak otomatik eklenir.")

# ====================================================================== 4) Yayinlar
with tabs[3]:
    st.markdown("#### Varyanta özgü yayınlar")
    explain("pubs")
    gl = X["nci_germline"]
    refs = gl[gl["c"] == C] if len(gl) else pd.DataFrame()
    if len(refs):
        st.markdown(f"**NCI TP53 germline veritabanı** – {refs['Family_ID'].nunique()} aile, {len(refs)} birey kaydı")
        st.dataframe(refs[[c for c in ["Family_ID", "Class", "Country", "Sex", "Short_topo", "Age_at_diagnosis", "PubMedLink"] if c in refs]]
                     .drop_duplicates().head(200), hide_index=True, width="stretch")
    if P1 and st.button("PubMed'de ara (canlı)"):
        from tp53.llm import pubmed_search
        try:
            pm = pubmed_search(P1, P3, C)
            st.session_state[f"pm|{C}"] = pm
        except Exception as ex:
            st.error(f"PubMed sorgusu başarısız: {ex}")
    pm = st.session_state.get(f"pm|{C}")
    if pm:
        st.markdown(f"Sorgu: `{pm['query']}` → **{pm['count']}** makale (ilk {len(pm['articles'])})")
        for a_ in pm["articles"]:
            st.markdown(f"- [{a_['title']}]({a_['url']}) — {a_['authors']} · *{a_['journal']}* {a_['year']}")
    if s(r.get("clinvar_id")):
        st.markdown(f"[ClinVar gönderimleri ve atıflar](https://www.ncbi.nlm.nih.gov/clinvar/variation/{r['clinvar_id']}/)")

# ====================================================================== 5) Gen degerlendirme
with tabs[4]:
    st.markdown("#### TP53 (tümör protein p53) — 17p13.1, NM_000546.6, P04637, 393 aa")
    st.markdown("Otozomal dominant **Li-Fraumeni sendromu** (MONDO:0018875). Hastalık mekanizması: fonksiyon kaybı ve "
                "missense varyantlarda dominant-negatif / fonksiyon kazanımı. ClinGen TP53 VCEP (affiliation 50013) "
                "gene özgü kuralları geçerlidir.")
    explain("gene_dist", "Bu sekme ne anlatıyor?")
    g = df.copy()
    g["sig"] = g.erepo_class.fillna(g.clinvar_germline).fillna("")
    g["grup"] = "VUS/çelişkili/yok"
    g.loc[g.sig.str.contains("athogenic") & ~g.sig.str.contains("onflict|enign"), "grup"] = "P/LP"
    g.loc[g.sig.str.contains("enign") & ~g.sig.str.contains("onflict|athogenic"), "grup"] = "B/LB"
    tab = g[g.clinvar_vcv.notna() | g.erepo_class.notna()].groupby(["consequence", "grup"]).size().rename("n").reset_index()
    st.altair_chart(alt.Chart(tab).mark_bar().encode(
        y=alt.Y("consequence:N", sort="-x", title=None), x=alt.X("n:Q", title="ClinVar/eRepo varyant sayısı"),
        color=alt.Color("grup:N", scale=alt.Scale(domain=["P/LP", "VUS/çelişkili/yok", "B/LB"], range=[RED, "#c9cdd6", GREEN])),
        tooltip=["consequence", "grup", "n"]).properties(height=28 * tab.consequence.nunique(),
                                                          title="Sonuç tipine göre sınıflama dağılımı")
        .configure_axisY(labelOverlap=False, labelLimit=200),
        width="stretch")
    explain("gene_hotspot", "Hotspot grafiği ne anlatıyor?")
    hs = pd.DataFrame([dict(kalıntı=v["residue"], tümör=v["tumorCount"]) for v in gene["hotspots"].values()])
    if len(hs):
        st.altair_chart(alt.Chart(hs.nlargest(25, "tümör")).mark_bar(color=TEAL).encode(
            x=alt.X("kalıntı:N", sort="-y", title=None), y=alt.Y("tümör:Q", title="Tümör sayısı"), tooltip=["kalıntı", "tümör"])
            .properties(height=220, title="cancerhotspots.org v2 – en sık TP53 hotspot kalıntıları"), width="stretch")
    feats = pd.DataFrame(gene["uniprot_features"])
    if len(feats):
        with st.expander("UniProt P04637 bölge, motif ve bağlanma yerleri"):
            if SHOW_HELP:
                st.info(HELP["gene_uniprot"])
            st.dataframe(feats[feats.type.isin(["Region", "Motif", "DNA binding", "Binding site", "Site"])],
                         hide_index=True, width="stretch")

# ====================================================================== 6) Iliskili hastaliklar
with tabs[5]:
    explain("dis_intro", "Bu sekme ne anlatıyor?")
    st.markdown("#### Bu varyant için bildirilen durumlar")
    tr = traits_of(r.get("clinvar_traits"))
    if tr:
        st.markdown(" ".join(f'<span class="pill">{t}</span>' for t in tr), unsafe_allow_html=True)
        st.caption("ClinVar'da laboratuvarların bu varyantı değerlendirdiği durumlar; varyantın bu hastalığa yol açtığı anlamına gelmez.")
    else:
        st.caption("ClinVar germline trait kaydı yok.")
    gl = X["nci_germline"]
    sub = gl[gl["c"] == C] if len(gl) else pd.DataFrame()
    if len(sub):
        c1, c2 = st.columns(2)
        with c1:
            cls = sub.drop_duplicates("Family_ID")["Class"].value_counts().rename_axis("Aile sınıfı").reset_index(name="aile")
            st.markdown("**NCI germline aile sınıfları**")
            st.dataframe(cls, hide_index=True, width="stretch")
        with c2:
            tp = sub[sub["Short_topo"].notna()]["Short_topo"].value_counts().head(15).rename_axis("Tümör").reset_index(name="n")
            st.markdown("**Taşıyıcılarda görülen tümörler (NCI)**")
            st.altair_chart(alt.Chart(tp).mark_bar(color=RED).encode(y=alt.Y("Tümör:N", sort="-x", title=None), x="n:Q")
                            .properties(height=22 * max(len(tp), 3)), width="stretch")
    st.markdown("#### TP53 ile ilişkili durumlar (gen düzeyi, ClinVar trait sıklığı)")
    allt = df["clinvar_traits"].dropna().map(traits_of).explode().dropna()
    st.dataframe(allt.value_counts().head(25).rename_axis("Durum").reset_index(name="TP53 varyant sayısı"),
                 hide_index=True, width="stretch")
    st.caption("Aynı hastalığın farklı yazımları (ör. 'Li-Fraumeni syndrome 1') tek satırda birleştirildi.")
    st.caption("Ana tablo: Li-Fraumeni sendromu (OMIM 151623, MONDO:0018875, AD). VCEP LFS kanser listesi ACMG sekmesindeki "
               "klinik formda yer alır.")

# ====================================================================== 7) Somatik
with tabs[6]:
    explain("som_intro", "Bu sekme ne anlatıyor?")
    c1, c2 = st.columns([1.2, 1])
    topo = X["nci_tumor_topography"]
    with c1:
        st.markdown(f"#### NCI TP53 somatik veritabanı – {f(r.get('nci_tumor_n'), 0)} tümör")
        t = topo[topo["c"] == C].head(20) if len(topo) else pd.DataFrame()
        if len(t):
            st.altair_chart(alt.Chart(t).mark_bar(color=NAVY).encode(y=alt.Y("Short_topo:N", sort="-x", title=None), x="n:Q",
                                                                       tooltip=["Short_topo", "n"]).properties(height=22 * len(t)),
                            width="stretch")
        st.caption(f"TCGA/ICGC/GENIE sayısı: {f(r.get('tcga_icgc_genie_count'), 0)} · COSMIC: {f(r.get('cosmic_id'), 0) if s(r.get('cosmic_id')) else '—'} · "
                   f"NCI hotspot: {s(r.get('nci_hotspot')) or '—'}")
    with c2:
        h = gene["hotspots"].get(str(int(float(r["codon"])))) if s(r.get("codon")) else None
        st.markdown("#### cancerhotspots.org v2")
        if h:
            st.markdown(f"Kalıntı **{h['residue']}**: {h['tumorCount']} tümör, q={h['qValue']}<br>"
                        f"Değişimler: {', '.join(f'{k}:{v}' for k, v in sorted(h['variants'].items(), key=lambda x: -x[1]))}",
                        unsafe_allow_html=True)
            comp = pd.DataFrame(list(h["composition"].items()), columns=["doku", "n"]).nlargest(12, "n")
            st.altair_chart(alt.Chart(comp).mark_bar(color=TEAL).encode(y=alt.Y("doku:N", sort="-x", title=None), x="n:Q")
                            .properties(height=20 * len(comp)), width="stretch")
        else:
            st.caption("Bu kalıntı anlamlı hotspot listesinde değil.")
    st.markdown("#### Klinik kanıt (CIViC, ClinVar somatik)")
    explain("som_civic", "CIViC tablosunun sütunları ne demek?")
    st.caption(f"ClinVar somatik klinik etki: {s(r.get('clinvar_somatic_tier')) or '—'} · Onkojenite: {s(r.get('clinvar_oncogenicity')) or '—'}")
    ce = X["civic_evidence"]
    if len(ce):
        scope_v = ce[ce.molecular_profile.str.fullmatch(f"TP53 {P1}", na=False)] if P1 else ce.iloc[0:0]
        gene_lvl = ce[ce.molecular_profile.str.contains(r"TP53 (?:Mutation|Deleterious|Loss|ALTERATION|Wildtype)", na=False, case=False)]
        st.markdown(f"**Varyant düzeyi:** {len(scope_v)} kanıt · **Gen düzeyi:** {len(gene_lvl)} kanıt")
        show = pd.concat([scope_v.assign(kapsam="Varyant"), gene_lvl.assign(kapsam="Gen")])
        if len(show):
            cols = ["kapsam", "evidence_level", "evidence_type", "significance", "disease", "therapies", "evidence_direction",
                    "molecular_profile", "evidence_statement", "citation", "evidence_civic_url"]
            st.dataframe(show[cols].sort_values(["kapsam", "evidence_level"], ascending=[False, True]), hide_index=True,
                         width="stretch", column_config={"evidence_civic_url": st.column_config.LinkColumn("CIViC")})
    st.caption("Tedavi bilgisi yalnızca literatür/kanıt özetidir; tedavi önerisi değildir.")

# ====================================================================== 8) LLM
with tabs[7]:
    from tp53 import llm
    explain("llm_intro", "Bu sekme ne yapıyor? (3 adım)")
    models = llm.ollama_models()
    if not models:
        st.warning("Ollama'ya ulaşılamadı (http://localhost:11434). `ollama serve` çalışıyor mu?")
    model = st.selectbox("Yerel model (Ollama)", models or [llm.DEFAULT_MODEL],
                         index=(models.index(llm.DEFAULT_MODEL) if llm.DEFAULT_MODEL in models else 0))
    st.markdown("Franklin, ClinVar gönderimlerindeki serbest metinden de novo, segregasyon ve vaka kanıtı çıkarır. "
                "Bu sekme aynı işi yerel bir LLM ile yapar; çıkan sayılar **öneridir**, onaylayınca klinik forma aktarılır.")
    if not s(r.get("clinvar_id")):
        st.info("Bu varyantın ClinVar kaydı yok.")
    else:
        if st.button("1) ClinVar gönderimlerini getir"):
            try:
                st.session_state[f"subs|{C}"] = llm.clinvar_submissions(r["clinvar_id"])
            except Exception as ex:
                st.error(f"ClinVar efetch başarısız: {ex}")
        subs = st.session_state.get(f"subs|{C}")
        if subs:
            sd = pd.DataFrame(subs)
            st.markdown(f"{len(sd)} gönderim, {(sd.comment.str.len() > 0).sum()} tanesinde serbest metin var.")
            st.dataframe(sd[["scv", "submitter", "classification", "date", "origin", "comment"]], hide_index=True,
                         width="stretch")
            if st.button("2) LLM ile klinik kanıt çıkar", disabled=not models):
                with st.spinner(f"{model} çalışıyor..."):
                    try:
                        st.session_state[f"llm|{C}"] = llm.extract_evidence(f"{C} ({P3})", subs, model)
                    except Exception as ex:
                        st.error(f"LLM hatası: {ex}")
        ev = st.session_state.get(f"llm|{C}")
        if ev:
            if "error" in ev:
                st.error(ev["error"])
            else:
                st.markdown(f"**Özet ({ev.get('_model')}):** {ev.get('summary_tr', '')}")
                sug = llm.to_clinical_suggestion(ev)
                st.json({k: v for k, v in ev.items() if not k.startswith("_")}, expanded=False)
                st.markdown("**Klinik form önerisi:** " + ", ".join(f"{k}={v}" for k, v in sug.items() if v))
                if st.button("3) Öneriyi klinik forma uygula (eRepo uzman kodları yerine)"):
                    dn = sug["denovo"]
                    st.session_state["pending_clinical"] = dict(
                        dn_strong_conf=sum(1 for g_, c_ in dn if g_ == LFS_STRONG and c_),
                        dn_strong_assumed=sum(1 for g_, c_ in dn if g_ == LFS_STRONG and not c_),
                        dn_mod_conf=sum(1 for g_, c_ in dn if g_ == LFS_MODERATE and c_),
                        dn_mod_assumed=sum(1 for g_, c_ in dn if g_ == LFS_MODERATE and not c_),
                        chompret=sug["chompret"], classic=sug["classic_lfs"], her2=sug["her2_lt40"], meioses=sug["meioses"],
                        families=sug["families"], vaf525=sug["vaf_5_25"], bs2=sug["bs2_females60"], use_expert=False)
                    st.rerun()
    st.divider()
    st.markdown("##### İsteğe bağlı: rapor paragrafı")
    st.caption("Yukarıdaki adımlardan bağımsızdır. ClinVar metnini okumaz; ACMG sekmesindeki son sınıfı ve karşılanan "
               "kriterleri bir genetik raporundaki yorum paragrafı gibi Türkçe yazar. Puanı değiştirmez.")
    if st.button("Rapor paragrafı yaz (LLM)", disabled=not models):
        with st.spinner("Yazılıyor..."):
            try:
                st.markdown(llm.explain(f"{C} ({P3})", clf.vcep(C, clinical_from_state()), model))
            except Exception as ex:
                st.error(f"LLM hatası: {ex}")

# ====================================================================== 9) Veri kaynaklari
with tabs[8]:
    explain("data_intro", "Bu sekme ne anlatıyor?")
    st.markdown("#### İndirilen veri setleri (data/raw/MANIFEST.json)")
    if MANIFEST:
        st.dataframe(pd.DataFrame([dict(Kaynak=k, Dosya=v["file"], KB=round(v["bytes"] / 1024, 1), İndirme=v["downloaded_utc"],
                                        Açıklama=v.get("note", ""), URL=v["url"]) for k, v in MANIFEST.items()]),
                     hide_index=True, width="stretch", column_config={"URL": st.column_config.LinkColumn()})
    st.markdown("#### Kurallar")
    st.markdown("- ClinGen TP53 VCEP spesifikasyonu v2.4 (CSpec GN009): https://cspec.clinicalgenome.org/cspec/ui/svi/doc/GN009\n"
                "- Fortuno et al., Genome Medicine 2025 (Tablo S1-S3): https://doi.org/10.1186/s13073-025-01536-3\n"
                "- Franklin yardım makaleleri (genel mod): https://help.genoox.com/en/collections/2537985")
    v = X["validation"]
    if len(v):
        cl = v[v.code == "SINIF"]
        st.markdown(f"#### Doğrulama: ClinGen uzman sınıflamalarıyla uyum **%{100 * (cl.bizim == cl.uzman).mean():.1f}** "
                    f"({(cl.bizim == cl.uzman).sum()}/{len(cl)})")
        st.dataframe(pd.crosstab(cl.uzman, cl.bizim), width="stretch")
        st.caption("Farkların çoğu uzman kayıtlarının eski spesifikasyon sürümleriyle (v2.0-2.2) yapılmasından kaynaklanır "
                   "(ör. BP4_Moderate ve Funk 2025 testi sonradan eklendi).")
