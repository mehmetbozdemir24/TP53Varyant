"""LLM destekli kanit cikarimi (yerel Ollama) ve canli literatur sorgulari.

Franklin, ClinVar gonderimlerindeki serbest metinden de novo, segregasyon, vaka ve fonksiyonel kanit cikarip
PS2/PS4/PP1/PS3 kurallarini uyguluyor. Bu modul ayni isi yerel bir LLM ile yapar:
  1) ClinVar VCV XML'inden her gonderimin (SCV) yorum metnini ceker,
  2) metni Ollama'daki modele (varsayilan gemma4:12b) yapilandirilmis JSON cikarimi icin gonderir,
  3) sonucu VCEP klinik formuna (acmg.Clinical) donusturulecek oneriler olarak dondurur.
LLM ciktisi her zaman OneRI'dir; kullanici arayuzde onaylamadan puana eklenmez.
"""
from __future__ import annotations

import json
import os
import re
import xml.etree.ElementTree as ET

import requests

OLLAMA = os.environ.get("OLLAMA_HOST", "http://localhost:11434")
DEFAULT_MODEL = os.environ.get("TP53_LLM_MODEL", "gemma4:12b")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UA = {"User-Agent": "TP53-Franklin research tool"}


# ------------------------------------------------------------------ ClinVar gonderim metinleri
def clinvar_submissions(variation_id: str | int, timeout=60):
    """ClinVar VCV XML -> [{scv, submitter, classification, date, comment, pmids, method, origin}]"""
    r = requests.get(EUTILS + "efetch.fcgi", params={"db": "clinvar", "rettype": "vcv", "is_variationid": "", "id": str(variation_id)},
                     headers=UA, timeout=timeout)
    r.raise_for_status()
    root = ET.fromstring(r.content)
    out = []
    for ca in root.iter("ClinicalAssertion"):
        acc = ca.find("ClinVarAccession")
        cls = ca.find(".//Classification")
        desc = None
        if cls is not None:
            for tag in ("GermlineClassification", "SomaticClinicalImpact", "OncogenicityClassification"):
                d = cls.find(tag)          # ElementTree: bos eleman False doner, 'or' zinciri kullanilmaz
                if d is not None:
                    desc = d.text
                    break
        comments = [c.text.strip() for c in ca.iter("Comment") if c.text and c.text.strip()]
        for at in ca.iter("Attribute"):
            if at.get("Type") in ("Description", "CompoundHeterozygote", "SingleHeterozygote", "Homozygote",
                                  "VariantAlleles", "IndividualsWithVariant") and at.text:
                comments.append(f"{at.get('Type')}: {at.text.strip()}")
        pmids = sorted({i.text for i in ca.iter("ID") if i.get("Source") == "PubMed" and i.text})
        methods = sorted({m.text for m in ca.iter("MethodType") if m.text})
        origins = sorted({o.text for o in ca.iter("Origin") if o.text})
        out.append(dict(scv=acc.get("Accession") if acc is not None else None,
                        submitter=acc.get("SubmitterName") if acc is not None else None,
                        date=(cls.get("DateLastEvaluated") if cls is not None else None),
                        classification=desc, comment=" ".join(comments)[:6000], pmids=pmids,
                        method=", ".join(methods), origin=", ".join(origins)))
    return out


# ------------------------------------------------------------------ PubMed
def pubmed_search(p1: str, p3: str | None = None, cdna: str | None = None, retmax=30, timeout=30):
    terms = [f'"{x}"' for x in (p1, p3[2:] if p3 and p3.startswith("p.") else p3, cdna) if x]
    q = f"TP53[tiab] AND ({' OR '.join(terms)})"
    r = requests.get(EUTILS + "esearch.fcgi", params={"db": "pubmed", "term": q, "retmax": retmax, "retmode": "json",
                                                       "sort": "relevance"}, headers=UA, timeout=timeout)
    r.raise_for_status()
    j = r.json()["esearchresult"]
    ids, count = j.get("idlist", []), int(j.get("count", 0))
    arts = []
    if ids:
        s = requests.get(EUTILS + "esummary.fcgi", params={"db": "pubmed", "id": ",".join(ids), "retmode": "json"},
                         headers=UA, timeout=timeout).json()["result"]
        for i in ids:
            d = s.get(i, {})
            arts.append(dict(pmid=i, title=d.get("title"), journal=d.get("fulljournalname") or d.get("source"),
                             year=(d.get("pubdate") or "")[:4], authors=", ".join(a["name"] for a in d.get("authors", [])[:3]),
                             url=f"https://pubmed.ncbi.nlm.nih.gov/{i}/"))
    return dict(query=q, count=count, articles=arts)


# ------------------------------------------------------------------ Ollama
def ollama_models(timeout=4):
    try:
        j = requests.get(f"{OLLAMA}/api/tags", timeout=timeout).json()
        return [m["name"] for m in j.get("models", [])]
    except Exception:
        return []


SCHEMA = {
    "type": "object",
    "properties": {
        "de_novo_cases": {"type": "array", "items": {"type": "object", "properties": {
            "cancer": {"type": "string"}, "age": {"type": ["number", "null"]},
            "parentage_confirmed": {"type": "boolean"}, "quote": {"type": "string"}}}},
        "probands_classic_lfs": {"type": "integer"},
        "probands_chompret": {"type": "integer"},
        "her2_breast_lt40": {"type": "integer"},
        "segregation_meioses": {"type": "integer"},
        "segregation_families": {"type": "integer"},
        "low_vaf_observations": {"type": "integer"},
        "healthy_female_carriers_over60": {"type": "integer"},
        "functional_studies": {"type": "array", "items": {"type": "object", "properties": {
            "assay": {"type": "string"}, "result": {"type": "string"}, "quote": {"type": "string"}}}},
        "lof_or_gof_notes": {"type": "string"},
        "evidence_quotes": {"type": "array", "items": {"type": "string"}},
        "summary_tr": {"type": "string"},
    },
    "required": ["probands_classic_lfs", "probands_chompret", "segregation_meioses", "summary_tr"],
}

PROMPT = """You are a clinical variant curation assistant applying the ClinGen TP53 VCEP v2.4 rules.
Below are ClinVar submitter comments for the germline TP53 variant {variant}.
Extract ONLY what is explicitly stated. Do not guess. Count unrelated probands, not families twice.
- de_novo_cases: probands where the variant arose de novo; parentage_confirmed=true only if both maternity and paternity confirmed.
- probands_classic_lfs / probands_chompret: number of unrelated probands described as meeting Classic LFS or (revised) Chompret criteria.
  If a text only says "individuals with Li-Fraumeni syndrome", count them as chompret unless "classic" is stated.
- her2_breast_lt40: unrelated probands with HER2-positive breast cancer diagnosed before age 40.
- segregation_meioses: informative meioses where the variant co-segregates with LFS cancers; segregation_families: number of families.
- low_vaf_observations: observations with variant allele fraction 5-35 % (possible clonal hematopoiesis).
- healthy_female_carriers_over60: unrelated women aged >= 60 without cancer carrying the variant.
- functional_studies: any functional assay results mentioned (transactivation, yeast, growth suppression, apoptosis, dominant negative).
- evidence_quotes: short verbatim quotes supporting each number.
- summary_tr: 3-5 sentence summary in Turkish.
Return JSON only.

TEXT:
{text}
"""


def extract_evidence(variant: str, submissions: list, model: str | None = None, timeout=900, max_chars=12000):
    """Gonderim metinlerinden yapilandirilmis klinik kanit cikarir (Ollama /api/chat, JSON sema ciktisi)."""
    model = model or DEFAULT_MODEL
    kw = re.compile(r"de novo|segregat|famil|proband|Li-?Fraumeni|Chompret|LFS|individual|patient|affected|"
                    r"carrier|mosaic|allele fraction|VAF|HER2|functional|transactivation|yeast|dominant", re.I)
    parts = []
    for s in submissions:
        c = (s.get("comment") or "").strip()
        if len(c) < 80 or re.fullmatch(r"(Description: )?not provided\.?", c, re.I):
            continue
        parts.append((len(kw.findall(c)), f"[{s.get('scv')} | {s.get('submitter')} | {s.get('classification')}]\n{c[:2500]}"))
    parts.sort(key=lambda x: -x[0])          # en bilgili gonderimler once
    text, used = "", 0
    for _, p in parts:
        if len(text) + len(p) > max_chars:
            break
        text += p + "\n\n"
        used += 1
    parts = parts[:used]
    if not text:
        return {"error": "Gonderimlerde serbest metin yok."}
    body = {"model": model, "stream": False, "think": False, "format": SCHEMA, "keep_alive": "10m",
            "options": {"temperature": 0, "num_ctx": 8192, "num_predict": 1500},
            "messages": [{"role": "user", "content": PROMPT.format(variant=variant, text=text)}]}
    r = requests.post(f"{OLLAMA}/api/chat", json=body, timeout=timeout)
    r.raise_for_status()
    content = r.json()["message"]["content"]
    try:
        out = json.loads(content)
    except json.JSONDecodeError:
        m = re.search(r"\{.*\}", content, re.S)
        out = json.loads(m.group(0)) if m else {"error": "JSON ayrisamadi", "raw": content[:2000]}
    out["_model"] = model
    out["_n_submissions_with_text"] = len(parts)
    return out


def to_clinical_suggestion(ev: dict):
    """LLM ciktisini acmg.Clinical alanlarina cevirir (kullanici onayi icin)."""
    from .acmg import LFS_MODERATE, LFS_STRONG
    strong_kw = ("breast", "choroid", "adrenocortical", "rhabdomyosarcoma", "osteosarcoma")
    dn = []
    for c in ev.get("de_novo_cases") or []:
        cancer = (c.get("cancer") or "").lower()
        age = c.get("age")
        strong = (("breast" in cancer and age is not None and age < 31) or "choroid" in cancer or
                  ("adrenocortical" in cancer and (age is None or age < 18)) or
                  (("rhabdo" in cancer or "osteosarcoma" in cancer) and (age is None or age < 46)))
        dn.append((LFS_STRONG if strong else LFS_MODERATE, bool(c.get("parentage_confirmed"))))
    return dict(denovo=dn, chompret=int(ev.get("probands_chompret") or 0), classic_lfs=int(ev.get("probands_classic_lfs") or 0),
                her2_lt40=int(ev.get("her2_breast_lt40") or 0), meioses=int(ev.get("segregation_meioses") or 0),
                families=max(1, int(ev.get("segregation_families") or 1)), vaf_5_25=int(ev.get("low_vaf_observations") or 0),
                bs2_females60=int(ev.get("healthy_female_carriers_over60") or 0))


def explain(variant: str, vcep_result: dict, model: str | None = None, timeout=180):
    """Siniflandirma sonucunu Turkce, Franklin'deki ozet paragrafina benzer sekilde anlatir."""
    model = model or DEFAULT_MODEL
    met = [e for e in vcep_result["evidence"] if e["met"] and e["points"] != 0]
    lines = "\n".join(f"- {e['label']} ({e['points']:+d}): {e['reason']}" for e in met)
    prompt = (f"TP53 {variant} varyantinin ClinGen TP53 VCEP v2.4 ile siniflandirmasi: {vcep_result['classification']} "
              f"(toplam {vcep_result['points']} puan).\nKarsilanan kriterler:\n{lines}\n\n"
              "Bunu bir klinik genetik raporundaki yorum paragrafi gibi, 5-7 cumlelik akici Turkce ile ozetle. "
              "Yeni bilgi uydurma, sadece verilen kanitlari kullan. Kalin yazi ve madde isareti kullanma.")
    r = requests.post(f"{OLLAMA}/api/chat", json={"model": model, "stream": False, "think": False,
                                                   "options": {"temperature": 0.2, "num_predict": 800},
                                                   "messages": [{"role": "user", "content": prompt}]}, timeout=timeout)
    r.raise_for_status()
    return r.json()["message"]["content"]
