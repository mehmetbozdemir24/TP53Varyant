"""TP53-Franklin: tum ham veri setlerini data/raw altina indirir.

Calistirma:  python scripts/download.py            (hepsi)
             python scripts/download.py gnomad     (tek kaynak)
Her indirme data/raw/MANIFEST.json dosyasina URL, tarih, boyut ve md5 ile yazilir.
"""
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
RAW.mkdir(parents=True, exist_ok=True)
MANIFEST = RAW / "MANIFEST.json"
H = {"User-Agent": "Mozilla/5.0 (TP53-Franklin research tool)"}
S = requests.Session()
S.headers.update(H)


def log(name, url, path, note=""):
    m = json.loads(MANIFEST.read_text(encoding="utf-8")) if MANIFEST.exists() else {}
    data = path.read_bytes()
    m[name] = {
        "url": url,
        "file": path.name,
        "bytes": len(data),
        "md5": hashlib.md5(data).hexdigest(),
        "downloaded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "note": note,
    }
    MANIFEST.write_text(json.dumps(m, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"  ok  {name:24s} {len(data)/1024:9.1f} KB -> {path.name}")


def get_file(name, url, fname, note=""):
    r = S.get(url, timeout=180)
    r.raise_for_status()
    p = RAW / fname
    p.write_bytes(r.content)
    log(name, url, p, note)


# ---------------------------------------------------------------- kaynaklar
def vcep_tables():
    """ClinGen TP53 VCEP v2 ek tablolari (Genome Medicine 2025, Tablo S1-S3)."""
    base = ("https://static-content.springer.com/esm/art%3A10.1186%2Fs13073-025-01536-3/"
            "MediaObjects/13073_2025_1536_MOESM{}_ESM.xlsx")
    notes = {1: "S1 splice +/-1,2 SpliceAI + PVS1 gucu",
             2: "S2 tum missense icin Align-GVGD, BayesDel, SpliceAI, PP3/BP4 kodu",
             3: "S3 tum missense icin Kato/Funk/Giacomelli/Kotler/Kawaguchi + PS3/BS3 kodu"}
    for i in (1, 2, 3):
        get_file(f"vcep_S{i}", base.format(i), f"vcep_TableS{i}.xlsx", notes[i])


def nci_tp53():
    """NIH NCI TP53 Database r21: somatik, germline, fonksiyonel (Kato vb.)."""
    for kind in ("TumorVariantDownload", "GermlineDownload", "FunctionDownload"):
        url = f"https://tp53.cancer.gov/static/data/{kind}_r21.csv"
        get_file(f"nci_{kind}", url, f"nci_{kind}_r21.csv", "NCI TP53 Database release 21")


def uniprot():
    get_file("uniprot_P04637", "https://rest.uniprot.org/uniprotkb/P04637.json",
             "uniprot_P04637.json", "UniProt P53_HUMAN: domainler, dogal varyantlar")
    url = "https://www.ebi.ac.uk/proteins/api/variation/P04637"
    r = S.get(url, headers={"Accept": "application/json"}, timeout=300)
    r.raise_for_status()
    p = RAW / "uniprot_variation_P04637.json"
    p.write_bytes(r.content)
    log("uniprot_variation", url, p, "EBI Proteins API varyasyon (ClinVar/UniProt/COSMIC kaynakli)")


def hotspots():
    for kind in ("single", "indel"):
        url = f"https://www.cancerhotspots.org/api/hotspots/{kind}"
        r = S.get(url, timeout=180)
        r.raise_for_status()
        rows = [x for x in r.json() if x.get("hugoSymbol") == "TP53"]
        p = RAW / f"cancerhotspots_{kind}_TP53.json"
        p.write_text(json.dumps(rows, indent=1), encoding="utf-8")
        log(f"hotspots_{kind}", url, p, f"cancerhotspots.org v2, TP53 satirlari ({len(rows)})")


def erepo():
    url = "https://erepo.clinicalgenome.org/evrepo/api/classifications?gene=TP53&format=json&matchLimit=5000"
    get_file("erepo_TP53", url, "erepo_TP53.json", "ClinGen eRepo TP53 VCEP uzman siniflamalari")


def myvariant():
    """dbNSFP (REVEL, AlphaMissense, BayesDel, VARITY, SIFT, PolyPhen ...) + CADD + dbscSNV + ClinVar.

    MyVariant.info hg19 tabanli; tum TP53 kodlayan SNV'leri scroll ile cekilir.
    """
    url = "https://myvariant.info/v1/query"
    fields = "dbnsfp,cadd.phred,cadd.consequence,dbscsnv,clinvar,gnomad_exome.af,gnomad_genome.af,dbsnp.rsid"
    out, params = [], {"q": "dbnsfp.genename:TP53", "fields": fields, "fetch_all": "true"}
    while True:
        r = S.get(url, params=params, timeout=300)
        r.raise_for_status()
        j = r.json()
        hits = j.get("hits", [])
        if not hits:
            break
        out.extend(hits)
        print(f"     myvariant {len(out)} kayit")
        params = {"scroll_id": j["_scroll_id"]}
        time.sleep(0.3)
    p = RAW / "myvariant_TP53_dbnsfp.json"
    p.write_text(json.dumps(out), encoding="utf-8")
    log("myvariant_dbnsfp", url + "?q=dbnsfp.genename:TP53", p,
        f"MyVariant.info dbNSFP+CADD+dbscSNV+ClinVar, {len(out)} varyant (hg19 _id)")


GNOMAD_Q = """
query($ds: DatasetId!) {
  region(chrom: "17", start: 7661000, stop: 7688500, reference_genome: GRCh38) {
    variants(dataset: $ds) {
      variant_id pos ref alt rsids consequence hgvsc hgvsp transcript_id flags
      exome  { ac an ac_hom ac_hemi filters populations { id ac an ac_hom } }
      genome { ac an ac_hom ac_hemi filters populations { id ac an ac_hom } }
    }
  }
}"""


def gnomad():
    url = "https://gnomad.broadinstitute.org/api"
    r = S.post(url, json={"query": GNOMAD_Q, "variables": {"ds": "gnomad_r4"}}, timeout=300)
    r.raise_for_status()
    j = r.json()
    if "errors" in j:
        raise RuntimeError(j["errors"][:2])
    v = j["data"]["region"]["variants"]   # TP53 tum gen bolgesi (UTR + intronlar dahil), hg38
    p = RAW / "gnomad_v4_TP53.json"
    p.write_text(json.dumps(v), encoding="utf-8")
    log("gnomad_v4", url, p, f"gnomAD v4.1 (gnomad_r4) TP53 varyantlari, alt populasyonlar dahil ({len(v)})")


def clinvar():
    """ClinVar: TP53 genindeki tum kayitlar (E-utilities esearch + esummary)."""
    eu = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    r = S.get(eu + "esearch.fcgi", params={"db": "clinvar", "term": "TP53[gene] AND single_gene[prop]",
                                            "retmax": 20000, "retmode": "json"}, timeout=120)
    r.raise_for_status()
    ids = r.json()["esearchresult"]["idlist"]
    print(f"     clinvar {len(ids)} kayit bulundu")
    docs = []
    for i in range(0, len(ids), 300):
        chunk = ids[i:i + 300]
        for attempt in range(4):
            try:
                rr = S.post(eu + "esummary.fcgi", data={"db": "clinvar", "id": ",".join(chunk),
                                                         "retmode": "json"}, timeout=180)
                rr.raise_for_status()
                res = rr.json()["result"]
                docs.extend(res[u] for u in res["uids"])
                break
            except Exception as e:  # NCBI gecici hatalari
                print(f"     tekrar {attempt+1}: {e}")
                time.sleep(3)
        time.sleep(0.4)
    print(f"     clinvar esummary {len(docs)}/{len(ids)}")
    p = RAW / "clinvar_TP53_esummary.json"
    p.write_text(json.dumps(docs), encoding="utf-8")
    log("clinvar_TP53", eu + "esummary.fcgi (db=clinvar, TP53[gene])", p, f"ClinVar ozetleri ({len(docs)})")


def civic():
    """CIViC gecelik TSV: klinik kanit ozetleri (somatik/terapotik), TP53 satirlari."""
    for kind in ("ClinicalEvidenceSummaries", "VariantSummaries"):
        url = f"https://civicdb.org/downloads/nightly/nightly-{kind}.tsv"
        r = S.get(url, timeout=300)
        r.raise_for_status()
        lines = r.text.splitlines()
        keep = [lines[0]] + [ln for ln in lines[1:]
                             if any(c == "TP53" or c.startswith("TP53 ") for c in ln.split("\t"))]
        p = RAW / f"civic_{kind}_TP53.tsv"
        p.write_text("\n".join(keep), encoding="utf-8")
        log(f"civic_{kind}", url, p, f"CIViC nightly, TP53 satirlari ({len(keep)-1})")


SOURCES = {"vcep": vcep_tables, "nci": nci_tp53, "uniprot": uniprot, "hotspots": hotspots,
           "erepo": erepo, "myvariant": myvariant, "gnomad": gnomad, "clinvar": clinvar, "civic": civic}

if __name__ == "__main__":
    todo = sys.argv[1:] or list(SOURCES)
    for k in todo:
        print(f"[{k}]", flush=True)
        try:
            SOURCES[k]()
        except Exception as e:
            print(f"  HATA {k}: {e.__class__.__name__}: {str(e)[:300]}")
    print("bitti ->", RAW)
