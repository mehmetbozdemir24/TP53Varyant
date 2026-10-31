"""ACMG/AMP siniflandirma motorlari.

1) VCEP modu  : ClinGen TP53 VCEP spesifikasyonu v2.4 (CSpec GN009), Tavtigian Bayes puan sistemi.
2) Franklin modu: Franklin (Genoox) yardim makalelerinde tarif edilen genel kurallar + Richards 2015 birlestirme.

Her iki motor da her kriter icin "karsilandi / karsilanmadi / uygulanmaz" kararini gerekceyle dondurur.
"""
from __future__ import annotations

import math
import re
from dataclasses import asdict, dataclass, field

import pandas as pd

from .common import EXONS, HOTSPOT_CODONS, grantham, is_nan

# ------------------------------------------------------------------ ortak yapilar
STRENGTH_POINTS = {"Very Strong": 8, "Strong": 4, "Moderate": 2, "Supporting": 1, "Stand Alone": 8,
                   "No Influence": 0}

CATEGORY = {  # Franklin'in ACMG sekmesindeki gruplar
    "BA1": "Population Data", "BS1": "Population Data", "BS2": "Population Data", "PM2": "Population Data",
    "PVS1": "Effect on Protein", "PS1": "Effect on Protein", "PM4": "Effect on Protein", "PM5": "Effect on Protein",
    "BP1": "Effect on Protein", "BP3": "Effect on Protein", "BP7": "Effect on Protein",
    "PM1": "Functional Data", "PP2": "Functional Data",
    "PP3": "In-silico Predictions", "BP4": "In-silico Predictions",
    "PS3": "Functional Studies", "BS3": "Functional Studies",
    "PP1": "Segregation Data", "BS4": "Segregation Data",
    "PS2": "De Novo Data", "PM6": "De Novo Data",
    "PS4": "Case Control Studies", "PM3": "Allelic Data", "BP2": "Allelic Data",
    "PP4": "Phenotype", "BP5": "Alternate Locus", "PP5": "Reputable Source Data", "BP6": "Reputable Source Data",
}
CATEGORY_ORDER = ["Functional Studies", "Segregation Data", "Case Control Studies", "Population Data",
                  "Effect on Protein", "Functional Data", "In-silico Predictions", "Reputable Source Data",
                  "De Novo Data", "Allelic Data", "Phenotype", "Alternate Locus"]
DESCRIPTION = {
    "PVS1": "Null varyant (nonsense, frameshift, kanonik +/-1,2 splice, baslangic kodonu), LOF hastalik mekanizmasi",
    "PS1": "Daha once patojenik kabul edilmis varyantla ayni amino asit degisimi",
    "PS2": "De novo (anne ve baba dogrulanmis)",
    "PS3": "Iyi kurulmus fonksiyonel calismalar zararli etki gosteriyor",
    "PS4": "Etkilenen bireylerde kontrollere gore anlamli artis / LFS kriterlerini karsilayan probandlar",
    "PM1": "Mutasyonel hotspot veya kritik fonksiyonel domain",
    "PM2": "Populasyon veritabanlarinda yok veya cok dusuk frekans",
    "PM3": "Resesif bozukluklarda trans konumda patojenik varyant",
    "PM4": "Tekrar olmayan bolgede in-frame indel / stop kaybi ile protein uzunluk degisimi",
    "PM5": "Bilinen patojenik varyantla ayni kalintida farkli missense degisim",
    "PM6": "Varsayilan de novo (dogrulamasiz)",
    "PP1": "Etkilenen aile bireylerinde kosegregasyon",
    "PP2": "Benign missense orani dusuk ve missense'in yaygin mekanizma oldugu gen",
    "PP3": "Hesaplamali kanitlar zararli etkiyi destekliyor",
    "PP4": "Fenotip/aile oykusu veya dusuk VAF (TP53: klonal hematopoez sinyali)",
    "PP5": "Saygin kaynak patojenik bildiriyor (kanit erisilemiyor)",
    "BA1": "Populasyon frekansi tek basina benign icin yeterince yuksek",
    "BS1": "Allel frekansi hastalik icin beklenenden yuksek",
    "BS2": "Erken tam penetrans beklenen durumda saglikli eriskinde gozlem",
    "BS3": "Iyi kurulmus fonksiyonel calismalar zararli etki gostermiyor",
    "BS4": "Etkilenen aile bireylerinde segregasyon yok",
    "BP1": "Yalnizca truncating varyantlarin hastalik yaptigi gende missense",
    "BP2": "Patojenik varyantla trans (dominant) veya cis konumda",
    "BP3": "Tekrar bolgesinde fonksiyonu bilinmeyen in-frame indel",
    "BP4": "Hesaplamali kanitlar etki olmadigini destekliyor",
    "BP5": "Alternatif molekuler nedeni olan vakada bulunmus",
    "BP6": "Saygin kaynak benign bildiriyor",
    "BP7": "Splice etkisi ongorulmeyen sinonim (veya derin intronik) varyant",
}


@dataclass
class Evidence:
    code: str                  # PS3
    met: bool
    strength: str | None = None  # Strong / Moderate / Supporting / Very Strong / Stand Alone / No Influence
    reason: str = ""
    source: str = ""
    auto: bool = True          # otomatik mi, klinik/kullanici girdisi mi
    applicable: bool = True    # TP53 VCEP'te uygulanabilir mi

    @property
    def label(self):
        if not self.met or self.strength in (None, self.default_strength()):
            return self.code
        return f"{self.code}_{self.strength.replace(' ', '')}"

    def default_strength(self):
        c = self.code
        return ("Very Strong" if c.startswith("PVS") else "Strong" if c.startswith(("PS", "BS")) else
                "Moderate" if c.startswith("PM") else "Supporting" if c.startswith(("PP", "BP")) else
                "Stand Alone" if c == "BA1" else None)

    @property
    def points(self):
        if not self.met or self.strength is None:
            return 0
        if self.code == "BA1":
            return -8
        p = STRENGTH_POINTS.get(self.strength, 0)
        return -p if self.code.startswith("B") else p

    @property
    def category(self):
        return CATEGORY.get(self.code, "Other")

    def to_dict(self):
        d = asdict(self)
        d.update(label=self.label, points=self.points, category=self.category,
                 description=DESCRIPTION.get(self.code, ""))
        return d


def _f(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return None


def _s(x):
    return None if is_nan(x) or x in ("", "None", "nan") else str(x)


VCEP_CLINICAL_CODES = ("PS2", "PS4", "PP1", "PP4", "BS2", "BS4", "BP2")


def points_to_class(points, ev):
    met = [e for e in ev if e.met and e.points != 0]
    if any(e.code == "BA1" for e in met):
        return "Benign"
    if points >= 10:
        return "Pathogenic"
    if points >= 6:
        return "Likely Pathogenic"
    if points <= -7:
        return "Benign"
    if points <= -2:
        return "Likely Benign"
    path = [e for e in met if e.points > 0]
    ben = [e for e in met if e.points < 0]
    if points == -1 and len(ben) >= 2 and [e.label for e in path] == ["PM2_Supporting"]:
        return "Likely Benign"
    return "Uncertain Significance"


# ================================================================== klinik puan tablolari (Genome Med 2025)
LFS_STRONG = "Guclu iliskili LFS kanseri"
LFS_MODERATE = "Orta iliskili LFS kanseri"
PS2_POINTS = {(LFS_STRONG, True): 4, (LFS_STRONG, False): 2, (LFS_MODERATE, True): 2, (LFS_MODERATE, False): 1}


def strength_from_points(pts, scale):
    """scale: [(alt_sinir, guc), ...] buyukten kucuge."""
    for lo, s in scale:
        if pts >= lo:
            return s
    return None


PS2_SCALE = [(8, "Very Strong"), (4, "Strong"), (2, "Moderate"), (1, "Supporting")]
PS4_SCALE = [(8, "Very Strong"), (4, "Strong"), (2, "Moderate"), (1, "Supporting")]


@dataclass
class Clinical:
    """Kullanicinin (veya LLM onerisinin) girdigi klinik kanitlar."""
    denovo: list = field(default_factory=list)   # [(LFS_STRONG|LFS_MODERATE, anne_baba_dogrulanmis:bool), ...]
    chompret: int = 0                             # Revised Chompret karsilayan akrabasiz proband
    classic_lfs: int = 0                          # Classic LFS karsilayan proband
    her2_lt40: int = 0                            # HER2+ meme ca < 40 yas proband
    meioses: int = 0                              # LFS kanseri ile kosegregasyon meiozu
    families: int = 1
    vaf_5_25: int = 0                             # PP4: VAF %5-25 bagimsiz gozlem
    vaf_25_35: int = 0                            # PP4: VAF %25-35 gozlem
    bs2_females60: int = 0                        # 60 yas kansersiz akrabasiz kadin tasiyici (tek kaynak)
    nonsegregation: bool = False                  # BS4
    trans_or_cis: bool = False                    # BP2
    use_expert: bool = True                       # eRepo uzman klinik kodlarini kullan

    def empty(self):
        return not (self.denovo or self.chompret or self.classic_lfs or self.her2_lt40 or self.meioses
                    or self.vaf_5_25 or self.vaf_25_35 or self.bs2_females60 or self.nonsegregation
                    or self.trans_or_cis)


def parse_code(label):
    """'PS4_Moderate' -> ('PS4','Moderate'); 'PS3' -> ('PS3', None); 'PS2_VeryStrong' -> ('PS2','Very Strong')."""
    base, _, mod = label.partition("_")
    m = mod.replace("_", "").lower()
    s = {"verystrong": "Very Strong", "strong": "Strong", "moderate": "Moderate", "supporting": "Supporting"}.get(m)
    return base, s


# ================================================================== siniflandirici
class TP53Classifier:
    def __init__(self, variants: pd.DataFrame):
        self.df = variants.set_index("c", drop=False) if "c" in variants.columns else variants
        er = self.df[self.df["erepo_class"].notna()]
        self.erepo = er
        cvp = self.df[self.df["clinvar_germline"].fillna("").str.contains("athogenic") &
                      (self.df["clinvar_stars"].fillna(0).astype(float) >= 1)]
        self.clinvar_pathogenic = cvp

    # -------------------------------------------------------------- arama
    def find(self, q: str):
        """'R175H', 'p.Arg175His', 'c.524G>A', 'TP53:c.524G>A', 'chr17-7578406-C-T' (hg19), '17-7675088-C-T' (hg38),
        'rs28934578' kabul eder. Eslesen varyant satirlarini dondurur."""
        from .common import norm_c, p3_to_p1
        q = q.strip().replace(" ", "")
        # 'TP53:', 'NM_000546.6:', 'NM_000546.5(TP53):' gibi onekleri at
        q = re.sub(r"^(?:NM_000546(?:\.\d+)?)?(?:\(?TP53\)?)?:?", "", q, flags=re.I)
        q = re.sub(r"^c\.([-*]?\d+)_\1(del|dup)[ACGT]?$", r"c.\1\2", q)    # c.215_215del -> c.215del
        if re.match(r"^p\.\(.*\)$", q, re.I):                            # p.(R175H) -> p.R175H
            q = "p." + q[3:-1]
        d = self.df
        if re.match(r"^rs\d+$", q, re.I):
            # ayni rsID birden fazla alleli kapsayabilir (rs1042522: P72R yaygin, P72H nadir) -> en sik ve kayitli olan once
            hit = d[d["rsid"].fillna("").str.lower() == q.lower()]
            score = (pd.to_numeric(hit["gnomad_af"], errors="coerce").fillna(0) * 1000
                     + hit["erepo_class"].notna().astype(int) * 10 + hit["clinvar_vcv"].notna().astype(int) * 5)
            return hit.loc[score.sort_values(ascending=False).index]
        m = re.match(r"^(?:chr)?17[-:_](\d+)[-:_]([ACGT]+)[-:_>]([ACGT]+)$", q, re.I)
        if m:
            pos, ref, alt = int(m.group(1)), m.group(2).upper(), m.group(3).upper()
            for col in ("hg19_pos", "hg38_pos"):
                hit = d[(pd.to_numeric(d[col], errors="coerce") == pos) & (d["g_ref"] == ref) & (d["g_alt"] == alt)]
                if len(hit):
                    return hit
            return d.iloc[0:0]
        if q.lower().startswith("c."):
            c = norm_c("c." + q[2:])
            return d[d["c"] == c]
        p1 = p3_to_p1(q if q.startswith("p.") else "p." + q) if re.search(r"[A-Za-z]{3}\d", q) else q.upper().replace("P.", "")
        p1 = p1.replace("X", "*").replace("TER", "*")
        if re.match(r"^(P\.)?\(?M(ET)?1\??\)?$", q.upper()):     # M1? / p.Met1? -> start kaybi
            return d[d["consequence"] == "start_lost"]
        ms = re.match(r"^([A-Z*])(\d+)\1$", p1)     # R213R -> R213= (sinonim)
        if ms:
            p1 = f"{ms.group(1)}{ms.group(2)}="
        hit = d[d["p1"] == p1]
        if len(hit) > 1:   # tek nukleotid degisimi ve en cok kanit tasiyan once
            score = (hit["c"].str.match(r"^c\.\d+[ACGT]>[ACGT]$").astype(int) * 100
                     + hit["erepo_class"].notna().astype(int) * 10 + hit["clinvar_vcv"].notna().astype(int) * 5
                     + pd.to_numeric(hit["nci_tumor_n"], errors="coerce").fillna(0).clip(upper=4))
            hit = hit.loc[score.sort_values(ascending=False).index]
        return hit

    # -------------------------------------------------------------- yardimcilar
    def _row(self, c):
        return self.df.loc[c] if isinstance(c, str) else c

    @staticmethod
    def _erepo_class_without_pm5(row):
        met = [x for x in (_s(row.get("erepo_met")) or "").split(";") if x]
        from .build import code_points
        pts = sum(code_points(x) for x in met if not x.startswith("PM5"))
        cls = "Pathogenic" if pts >= 10 else "Likely Pathogenic" if pts >= 6 else "Other"
        clinical = any(x.split("_")[0] in ("PS2", "PS4", "PP1") for x in met)
        return cls, clinical

    # ============================================================== VCEP v2.4
    def vcep(self, c, clinical: Clinical | None = None):
        r = self._row(c)
        clinical = clinical or Clinical()
        csq = _s(r.get("consequence")) or "other"
        p1 = _s(r.get("p1"))
        codon = _f(r.get("codon"))
        spliceai = _f(r.get("spliceai"))
        ev: dict[str, Evidence] = {}
        notes = []
        src_cspec = "ClinGen TP53 VCEP CSpec v2.4"

        # ---------------- PVS1
        pvs1 = Evidence("PVS1", False, reason="Null varyant degil.", source=src_cspec)
        ptc = r.get("ptc_codon")
        ptc_i = int(float(ptc)) if _f(ptc) is not None else None
        if csq in ("nonsense", "frameshift"):
            if ptc_i is not None:
                if ptc_i < 351:
                    pvs1 = Evidence("PVS1", True, "Very Strong", f"Erken stop kodonu p.{ptc_i} (< p.Lys351): NMD beklenir.", src_cspec)
                elif ptc_i <= 355:
                    pvs1 = Evidence("PVS1", True, "Strong", f"Erken stop p.{ptc_i} (p.Lys351-Ala355): NMD beklenmez.", src_cspec)
                else:
                    pvs1 = Evidence("PVS1", True, "Moderate", f"Erken stop p.{ptc_i} (p.Gly356-Asp393): NMD beklenmez.", src_cspec)
            elif ptc == "beyond":
                pvs1 = Evidence("PVS1", True, "Moderate", "Cerceve kaymasi dogal stop kodonunun otesinde sonlaniyor.", src_cspec)
            elif codon is not None and codon < 330:
                pvs1 = Evidence("PVS1", True, "Very Strong",
                                f"Cerceve kaymasi p.{int(codon)}'de basliyor; yeni stop konumu hesaplanamadi, "
                                "p.Lys351 oncesi varsayildi.", src_cspec)
            else:
                pvs1.reason = "Cerceve kaymasinin stop konumu belirlenemedi; elle degerlendirin."
        elif csq == "stop_lost" and re.search(r"(del|dup|ins)", r["c"]):
            pvs1 = Evidence("PVS1", True, "Moderate", "Cerceve kaymasi dogal stop kodonunu asiyor (uzama): VCEP PVS1_Moderate.",
                            src_cspec)
        elif csq == "start_lost":
            pvs1 = Evidence("PVS1", True, "Very Strong", "Baslangic kodonu varyanti (VCEP: PVS1 uygulanabilir).", src_cspec)
        elif csq == "splice_canonical":
            s1 = _s(r.get("s1_pvs1"))
            if s1:
                if "N/A" in s1:
                    pvs1 = Evidence("PVS1", False, None, f"Tablo S1: {s1} ({_s(r.get('s1_rationale')) or ''}).",
                                    "VCEP Tablo S1")
                else:
                    st = "Moderate" if "Moderate" in s1 else "Strong" if "Strong" in s1 else "Very Strong"
                    pvs1 = Evidence("PVS1", True, st, f"Tablo S1: {s1}. {_s(r.get('s1_rationale')) or ''}",
                                    "VCEP Tablo S1 (SpliceAI + RNA verisi)")
            else:
                pvs1 = Evidence("PVS1", True, "Very Strong",
                                "Kanonik +/-1,2 splice varyanti; Tablo S1'de yok, varsayilan tam guc.", src_cspec)
        ev["PVS1"] = pvs1
        pvs1_full = pvs1.met and pvs1.strength == "Very Strong"

        # ---------------- PP3 / BP4 / BP7  (in-silico)
        pp3 = Evidence("PP3", False, reason="Hesaplamali kanit esigi saglanmadi.", source="VCEP Tablo S2 / SpliceAI")
        bp4 = Evidence("BP4", False, reason="Hesaplamali kanit esigi saglanmadi.", source="VCEP Tablo S2 / SpliceAI")
        bp7 = Evidence("BP7", False, reason="Sinonim/derin intronik degil veya splice etkisi dislanamadi.", source=src_cspec)
        pp3_splice = False
        if pvs1.met:
            pp3.reason = "PVS1 uygulandigi icin PP3 kullanilmaz."
        elif csq == "missense":
            code = _s(r.get("s2_code"))
            agv = _s(r.get("agvgd")); bd = _f(r.get("bayesdel_vcep"))
            if code is None and bd is not None and agv:
                cls = int(re.sub(r"\D", "", agv) or 0)
                code = ("PP3_moderate" if cls >= 65 and bd >= 0.16 else "PP3" if 25 <= cls <= 55 and bd >= 0.16 else
                        "BP4_moderate" if bd <= -0.008 and cls < 65 else "BP4" if -0.008 < bd < 0.16 and cls < 65 else
                        "No evidence")
            txt = f"Align-GVGD {agv}, BayesDel {bd:.3f}" if bd is not None else f"Align-GVGD {agv}"
            if code and code.lower().startswith("pp3"):
                pp3 = Evidence("PP3", True, "Moderate" if "moderate" in code.lower() else "Supporting",
                               f"{txt} -> {code}.", "VCEP Tablo S2 (Align-GVGD + BayesDel)")
            elif code and code.lower().startswith("bp4"):
                if spliceai is not None and spliceai >= 0.2:
                    bp4.reason = f"{txt} BP4'u destekliyor ama SpliceAI {spliceai:.2f} >= 0.2 (splice etkisi)."
                else:
                    bp4 = Evidence("BP4", True, "Moderate" if "moderate" in code.lower() else "Supporting",
                                   f"{txt}, SpliceAI {spliceai if spliceai is not None else 'yok'} -> {code}.",
                                   "VCEP Tablo S2 (Align-GVGD + BayesDel)")
            else:
                pp3.reason = bp4.reason = f"{txt}: Tablo S2'de kanit yok."
            if spliceai is not None and spliceai >= 0.2:
                pp3_splice = True
                if not pp3.met:
                    pp3 = Evidence("PP3", True, "Supporting", f"SpliceAI {spliceai:.2f} >= 0.2: splice etkisi ongoruluyor.",
                                   "SpliceAI")
                else:
                    pp3.reason += f" Ayrica SpliceAI {spliceai:.2f} >= 0.2."
        elif csq in ("synonymous", "intronic", "splice_region", "5_utr", "3_utr"):
            if spliceai is None:
                pp3.reason = bp4.reason = "SpliceAI skoru bulunamadi."
            elif spliceai >= 0.2:
                pp3 = Evidence("PP3", True, "Supporting", f"SpliceAI {spliceai:.2f} >= 0.2.", "SpliceAI")
                pp3_splice = True
            elif spliceai <= 0.1 and csq in ("synonymous", "intronic", "splice_region"):
                bp4 = Evidence("BP4", True, "Supporting", f"SpliceAI {spliceai:.2f} <= 0.1.", "SpliceAI")
                cpos, off = _f(r.get("c_pos")), _f(r.get("c_offset"))
                core = False
                if csq == "synonymous" and cpos is not None:
                    for ex, a, b in EXONS:
                        if a <= cpos <= b:
                            core = cpos >= b - 2 or cpos == a
                    ok = not core
                    why = "sinonim, ekson cekirdek motifi disinda"
                else:
                    ok = off is not None and (off >= 7 or off <= -21)
                    why = f"intronik ofset {int(off) if off is not None else '?'} (+7 / -21 otesinde)"
                if ok:
                    bp7 = Evidence("BP7", True, "Supporting", f"{why}; SpliceAI {spliceai:.2f} <= 0.1.", src_cspec)
                else:
                    bp7.reason = "Ekson cekirdek motifinde veya +7/-21 icinde; BP7 uygulanmaz."
            else:
                pp3.reason = bp4.reason = f"SpliceAI {spliceai:.2f} ara bolgede (0.1-0.2)."
        elif csq == "inframe_deletion":
            pp3.reason = bp4.reason = "Tek amino asit delesyonu icin VCEP BayesDel tablosu gerekir (bu surumde yok)."
        else:
            pp3.reason = bp4.reason = f"{csq} icin in-silico kural tanimli degil."
        ev["PP3"], ev["BP4"], ev["BP7"] = pp3, bp4, bp7

        # ---------------- PS3 / BS3
        ps3 = Evidence("PS3", False, reason="Uygun fonksiyonel veri yok.", source="VCEP Tablo S3")
        bs3 = Evidence("BS3", False, reason="Uygun fonksiyonel veri yok.", source="VCEP Tablo S3")
        code = _s(r.get("s3_code"))
        assays = ", ".join(f"{k} {_s(r.get(k.lower()))}" for k in ("Kato", "Funk", "Giacomelli", "Kotler", "Kawaguchi")
                           if _s(r.get(k.lower())) and _s(r.get(k.lower())) != "NA")
        if csq == "missense" and code and code != "No evidence":
            if pp3_splice:
                ps3.reason = bs3.reason = (f"{code} onerilirdi, fakat PP3 SpliceAI'ye dayandigi icin protein duzeyi "
                                           "testleri kullanilmaz (VCEP uyarisi).")
            elif code.startswith("PS3"):
                st = "Moderate" if "Moderate" in code else "Supporting" if "Supporting" in code else "Strong"
                if pvs1_full:
                    ps3.reason = "PVS1 tam guc uygulandi; PS3 kullanilmaz."
                else:
                    if pvs1.met and pvs1.strength == "Strong" and st == "Strong":
                        st = "Moderate"
                    ps3 = Evidence("PS3", True, st, f"{assays} -> {code}.", "VCEP Tablo S3 (Kato/Funk/Giacomelli/Kotler/Kawaguchi)")
            elif code.startswith("BS3"):
                st = "Supporting" if "Supporting" in code else "Strong"
                bs3 = Evidence("BS3", True, st, f"{assays} -> {code}.", "VCEP Tablo S3")
        elif csq == "missense":
            ps3.reason = bs3.reason = f"Tablo S3: {assays or 'veri yok'} -> kanit yok."
        elif pvs1_full:
            ps3.reason = "PVS1 tam guc uygulandi; PS3 kullanilmaz."
        ev["PS3"], ev["BS3"] = ps3, bs3

        # ---------------- PM1
        pm1 = Evidence("PM1", False, reason="Missense degil veya hotspot degil.", source="VCEP CSpec + cancerhotspots.org v2")
        if csq == "missense":
            n = _f(r.get("hotspot_aa_count")) or 0
            if codon is not None and int(codon) in HOTSPOT_CODONS:
                pm1 = Evidence("PM1", True, "Moderate", f"Kodon {int(codon)} VCEP sabit hotspot listesinde "
                               "(175, 245, 248, 249, 273, 282).", "ClinGen TP53 VCEP CSpec v2.4")
            elif n >= 10:
                pm1 = Evidence("PM1", True, "Moderate", f"cancerhotspots.org'da ayni degisim {int(n)} somatik gozlem (>= 10).",
                               "cancerhotspots.org v2")
            elif n >= 2:
                pm1 = Evidence("PM1", True, "Supporting", f"cancerhotspots.org'da ayni degisim {int(n)} somatik gozlem (2-9).",
                               "cancerhotspots.org v2")
            else:
                k = _f(r.get("tcga_icgc_genie_count"))
                pm1.reason = (f"Hotspot kodonu degil; cancerhotspots.org (anlamli hotspot listesi) gozlem sayisi {int(n)}."
                              + (f" Bilgi: NCI TCGA/ICGC/GENIE sayisi {int(k)}; kalinti listede yoksa "
                                 "cancerhotspots.org sitesinden elle kontrol edin." if k else ""))
        ev["PM1"] = pm1

        # ---------------- PM2 / BS1 / BA1
        ac, an = _f(r.get("gnomad_ac")), _f(r.get("gnomad_an"))
        af = _f(r.get("gnomad_af")) or 0.0
        faf = _f(r.get("gnomad_faf95_popmax")) or 0.0
        fpop = _s(r.get("gnomad_faf95_pop"))
        mpop = _f(r.get("gnomad_max_pop_af_multi")) or 0.0
        src_g = "gnomAD v4.1 (exome+genome)"
        pm2 = Evidence("PM2", False, reason="", source=src_g)
        bs1 = Evidence("BS1", False, reason=f"FAF95 popmax {faf:.2e} < 0.0003.", source=src_g)
        ba1 = Evidence("BA1", False, reason=f"FAF95 popmax {faf:.2e} < 0.001.", source=src_g)
        if an is None or not ac:
            pm2 = Evidence("PM2", True, "Supporting", "gnomAD v4.1'de gozlenmemis.", src_g)
        elif af < 0.00003 and mpop < 0.00004:
            pm2 = Evidence("PM2", True, "Supporting",
                           f"gnomAD AF {af:.2e} (AC {int(ac)}/{int(an)}) < 0.00003; coklu alel goren gruplarda max AF {mpop:.2e} < 0.00004.",
                           src_g)
        else:
            pm2.reason = (f"gnomAD AF {af:.2e} (AC {int(ac)}) esik ustunde. VCEP: dusuk VAF (< 0.35) alleller CHIP olabilir; "
                          "gnomAD allel dengesi histogramindan duzeltme yapilabilir.")
        if faf >= 0.001:
            ba1 = Evidence("BA1", True, "Stand Alone", f"FAF95 {faf:.2e} >= 0.001 ({fpop}).", src_g)
        elif faf >= 0.0003:
            bs1 = Evidence("BS1", True, "Strong", f"FAF95 {faf:.2e} 0.0003-0.001 araliginda ({fpop}).", src_g)
        ev["PM2"], ev["BS1"], ev["BA1"] = pm2, bs1, ba1

        # ---------------- PS1 / PM5 (eRepo uzman siniflamalarina gore)
        ps1 = Evidence("PS1", False, reason="Ayni amino asit degisimine sahip uzman siniflamali P/LP varyant yok.",
                       source="ClinGen eRepo (TP53 VCEP)")
        pm5 = Evidence("PM5", False, reason="Ayni kodonda uzman siniflamali P/LP farkli missense yok.",
                       source="ClinGen eRepo (TP53 VCEP)")
        if csq == "missense" and p1:
            ok_splice = spliceai is None or spliceai < 0.2
            same = self.erepo[(self.erepo["p1"] == p1) & (self.erepo.index != r["c"])]
            if len(same) and ok_splice:
                cls = set(same["erepo_class"])
                if "Pathogenic" in cls:
                    ps1 = Evidence("PS1", True, "Strong", f"Ayni degisim ({', '.join(same.index)}) uzman paneli: Pathogenic.",
                                   "ClinGen eRepo")
                elif "Likely Pathogenic" in cls:
                    ps1 = Evidence("PS1", True, "Moderate", f"Ayni degisim ({', '.join(same.index)}) uzman paneli: Likely Pathogenic.",
                                   "ClinGen eRepo")
            g_self = _f(r.get("grantham"))
            other = self.erepo[(self.erepo["codon"].astype(float) == codon) & (self.erepo["p1"] != p1) &
                               (self.erepo["consequence"] == "missense")] if codon is not None else self.erepo.iloc[0:0]
            n_p, n_lp, used, skipped = 0, 0, [], []
            for oc, orow in other.iterrows():
                cls, clin = self._erepo_class_without_pm5(orow)
                g_o = _f(orow.get("grantham"))
                if g_self is not None and g_o is not None and g_self < g_o:
                    skipped.append(f"{orow['p1']} (Grantham {int(g_o)} > {int(g_self)})")
                    continue
                if cls == "Pathogenic":
                    n_p += 1; used.append(f"{orow['p1']} P")
                elif cls == "Likely Pathogenic" and clin:
                    n_lp += 1; used.append(f"{orow['p1']} LP")
            if not ok_splice:
                pm5.reason = f"SpliceAI {spliceai:.2f} >= 0.2; PM5 icin splice etkisi dislanmali."
            elif n_p >= 2:
                pm5 = Evidence("PM5", True, "Strong", f"Ayni kodonda PM5'siz Pathogenic: {', '.join(used)}.", "ClinGen eRepo")
            elif n_p == 1:
                pm5 = Evidence("PM5", True, "Moderate", f"Ayni kodonda PM5'siz Pathogenic: {', '.join(used)}.", "ClinGen eRepo")
            elif n_lp >= 1:
                pm5 = Evidence("PM5", True, "Supporting", f"Ayni kodonda klinik kanitli Likely Pathogenic: {', '.join(used)}.",
                               "ClinGen eRepo")
            if skipped and not pm5.met:
                pm5.reason += f" Grantham kosulu saglanmayanlar: {', '.join(skipped)}."
            if ps1.met and pm5.met and ps1.points + pm5.points > 4:
                pm5.strength = {4: "Supporting", 2: "Moderate"}.get(4 - ps1.points) if ps1.points < 4 else None
                if pm5.strength is None:
                    pm5.met = False
                    pm5.reason = "PS1 + PM5 toplami 4 puani gecemez; PM5 dusuruldu."
        ev["PS1"], ev["PM5"] = ps1, pm5

        # ---------------- klinik kriterler
        self._clinical(r, clinical, ev, notes)

        # ---------------- TP53'te uygulanmayanlar
        for code in ("PM3", "PM4", "PM6", "PP2", "PP5", "BP1", "BP3", "BP5", "BP6"):
            why = {"PM6": "VCEP: PM6 kullanilmaz, tum de novo kanit PS2'de.",
                   "PP5": "ClinGen SVI: PP5 kaldirildi.", "BP6": "ClinGen SVI: BP6 kaldirildi."}.get(
                code, "TP53 VCEP v2.4: bu kriter TP53 icin uygulanmaz.")
            ev[code] = Evidence(code, False, reason=why, source=src_cspec, applicable=False)

        evs = list(ev.values())
        pts = sum(e.points for e in evs)
        auto_pts = sum(e.points for e in evs if e.auto)
        res = dict(points=pts, classification=points_to_class(pts, evs),
                   auto_points=auto_pts, auto_classification=points_to_class(auto_pts, [e for e in evs if e.auto]),
                   met=[e.label for e in evs if e.met and e.points != 0],
                   evidence=[e.to_dict() for e in evs], notes=notes, engine="ClinGen TP53 VCEP v2.4")
        return res

    def _clinical(self, r, cl: Clinical, ev, notes):
        src = "Kullanici girdisi"
        expert = {}
        if cl.use_expert and _s(r.get("erepo_met")):
            for lab in _s(r.get("erepo_met")).split(";"):
                base, st = parse_code(lab)
                if base in VCEP_CLINICAL_CODES:
                    expert[base] = st or Evidence(base, True).default_strength()
        ba_bs1 = ev["BA1"].met or ev["BS1"].met

        def put(code, strength, reason, source):
            ev[code] = Evidence(code, strength is not None, strength, reason, source, auto=False)

        # PS2
        pts = sum(PS2_POINTS.get((g, bool(conf)), 0) for g, conf in cl.denovo)
        if pts:
            put("PS2", strength_from_points(pts, PS2_SCALE), f"{len(cl.denovo)} de novo proband, toplam {pts} puan.", src)
        elif "PS2" in expert:
            put("PS2", expert["PS2"], "ClinGen uzman paneli (eRepo) klinik kaniti.", "ClinGen eRepo")
        else:
            ev["PS2"] = Evidence("PS2", False, reason="De novo verisi girilmedi.", source="Klinik veri", auto=False)
        # PS4
        pts = 0.5 * cl.chompret + 1.0 * cl.classic_lfs + 0.5 * cl.her2_lt40
        if pts and not ev["PM2"].met:
            ev["PS4"] = Evidence("PS4", False, reason=f"{pts} puan, fakat PM2_Supporting karsilanmadigi icin PS4 uygulanmaz.",
                                 source=src, auto=False)
        elif pts:
            put("PS4", strength_from_points(pts, PS4_SCALE),
                f"Chompret {cl.chompret} x0.5 + Classic LFS {cl.classic_lfs} x1 + HER2+ <40y {cl.her2_lt40} x0.5 = {pts} puan.", src)
        elif "PS4" in expert:
            put("PS4", expert["PS4"], "ClinGen uzman paneli (eRepo) klinik kaniti.", "ClinGen eRepo")
        else:
            ev["PS4"] = Evidence("PS4", False, reason="Proband/fenotip verisi girilmedi.", source="Klinik veri", auto=False)
        # PP1
        if cl.meioses and not ba_bs1:
            st = ("Strong" if cl.meioses >= 7 and cl.families > 1 else "Moderate" if cl.meioses >= 5 else
                  "Supporting" if cl.meioses >= 3 else None)
            put("PP1", st, f"{cl.meioses} meioz, {cl.families} aile.", src)
        elif "PP1" in expert:
            put("PP1", expert["PP1"], "ClinGen uzman paneli (eRepo) klinik kaniti.", "ClinGen eRepo")
        else:
            ev["PP1"] = Evidence("PP1", False, reason="Kosegregasyon verisi girilmedi." if not ba_bs1 else "BA1/BS1 varken PP1 uygulanmaz.",
                                 source="Klinik veri", auto=False)
        # PP4 (dusuk VAF)
        if (cl.vaf_5_25 or cl.vaf_25_35) and not ba_bs1:
            st = "Moderate" if cl.vaf_5_25 >= 2 else "Supporting"
            put("PP4", st, f"VAF %5-25: {cl.vaf_5_25}, VAF %25-35: {cl.vaf_25_35} gozlem (multigen panel).", src)
        elif "PP4" in expert:
            put("PP4", expert["PP4"], "ClinGen uzman paneli (eRepo) klinik kaniti.", "ClinGen eRepo")
        else:
            ev["PP4"] = Evidence("PP4", False, reason="Dusuk VAF verisi girilmedi.", source="Klinik veri", auto=False)
        # BS2
        n = cl.bs2_females60
        if n >= 2:
            put("BS2", "Strong" if n >= 8 else "Moderate" if n >= 4 else "Supporting",
                f"{n} akrabasiz, 60 yas ustu kansersiz kadin tasiyici (tek kaynak).", src)
        elif "BS2" in expert:
            put("BS2", expert["BS2"], "ClinGen uzman paneli (eRepo) klinik kaniti.", "ClinGen eRepo")
        else:
            ev["BS2"] = Evidence("BS2", False, reason="Saglikli tasiyici verisi girilmedi.", source="Klinik veri", auto=False)
        # BS4 / BP2
        if cl.nonsegregation:
            put("BS4", "Strong", "Etkilenen aile bireylerinde segregasyon yok.", src)
        elif "BS4" in expert:
            put("BS4", expert["BS4"], "ClinGen uzman paneli (eRepo) klinik kaniti.", "ClinGen eRepo")
        else:
            ev["BS4"] = Evidence("BS4", False, reason="Segregasyon yoklugu bildirilmedi.", source="Klinik veri", auto=False)
        if cl.trans_or_cis:
            put("BP2", "Supporting", "Baska P/LP TP53 varyantiyla trans (veya fazsiz >= 3 gozlem). Guc CSpec'te net degil.", src)
        elif "BP2" in expert:
            put("BP2", expert["BP2"], "ClinGen uzman paneli (eRepo) klinik kaniti.", "ClinGen eRepo")
        else:
            ev["BP2"] = Evidence("BP2", False, reason="Allelik veri girilmedi.", source="Klinik veri", auto=False)
        if expert:
            notes.append("Klinik kriterler ClinGen eRepo uzman kaydindan alindi: " +
                         ", ".join(f"{k}_{v.replace(' ', '')}" for k, v in expert.items()))

    # ============================================================== Franklin modu (genel ACMG)
    def franklin(self, c):
        """Franklin yardim makalelerindeki genel kurallar (help.genoox.com) + Richards 2015 birlestirme.
        Franklin'in ozel ML bilesenleri (PM1 hotspot modeli, ClinVar metin madenciligi) burada aciklanan
        yaklasik yontemlerle taklit edilir; birebir ayni sonuc garanti edilmez."""
        r = self._row(c)
        csq = _s(r.get("consequence")) or "other"
        p1 = _s(r.get("p1"))
        codon = _f(r.get("codon"))
        ev = {}
        H = "Franklin yardim merkezi"
        # PVS1
        v = self.vcep(c, Clinical(use_expert=False))
        vp = next(e for e in v["evidence"] if e["code"] == "PVS1")
        ev["PVS1"] = Evidence("PVS1", vp["met"], vp["strength"], vp["reason"] + " (ClinGen SVI PVS1 karar agaci)", H)
        # PM2 / BS1 / BA1 : Franklin TP53 esikleri (sayfada gosterilen): PM2 %0.05, BS1 %0.1, BA1 %5
        af = _f(r.get("gnomad_af")) or 0.0
        maxpop = max([d["af"] for d in __import__("json").loads(_s(r.get("gnomad_pops")) or "{}").values()] or [0.0])
        ev["PM2"] = Evidence("PM2", maxpop < 0.0005, "Moderate" if maxpop < 0.0005 else None,
                             f"gnomAD max alt populasyon AF {maxpop:.2e} {'<' if maxpop < 0.0005 else '>='} %0.05 (TP53 gen esigi).", H)
        ev["BA1"] = Evidence("BA1", maxpop >= 0.05, "Stand Alone" if maxpop >= 0.05 else None, f"Max AF {maxpop:.2e} vs %5.", H)
        ev["BS1"] = Evidence("BS1", 0.001 <= maxpop < 0.05, "Strong" if 0.001 <= maxpop < 0.05 else None,
                             f"Max AF {maxpop:.2e} vs %0.1 (TP53 gen esigi).", H)
        # PS1 / PM5 / PP5 / BP6 : ClinVar (>= 1 yildiz)
        cv = self.df[self.df["clinvar_stars"].fillna(0).astype(float) >= 1]
        sig = cv["clinvar_germline"].fillna("")
        ps1 = Evidence("PS1", False, reason="ClinVar'da ayni amino asit degisimli P/LP baska varyant yok.", source="ClinVar")
        pm5 = Evidence("PM5", False, reason="ClinVar'da ayni kodonda P/LP farkli missense yok.", source="ClinVar")
        if csq == "missense" and p1:
            same = cv[(cv["p1"] == p1) & (cv["c"] != r["c"]) & sig.str.contains("athogenic") & ~sig.str.contains("onflict")]
            if len(same):
                st = "Strong" if same["clinvar_germline"].str.startswith("Pathogenic").any() else "Moderate"
                ps1 = Evidence("PS1", True, st, f"ClinVar: {', '.join(same['c'])} ({', '.join(same['clinvar_germline'])}).", "ClinVar")
            oth = cv[(cv["codon"].astype(float) == codon) & (cv["p1"] != p1) & (cv["consequence"] == "missense")
                     & sig.str.contains("athogenic") & ~sig.str.contains("onflict")]
            if len(oth):
                st = "Moderate" if oth["clinvar_germline"].str.startswith("Pathogenic").any() else "Supporting"
                pm5 = Evidence("PM5", True, st, f"ClinVar ayni kodon: {', '.join(oth['p1'].unique()[:6])}.", "ClinVar")
        ev["PS1"], ev["PM5"] = ps1, pm5
        g = _s(r.get("clinvar_germline")) or ""
        ev["PP5"] = Evidence("PP5", "athogenic" in g and "onflict" not in g, "No Influence" if "athogenic" in g else None,
                             f"ClinVar: {g or 'kayit yok'} ({_s(r.get('clinvar_review')) or '-'}). Franklin cift sayimi onlemek icin "
                             "PP5'i 'No influence' olarak gosterir.", "ClinVar")
        # Franklin BP6'yi hala uyguluyor ve Strong sayiyor (V73M, R110H, R290H sayfalarinda goruldu, 31.10.2026)
        bp6 = "enign" in g and "onflict" not in g and "athogenic" not in g
        ev["BP6"] = Evidence("BP6", bp6, "Strong" if bp6 else None,
                             f"ClinVar: {g or 'kayit yok'} ({_s(r.get('clinvar_review')) or '-'}). Franklin BP6'yi Strong uygular.",
                             "ClinVar")
        # PM1 : Franklin hotspot yaklasimi -> iki benign missense arasindaki P missense yogunlugu
        pm1 = Evidence("PM1", False, reason="Missense degil.", source="Franklin-benzeri hotspot penceresi (ClinVar)")
        if csq == "missense" and codon is not None:
            mis = cv[cv["consequence"] == "missense"].copy()
            mis["cod"] = mis["codon"].astype(float)
            ben = mis[mis["clinvar_germline"].fillna("").str.contains("enign")]["cod"]
            pat = mis[mis["clinvar_germline"].fillna("").str.contains("athogenic") &
                      ~mis["clinvar_germline"].fillna("").str.contains("onflict")]["cod"]
            lo = ben[ben < codon].max() if (ben < codon).any() else 1
            hi = ben[ben > codon].min() if (ben > codon).any() else 393
            if (ben == codon).any():
                pm1.reason = "Ayni kodonda benign missense var."
            else:
                inside = pat[(pat > lo) & (pat < hi)]
                width = max(hi - lo - 1, 1)
                dens = inside.nunique() / width
                n = len(inside)
                st = "Strong" if n >= 30 and dens >= 0.5 else "Moderate" if n >= 10 and dens >= 0.3 else "Supporting" if n >= 5 and dens >= 0.2 else None
                pm1 = Evidence("PM1", st is not None, st, f"Benign missense sinirlari kodon {int(lo)}-{int(hi)}; aradaki P/LP missense "
                               f"{n} varyant, {inside.nunique()} kodon (yogunluk {dens:.2f}). Esikler bu araca ozgu.",
                               "Franklin-benzeri hotspot penceresi (ClinVar)")
        ev["PM1"] = pm1
        # PP2 / BP1
        allmis = cv[cv["consequence"] == "missense"]["clinvar_germline"].fillna("")
        npm = (allmis.str.contains("athogenic") & ~allmis.str.contains("onflict")).sum()
        nbm = allmis.str.contains("enign").sum()
        pp2 = csq == "missense" and npm >= 8 and npm >= 3 * nbm
        ev["PP2"] = Evidence("PP2", pp2, "Supporting" if pp2 else None,
                             f"Gende ClinVar P/LP missense {npm}, B/LB missense {nbm} (kural: P >= 3xB ve P >= 8). "
                             "Not: Franklin R175H icin PP2'yi karsilanmamis gosteriyor.", H)
        ev["BP1"] = Evidence("BP1", False, reason="TP53'te missense patojenik varyantlar var.", source=H)
        # PP3 / BP4 : REVEL (+ SpliceAI / dbscSNV)
        rv = _f(r.get("revel")); sa = _f(r.get("spliceai")); ada = _f(r.get("dbscsnv_ada"))
        pp3 = Evidence("PP3", False, reason="", source="REVEL / SpliceAI / dbscSNV")
        bp4 = Evidence("BP4", False, reason="", source="REVEL / SpliceAI / dbscSNV")
        if ev["PVS1"].met:
            pp3.reason = bp4.reason = "Null varyant; in-silico kurallari uygulanmaz."
        elif csq == "missense" and rv is not None:
            st = "Strong" if rv >= 0.932 else "Moderate" if rv >= 0.773 else "Supporting" if rv >= 0.644 else None
            if st:
                if pm1.met and st == "Strong":
                    st = "Moderate"
                elif pm1.met and st == "Moderate":
                    st = "Supporting"
                pp3 = Evidence("PP3", True, st, f"REVEL {rv:.3f}" + (" (PM1 ile birlikte guc dusuruldu)" if pm1.met else ""),
                               "REVEL (ClinGen SVI 2022 esikleri)")
            elif rv < 0.15:
                bp4 = Evidence("BP4", True, "Supporting", f"REVEL {rv:.3f} < 0.15.", "REVEL")
            else:
                pp3.reason = bp4.reason = f"REVEL {rv:.3f} ara bolgede."
        if not pp3.met and not ev["PVS1"].met and ((sa is not None and sa > 0.5) or (ada is not None and ada > 0.7)):
            pp3 = Evidence("PP3", True, "Supporting", f"Splice: SpliceAI {sa}, dbscSNV Ada {ada}.", "SpliceAI / dbscSNV")
        elif not bp4.met and csq != "missense" and sa is not None and sa < 0.1 and (ada is None or ada < 0.15) and not ev["PVS1"].met:
            bp4 = Evidence("BP4", True, "Supporting", f"SpliceAI {sa} < 0.1.", "SpliceAI")
        ev["PP3"], ev["BP4"] = pp3, bp4
        ev["BP7"] = Evidence("BP7", csq == "synonymous" and sa is not None and sa < 0.1,
                             "Supporting" if csq == "synonymous" and sa is not None and sa < 0.1 else None,
                             "Sinonim ve splice etkisi ongorulmuyor." if csq == "synonymous" else "Sinonim degil.", "SpliceAI")
        # PS3 / PS4 / PP1 : Franklin bunlari ClinVar gonderim metninden cikarir; burada NCI verisiyle yaklasik
        # PS3 / BS3: Franklin'in TP53 sayfalarindaki guc, VCEP Tablo S3 koduyla ayni (V73M BS3, R110H BS3_Supporting,
        # R290H BS3). Once S3 kodu, yoksa Kato 2003 sinifi kullanilir.
        kato = _s(r.get("kato")) or _s(r.get("kato_nci"))
        s3 = _s(r.get("s3_code")) or ""
        m3 = re.match(r"^(PS3|BS3)(?:_(Very ?Strong|Strong|Moderate|Supporting))?$", s3.strip(), re.I)
        ps3 = Evidence("PS3", False, reason=f"VCEP S3: {s3 or 'yok'}; Kato: {kato or 'yok'}.", source="VCEP Tablo S3 / Kato 2003")
        bs3 = Evidence("BS3", False, reason=ps3.reason, source="VCEP Tablo S3 / Kato 2003")
        if m3:
            code3 = m3.group(1).upper()
            st3 = (m3.group(2) or "Strong").lower().replace(" ", "")
            st3 = {"verystrong": "Very Strong", "strong": "Strong", "moderate": "Moderate", "supporting": "Supporting"}[st3]
            e3 = Evidence(code3, True, st3, f"VCEP Tablo S3 fonksiyonel kod: {s3}.", "VCEP Tablo S3")
            ps3, bs3 = (e3, bs3) if code3 == "PS3" else (ps3, e3)
        elif csq == "missense" and kato:
            kl = kato.lower()
            if kl.startswith("non"):
                ps3 = Evidence("PS3", True, "Strong", "Kato 2003 maya transaktivasyonu: fonksiyonsuz.", "NCI TP53 (Kato 2003)")
            elif kl.startswith("partial"):
                ps3 = Evidence("PS3", True, "Supporting", "Kato 2003: kismen fonksiyonel.", "NCI TP53 (Kato 2003)")
            elif kl.startswith("functional") or "supertrans" in kl:
                bs3 = Evidence("BS3", True, "Strong", f"Kato 2003: {kato}.", "NCI TP53 (Kato 2003)")
        ev["PS3"], ev["BS3"] = ps3, bs3
        # BS2: Franklin TP53'te uzman paneli kaydindaki BS2'yi gosteriyor; burada eRepo kodu kullanilir
        bs2c = next((x for x in (_s(r.get("erepo_met")) or "").split(";") if x.strip().startswith("BS2")), None)
        if bs2c:
            st2 = bs2c.strip().split("_")[1].title() if "_" in bs2c else "Strong"
            ev["BS2"] = Evidence("BS2", True, st2, f"ClinGen uzman kaydi: {bs2c.strip()} (60 yas ustu kansersiz tasiyicilar).",
                                 "ClinGen eRepo")
        # PS4: yalnizca Li-Fraumeni fenotipli aileler (LFS, LFL, Chompret); populasyonda sik varyantlarda uygulanmaz
        fam = _f(r.get("nci_germ_lfs_families")) or 0
        common = ev["BA1"].met or ev["BS1"].met
        ps4_st = None if common else "Moderate" if fam >= 5 else "Supporting" if fam >= 2 else None
        ev["PS4"] = Evidence("PS4", ps4_st is not None, ps4_st,
                             f"NCI germline veritabaninda Li-Fraumeni fenotipli (LFS/LFL/Chompret) {int(fam)} aile"
                             + ("; varyant populasyonda sik oldugu icin uygulanmadi." if common and fam >= 2 else ".")
                             + " Franklin bunu ClinVar metninden cikarir.", "NCI TP53 germline - yaklasik")
        ev["PP1"] = Evidence("PP1", False, reason="Kosegregasyon verisi otomatik cikarilamadi (LLM sekmesi ile ClinVar metninden bakin).",
                             source=H)
        for k in ("PS2", "PM3", "PM4", "PM6", "PP4", "BS2", "BS3", "BS4", "BP2", "BP3", "BP5"):
            ev.setdefault(k, Evidence(k, False, reason="Otomatik degerlendirilmedi (kullanici kaniti gerekir).", source=H, auto=False))
        evs = list(ev.values())
        # Franklin catisan kanitlari puanla birlestiriyor (R110H: PM2+PM5 ile BP6+BS3_Sup+BS2_Sup -> Likely Benign);
        # bu, Tavtigian 2020 Bayes puan sistemiyle uyumlu. Richards 2015 tablosu sadece karsilastirma icin dondurulur.
        pts = sum(e.points for e in evs if e.met)
        cls = points_to_class(pts, evs)
        return dict(classification=cls, points=pts, richards=richards_combine(evs), met=[e.label for e in evs if e.met],
                    evidence=[e.to_dict() for e in evs], engine="Franklin-benzeri genel ACMG (puan sistemi)")


def richards_combine(evs):
    """Richards et al. 2015 Tablo 5 birlestirme kurallari (No Influence sayilmaz)."""
    met = [e for e in evs if e.met and e.strength not in (None, "No Influence")]
    def n(prefix_strength):
        return sum(1 for e in met if e.code[0] == prefix_strength[0] and e.strength == prefix_strength[1])
    pvs = sum(1 for e in met if e.code.startswith("P") and e.strength == "Very Strong")
    ps = sum(1 for e in met if e.code.startswith("P") and e.strength == "Strong")
    pm = sum(1 for e in met if e.code.startswith("P") and e.strength == "Moderate")
    pp = sum(1 for e in met if e.code.startswith("P") and e.strength == "Supporting")
    ba = any(e.code == "BA1" for e in met)
    bs = sum(1 for e in met if e.code.startswith("B") and e.strength == "Strong")
    bp = sum(1 for e in met if e.code.startswith("B") and e.strength == "Supporting")
    path = ((pvs >= 1 and (ps >= 1 or pm >= 2 or (pm == 1 and pp == 1) or pp >= 2)) or ps >= 2 or
            (ps == 1 and (pm >= 3 or (pm == 2 and pp >= 2) or (pm == 1 and pp >= 4))))
    lpath = ((pvs >= 1 and pm == 1) or (ps == 1 and 1 <= pm <= 2) or (ps == 1 and pp >= 2) or pm >= 3 or
             (pm == 2 and pp >= 2) or (pm == 1 and pp >= 4))
    ben = ba or bs >= 2
    lben = (bs == 1 and bp == 1) or bp >= 2
    if (path or lpath) and (ben or lben):
        return "Uncertain Significance"
    if path:
        return "Pathogenic"
    if lpath:
        return "Likely Pathogenic"
    if ben:
        return "Benign"
    if lben:
        return "Likely Benign"
    return "Uncertain Significance"
