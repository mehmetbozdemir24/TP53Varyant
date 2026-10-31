"""Ortak sabitler ve yardimci fonksiyonlar."""
from __future__ import annotations

import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
MODELS = ROOT / "models"
PROC.mkdir(parents=True, exist_ok=True)
MODELS.mkdir(parents=True, exist_ok=True)

TRANSCRIPT = "NM_000546.6"
ENST = "ENST00000269305"
UNIPROT = "P04637"

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
       "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
       "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*", "Sec": "U"}
AA1 = {v: k for k, v in AA3.items()}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

# ClinGen TP53 VCEP domain tanimlari (CSpec GN009, PVS1 bolumu)
DOMAINS = [
    ("TAD1", 17, 25), ("TAD2", 48, 56), ("Prolin bolgesi", 64, 92), ("DNA baglanma domaini", 100, 292),
    ("Mentese (hinge)", 293, 324), ("Oligomerizasyon domaini", 325, 356), ("C-terminal (bazik) domain", 368, 387),
]
HOTSPOT_CODONS = {175, 245, 248, 249, 273, 282}  # VCEP PM1 sabit hotspot kodonlari

# NM_000546 kodlayan ekson sinirlari (c. koordinat); ekson 1 kodlamayan
EXONS = [(2, 1, 74), (3, 75, 96), (4, 97, 375), (5, 376, 559), (6, 560, 672), (7, 673, 782),
         (8, 783, 919), (9, 920, 993), (10, 994, 1100), (11, 1101, 1182)]


def domain_of(codon):
    if codon is None or (isinstance(codon, float) and math.isnan(codon)):
        return None
    for name, a, b in DOMAINS:
        if a <= codon <= b:
            return name
    return "Diger"


def exon_of_cpos(cpos: int):
    for ex, a, b in EXONS:
        if a <= cpos <= b:
            return ex
    return None


# ------------------------------------------------------------------ HGVS
def norm_c(s):
    """'NM_000546.5:c.215delC' -> 'c.215del'. Bos/uyumsuzsa None."""
    if s is None or (isinstance(s, float) and math.isnan(s)):
        return None
    s = str(s).strip()
    if ":" in s:
        s = s.split(":", 1)[1]
    s = s.split(" ")[0]
    if not s.startswith("c."):
        return None
    s = re.sub(r"(del|dup)[ACGTN]+$", r"\1", s)        # delC -> del, dupA -> dup
    s = re.sub(r"del[ACGTN]+ins", "delins", s)          # delCinsA -> delinsA
    return s


def p3_to_p1(p):
    """'p.Arg175His' -> 'R175H'; 'p.Arg213Ter'/'p.Arg213*' -> 'R213*'; fs, del vb. basitce kisaltilir."""
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return None
    p = str(p).strip().replace("(", "").replace(")", "")
    if p.startswith("p."):
        p = p[2:]
    out = p
    for k, v in AA3.items():
        out = out.replace(k, v)
    return out


def p1_to_p3(p1):
    if not p1:
        return None
    m = re.match(r"^([A-Z\*])(\d+)([A-Z\*=]|fs.*|del.*|dup.*|ins.*|delins.*)$", p1)
    if not m:
        return "p." + p1
    a, n, b = m.groups()
    a3 = AA1.get(a, a)
    if len(b) == 1 and b not in "=":
        b = AA1.get(b, b)
    return f"p.{a3}{n}{b}"


MISSENSE_RE = re.compile(r"^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])$")


def parse_cpos(c):
    """'c.524G>A' -> (524, 0); 'c.375+1G>A' -> (375, 1); 'c.376-2A>G' -> (376, -2). UTR icin None."""
    if not c:
        return None
    m = re.match(r"^c\.(\d+)([+-]\d+)?", c)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2) or 0)


def consequence(c, p1):
    """c. ve p.(tek harf) uzerinden sonuc tipi."""
    c = c if isinstance(c, str) else None
    p1 = p1 if isinstance(p1, str) else None
    pos = parse_cpos(c) if c else None
    m = re.match(r"^c\.[-*]?\d+([+-]\d+)", c or "")
    if m and "_" not in c.split("+")[0][2:] or (m and re.match(r"^c\.[-*]?\d+[+-]\d+[ACGT]>", c)):
        off = int(m.group(1))
        return "splice_canonical" if abs(off) <= 2 else ("splice_region" if abs(off) <= 8 else "intronic")
    if c and c.startswith("c.-"):
        return "5_utr"
    if c and c.startswith("c.*"):
        return "3_utr"
    if pos and pos[1] != 0:
        return "splice_canonical" if abs(pos[1]) <= 2 else ("splice_region" if abs(pos[1]) <= 8 else "intronic")
    if p1:
        if re.match(r"^M1(\D|$)", p1) and "fs" not in p1:
            return "start_lost"
        if "fs" in p1:
            return "frameshift"
        if p1.endswith("*") and not p1.startswith("*"):
            return "nonsense"
        if p1.startswith("*"):
            return "stop_lost"
        if p1.endswith("=") or (MISSENSE_RE.match(p1) and p1[0] == p1[-1]):
            return "synonymous"
        if MISSENSE_RE.match(p1):
            return "missense"
        if "delins" in p1:
            return "inframe_delins"
        if "del" in p1:
            return "inframe_deletion"
        if "dup" in p1 or "ins" in p1:
            return "inframe_insertion"
    if c and ("del" in c or "dup" in c or "ins" in c):
        return "indel_other"
    return "other"


# ------------------------------------------------------------------ Grantham
_G = """S R 110|S L 145|S P 74|S T 58|S A 99|S V 124|S G 56|S I 142|S F 155|S Y 144|S C 112|S H 89|S Q 68|S N 46|S K 121|S D 65|S E 80|S M 135|S W 177|R L 102|R P 103|R T 71|R A 112|R V 96|R G 125|R I 97|R F 97|R Y 77|R C 180|R H 29|R Q 43|R N 86|R K 26|R D 96|R E 54|R M 91|R W 101|L P 98|L T 92|L A 96|L V 32|L G 138|L I 5|L F 22|L Y 36|L C 198|L H 99|L Q 113|L N 153|L K 107|L D 172|L E 138|L M 15|L W 61|P T 38|P A 27|P V 68|P G 42|P I 95|P F 114|P Y 110|P C 169|P H 77|P Q 76|P N 91|P K 103|P D 108|P E 93|P M 87|P W 147|T A 58|T V 69|T G 59|T I 89|T F 103|T Y 92|T C 149|T H 47|T Q 42|T N 65|T K 78|T D 85|T E 65|T M 81|T W 128|A V 64|A G 60|A I 94|A F 113|A Y 112|A C 195|A H 86|A Q 91|A N 111|A K 106|A D 126|A E 107|A M 84|A W 148|V G 109|V I 29|V F 50|V Y 55|V C 192|V H 84|V Q 96|V N 133|V K 97|V D 152|V E 121|V M 21|V W 88|G I 135|G F 153|G Y 147|G C 159|G H 98|G Q 87|G N 80|G K 127|G D 94|G E 98|G M 127|G W 184|I F 21|I Y 33|I C 198|I H 94|I Q 109|I N 149|I K 102|I D 168|I E 134|I M 10|I W 61|F Y 22|F C 205|F H 100|F Q 116|F N 158|F K 102|F D 177|F E 140|F M 28|F W 40|Y C 194|Y H 83|Y Q 99|Y N 143|Y K 85|Y D 160|Y E 122|Y M 36|Y W 37|C H 174|C Q 154|C N 139|C K 202|C D 154|C E 170|C M 196|C W 215|H Q 24|H N 68|H K 32|H D 81|H E 40|H M 87|H W 115|Q N 46|Q K 53|Q D 61|Q E 29|Q M 101|Q W 130|N K 94|N D 23|N E 42|N M 142|N W 174|K D 101|K E 56|K M 95|K W 110|D E 45|D M 160|D W 181|E M 126|E W 152|M W 67"""
GRANTHAM = {}
for _t in _G.split("|"):
    _a, _b, _v = _t.split()
    GRANTHAM[(_a, _b)] = GRANTHAM[(_b, _a)] = int(_v)


def grantham(a, b):
    if a == b:
        return 0
    return GRANTHAM.get((a, b))


# Kyte-Doolittle hidrofobisite, hacim (A^3), yuk (pH7)
HYDRO = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5, L=3.8, K=-3.9,
             M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
VOLUME = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4, G=60.1, H=153.2, I=166.7,
              L=166.7, K=168.6, M=162.9, F=189.9, P=112.7, S=89.0, T=116.1, W=227.8, Y=193.6, V=140.0)
CHARGE = dict(R=1, K=1, H=0.1, D=-1, E=-1)


def is_nan(x):
    return x is None or (isinstance(x, float) and math.isnan(x))


CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for _x in _b:
    for _y in _b:
        for _z in _b:
            CODON[_x + _y + _z] = _aa[_i]
            _i += 1


def translate(seq):
    return "".join(CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


HG19_HG38_OFFSET = 96682   # TP53 bolgesinde hg38 = hg19 + 96682 (bitisik liftover)
