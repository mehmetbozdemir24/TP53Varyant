"""Uctan uca arayuz testi (Streamlit AppTest): farkli turde varyantlari girer, tum sekmeleri cizdirir,
hata/uyari olup olmadigini ve temel sonuclari raporlar.

    python scripts/test_app.py
"""
from __future__ import annotations

import re
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from streamlit.testing.v1 import AppTest  # noqa: E402

CASES = [
    # (girdi, beklenen sonuc tipi veya None, aciklama)
    ("R175H", "missense", "hotspot missense, kisa protein yazimi"),
    ("p.Arg248Gln", "missense", "uzun protein yazimi"),
    ("TP53:c.743G>A", "missense", "Franklin baslik yazimi"),
    ("chr17-7577538-C-T", "missense", "hg19 koordinat"),
    ("17-7674220-C-T", "missense", "hg38 koordinat"),
    ("rs1042522", "missense", "rsID (P72R)"),
    ("  r175h ", "missense", "kucuk harf + bosluk"),
    ("P72R", "missense", "yaygin polimorfizm"),
    ("R337H", "missense", "dusuk penetranslı kurucu varyant"),
    ("V217M", "missense", "celiskili ClinVar"),
    ("R213*", "nonsense", "nonsense"),
    ("p.Arg213Ter", "nonsense", "nonsense uzun yazim"),
    ("c.215del", "frameshift", "frameshift"),
    ("c.1179dup", None, "stop uzamasi"),
    ("c.375+1G>A", "splice_canonical", "kanonik splice"),
    ("c.375+5G>A", None, "splice bolgesi"),
    ("c.672+62A>G", None, "intronik"),
    ("c.*1175A>C", "3_utr", "3' UTR benign"),
    ("M1I", None, "start kaybi"),
    ("R213R", None, "sinonim"),
    ("G245S", "missense", "hotspot"),
    ("c.1101-2A>G", None, "splice akseptor"),
    ("R999H", "YOK", "gecersiz kodon"),
    ("BRCA1", "YOK", "anlamsiz girdi"),
    ("", "YOK", "bos girdi"),
]


def run_case(q):
    at = AppTest.from_file(str(ROOT / "app.py"), default_timeout=120)
    at.session_state["q"] = q
    at.run()
    return at


def summarize(at):
    errs = [e.value for e in at.error]
    excs = [str(e.value)[:300] for e in at.exception]
    head = next((m.value for m in at.markdown if 'class="vhead"' in m.value), "")
    metrics = {m.label: m.value for m in at.metric}
    return errs, excs, head, metrics


def main():
    bad = 0
    t0 = time.time()
    for q, exp, note in CASES:
        at = run_case(q)
        errs, excs, head, met = summarize(at)
        csq = met.get("Sonuç tipi")
        ok = not excs
        if exp == "YOK":
            ok = ok and any("bulunamad" in e for e in errs)
        elif exp:
            ok = ok and csq == exp and not errs
        else:
            ok = ok and not errs and bool(head)
        cls = [m.value for m in at.markdown if "Pathogenic" in m.value or "Benign" in m.value or "Uncertain" in m.value]
        flag = "OK " if ok else "HATA"
        bad += not ok
        title = head.split("<h1>")[1].split("</h1>")[0] if "<h1>" in head else "-"
        print(f"[{flag}] {q!r:24} {note:28} -> {title:22} csq={csq} ML={met.get('TP53 ML skoru')}")
        for e in excs:
            print("       ISTISNA:", e)
        for e in errs:
            if exp != "YOK":
                print("       HATA MESAJI:", e)
    print(f"\n{len(CASES)} vaka, {bad} hata, {time.time() - t0:.0f} sn")

    # klinik form: LLM onerisinin forma aktarilmasi (pending_clinical) puani degistirmeli
    at = run_case("R175H")
    before = [m.value for m in at.markdown if "Önerilen sınıflandırma" in m.value][:1]
    at.session_state["pending_clinical"] = dict(chompret=3, classic=1, meioses=5, families=1, vaf525=1, use_expert=False)
    at.run()
    after = [m.value for m in at.markdown if "Önerilen sınıflandırma" in m.value][:1]
    print("\nForma uygula testi (R175H):")
    pts = lambda h: re.sub(r"<[^>]+>", " ", h[0])[:200] if h else "-"
    print("  once :", " ".join(pts(before).split()))
    print("  sonra:", " ".join(pts(after).split()))
    print("  istisna:", [str(e.value)[:200] for e in at.exception] or "yok")
    return bad


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
