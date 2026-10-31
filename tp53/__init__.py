"""TP53-Franklin: TP53'e odakli, Franklin (Genoox) benzeri varyant yorumlama araci.

Moduller
--------
common  : yollar, HGVS yardimcilari, amino asit tablolari, Grantham, domainler
build   : ham veri setlerinden birlesik varyant tablosu (data/processed/variants.parquet)
acmg    : ClinGen TP53 VCEP v2.4 kurallari (Bayes puan) + Franklin tarzi genel ACMG modu
ml      : TP53'e ozel "Aggregated Prediction" meta-modeli ve fonksiyon/hipomorf modeli
llm     : ClinVar gonderim metinlerinden yerel LLM (Ollama) ile klinik kanit cikarimi
"""
__version__ = "0.1.0"
