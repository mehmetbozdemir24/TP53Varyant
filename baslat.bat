@echo off
cd /d "%~dp0"
if not exist data\processed\variants.parquet (
    echo Islenmis veri yok, hat calistiriliyor...
    python run_pipeline.py
)
python -m streamlit run app.py
