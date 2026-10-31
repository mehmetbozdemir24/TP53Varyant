import importlib, json, urllib.request
for m in ["streamlit", "altair", "plotly", "scipy", "pyarrow", "sklearn", "torch", "transformers", "xgboost", "lightgbm", "openpyxl", "requests"]:
    try:
        mod = importlib.import_module(m)
        extra = ""
        if m == "torch":
            extra = f" cuda={mod.cuda.is_available()}" + (f" {mod.cuda.get_device_name(0)}" if mod.cuda.is_available() else "")
        print(f"{m:13s} {getattr(mod, '__version__', '?')}{extra}", flush=True)
    except Exception as e:
        print(f"{m:13s} YOK ({e.__class__.__name__})", flush=True)
try:
    j = json.load(urllib.request.urlopen("http://localhost:11434/api/tags", timeout=5))
    print("ollama:", [x["name"] for x in j.get("models", [])])
except Exception as e:
    print("ollama: calismiyor", e.__class__.__name__)
