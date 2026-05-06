import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve
import json

results = {}

for model in ["nt100m", "nt500m"]:
    path = f"results/perturbation_{model}.csv"
    try:
        df = pd.read_csv(path)
    except FileNotFoundError:
        print(f"Skipping {model} — results not found")
        continue

    t1 = df[df["perturbation_type"] == "Type1"]
    t2 = df[df["perturbation_type"] == "Type2"]
    t3 = df[df["perturbation_type"] == "Type3"]

    FFR = t1["flipped"].mean()
    SFR = t3["flipped"].mean()
    BSR = FFR / SFR if SFR > 0 else float('inf')

    print(f"\n=== {model.upper()} ===")
    print(f"FFR (consensus flip rate):     {FFR:.3f}")
    print(f"SFR (silent flip rate):        {SFR:.3f}")
    print(f"ESE flip rate (Type2):         {t2['flipped'].mean():.3f}")
    print(f"BSR (FFR/SFR):                 {BSR:.2f}")

    results[model] = {
        "FFR": round(FFR, 4),
        "SFR": round(SFR, 4),
        "ESR": round(t2["flipped"].mean(), 4),
        "BSR": round(BSR, 4)
    }

    # ClinVar AUROC
    try:
        cv   = pd.read_csv("data/processed/clinvar_gencode_merged.csv")
        t1_m = pd.merge(cv, t1, on=["chrom","pos"], how="inner")
        if len(t1_m) > 10:
            auroc = roc_auc_score(t1_m["label"], t1_m["conf_change"])
            auprc = average_precision_score(t1_m["label"], t1_m["conf_change"])
            print(f"ClinVar AUROC:                 {auroc:.3f}")
            print(f"ClinVar AUPRC:                 {auprc:.3f}")
            results[model]["AUROC"] = round(auroc, 4)
            results[model]["AUPRC"] = round(auprc, 4)
    except Exception as e:
        print(f"ClinVar eval failed: {e}")

with open("results/metrics_summary.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nSaved to results/metrics_summary.json")
