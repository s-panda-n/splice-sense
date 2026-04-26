import torch
import pandas as pd
import numpy as np
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from perturbation.mutagenesis import get_all_mutations
from tqdm import tqdm

MODEL_PATH    = "models/checkpoints/dnabert2-splice"
BASELINE_PATH = "results/baseline_predictions.csv"
OUT_PATH      = "results/perturbation_results.csv"
BATCH_SIZE    = 64

print("Loading model...")
tokenizer = AutoTokenizer.from_pretrained(MODEL_PATH, trust_remote_code=True)
model     = AutoModelForSequenceClassification.from_pretrained(
    MODEL_PATH, trust_remote_code=True
)
model.eval()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

print("Loading baseline predictions...")
df = pd.read_csv(BASELINE_PATH)

def predict_single(sequence):
    inputs = tokenizer(
        sequence,
        return_tensors="pt",
        truncation=True,
        padding="max_length",
        max_length=512
    ).to(device)
    with torch.no_grad():
        logits = model(**inputs).logits
    probs = torch.softmax(logits, dim=-1).cpu().numpy()[0]
    pred  = int(np.argmax(probs))
    conf  = float(np.max(probs))
    return pred, conf

records = []
for idx, row in tqdm(df.iterrows(), total=len(df)):
    seq           = row["sequence"]
    site_type     = row["site_type"]
    baseline_pred = int(row["baseline_pred"])
    baseline_conf = float(row["baseline_confidence"])

    mutations = get_all_mutations(seq, site_type)

    for mut in mutations:
        perturbed_pred, perturbed_conf = predict_single(mut["mutated_seq"])
        records.append({
            "seq_id":           idx,
            "chrom":            row["chrom"],
            "pos":              row["pos"],
            "site_type":        site_type,
            "perturbation_type": mut["type"],
            "mut_position":     mut["position"],
            "ref":              mut["ref"],
            "alt":              mut["alt"],
            "baseline_pred":    baseline_pred,
            "baseline_conf":    baseline_conf,
            "perturbed_pred":   perturbed_pred,
            "perturbed_conf":   perturbed_conf,
            "flipped":          baseline_pred != perturbed_pred,
            "conf_change":      baseline_conf - perturbed_conf
        })

results = pd.DataFrame(records)
results.to_csv(OUT_PATH, index=False)
print(f"\nSaved {len(results)} perturbation results to {OUT_PATH}")

# Quick summary
for ptype in ["Type1", "Type2", "Type3"]:
    subset = results[results["perturbation_type"] == ptype]
    ffr = subset["flipped"].mean() * 100
    print(f"{ptype} flip rate: {ffr:.1f}%")
