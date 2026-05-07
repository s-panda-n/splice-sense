import torch
import pandas as pd
import numpy as np
import argparse
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from perturbation.mutagenesis import get_all_mutations
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("--model", choices=["nt100m", "nt500m"], required=True)
args = parser.parse_args()

MODEL_PATHS = {
    "nt100m": "models/checkpoints/nt-100m-splice",
    "nt500m": "models/checkpoints/nt-500m-splice"
}

MODEL_PATH    = MODEL_PATHS[args.model]
BASELINE_PATH = f"results/baseline_{args.model}.csv"
OUT_PATH      = f"results/perturbation_{args.model}.csv"

print(f"Loading {args.model} from {MODEL_PATH}...")
tokenizer = AutoTokenizer.from_pretrained(MODEL_PATH, trust_remote_code=True)
model     = AutoModelForSequenceClassification.from_pretrained(
    MODEL_PATH, trust_remote_code=True
)
model.eval()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)
print(f"Running on: {device}")

df = pd.read_csv(BASELINE_PATH)
print(f"Running perturbations on {len(df):,} sequences...")

def predict(sequence):
    inputs = tokenizer(
        sequence, return_tensors="pt",
        truncation=True, padding="max_length", max_length=512
    ).to(device)
    with torch.no_grad():
        logits = model(**inputs).logits
    probs = torch.softmax(logits, dim=-1).cpu().numpy()[0]
    return int(np.argmax(probs)), float(np.max(probs))

records = []
for idx, row in tqdm(df.iterrows(), total=len(df)):
    seq           = row["sequence"]
    site_type     = row["site_type"]
    baseline_pred = int(row["baseline_pred"])
    baseline_conf = float(row["baseline_confidence"])

    for mut in get_all_mutations(seq, site_type):
        pred, conf = predict(mut["mutated_seq"])
        records.append({
            "seq_id":            idx,
            "chrom":             row["chrom"],
            "pos":               row["pos"],
            "site_type":         site_type,
            "perturbation_type": mut["type"],
            "mut_position":      mut["position"],
            "ref":               mut["ref"],
            "alt":               mut["alt"],
            "baseline_pred":     baseline_pred,
            "baseline_conf":     baseline_conf,
            "perturbed_pred":    pred,
            "perturbed_conf":    conf,
            "flipped":           baseline_pred != pred,
            "conf_change":       baseline_conf - conf
        })

results = pd.DataFrame(records)
results.to_csv(OUT_PATH, index=False)
print(f"\nSaved {len(results):,} rows to {OUT_PATH}")
for ptype in ["Type1", "Type2", "Type3"]:
    flip = results[results["perturbation_type"]==ptype]["flipped"].mean()*100
    print(f"{ptype} flip rate: {flip:.1f}%")
