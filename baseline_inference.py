import torch
import pandas as pd
import numpy as np
import argparse
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("--model", choices=["nt100m", "nt500m"], required=True)
args = parser.parse_args()

MODEL_PATHS = {
    "nt100m": "models/checkpoints/nt-100m-splice",
    "nt500m": "models/checkpoints/nt-500m-splice"
}

MODEL_PATH = MODEL_PATHS[args.model]
OUT_PATH   = f"results/baseline_{args.model}.csv"
BATCH_SIZE = 64

print(f"Loading {args.model} from {MODEL_PATH}...")
tokenizer = AutoTokenizer.from_pretrained(MODEL_PATH, trust_remote_code=True)
model     = AutoModelForSequenceClassification.from_pretrained(
    MODEL_PATH, trust_remote_code=True
)
model.eval()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)
print(f"Running on: {device}")

df = pd.read_csv("data/processed/splice_sites_10k.csv")
print(f"Running baseline inference on {len(df):,} sequences...")

all_preds, all_confs = [], []

for i in tqdm(range(0, len(df), BATCH_SIZE)):
    batch  = df["sequence"].iloc[i:i+BATCH_SIZE].tolist()
    inputs = tokenizer(
        batch, return_tensors="pt",
        truncation=True, padding="max_length", max_length=512
    ).to(device)
    with torch.no_grad():
        logits = model(**inputs).logits
    probs = torch.softmax(logits, dim=-1).cpu().numpy()
    preds = np.argmax(probs, axis=-1)
    confs = np.max(probs, axis=-1)
    all_preds.extend(preds.tolist())
    all_confs.extend(confs.tolist())

df["baseline_pred"]       = all_preds
df["baseline_confidence"] = all_confs
df.to_csv(OUT_PATH, index=False)

print(f"\nSaved to {OUT_PATH}")
print(f"Prediction distribution:")
print(df["baseline_pred"].value_counts())
# 0=neither, 1=acceptor, 2=donor
print(f"\nMean confidence: {df['baseline_confidence'].mean():.3f}")
