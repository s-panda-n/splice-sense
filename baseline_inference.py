import torch
import pandas as pd
import numpy as np
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from tqdm import tqdm

MODEL_PATH = "models/checkpoints/dnabert2-splice"
DATA_PATH  = "data/processed/splice_sites_10k.csv"
OUT_PATH   = "results/baseline_predictions.csv"
BATCH_SIZE = 64

print("Loading model...")
tokenizer = AutoTokenizer.from_pretrained(MODEL_PATH, trust_remote_code=True)
model     = AutoModelForSequenceClassification.from_pretrained(
    MODEL_PATH, trust_remote_code=True
)
model.eval()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)
print(f"Running on: {device}")

print("Loading sequences...")
df = pd.read_csv(DATA_PATH)
print(f"Sequences: {len(df)}")

all_preds   = []
all_confs   = []
all_logits  = []

for i in tqdm(range(0, len(df), BATCH_SIZE)):
    batch_seqs = df["sequence"].iloc[i:i+BATCH_SIZE].tolist()
    inputs = tokenizer(
        batch_seqs,
        return_tensors="pt",
        truncation=True,
        padding="max_length",
        max_length=512
    ).to(device)

    with torch.no_grad():
        outputs = model(**inputs)

    logits = outputs.logits.cpu().numpy()
    probs  = torch.softmax(outputs.logits, dim=-1).cpu().numpy()
    preds  = np.argmax(logits, axis=-1)
    confs  = np.max(probs, axis=-1)

    all_preds.extend(preds.tolist())
    all_confs.extend(confs.tolist())
    all_logits.extend(logits.tolist())

df["baseline_pred"]       = all_preds
df["baseline_confidence"] = all_confs
df["baseline_logits"]     = all_logits

df.to_csv(OUT_PATH, index=False)
print(f"\nSaved to {OUT_PATH}")
print(f"\nPrediction distribution:")
print(df["baseline_pred"].value_counts())
# 0 = neither, 1 = acceptor, 2 = donor
