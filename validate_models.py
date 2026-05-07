import torch
import pandas as pd
import numpy as np
import argparse
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from sklearn.metrics import accuracy_score, classification_report
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("--model", choices=["nt100m", "nt500m"], required=True)
args = parser.parse_args()

MODEL_PATHS = {
    "nt100m": "models/checkpoints/nt-100m-splice",
    "nt500m": "models/checkpoints/nt-500m-splice"
}
MODEL_PATH = MODEL_PATHS[args.model]
BATCH_SIZE = 64

# Label mapping used during training
# donor=1, acceptor=2, neither=0
LABEL_MAP = {"donor": 1, "acceptor": 2}

print(f"Loading {args.model}...")
tokenizer = AutoTokenizer.from_pretrained(MODEL_PATH, trust_remote_code=True)
model     = AutoModelForSequenceClassification.from_pretrained(
    MODEL_PATH, trust_remote_code=True
)
model.eval()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

df     = pd.read_csv("data/processed/splice_sites_10k.csv")
sample = df.sample(500, random_state=42)

preds, trues = [], []
for i in tqdm(range(0, len(sample), BATCH_SIZE)):
    batch  = sample.iloc[i:i+BATCH_SIZE]
    inputs = tokenizer(
        batch["sequence"].tolist(), return_tensors="pt",
        truncation=True, padding="max_length", max_length=512
    ).to(device)
    with torch.no_grad():
        logits = model(**inputs).logits
    pred = torch.argmax(logits, dim=-1).cpu().numpy()
    true = batch["site_type"].map(LABEL_MAP).values
    preds.extend(pred.tolist())
    trues.extend(true.tolist())

print(f"\nAccuracy: {accuracy_score(trues, preds):.4f}")
print(classification_report(trues, preds,
      target_names=["neither", "donor", "acceptor"],
      labels=[0, 1, 2]))
