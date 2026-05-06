import torch
import pandas as pd
import numpy as np
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from sklearn.metrics import accuracy_score, classification_report
from tqdm import tqdm

MODEL_PATH = "models/checkpoints/nt-100m-splice"
BATCH_SIZE = 64

tokenizer = AutoTokenizer.from_pretrained(MODEL_PATH, trust_remote_code=True)
model     = AutoModelForSequenceClassification.from_pretrained(
    MODEL_PATH, trust_remote_code=True
)
model.eval()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

df = pd.read_csv("data/processed/splice_sites_10k.csv")
sample = df.sample(500, random_state=42)

preds, trues = [], []
for i in tqdm(range(0, len(sample), BATCH_SIZE)):
    batch = sample.iloc[i:i+BATCH_SIZE]
    inputs = tokenizer(
        batch["sequence"].tolist(), return_tensors="pt",
        truncation=True, padding="max_length", max_length=512
    ).to(device)
    with torch.no_grad():
        logits = model(**inputs).logits
    pred = torch.argmax(logits, dim=-1).cpu().numpy()
    true = batch["site_type"].map({"donor": 2, "acceptor": 1}).values
    preds.extend(pred)
    trues.extend(true)

print(f"\nAccuracy: {accuracy_score(trues, preds):.4f}")
print(classification_report(trues, preds,
      target_names=["neither", "acceptor", "donor"]))
