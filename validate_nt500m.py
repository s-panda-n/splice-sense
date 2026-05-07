import torch
import pandas as pd
import numpy as np
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from peft import PeftModel
from sklearn.metrics import accuracy_score, classification_report
from tqdm import tqdm

BATCH_SIZE = 32
LABEL_MAP  = {"donor": 1, "acceptor": 2}

print("Loading NT-500M with LoRA...")
base_model = AutoModelForSequenceClassification.from_pretrained(
    "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species",
    trust_remote_code=True,
    num_labels=3
)
model     = PeftModel.from_pretrained(base_model, "models/checkpoints/nt-500m-splice")
tokenizer = AutoTokenizer.from_pretrained(
    "models/checkpoints/nt-500m-splice", trust_remote_code=True
)
model.eval()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)
print(f"Model loaded on {device}")

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
