import torch
import numpy as np
import pandas as pd
from transformers import AutoTokenizer, BertForSequenceClassification, BertConfig

print("=== Environment Check ===")
print(f"PyTorch: {torch.__version__}")
print(f"Transformers: {__import__('transformers').__version__}")
print(f"NumPy: {np.__version__}")
print(f"CUDA available: {torch.cuda.is_available()}")

print("\n=== Loading DNABERT-2 ===")
tokenizer = AutoTokenizer.from_pretrained(
    "zhihan1996/DNABERT-2-117M", trust_remote_code=True
)
print("Tokenizer loaded OK")

model = BertForSequenceClassification.from_pretrained(
    "zhihan1996/DNABERT-2-117M",
    num_labels=3,
    trust_remote_code=True,
    ignore_mismatched_sizes=True
)
print("Model loaded OK")

print("\n=== Test Inference ===")
test_seq = "ACGTACGT" * 50
inputs = tokenizer(
    test_seq, return_tensors="pt",
    truncation=True, padding="max_length", max_length=512
)
with torch.no_grad():
    logits = model(**inputs).logits

print(f"Output shape: {logits.shape}")
print(f"Logits: {logits}")

print("\n=== Loading Data ===")
df = pd.read_csv("data/processed/splice_sites_10k.csv")
print(f"Splice sites loaded: {df.shape}")

print("\n=== All checks passed ===")
