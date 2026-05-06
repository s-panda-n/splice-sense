import torch
import pandas as pd
import numpy as np
from transformers import (AutoTokenizer, AutoModelForSequenceClassification,
                           TrainingArguments, Trainer)
from torch.utils.data import Dataset
from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import train_test_split

MODEL_NAME = "InstaDeepAI/nucleotide-transformer-v2-100m-multi-species"

print("Loading data...")
df = pd.read_csv("data/processed/splice_sites_all.csv")
print(f"Total canonical splice sites: {len(df):,}")

df["label"] = df["site_type"].map({"donor": 1, "acceptor": 2})

print("Generating negative examples...")
negatives = df.sample(len(df) // 2, random_state=42).copy()
negatives["sequence"] = negatives["sequence"].apply(
    lambda s: ''.join(np.random.permutation(list(s)))
)
negatives["label"] = 0

all_data = pd.concat([df[["sequence","label"]], negatives[["sequence","label"]]])
all_data = all_data.sample(frac=1, random_state=42).reset_index(drop=True)

train_df, test_df = train_test_split(all_data, test_size=0.1, random_state=42,
                                      stratify=all_data["label"])
train_sample = train_df.sample(30000, random_state=42)
test_sample  = test_df.sample(5000,  random_state=42)
print(f"Train: {len(train_sample):,} | Test: {len(test_sample):,}")

class SpliceDataset(Dataset):
    def __init__(self, sequences, labels, tokenizer, max_length=512):
        self.encodings = tokenizer(
            sequences.tolist(),
            truncation=True,
            padding="max_length",
            max_length=max_length
        )
        self.labels = labels.tolist()

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        item = {k: torch.tensor(v[idx]) for k, v in self.encodings.items()}
        item["labels"] = torch.tensor(self.labels[idx])
        return item

print(f"Loading {MODEL_NAME}...")
tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME, trust_remote_code=True)
model     = AutoModelForSequenceClassification.from_pretrained(
    MODEL_NAME, trust_remote_code=True, num_labels=3
)
print("Model loaded OK")

print("Tokenizing...")
train_dataset = SpliceDataset(train_sample["sequence"], train_sample["label"], tokenizer)
test_dataset  = SpliceDataset(test_sample["sequence"],  test_sample["label"],  tokenizer)

def compute_metrics(eval_pred):
    logits, labels = eval_pred
    preds = np.argmax(logits, axis=-1)
    return {
        "accuracy": accuracy_score(labels, preds),
        "f1_macro": f1_score(labels, preds, average="macro")
    }

args = TrainingArguments(
    output_dir="models/checkpoints/nt-100m-splice",
    num_train_epochs=3,
    per_device_train_batch_size=32,
    per_device_eval_batch_size=64,
    learning_rate=3e-5,
    evaluation_strategy="epoch",
    save_strategy="epoch",
    load_best_model_at_end=True,
    metric_for_best_model="accuracy",
    logging_steps=50,
    report_to="none"
)

trainer = Trainer(
    model=model,
    args=args,
    train_dataset=train_dataset,
    eval_dataset=test_dataset,
    compute_metrics=compute_metrics
)

print("Starting fine-tuning...")
trainer.train()
trainer.save_model("models/checkpoints/nt-100m-splice")
tokenizer.save_pretrained("models/checkpoints/nt-100m-splice")
print("Done — saved to models/checkpoints/nt-100m-splice")
