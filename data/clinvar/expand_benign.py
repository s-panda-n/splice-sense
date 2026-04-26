import pandas as pd

# Load current filtered ClinVar
cv = pd.read_csv("data/processed/clinvar_splice_filtered.csv")

print(f"Current dataset: {len(cv)} variants")
print(cv["label"].value_counts())

# Separate pathogenic and benign
pathogenic = cv[cv["label"] == 1].copy()
benign     = cv[cv["label"] == 0].copy()

print(f"\nPathogenic: {len(pathogenic)}")
print(f"Benign:     {len(benign)}")

# Strategy: undersample pathogenic to 2000, keep all benign
# Then note in paper that class imbalance exists and AUPRC is primary metric
pathogenic_sample = pathogenic.sample(2000, random_state=42)

# Also grab "Uncertain significance" variants as potential negatives?
# No — too noisy. Just use what we have and note limitation.

balanced = pd.concat([pathogenic_sample, benign], ignore_index=True)
balanced = balanced.sample(frac=1, random_state=42)  # shuffle

print(f"\nBalanced dataset: {len(balanced)} variants")
print(balanced["label"].value_counts())

balanced.to_csv("data/processed/clinvar_balanced.csv", index=False)
print("Saved to data/processed/clinvar_balanced.csv")

# Also save full unbalanced for AUROC calculation
print("\nKeeping full dataset for AUROC — will use AUPRC as primary metric")
print("Note this imbalance in Methods section of paper")
