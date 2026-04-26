import pandas as pd
import re

print("Loading ClinVar...")
df = pd.read_csv("data/clinvar/variant_summary.txt", sep="\t", low_memory=False)

def is_splice_site(name):
    if pd.isna(name):
        return False
    return bool(re.search(r'c\.\d+[+-][12][ACGT]>', str(name)))

print("Filtering...")
mask = (
    df["Name"].apply(is_splice_site) &
    (df["Type"] == "single nucleotide variant") &
    (df["Assembly"] == "GRCh38") &
    (df["ClinicalSignificance"].isin([
        "Pathogenic",
        "Likely pathogenic",
        "Benign",
        "Likely benign"
    ]))
)

filtered = df[mask].copy()

# Simplify labels — binary classification
filtered["label"] = filtered["ClinicalSignificance"].map({
    "Pathogenic":        1,
    "Likely pathogenic": 1,
    "Benign":            0,
    "Likely benign":     0
})

print(f"\nFinal dataset: {len(filtered)} variants")
print(filtered["ClinicalSignificance"].value_counts())
print(f"\nBinary label balance:")
print(filtered["label"].value_counts())

filtered.to_csv("data/processed/clinvar_splice_filtered.csv", index=False)
print("\nSaved to data/processed/clinvar_splice_filtered.csv")
