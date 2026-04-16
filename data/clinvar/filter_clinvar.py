# data/clinvar/filter_clinvar.py
import pandas as pd

df = pd.read_csv("data/clinvar/variant_summary.txt", sep="\t", low_memory=False)

filtered = df[
    (df["Type"] == "single nucleotide variant") &
    (df["Assembly"] == "GRCh38") &
    (df["ClinicalSignificance"].isin(["Pathogenic", "Benign"])) &
    (df["Name"].str.contains("splice", case=False, na=False))
]

print(f"Filtered to {len(filtered):,} splice site variants")
print(filtered[["GeneSymbol", "ClinicalSignificance", "Chromosome", "Start", "ReferenceAllele", "AlternateAllele"]].head(10))

filtered.to_csv("data/processed/clinvar_splice_filtered.csv", index=False)
print("Saved to data/processed/clinvar_splice_filtered.csv")