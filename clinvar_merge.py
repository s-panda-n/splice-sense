import pandas as pd

print("Loading data...")
# Use ALL splice sites not just 10k sample
splice = pd.read_csv("data/processed/splice_sites_all.csv")
cv     = pd.read_csv("data/processed/clinvar_balanced.csv")

print(f"Splice sites: {len(splice):,}")
print(f"ClinVar variants: {len(cv):,}")

# Fix chromosome format
cv["Chromosome"] = "chr" + cv["Chromosome"].astype(str)

# Exact match first
merged = pd.merge(
    splice, cv,
    left_on=["chrom", "pos"],
    right_on=["Chromosome", "Start"],
    how="inner"
)
print(f"Exact match overlaps: {len(merged)}")

# If still low, try ±2bp
if len(merged) < 100:
    print("Trying ±2bp window...")
    splice["pos_min"] = splice["pos"] - 2
    splice["pos_max"] = splice["pos"] + 2

    records = []
    cv_grouped = cv.groupby("Chromosome")

    for chrom, cv_chrom in cv_grouped:
        splice_chrom = splice[splice["chrom"] == chrom]
        if len(splice_chrom) == 0:
            continue
        for _, cv_row in cv_chrom.iterrows():
            matches = splice_chrom[
                (splice_chrom["pos_min"] <= cv_row["Start"]) &
                (splice_chrom["pos_max"] >= cv_row["Start"])
            ]
            for _, sp_row in matches.iterrows():
                records.append({**sp_row.to_dict(), **cv_row.to_dict()})

    merged = pd.DataFrame(records)
    print(f"Window match overlaps: {len(merged)}")

print(f"\nFinal overlaps: {len(merged)}")
print(f"Label balance:\n{merged['label'].value_counts()}")
print(merged[["chrom","pos","site_type","ClinicalSignificance","label"]].head(10))

merged.to_csv("data/processed/clinvar_gencode_merged.csv", index=False)
print("Saved.")
