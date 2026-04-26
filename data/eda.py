import pandas as pd

# ── Splice sites ──────────────────────────
df = pd.read_csv("data/processed/splice_sites_10k.csv")
print("=== SPLICE SITES ===")
print(f"Shape: {df.shape}")
print(f"\nSite type balance:\n{df['site_type'].value_counts()}")
print(f"\nStrand balance:\n{df['strand'].value_counts()}")
print(f"\nSequence lengths: {df['sequence'].str.len().unique()}")
print(f"\nTop chromosomes:\n{df['chrom'].value_counts().head(5)}")

# GT/AG consensus check
donors    = df[df["site_type"] == "donor"]
acceptors = df[df["site_type"] == "acceptor"]
gt_pct = (donors["sequence"].str[200:202] == "GT").mean() * 100
ag_pct = (acceptors["sequence"].str[198:200] == "AG").mean() * 100
print(f"\nGT consensus at donor sites: {gt_pct:.1f}%")
print(f"AG consensus at acceptor sites: {ag_pct:.1f}%")

# ── ClinVar ───────────────────────────────
cv = pd.read_csv("data/processed/clinvar_splice_filtered.csv")
print("\n=== CLINVAR ===")
print(f"Shape: {cv.shape}")
print(f"\nClass balance:\n{cv['ClinicalSignificance'].value_counts()}")
print(f"\nTop genes:\n{cv['GeneSymbol'].value_counts().head(10)}")
print(f"\nMissing values:\n{cv.isnull().sum()[cv.isnull().sum() > 0]}")