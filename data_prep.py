# data_prep.py
import gtfparse
from pyfaidx import Fasta
import pandas as pd
from tqdm import tqdm

GTF_PATH   = "data/raw/gencode.v47.basic.annotation.gtf.gz"
FASTA_PATH = "data/raw/GRCh38.primary_assembly.genome.fa"
OUT_PATH   = "data/processed/splice_sites.csv"
WINDOW     = 200

print("Loading GTF...")
gtf = gtfparse.read_gtf(GTF_PATH)

# Convert to pandas first, then filter
gtf = gtf.to_pandas()
exons = gtf[gtf["feature"] == "exon"].copy()
print(f"Found {len(exons):,} exon records")

print("Loading genome...")
genome = Fasta(FASTA_PATH)

records = []
for _, exon in tqdm(exons.iterrows(), total=len(exons)):
    chrom  = exon["seqname"]
    strand = exon["strand"]
    start  = int(exon["start"])
    end    = int(exon["end"])

    if chrom not in genome:
        continue

    for site_pos, site_type in [(end, "donor"), (start, "acceptor")]:
        w_start = max(0, site_pos - WINDOW)
        w_end   = site_pos + WINDOW
        seq     = str(genome[chrom][w_start:w_end])

        if len(seq) == WINDOW * 2:
            records.append({
                "chrom":     chrom,
                "pos":       site_pos,
                "strand":    strand,
                "site_type": site_type,
                "sequence":  seq
            })

df = pd.DataFrame(records).drop_duplicates(subset=["chrom", "pos", "site_type"])
print(f"Extracted {len(df):,} splice sites")

df.to_csv(OUT_PATH, index=False)
df.sample(10000, random_state=42).to_csv("data/processed/splice_sites_10k.csv", index=False)
print("Done — saved splice_sites.csv and splice_sites_10k.csv")
