import gtfparse
from pyfaidx import Fasta
import pandas as pd
from tqdm import tqdm

# Change these lines in data/data_prep_v2.py
GTF_PATH   = "data/raw/gencode.v47.basic.annotation.gtf.gz"
FASTA_PATH = "data/raw/GRCh38.primary_assembly.genome.fa"
WINDOW     = 200

def reverse_complement(seq):
    comp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return ''.join(comp.get(b,'N') for b in reversed(seq.upper()))

print("Loading GTF...")
gtf   = gtfparse.read_gtf(GTF_PATH).to_pandas()
exons = gtf[gtf["feature"] == "exon"].copy()
print(f"Found {len(exons):,} exon records")

print("Loading genome...")
genome = Fasta(FASTA_PATH)

records = []
for _, exon in tqdm(exons.iterrows(), total=len(exons)):
    chrom  = exon["seqname"]
    strand = exon["strand"]
    start  = int(exon["start"])  # 1-based
    end    = int(exon["end"])    # 1-based

    if chrom not in genome:
        continue

    if strand == "+":
        # ── Donor: GT immediately after exon end ──────────────
        # pyfaidx is 0-based: genome[chrom][end] = first intron base
        d_center = end  # 0-based position of G in GT
        d_start  = max(0, d_center - WINDOW)
        d_end    = d_center + WINDOW
        d_seq    = str(genome[chrom][d_start:d_end]).upper()

        # ── Acceptor: AG ends just before exon start ──────────
        # pyfaidx: genome[chrom][start-2:start] = AG (last 2 intron bases)
        a_center = start - 2  # 0-based position of A in AG
        a_start  = max(0, a_center - WINDOW)
        a_end    = a_center + WINDOW
        a_seq    = str(genome[chrom][a_start:a_end]).upper()

        if len(d_seq) == WINDOW * 2:
            records.append({
                "chrom": chrom, "pos": d_center,
                "strand": strand, "site_type": "donor",
                "sequence": d_seq
            })
        if len(a_seq) == WINDOW * 2:
            records.append({
                "chrom": chrom, "pos": a_center,
                "strand": strand, "site_type": "acceptor",
                "sequence": a_seq
            })

    else:  # minus strand
        # ── Donor on minus: AC on forward = GT on minus ───────
        # AC is just before exon start: genome[chrom][start-3:start-1]
        # Center on A (first base of AC): start-3 (0-based)
        d_center = start - 3
        d_start  = max(0, d_center - WINDOW)
        d_end    = d_center + WINDOW
        d_seq    = str(genome[chrom][d_start:d_end]).upper()
        d_seq    = reverse_complement(d_seq)

        # ── Acceptor on minus: CT on forward = AG on minus ────
        # CT is just after exon end: genome[chrom][end:end+2]
        # Center on C: end (0-based)
        a_center = end
        a_start  = max(0, a_center - WINDOW)
        a_end    = a_center + WINDOW
        a_seq    = str(genome[chrom][a_start:a_end]).upper()
        a_seq    = reverse_complement(a_seq)

        if len(d_seq) == WINDOW * 2:
            records.append({
                "chrom": chrom, "pos": d_center,
                "strand": strand, "site_type": "donor",
                "sequence": d_seq
            })
        if len(a_seq) == WINDOW * 2:
            records.append({
                "chrom": chrom, "pos": a_center,
                "strand": strand, "site_type": "acceptor",
                "sequence": a_seq
            })

df = pd.DataFrame(records).drop_duplicates(subset=["chrom", "pos", "site_type"])
print(f"\nTotal extracted: {len(df):,} splice sites")

# ── Filter to canonical GT-AG only ────────────────────────────
donors    = df[df["site_type"] == "donor"]
acceptors = df[df["site_type"] == "acceptor"]

canonical_donors    = donors[donors["sequence"].str[WINDOW:WINDOW+2] == "GT"]
canonical_acceptors = acceptors[acceptors["sequence"].str[WINDOW-2:WINDOW] == "AG"]

df_canonical = pd.concat([canonical_donors, canonical_acceptors])
print(f"After canonical filter: {len(df_canonical):,} splice sites")

gt_pct = (canonical_donors["sequence"].str[200:202] == "GT").mean() * 100
ag_pct = (canonical_acceptors["sequence"].str[198:200] == "AG").mean() * 100
print(f"GT at donor position:    {gt_pct:.1f}%")
print(f"AG at acceptor position: {ag_pct:.1f}%")

# Save
df_canonical.to_csv("data/processed/splice_sites_all.csv", index=False)
df_canonical.sample(10000, random_state=42).to_csv(
    "data/processed/splice_sites_10k.csv", index=False
)
print("\nSaved splice_sites_all.csv and splice_sites_10k.csv")
