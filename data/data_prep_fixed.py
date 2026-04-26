import gtfparse
from pyfaidx import Fasta
import pandas as pd
from tqdm import tqdm

GTF_PATH   = "data/raw/gencode.v47.basic.annotation.gtf.gz"
FASTA_PATH = "data/raw/GRCh38.primary_assembly.genome.fa"
WINDOW     = 200

def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))

print("Loading GTF...")
gtf = gtfparse.read_gtf(GTF_PATH).to_pandas()
exons = gtf[gtf["feature"] == "exon"].copy()
print(f"Found {len(exons):,} exon records")

print("Loading genome...")
genome = Fasta(FASTA_PATH)

records = []
for _, exon in tqdm(exons.iterrows(), total=len(exons)):
    chrom  = exon["seqname"]
    strand = exon["strand"]
    start  = int(exon["start"])  # 1-based, start of exon
    end    = int(exon["end"])    # 1-based, end of exon

    if chrom not in genome:
        continue

    if strand == "+":
        # Donor site: end of exon, GT follows immediately after
        # Window: 200bp of exon + 200bp of intron
        donor_pos    = end
        acceptor_pos = start - 1  # base before exon start = last intron base

        # Donor: exon context [end-200 : end] + intron context [end : end+200]
        d_start = max(0, donor_pos - WINDOW)
        d_end   = donor_pos + WINDOW
        d_seq   = str(genome[chrom][d_start:d_end])

        # Acceptor: intron context [start-200 : start] + exon context [start : start+200]  
        a_start = max(0, (start - 1) - WINDOW)
        a_end   = (start - 1) + WINDOW
        a_seq   = str(genome[chrom][a_start:a_end])

        if len(d_seq) == WINDOW * 2:
            records.append({
                "chrom": chrom, "pos": donor_pos,
                "strand": strand, "site_type": "donor", "sequence": d_seq
            })
        if len(a_seq) == WINDOW * 2:
            records.append({
                "chrom": chrom, "pos": start - 1,
                "strand": strand, "site_type": "acceptor", "sequence": a_seq
            })

    else:  # minus strand
        # On minus strand, exon start/end are flipped relative to transcription
        # Donor site: start of exon (lowest coordinate) — GT on minus strand = AC on forward
        # Acceptor site: end of exon (highest coordinate) — AG on minus strand = CT on forward

        # Donor: centered on exon start, reverse complemented
        d_start = max(0, start - 1 - WINDOW)
        d_end   = start - 1 + WINDOW
        d_seq   = str(genome[chrom][d_start:d_end])
        d_seq   = reverse_complement(d_seq)

        # Acceptor: centered on exon end, reverse complemented
        a_start = max(0, end - WINDOW)
        a_end   = end + WINDOW
        a_seq   = str(genome[chrom][a_start:a_end])
        a_seq   = reverse_complement(a_seq)

        if len(d_seq) == WINDOW * 2:
            records.append({
                "chrom": chrom, "pos": start,
                "strand": strand, "site_type": "donor", "sequence": d_seq
            })
        if len(a_seq) == WINDOW * 2:
            records.append({
                "chrom": chrom, "pos": end,
                "strand": strand, "site_type": "acceptor", "sequence": a_seq
            })

df = pd.DataFrame(records).drop_duplicates(subset=["chrom", "pos", "site_type"])
print(f"\nExtracted {len(df):,} splice sites")

# Validate before saving
donors    = df[df["site_type"] == "donor"]
acceptors = df[df["site_type"] == "acceptor"]
gt_pct = (donors["sequence"].str[200:202] == "GT").mean() * 100
ag_pct = (acceptors["sequence"].str[198:200] == "AG").mean() * 100
print(f"GT at donor position:    {gt_pct:.1f}%  (expect >85%)")
print(f"AG at acceptor position: {ag_pct:.1f}%  (expect >85%)")

if gt_pct > 80 and ag_pct > 80:
    df.to_csv("data/processed/splice_sites_all.csv", index=False)
    df.sample(10000, random_state=42).to_csv("data/processed/splice_sites_10k.csv", index=False)
    print("\nSaved splice_sites_all.csv and splice_sites_10k.csv")
else:
    print("\nWARNING: Consensus too low — check extraction logic before saving")
