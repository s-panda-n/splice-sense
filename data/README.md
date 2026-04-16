# Data Download Instructions

Run these commands on Greene HPC inside a tmux session.

## 1. GENCODE v47

```bash
cd data/raw

# GTF annotation (~50MB)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz

# Reference genome FASTA (~3GB) 
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz

# Decompress and index
gunzip GRCh38.primary_assembly.genome.fa.gz
python -c "from pyfaidx import Fasta; Fasta('GRCh38.primary_assembly.genome.fa')"
```

## 2. ClinVar

```bash
cd data/clinvar
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
gunzip variant_summary.txt.gz
```
