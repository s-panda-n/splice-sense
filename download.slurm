#!/bin/bash
#SBATCH --job-name=download_data
#SBATCH --account=bmsc_ga_4551-2026sp
#SBATCH --partition=n2c48m24
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=slurm/logs/download.out

cd /scratch/spp9400/splice-sense/data/raw

# GTF
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz

# FASTA
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz

# Decompress and index
gunzip GRCh38.primary_assembly.genome.fa.gz
source /share/apps/pyenv/py3.9/etc/profile.d/conda.sh
conda activate genomics
python -c "from pyfaidx import Fasta; Fasta('GRCh38.primary_assembly.genome.fa'); print('indexed')"

echo "Done"
