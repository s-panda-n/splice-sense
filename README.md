# genomic-lm-splice-perturbation

**How Brittle Are Genomic Language Models?**  
A biologically-grounded perturbation sensitivity study on splice site prediction.

**Authors:** Amy Kim · Spandan Patil  
**Course:** AI in Genomics — NYU Courant Institute, Spring 2026

---

## Research Question
Do genomic language models (DNABERT-2, Nucleotide Transformer, Evo) change their
splice site predictions for the right biological reasons — flipping on disease-causing
mutations and staying stable on harmless ones?

## Repository Structure
data/               # processed CSVs (raw files not committed)
models/configs/     # training configuration files
perturbation/       # mutagenesis engine and perturbation type selectors
analysis/           # PSS, FFR/SFR/BSR metrics, saliency maps
figures/            # plotting scripts and output PNGs
notebooks/          # exploratory analysis
slurm/              # HPC job scripts
results/            # experiment output CSVs
## Setup
```bash
conda create -n genomics python=3.10 -y
conda activate genomics
pip install -r requirements.txt
```

## Data
- GENCODE v47 (GTF + GRCh38 FASTA) — download instructions in `data/README.md`
- ClinVar variant summary — download instructions in `data/README.md`
- Raw files are not committed to this repo due to size

## Reproducing Results
```bash
bash run_all.sh
```
