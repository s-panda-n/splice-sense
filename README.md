```markdown
# splice-sense

**How Brittle Are Genomic Language Models?**  
A biologically-grounded perturbation sensitivity study on splice site prediction.

**Authors:** Amy Kim · Spandan Patil  
**Course:** AI in Genomics · NYU Courant Institute · Spring 2026

---

## Research Question

Do genomic language models change their splice site predictions for the **right biological reasons** — flipping on mutations that destroy the GT/AG consensus signal, while remaining stable on biologically harmless changes?

We evaluate two Nucleotide Transformer models (100M and 500M parameters) and validate findings against real human variants from ClinVar with known pathogenicity labels.

---

## Key Results

| Metric | NT-100M | NT-500M |
|---|---|---|
| Fine-tuning accuracy | 99.2% | 99.1% |
| FFR (consensus flip rate) | 19.8% | 18.6% |
| SFR (silent flip rate) | 0.2% | 0.2% |
| BSR (FFR / SFR) | 106.0 | 100.5 |
| ClinVar AUROC | 0.587 | 0.841 |
| ClinVar AUPRC | 0.960 | 0.988 |
| Mean confidence change (Type 1) | 0.032 | 0.062 |

---

## Repository Structure

```
splice-sense/
├── data/
│   ├── raw/                          # GENCODE files (not committed — see Step 1)
│   ├── processed/                    # Processed CSVs (not committed — see Step 2)
│   └── clinvar/                      # ClinVar raw download + filter scripts
│       ├── filter_clinvar.py         # Filter ClinVar to splice site variants
│       └── expand_benign.py          # Balance pathogenic/benign classes
├── models/
│   └── checkpoints/                  # Fine-tuned weights (not committed — see Step 3)
├── perturbation/
│   ├── __init__.py
│   └── mutagenesis.py                # Type 1/2/3 perturbation engine
├── analysis/
│   └── metrics.py                    # FFR, SFR, BSR, AUROC/AUPRC computation
├── figures/
│   ├── plot_all.py                   # Generates all 4 paper figures
│   └── outputs/                      # Generated figure PNGs
├── results/
│   ├── baseline_nt100m.csv           # Baseline predictions — NT-100M
│   ├── baseline_nt500m.csv           # Baseline predictions — NT-500M
│   ├── perturbation_nt100m.csv       # Perturbation results — NT-100M (360k rows)
│   ├── perturbation_nt500m.csv       # Perturbation results — NT-500M (360k rows)
│   ├── metrics_summary.json          # FFR, SFR, BSR, AUROC per model
│   └── analysis/                     # Extended analysis outputs
├── slurm/                            # NYU Greene HPC job scripts
├── data_prep_v2.py                   # Extract splice sites from GENCODE
├── clinvar_merge.py                  # Merge ClinVar variants with GENCODE sequences
├── finetune_nt100m.py                # Fine-tune NT-100M on splice site classification
├── finetune_nt500m.py                # Fine-tune NT-500M with LoRA
├── load_model.py                     # Unified model loader (handles LoRA for NT-500M)
├── baseline_inference.py             # Baseline predictions on 10k sequences
├── run_perturbation.py               # Full perturbation experiment
├── validate_models.py                # Validate NT-100M accuracy
├── validate_nt500m.py                # Validate NT-500M (LoRA-specific loader)
└── requirements.txt
```

---

## Environment Setup

Python 3.10.14 is required. Other versions will cause compatibility issues with the model libraries.

```bash
conda create --prefix /path/to/envs/genomics python=3.10.14 -y
conda activate /path/to/envs/genomics

pip install \
    torch==2.1.0 \
    transformers==4.40.0 \
    tokenizers==0.19.1 \
    datasets==2.19.0 \
    peft==0.9.0 \
    accelerate==0.29.3 \
    scikit-learn==1.4.2 \
    "numpy==1.26.4" \
    pandas==2.2.2 \
    einops==0.7.0 \
    pyfaidx==0.7.2 \
    gtfparse==2.1.0 \
    matplotlib==3.8.4 \
    seaborn==0.13.2 \
    tqdm==4.66.2
```

Verify the install:
```bash
python -c "import torch, transformers; print(torch.__version__, transformers.__version__)"
# Expected: 2.1.0+cu121  4.40.0
```

> **Note:** All experiments were run on NYU Greene HPC with NVIDIA A100 40GB GPUs. A GPU with at least 16GB VRAM is required for fine-tuning. Inference-only steps (Steps 4-5) require ~8GB VRAM.

---

## Step 1 — Download Raw Data

All data is publicly available with no login required.

**GENCODE v47:**
```bash
mkdir -p data/raw && cd data/raw

# GTF annotation (~50MB)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz

# Reference genome FASTA (~3GB, takes 20-40 mins)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

# Index the FASTA
python -c "from pyfaidx import Fasta; Fasta('GRCh38.primary_assembly.genome.fa'); print('Indexed OK')"
cd ../..
```

**ClinVar:**
```bash
mkdir -p data/clinvar && cd data/clinvar
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
gunzip variant_summary.txt.gz
cd ../..
```

---

## Step 2 — Generate Processed Data

```bash
# Extract canonical GT-AG splice sites from GENCODE
# Runtime: ~30-45 mins | RAM: 32GB+
python data_prep_v2.py
# Output: data/processed/splice_sites_all.csv  (331,332 sequences)
#         data/processed/splice_sites_10k.csv  (10,000 sample for experiments)

# Filter ClinVar to canonical splice site variants
python data/clinvar/filter_clinvar.py
# Output: data/processed/clinvar_splice_filtered.csv (36,347 variants)

# Balance pathogenic/benign classes
python data/clinvar/expand_benign.py
# Output: data/processed/clinvar_balanced.csv (2,250 variants)

# Merge ClinVar with GENCODE splice sites on genomic coordinates
python clinvar_merge.py
# Output: data/processed/clinvar_gencode_merged.csv (1,146 overlapping variants)
```

---

## Step 3 — Fine-Tune Models

Both models are fine-tuned on splice site classification: donor (label=1), acceptor (label=2), neither (label=0). Training data is 30,000 sequences sampled from `splice_sites_all.csv`. The 10,000 experiment sequences in `splice_sites_10k.csv` are held out and never seen during training.

**NT-100M — full fine-tuning:**
```bash
# Runtime: ~45 mins on A100
python finetune_nt100m.py
# Output: models/checkpoints/nt-100m-splice/
# Expected test accuracy: ~99%
```

**NT-500M — LoRA fine-tuning (r=8):**
```bash
# Runtime: ~2 hours on A100
python finetune_nt500m.py
# Output: models/checkpoints/nt-500m-splice/
# Expected test accuracy: ~99%
```

**Validate both models before proceeding:**
```bash
python validate_models.py --model nt100m    # NT-100M
python validate_nt500m.py                   # NT-500M (LoRA requires different loader)
# Both should show >95% accuracy
```

---

## Step 4 — Baseline Inference

Run each fine-tuned model on the 10,000 unperturbed sequences to record baseline predictions and confidence scores.

```bash
python baseline_inference.py --model nt100m
# Output: results/baseline_nt100m.csv | Runtime: ~15 mins

python baseline_inference.py --model nt500m
# Output: results/baseline_nt500m.csv | Runtime: ~20 mins
```

---

## Step 5 — Perturbation Experiment

For each of the 10,000 sequences, 36 mutations are applied and the model is run on each mutant. This produces 360,000 inference calls per model.

```bash
python run_perturbation.py --model nt100m
# Output: results/perturbation_nt100m.csv | Runtime: 4-6 hrs on A100

python run_perturbation.py --model nt500m
# Output: results/perturbation_nt500m.csv | Runtime: 6-8 hrs on A100
```

**Perturbation types:**

| Type | Positions | Biological meaning | Expected model behavior |
|---|---|---|---|
| Type 1 | 200-201 (donor) / 198-199 (acceptor) | Consensus GT/AG — defines the splice site | Should flip |
| Type 2 | 170-199 | Exonic splicing enhancer region | May flip |
| Type 3 | 210-380 | Deep intronic — biologically silent | Should not flip |

---

## Step 6 — Metrics and Figures

```bash
# Compute FFR, SFR, BSR, AUROC, AUPRC
python analysis/metrics.py
# Output: results/metrics_summary.json

# Generate all 4 paper figures
python figures/plot_all.py
# Output: figures/outputs/fig1_pss_heatmap.png
#         figures/outputs/fig2_ffr_sfr.png
#         figures/outputs/fig3_roc.png
#         figures/outputs/fig4_clinvar_violin.png
```

---

## Metric Definitions

| Metric | Definition |
|---|---|
| **FFR** (Functional Flip Rate) | Fraction of Type 1 mutations that cause a prediction flip. Higher = more sensitive at critical positions. |
| **SFR** (Silent Flip Rate) | Fraction of Type 3 mutations that cause a prediction flip. Lower = more stable on biologically meaningless changes. |
| **BSR** (Biological Specificity Ratio) | FFR / SFR. A BSR of 100 means the model is 100x more likely to flip on a splice-destroying mutation than a silent one. |
| **AUROC** | Area under the ROC curve for ClinVar pathogenicity prediction using confidence change as score. |
| **AUPRC** | Area under the precision-recall curve. Primary metric due to class imbalance in ClinVar. |

---

## Reproducibility Notes

- **Random seeds:** All scripts use `random_state=42`
- **Processed data:** Not committed due to size (134MB+). Regenerate using Steps 2-3 above
- **Model checkpoints:** Not committed due to size. Regenerate using Step 3
- **Total GPU time:** ~12-15 hours across all steps
- **Hardware:** NVIDIA A100 40GB. Results may vary slightly on different hardware

---

## References

- Dalla-Torre et al. (2023). The Nucleotide Transformer. *bioRxiv*
- Landrum et al. (2018). ClinVar. *Nucleic Acids Research*
- Frankish et al. (2023). GENCODE. *Nucleic Acids Research*
- Ribeiro et al. (2020). CheckList: Beyond Accuracy. *ACL 2020*
- Grešová et al. (2023). Genomic Benchmarks. *BMC Genomic Data*
```
