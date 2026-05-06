#!/bin/bash
#SBATCH --job-name=baseline_nt100m
#SBATCH --account=bmsc_ga_4551-2026sp
#SBATCH --partition=c12m85-a100-1
#SBATCH --gres=gpu:1
#SBATCH --mem=60GB
#SBATCH --time=01:00:00
#SBATCH --output=slurm/logs/baseline_nt100m.out

source /home/spp9400/miniforge3/etc/profile.d/conda.sh
conda activate /home/spp9400/.conda/envs/genomics
cd /scratch/spp9400/splice-sense
python baseline_inference.py --model nt100m
