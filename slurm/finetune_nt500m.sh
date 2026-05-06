#!/bin/bash
#SBATCH --job-name=finetune_nt500m
#SBATCH --account=bmsc_ga_4551-2026sp
#SBATCH --partition=c12m85-a100-1
#SBATCH --gres=gpu:1
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=slurm/logs/finetune_nt500m.out

source /home/spp9400/miniforge3/etc/profile.d/conda.sh
conda activate /home/spp9400/.conda/envs/genomics

cd /scratch/spp9400/splice-sense
which python
python --version
python -c "import torch; print('torch:', torch.__version__)"
python finetune_nt500m.py
