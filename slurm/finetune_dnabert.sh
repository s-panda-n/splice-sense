#!/bin/bash
#SBATCH --job-name=finetune_dnabert
#SBATCH --account=bmsc_ga_4551-2026sp
#SBATCH --partition=c12m85-a100-1
#SBATCH --gres=gpu:1
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#SBATCH --output=slurm/logs/finetune_dnabert.out

source /home/spp9400/miniforge3/etc/profile.d/conda.sh
conda activate /home/spp9400/miniforge3/envs/splice

cd /scratch/spp9400/splice-sense
which python
python --version
python -c "import torch; print('torch version:', torch.__version__)"
python finetune_dnabert.py
