#!/bin/bash
#SBATCH --job-name=data_prep_v2
#SBATCH --account=bmsc_ga_4551-2026sp
#SBATCH --partition=c12m85-a100-1
#SBATCH --gres=gpu:1
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=slurm/logs/data_prep_v2.out

source /home/spp9400/miniforge3/etc/profile.d/conda.sh
conda activate splice
cd /scratch/spp9400/splice-sense
python data_prep_v2.py
