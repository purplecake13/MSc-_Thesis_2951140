#!/bin/bash
#SBATCH --account=elhamsak-pd-thesis
#SBATCH --job-name=cx_02_combine
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/cx/cx_02_combine_platforms_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/cx/cx_02_combine_platforms_%j.err

module purge
module load bear-apps/2023a
module load Miniforge3
module load Python/3.11.3-GCCcore-12.3.0

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /rds/projects/e/elhamsak-pd-thesis/envs/py_analysis

# source activate base  # adjust to whichever env has pandas -- same as cx_01

mkdir -p /rds/projects/e/elhamsak-pd-thesis/logs/cx

# python /rds/projects/e/elhamsak-pd-thesis/scripts/cx_02_combine_platforms.py
/rds/projects/e/elhamsak-pd-thesis/envs/py_analysis/bin/python /rds/projects/e/elhamsak-pd-thesis/scripts/cx_02_combine_platforms.py
