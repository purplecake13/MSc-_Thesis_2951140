#!/bin/bash
#SBATCH --job-name=pd_merge_genus
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/10_pd_merge/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/10_pd_merge/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

set -euo pipefail

mkdir -p /rds/projects/e/elhamsak-pd-thesis/logs/10_pd_merge

module purge
module load bear-apps/2023a
module load Python/3.11.3-GCCcore-12.3.0

python3 /rds/projects/e/elhamsak-pd-thesis/scripts/pd_04_merge_genus_level.py
