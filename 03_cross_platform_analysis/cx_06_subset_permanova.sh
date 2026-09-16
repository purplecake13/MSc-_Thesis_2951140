#!/bin/bash
#SBATCH --job-name=cx06_subset_permanova
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/cx06_subset_permanova_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/cx06_subset_permanova_%j.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --account=elhamsak-pd-thesis

set -euo pipefail

module purge
module load bear-apps/2024a
module load R/4.5.0-gfbf-2024a

mkdir -p /rds/projects/e/elhamsak-pd-thesis/logs

Rscript /rds/projects/e/elhamsak-pd-thesis/scripts/cx_06_subset_permanova.R
