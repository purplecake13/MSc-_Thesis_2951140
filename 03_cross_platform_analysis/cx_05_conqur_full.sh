#!/bin/bash
#SBATCH --job-name=cx05_conqur_full
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/cx05_conqur/cx05_conqur_full_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/cx05_conqur/cx05_conqur_full_%j.err
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --account=elhamsak-pd-thesis

set -euo pipefail

module purge
module load bear-apps/2024a
module load R/4.5.0-gfbf-2024a  # ConQuR is already installed under this version

mkdir -p /rds/projects/e/elhamsak-pd-thesis/logs/cx05_conqur

Rscript /rds/projects/e/elhamsak-pd-thesis/scripts/cx_05_conqur_full.R
