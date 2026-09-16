#!/bin/bash
#SBATCH --job-name=filter_normalize
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/cx_01/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/cx_01/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

set -euo pipefail

mkdir -p /rds/projects/e/elhamsak-pd-thesis/logs/10_merge

module purge
module load bear-apps/2023a
module load Python/3.11.3-GCCcore-12.3.0
which python; python -c "import sys; print(sys.executable)"

python3 /rds/projects/e/elhamsak-pd-thesis/scripts/cx_01_filter_normalize.py
