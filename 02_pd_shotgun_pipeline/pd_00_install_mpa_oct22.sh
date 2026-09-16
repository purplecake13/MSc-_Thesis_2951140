#!/bin/bash
#SBATCH --job-name=mpa_oct22_install
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/08_metaphlan4/oct22_install_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/08_metaphlan4/oct22_install_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00

module purge
module load bear-apps/2023a
module load MetaPhlAn/4.1.1-foss-2023a

metaphlan --install --index mpa_vOct22_CHOCOPhlAnSGB_202403 --bowtie2db /rds/projects/e/elhamsak-ad-thesis/reference_genomes/MetaPhlAn4

echo "Oct22 install complete: $(date)"
