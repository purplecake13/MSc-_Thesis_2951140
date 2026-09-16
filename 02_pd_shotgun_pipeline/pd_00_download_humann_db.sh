#!/bin/bash
#SBATCH --job-name=humann_db_download
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/09_humann3/db_download_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/09_humann3/db_download_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00

module purge
module load bear-apps/2023a
module load humann/4.0.0a1-foss-2023a

DB_DIR=/rds/projects/e/elhamsak-pd-thesis/reference_genomes/humann3_db
mkdir -p "$DB_DIR"

# chocophlan and utility_mapping already downloaded successfully (job 52034379)
# humann_databases --download chocophlan full "$DB_DIR" --update-config yes
# humann_databases --download utility_mapping full "$DB_DIR" --update-config yes

humann_databases --download uniref uniref90_ec_filtered_diamond "$DB_DIR" --update-config yes
humann_databases --download uniref uniref90_diamond "$DB_DIR"

echo "DB download complete: $(date)"
