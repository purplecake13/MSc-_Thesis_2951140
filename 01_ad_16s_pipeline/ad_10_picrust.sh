#!/bin/bash
#SBATCH --job-name=picrust2_ad
#SBATCH --array=0-7%4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/picrust2/picrust2_ad_%A_%a.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/picrust2/picrust2_ad_%A_%a.err

module purge
module load bear-apps/2023a
module load Miniforge3
source activate picrust2

STUDIES=(binyinli2019 cirstea2022 ling2021 liu2019 ueda2021 yamashiro2024 yildirim2022 zhuang2018)
STUDY=${STUDIES[$SLURM_ARRAY_TASK_ID]}
BASE=/rds/projects/e/elhamsak-ad-thesis/ad_qiime2/$STUDY
OUTDIR=/rds/projects/e/elhamsak-ad-thesis/ad_picrust2/$STUDY

# clear output dir manually before rerunning, or the pipeline will error out.
if [ -d "$OUTDIR" ]; then
  echo "Removing existing output dir: $OUTDIR"
  rm -rf "$OUTDIR"
fi

picrust2_pipeline.py \
  -s "$BASE/picrust2_input/repseqs_export/dna-sequences.fasta" \
  -i "$BASE/picrust2_input/table_export/feature-table.biom" \
  -o "/rds/projects/e/elhamsak-ad-thesis/picrust2_ad_output/$STUDY" \
  -p "$SLURM_CPUS_PER_TASK" \
  --stratified \
  --verbose
  

echo "PICRUSt2 completed for $STUDY: $(date)"
