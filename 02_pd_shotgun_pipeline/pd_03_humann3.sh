#!/bin/bash
#SBATCH --job-name=humann3
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/09_humann3/%x_%A_%a.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/09_humann3/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

set -euo pipefail

# ---- argument parsing ----
HOST_REMOVED_DIR=""
OUTPUT_DIR=""
SAMPLE_LIST=""
TAXPROFILE_DIR=""
NUCLEOTIDE_DB=""
PROTEIN_DB=""
UTILITY_DB=""
R1_SUFFIX="_host_removed_1.fastq.gz"
R2_SUFFIX="_host_removed_2.fastq.gz"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --host-removed-dir) HOST_REMOVED_DIR="$2"; shift 2 ;;
    --output-dir)       OUTPUT_DIR="$2"; shift 2 ;;
    --sample-list)      SAMPLE_LIST="$2"; shift 2 ;;
    --taxprofile-dir)   TAXPROFILE_DIR="$2"; shift 2 ;;
    --nucleotide-db)    NUCLEOTIDE_DB="$2"; shift 2 ;;
    --protein-db)       PROTEIN_DB="$2"; shift 2 ;;
    --utility-db)       UTILITY_DB="$2"; shift 2 ;;
    --r1-suffix)        R1_SUFFIX="$2"; shift 2 ;;
    --r2-suffix)        R2_SUFFIX="$2"; shift 2 ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done
# NOTE: no trailing bare 'shift' here -- this exact bug (redundant shift after
# the case block) silently mis-assigned flags in 06_fastp.sh earlier in this
# project. Confirmed absent here.

for req in HOST_REMOVED_DIR OUTPUT_DIR SAMPLE_LIST TAXPROFILE_DIR NUCLEOTIDE_DB PROTEIN_DB UTILITY_DB; do
  if [ -z "${!req}" ]; then
    echo "ERROR: missing required argument --${req,,}"
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
if [ -z "$SAMPLE" ]; then
  echo "ERROR: no sample found at line ${SLURM_ARRAY_TASK_ID} of ${SAMPLE_LIST}"
  exit 1
fi

R1="${HOST_REMOVED_DIR}/${SAMPLE}${R1_SUFFIX}"
R2="${HOST_REMOVED_DIR}/${SAMPLE}${R2_SUFFIX}"
TAXPROFILE="${TAXPROFILE_DIR}/${SAMPLE}_metaphlan_profile.txt"

for f in "$R1" "$R2" "$TAXPROFILE"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: expected input file not found: $f"
    exit 1
  fi
done

# Final output check -- HUMAnN doesn't fail loudly on a pre-existing output dir
# the way MetaPhlAn does on a stale .bz2, it may just skip/resume oddly.
# Force a clean slate for this sample explicitly rather than trusting resume behaviour.
if [ -f "${OUTPUT_DIR}/${SAMPLE}_concat_2_genefamilies.tsv" ]; then
  echo "Final output already exists for ${SAMPLE} -- skipping. Delete it manually to force a re-run."
  exit 0
fi
rm -rf "${OUTPUT_DIR}/${SAMPLE}_concat_humann_temp"

module purge
module load bear-apps/2023a
module load humann/4.0.0a1-foss-2023a
module load Bowtie2/2.5.4-GCC-12.3.0

# Build concatenated input just-in-time in a scratch location, run, then delete.
# This avoids ever holding a second permanent full-size copy of the dataset.
SCRATCH_DIR="${OUTPUT_DIR}/_scratch"
mkdir -p "$SCRATCH_DIR"
CONCAT="${SCRATCH_DIR}/${SAMPLE}_concat.fastq.gz"

cat "$R1" "$R2" > "$CONCAT"

# gzip integrity check on the concat file before committing compute time to it
if ! gzip -t "$CONCAT" 2>/dev/null; then
  echo "ERROR: concatenated input failed gzip integrity check: $CONCAT"
  rm -f "$CONCAT"
  exit 1
fi

set +e
humann \
  --input "$CONCAT" \
  --output "$OUTPUT_DIR" \
  --output-basename "${SAMPLE}_concat" \
  --taxonomic-profile "$TAXPROFILE" \
  --nucleotide-database "$NUCLEOTIDE_DB" \
  --protein-database "$PROTEIN_DB" \
  --utility-database "$UTILITY_DB" \
  --threads "${SLURM_CPUS_PER_TASK}"

HUMANN_EXIT=$?
set -e

# Clean up the temp concat file and humann's own intermediate temp dir regardless of outcome
rm -f "$CONCAT"
rm -rf "${OUTPUT_DIR}/${SAMPLE}_concat_humann_temp"

if [ $HUMANN_EXIT -ne 0 ]; then
  echo "ERROR: humann exited non-zero for ${SAMPLE}"
  exit $HUMANN_EXIT
fi

echo "Completed ${SAMPLE}"
