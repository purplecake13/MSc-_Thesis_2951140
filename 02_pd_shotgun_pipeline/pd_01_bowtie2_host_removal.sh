#!/bin/bash
#SBATCH --job-name=bowtie2_host_removal
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
# Log paths are set dynamically below to allow for timestamped subdirectories.
# =============================================================================
# 07_bowtie2_host_removal.sh — reusable Bowtie2 human-read removal, SLURM array job
#
# Aligns fastp-trimmed paired-end shotgun reads against a host (human) reference
# index and retains only reads that do NOT align (i.e. non-host / microbial reads).
#
# Usage (sbatch flags override SBATCH headers above — no editing needed):
#   sbatch --job-name=bt2_bedarf --array=1-59 \
#     07_bowtie2_host_removal.sh \
#       --data-dir /rds/projects/.../pd_trimmed/bedarf2017 \
#       --output-dir /rds/projects/.../pd_host_removed/bedarf2017 \
#       --sample-list /rds/projects/.../pd_raw_data/bedarf2017/bedarf2017_samples.txt \
#       --index /rds/projects/e/elhamsak-ad-thesis/reference_genomes/bowtie2_hg38/GRCh38_noalt_as
#
# sample-list = one sample basename per line (same convention as 06_fastp.sh)
# array index N picks line N from --sample-list (1-indexed, matches SLURM_ARRAY_TASK_ID)
#
# Expects input files named: <data-dir>/<sample>_1.trimmed.fastq.gz
#                             <data-dir>/<sample>_2.trimmed.fastq.gz
# (i.e. the default output naming from 06_fastp.sh)
# =============================================================================

set -euo pipefail

# ---- defaults -----------------------------------------------------------
THREADS=8
INPUT_SUFFIX_R1="_1.trimmed.fastq.gz"
INPUT_SUFFIX_R2="_2.trimmed.fastq.gz"

# ---- parse args -----------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir) DATA_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --sample-list) SAMPLE_LIST="$2"; shift 2 ;;
    --index) BT2_INDEX="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --r1-suffix) INPUT_SUFFIX_R1="$2"; shift 2 ;;
    --r2-suffix) INPUT_SUFFIX_R2="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

# ---- validate required args ------------------------------------------------
: "${DATA_DIR:?--data-dir is required}"
: "${OUTPUT_DIR:?--output-dir is required}"
: "${SAMPLE_LIST:?--sample-list is required}"
: "${BT2_INDEX:?--index is required (path prefix to .bt2 index files, without extension)}"

# Sanity check the index actually exists before launching the array
if [[ ! -f "${BT2_INDEX}.1.bt2" ]]; then
  echo "ERROR: Bowtie2 index not found at prefix: ${BT2_INDEX}" >&2
  echo "  expected: ${BT2_INDEX}.1.bt2 (and siblings)" >&2
  exit 1
fi

# ---- dynamic log directory setup -------------------------------------------
TIMESTAMP=$(date +"%Y%m%d_%H%M")
LOG_DIR="/rds/projects/e/elhamsak-pd-thesis/logs/07_bowtie2/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}"
mkdir -p "$OUTPUT_DIR" "$OUTPUT_DIR/logs" "$LOG_DIR"

# Redirect stdout and stderr to per-task log files within the timestamped folder
exec > "${LOG_DIR}/task_${SLURM_ARRAY_TASK_ID}.out" \
     2> "${LOG_DIR}/task_${SLURM_ARRAY_TASK_ID}.err"

# ---- resolve this array task's sample --------------------------------------
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
if [[ -z "$SAMPLE" ]]; then
  echo "ERROR: no sample found at line ${SLURM_ARRAY_TASK_ID} of ${SAMPLE_LIST}" >&2
  exit 1
fi

R1="${DATA_DIR}/${SAMPLE}${INPUT_SUFFIX_R1}"
R2="${DATA_DIR}/${SAMPLE}${INPUT_SUFFIX_R2}"

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
  echo "ERROR: missing input files for sample ${SAMPLE}" >&2
  echo "  expected: $R1" >&2
  echo "  expected: $R2" >&2
  exit 1
fi

echo "=== Bowtie2 host removal: ${SAMPLE} ==="
echo "Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Log dir: ${LOG_DIR}"
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Index: ${BT2_INDEX}"

# ---- module setup -----------------------------------------------------------
echo "Loading Bowtie2 (2024a)..."
module purge
module load bear-apps/2024a
module load Bowtie2/2.5.4-GCC-13.3.0

# ---- output paths -----------------------------------------------------------
# --un-conc-gz writes reads where NEITHER mate aligned concordantly to the host
# genome — these are the retained, non-host (microbial) reads we want to keep.
# %-> Bowtie2 substitutes 1/2 for R1/R2 in the output filenames.
UNALIGNED_PATTERN="${OUTPUT_DIR}/${SAMPLE}_host_removed_%.fastq.gz"
SAM_DISCARD="/dev/null"   # we don't need the alignment itself, only the unaligned reads
SUMMARY_FILE="${OUTPUT_DIR}/${SAMPLE}_bowtie2_host_removed_summary.txt"

echo "Running Bowtie2..."
bowtie2 \
  -x "$BT2_INDEX" \
  -1 "$R1" -2 "$R2" \
  --un-conc-gz "$UNALIGNED_PATTERN" \
  --no-unal \
  -S "$SAM_DISCARD" \
  -p "$THREADS" \
  2> "$SUMMARY_FILE"

echo "--- Bowtie2 alignment summary ---"
cat "$SUMMARY_FILE"
echo "----------------------------------"

# Bowtie2 with --un-conc-gz and a pattern ending in "_%.fastq.gz" produces
# files literally named "..._host_removed_1.fastq.gz" / "..._host_removed_2.fastq.gz"
FINAL_R1="${OUTPUT_DIR}/${SAMPLE}_host_removed_1.fastq.gz"
FINAL_R2="${OUTPUT_DIR}/${SAMPLE}_host_removed_2.fastq.gz"

if [[ ! -f "$FINAL_R1" || ! -f "$FINAL_R2" ]]; then
  echo "ERROR: expected host_removed output files not found for ${SAMPLE}" >&2
  echo "  expected: $FINAL_R1" >&2
  echo "  expected: $FINAL_R2" >&2
  exit 1
fi

echo "=== Done: ${SAMPLE} ==="
echo "Output: ${FINAL_R1}"
echo "Output: ${FINAL_R2}"
echo "Summary: ${SUMMARY_FILE}"

# ---- sanity flag on host removal rate ---------------------------------------
# Typical expectation for gut/stool shotgun metagenomics: 1-10% of reads align
# to host. Parse the overall alignment rate from the summary for a quick flag.
OVERALL_RATE=$(grep "overall alignment rate" "$SUMMARY_FILE" | grep -oE '[0-9]+\.[0-9]+%' | head -1)
echo "Host alignment rate (reads removed): ${OVERALL_RATE:-unknown}"
echo "NOTE: expected range for gut/stool samples is ~1-10%. Values well outside"
echo "this range may indicate contamination, mislabelled sample type, or an"
echo "index/species mismatch — flag for manual review if so."