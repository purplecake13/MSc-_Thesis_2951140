#!/bin/bash
#SBATCH --job-name=metaphlan4_run
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
# NOTE: log paths set dynamically below, same pattern as 07_bowtie2_host_removal.sh
# =============================================================================
# 08_metaphlan4.sh — reusable MetaPhlAn4 taxonomic profiling, SLURM array job
#
# Runs MetaPhlAn4 directly on paired-end FASTQ input. Designed for two cases:
#   (a) Datasets confirmed ALREADY host-depleted at deposition (e.g. Bedarf,
#       Boktor RUMC, Boktor TBC) -- point --data-dir at fastp-trimmed reads
#       directly, skip Bowtie2 entirely.
#   (b) Datasets requiring genuine host removal (e.g. Clasen) -- point
#       --data-dir at Bowtie2 host-removed output instead, once that step
#       is complete.
#
# Usage (sbatch flags override SBATCH headers above — no editing needed):
#   sbatch --job-name=mpa4_boktor_tbc --array=1-41 \
#     08_metaphlan4.sh \
#       --data-dir /rds/projects/.../pd_trimmed/boktor2023/tbc \
#       --output-dir /rds/projects/.../pd_metaphlan4/boktor2023/tbc \
#       --sample-list /rds/projects/.../pd_raw_data/boktor2023/tbc/boktor2023tbc_samples.txt \
#       --bowtie2db /rds/projects/e/elhamsak-ad-thesis/reference_genomes/MetaPhlAn4 \
#       --mpa-index mpa_vJan25_CHOCOPhlAnSGB_202503
#
# sample-list = one sample basename per line (same convention as 06_fastp.sh / 07_bowtie2_host_removal.sh)
# array index N picks line N from --sample-list (1-indexed, matches SLURM_ARRAY_TASK_ID)
#
# Default expects input files named: <data-dir>/<sample>_1.trimmed.fastq.gz
#                                     <data-dir>/<sample>_2.trimmed.fastq.gz
# (i.e. fastp output naming). Override --r1-suffix/--r2-suffix for host-removed
# input instead, e.g. --r1-suffix _host_removed_1.fastq.gz
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
    --bowtie2db) BOWTIE2DB="$2"; shift 2 ;;
    --mpa-index) MPA_INDEX="$2"; shift 2 ;;
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
: "${BOWTIE2DB:?--bowtie2db is required (path to MetaPhlAn4 bowtie2db directory)}"
: "${MPA_INDEX:?--mpa-index is required (e.g. mpa_vJan25_CHOCOPhlAnSGB_202503)}"

if [[ ! -d "${BOWTIE2DB}" ]]; then
  echo "ERROR: MetaPhlAn4 database directory not found: ${BOWTIE2DB}" >&2
  exit 1
fi

if [[ ! -f "${SAMPLE_LIST}" ]]; then
  echo "ERROR: sample list not found: ${SAMPLE_LIST}" >&2
  exit 1
fi

# ---- dynamic log directory setup -------------------------------------------
LOG_DIR="/rds/projects/e/elhamsak-pd-thesis/logs/08_metaphlan4/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}"
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

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

echo "=== MetaPhlAn4: ${SAMPLE} ==="
echo "Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Log dir: ${LOG_DIR}"
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Bowtie2DB: ${BOWTIE2DB}"
echo "MPA index: ${MPA_INDEX}"

# ---- module setup -----------------------------------------------------------
echo "Loading MetaPhlAn4..."
module purge
module load bear-apps/2023a
module load MetaPhlAn/4.1.1-foss-2023a

# ---- output paths -----------------------------------------------------------
BOWTIE2OUT="${OUTPUT_DIR}/${SAMPLE}_bowtie2.bz2"
PROFILE_OUT="${OUTPUT_DIR}/${SAMPLE}_metaphlan_profile.txt"

echo "Running MetaPhlAn4..."

metaphlan \
  "${R1},${R2}" \
  --bowtie2db "${BOWTIE2DB}" \
  --index "${MPA_INDEX}" \
  --input_type fastq \
  -t rel_ab_w_read_stats \
  --nproc "${THREADS}" \
  --bowtie2out "${BOWTIE2OUT}" \
  -o "${PROFILE_OUT}"
  
# ---- verify output actually has content, not just an empty file -----------
if [[ ! -s "${PROFILE_OUT}" ]]; then
  echo "ERROR: MetaPhlAn4 profile output missing or empty for ${SAMPLE}" >&2
  echo "  expected: ${PROFILE_OUT}" >&2
  exit 1
fi

# Quick sanity check: a real profile should have more than just header/comment lines
N_TAXA_LINES=$(grep -vc "^#" "${PROFILE_OUT}" || true)
if [[ "${N_TAXA_LINES}" -eq 0 ]]; then
  echo "ERROR: MetaPhlAn4 profile for ${SAMPLE} contains no taxon rows (only comments/header)" >&2
  echo "  this likely means classification failed silently -- check ${BOWTIE2OUT} and rerun" >&2
  exit 1
fi

echo "=== Done: ${SAMPLE} ==="
echo "Profile: ${PROFILE_OUT} (${N_TAXA_LINES} taxon rows)"
echo "Bowtie2 intermediate: ${BOWTIE2OUT}"