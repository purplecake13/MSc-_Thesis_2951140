#!/bin/bash
#SBATCH --job-name=fastp_merge
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/08_ad_qiime2/fastp_merge/%x_%A_%a.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/08_ad_qiime2/fastp_merge/%x_%A_%a.err
# =============================================================================
# /rds/projects/e/elhamsak-pd-thesis/scripts/ad_02a_fastp_merge.sh — fastp merge mode (-m), per-sample array job
#
# Third candidate merging approach for datasets where DADA2's overlap-based
# merging underperforms (e.g. Tran2019: V3-V4 amplicon, insufficient R1/R2
# overlap, DADA2 merge capped ~45-55%). fastp's merge mode is STILL an
# overlap-based approach (same family as DADA2/VSEARCH, not architecturally
# different like concatenation) — it may or may not outperform either,
# this is an empirical comparison, not a guaranteed fix.
#
# Runs as a SLURM array job, one task per sample, using a sample list file
# (one SRR accession per line) to index into the array.
#
# Usage:
#   sbatch --array=1-56 \
#     --job-name=fastp_merge_tran2019 \
#     11_fastp_merge_array.sh \
#       --sample-list /rds/projects/e/elhamsak-ad-thesis/ad_raw_data/tran2019/tran2019_samples.txt \
#       --data-dir /rds/projects/e/elhamsak-ad-thesis/ad_trimmed/tran2019 \
#       --output-dir /rds/projects/e/elhamsak-ad-thesis/ad_qiime2/tran2019/fastp_merge
#
# --array=1-56 MUST match the line count of --sample-list exactly (56 for
# the corrected human-only Tran2019 sample list). Check with:
#   wc -l tran2019_samples.txt
# before submitting — a mismatched array range will either skip samples
# or fail on out-of-range line lookups.
#
# After all array tasks complete, aggregate per-sample JSON reports with
# 11b_fastp_merge_summary.py (or MultiQC) to get one merge-rate table
# comparable against the DADA2 and VSEARCH-merge baselines.
# =============================================================================

set -euo pipefail

# ---- ensure log directory exists (in case it wasn't created before sbatch) --
mkdir -p /rds/projects/e/elhamsak-pd-thesis/logs/08_ad_qiime2/fastp_merge

# ---- defaults -----------------------------------------------------------
R1_SUFFIX="_1.trimmed.fastq.gz"
R2_SUFFIX="_2.trimmed.fastq.gz"
THREADS=4

# ---- parse args -----------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-list) SAMPLE_LIST="$2"; shift 2 ;;
    --data-dir) DATA_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --r1-suffix) R1_SUFFIX="$2"; shift 2 ;;
    --r2-suffix) R2_SUFFIX="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

# ---- validate required args ------------------------------------------------
: "${SAMPLE_LIST:?--sample-list is required}"
: "${DATA_DIR:?--data-dir is required}"
: "${OUTPUT_DIR:?--output-dir is required}"
: "${SLURM_ARRAY_TASK_ID:?This script must be run as a SLURM array job (sbatch --array=...)}"

if [[ ! -f "$SAMPLE_LIST" ]]; then
  echo "ERROR: sample list not found: $SAMPLE_LIST" >&2
  exit 1
fi

# ---- look up this array task's sample -------------------------------------
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

if [[ -z "$SAMPLE" ]]; then
  echo "ERROR: no sample found at line ${SLURM_ARRAY_TASK_ID} of ${SAMPLE_LIST}" >&2
  echo "Check that --array range matches the sample list line count (wc -l ${SAMPLE_LIST})" >&2
  exit 1
fi

R1="${DATA_DIR}/${SAMPLE}${R1_SUFFIX}"
R2="${DATA_DIR}/${SAMPLE}${R2_SUFFIX}"

if [[ ! -f "$R1" ]]; then
  echo "ERROR: R1 file not found: $R1" >&2
  exit 1
fi
if [[ ! -f "$R2" ]]; then
  echo "ERROR: R2 file not found: $R2" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

MERGED_OUT="${OUTPUT_DIR}/${SAMPLE}_merged.fastq.gz"
UNMERGED_R1="${OUTPUT_DIR}/${SAMPLE}_unmerged_R1.fastq.gz"
UNMERGED_R2="${OUTPUT_DIR}/${SAMPLE}_unmerged_R2.fastq.gz"
JSON_REPORT="${OUTPUT_DIR}/${SAMPLE}_fastp_merge.json"
HTML_REPORT="${OUTPUT_DIR}/${SAMPLE}_fastp_merge.html"

echo "=== fastp merge mode ==="
echo "Array task:   ${SLURM_ARRAY_TASK_ID}"
echo "Sample:       $SAMPLE"
echo "R1 input:     $R1"
echo "R2 input:     $R2"
echo "Threads:      $THREADS"
echo "Output dir:   $OUTPUT_DIR"
echo "  merged:     $MERGED_OUT"
echo "  unmerged:   $UNMERGED_R1 / $UNMERGED_R2"
echo "  report:     $JSON_REPORT"
echo ""

# ---- module setup -----------------------------------------------------------
echo "Loading fastp..."
module purge
module load bear-apps/2023a
module load fastp/0.24.0-GCC-12.3.0

# ---- run fastp merge ---------------------------------------------------------
echo "=== Running fastp -m (merge mode) ==="
fastp \
  -i "$R1" \
  -I "$R2" \
  -m \
  --merged_out "$MERGED_OUT" \
  --out1 "$UNMERGED_R1" \
  --out2 "$UNMERGED_R2" \
  -j "$JSON_REPORT" \
  -h "$HTML_REPORT" \
  -w "$THREADS"

echo ""
echo "=== fastp merge complete for ${SAMPLE} ==="
echo "Merged reads written to:   $MERGED_OUT"
echo "JSON report:               $JSON_REPORT"
echo ""
echo "NOTE: merge rate for this sample = merged_reads / total_reads, found in"
echo "the JSON report under read1_after_filtering / merging summary fields."
echo "After all array tasks complete, aggregate across samples (e.g. with"
echo "MultiQC, which natively parses fastp JSON, or a small summary script)"
echo "to get one comparable merge-rate table against the DADA2 and"
echo "VSEARCH-merge baselines."