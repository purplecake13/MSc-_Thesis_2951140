#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
# Log paths are set dynamically below to allow for timestamped subdirectories.
# =============================================================================
# 06_fastp.sh — reusable fastp paired-end trimming, SLURM array job
#
# Usage (sbatch flags override SBATCH headers above — no editing needed):
#   sbatch --job-name=fastp_bedarf --array=1-59 \
#     06_fastp.sh \
#       --data-dir /rds/projects/.../raw_data/PD_shotgun/Bedarf_2017 \
#       --output-dir /rds/projects/.../processed/fastp/Bedarf_2017 \
#       --sample-list /rds/projects/.../metadata/bedarf2017_samples.txt \
#       --platform shotgun \
#       --min-len-absolute 45
#       (use --min-len-absolute instead of --read-length when a dataset has
#        heterogeneous native read lengths across runs/batches — e.g. Bedarf,
#        which mixes ~90bp and ~140bp runs. A single read-length-derived
#        fraction cannot be correct for both batches at once.)
#
#   sbatch --job-name=fastp_tran --array=1-176 \
#     06_fastp.sh \
#       --data-dir /rds/projects/.../raw_data/AD_16S/Tran_2019 \
#       --output-dir /rds/projects/.../processed/fastp/Tran_2019 \
#       --sample-list /rds/projects/.../metadata/tran2019_samples.txt \
#       --platform amplicon \
#       --read-length 250 \
#       --trim-primers \
#       --fwd-primer CCTACGGGNGGCWGCAG \
#       --rev-primer GACTACHVGGGTATCTAATCC
#
#   sbatch --job-name=fastp_yildirim2022_sw --array=1-98 \
#     06_fastp.sh \
#       --data-dir /rds/projects/.../ad_raw_data/no_metadata/yildirim2022 \
#       --output-dir /rds/projects/.../ad_trimmed/yildirim2022 \
#       --sample-list /rds/projects/.../yildirim2022_samples.txt \
#       --platform amplicon --trim-primers \
#       --fwd-primer CCTACGGGNGGCWGCAG --rev-primer GACTACHVGGGTATCTAATCC \
#       --sliding-window --cut-window-size 4 --cut-mean-quality 20 \
#       --min-len-absolute 100
#       (--sliding-window enables fastp's --cut_right positional quality
#        trimming, in place of/alongside the length-floor filter. Use this
#        instead of tuning READ_LENGTH/MIN_LEN_FRACTION when reads have
#        heterogeneous quality dropoff rather than heterogeneous native
#        length. Applied uniformly across all AD 16S datasets per
#        supervisor decision, June 2026.)
#
# sample-list = one sample basename per line (no _R1/_R2 suffix, no extension)
# array index N picks line N from --sample-list (1-indexed, matches SLURM_ARRAY_TASK_ID)
# =============================================================================

set -euo pipefail

# ---- defaults -----------------------------------------------------------
PLATFORM="shotgun"          # shotgun | amplicon
READ_LENGTH=150
MIN_LEN_FRACTION=0.75
MIN_LEN_ABSOLUTE=""          # if set, overrides READ_LENGTH x MIN_LEN_FRACTION entirely
TRIM_PRIMERS=false
FWD_PRIMER=""
REV_PRIMER=""
TRIM_POLY_G=false
SLIDING_WINDOW=false
CUT_WINDOW_SIZE=4
CUT_MEAN_QUALITY=20
QUAL_THRESHOLD=15
THREADS=4

# ---- parse args -----------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir) DATA_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --sample-list) SAMPLE_LIST="$2"; shift 2 ;;
    --platform) PLATFORM="$2"; shift 2 ;;
    --read-length) READ_LENGTH="$2"; shift 2 ;;
    --min-len-fraction) MIN_LEN_FRACTION="$2"; shift 2 ;;
    --min-len-absolute) MIN_LEN_ABSOLUTE="$2"; shift 2 ;;
    --trim-primers) TRIM_PRIMERS=true; shift 1 ;;
    --fwd-primer) FWD_PRIMER="$2"; shift 2 ;;
    --rev-primer) REV_PRIMER="$2"; shift 2 ;;
    --trim-poly-g) TRIM_POLY_G=true; shift 1 ;;
    --sliding-window) SLIDING_WINDOW=true; shift 1 ;;
    --cut-window-size) CUT_WINDOW_SIZE="$2"; shift 2 ;;
    --cut-mean-quality) CUT_MEAN_QUALITY="$2"; shift 2 ;;
    --qual-threshold) QUAL_THRESHOLD="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

# ---- validate required args ------------------------------------------------
: "${DATA_DIR:?--data-dir is required}"
: "${OUTPUT_DIR:?--output-dir is required}"
: "${SAMPLE_LIST:?--sample-list is required}"

if [[ "$TRIM_PRIMERS" == true && ( -z "$FWD_PRIMER" || -z "$REV_PRIMER" ) ]]; then
  echo "ERROR: --trim-primers requires --fwd-primer and --rev-primer" >&2
  exit 1
fi

# ---- dynamic log directory setup -------------------------------------------
TIMESTAMP=$(date +"%Y%m%d_%H%M")
LOG_DIR="/rds/projects/e/elhamsak-pd-thesis/logs/06_fastp/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}"
mkdir -p "$OUTPUT_DIR" "$OUTPUT_DIR/reports" "$LOG_DIR"

# Redirect stdout and stderr to per-task log files within the timestamped folder
exec > "${LOG_DIR}/task_${SLURM_ARRAY_TASK_ID}.out" \
     2> "${LOG_DIR}/task_${SLURM_ARRAY_TASK_ID}.err"

# ---- resolve this array task's sample --------------------------------------
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
if [[ -z "$SAMPLE" ]]; then
  echo "ERROR: no sample found at line ${SLURM_ARRAY_TASK_ID} of ${SAMPLE_LIST}" >&2
  exit 1
fi

R1="${DATA_DIR}/${SAMPLE}_1.fastq.gz"
R2="${DATA_DIR}/${SAMPLE}_2.fastq.gz"

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
  echo "ERROR: missing input files for sample ${SAMPLE}" >&2
  echo "  expected: $R1" >&2
  echo "  expected: $R2" >&2
  exit 1
fi

echo "=== fastp: ${SAMPLE} (platform=${PLATFORM}) ==="
echo "Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Log dir: ${LOG_DIR}"
echo "R1: ${R1}"
echo "R2: ${R2}"

# ---- compute min length -----------------------------------------------------
# If --min-len-absolute is set, it takes priority over the read-length-derived
# fraction. Use this for datasets with heterogeneous native read lengths across
# batches/runs, where no single READ_LENGTH x MIN_LEN_FRACTION is valid for all
# reads (e.g. Bedarf 2017's mixed ~90bp/~140bp runs).

if [[ -n "$MIN_LEN_ABSOLUTE" ]]; then
  MIN_LEN="$MIN_LEN_ABSOLUTE"
  echo "Minimum length after trimming: ${MIN_LEN} (absolute override, not derived from read length)"
else
  MIN_LEN=$(awk -v rl="$READ_LENGTH" -v frac="$MIN_LEN_FRACTION" 'BEGIN { printf "%d", rl * frac }')
  echo "Minimum length after trimming: ${MIN_LEN} (= ${MIN_LEN_FRACTION} x ${READ_LENGTH}bp read length)"
fi

# ---- primer trimming step (amplicon only, optional) ------------------------
WORK_R1="$R1"
WORK_R2="$R2"

if [[ "$TRIM_PRIMERS" == true ]]; then
  PRIMER_TRIMMED_DIR="${OUTPUT_DIR}/primer_trimmed"
  mkdir -p "$PRIMER_TRIMMED_DIR"
  PT_R1="${PRIMER_TRIMMED_DIR}/${SAMPLE}_1.primertrimmed.fastq.gz"
  PT_R2="${PRIMER_TRIMMED_DIR}/${SAMPLE}_2.primertrimmed.fastq.gz"

  echo "Loading Cutadapt (2024a)..."
  module purge
  module load bear-apps/2024a
  module load cutadapt/5.1-GCCcore-13.3.0
  
  echo "Trimming primers with cutadapt (fwd=${FWD_PRIMER}, rev=${REV_PRIMER})"
  cutadapt \
    -g "$FWD_PRIMER" \
    -G "$REV_PRIMER" \
    --discard-untrimmed \
    -o "$PT_R1" -p "$PT_R2" \
    "$R1" "$R2" \
    > "${OUTPUT_DIR}/reports/${SAMPLE}_cutadapt.log" 2>&1

  WORK_R1="$PT_R1"
  WORK_R2="$PT_R2"
fi

# ---- build fastp flags -------------------------------------------------------
FASTP_EXTRA_FLAGS=()
if [[ "$TRIM_POLY_G" == true ]]; then
  FASTP_EXTRA_FLAGS+=(--trim_poly_g)
fi

if [[ "$SLIDING_WINDOW" == true ]]; then
  FASTP_EXTRA_FLAGS+=(--cut_right --cut_window_size "$CUT_WINDOW_SIZE" --cut_mean_quality "$CUT_MEAN_QUALITY")
  echo "Sliding-window trimming enabled: window=${CUT_WINDOW_SIZE}bp, mean_qual>=${CUT_MEAN_QUALITY}"
fi

OUT_R1="${OUTPUT_DIR}/${SAMPLE}_1.trimmed.fastq.gz"
OUT_R2="${OUTPUT_DIR}/${SAMPLE}_2.trimmed.fastq.gz"
JSON_REPORT="${OUTPUT_DIR}/reports/${SAMPLE}_fastp.json"
HTML_REPORT="${OUTPUT_DIR}/reports/${SAMPLE}_fastp.html"

echo "Loading fastp (2023a)..."
module purge
module load bear-apps/2023a
module load fastp/0.24.0-GCC-12.3.0

echo "Running fastp..."
fastp \
  -i "$WORK_R1" -I "$WORK_R2" \
  -o "$OUT_R1" -O "$OUT_R2" \
  --detect_adapter_for_pe \
  -q "$QUAL_THRESHOLD" \
  -l "$MIN_LEN" \
  --thread "$THREADS" \
  "${FASTP_EXTRA_FLAGS[@]}" \
  -j "$JSON_REPORT" \
  -h "$HTML_REPORT"

echo "=== Done: ${SAMPLE} ==="
echo "Output: ${OUT_R1}"
echo "Output: ${OUT_R2}"
echo "Report: ${HTML_REPORT}"