#!/bin/bash
#SBATCH --job-name=dada2_denoise
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/08_ad_qiime2/dada2_denoise/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/08_ad_qiime2/dada2_denoise/%x_%j.err
# =============================================================================
# /rds/projects/e/elhamsak-pd-thesis/scripts/ad_04_dada2_denoise.sh — reusable QIIME2 DADA2 denoise wrapper (paired or single)
#
# Runs qiime dada2 denoise-paired OR denoise-single on a demultiplexed QIIME2
# artifact and produces a feature table, representative sequences, and
# denoising stats.
#
# NOT a SLURM array job (DADA2 processes all samples in the demux artifact
# as one unit) — submit as a single job per dataset.
#
# --mode paired (default): standard workflow, e.g. Cirstea/Yamashiro-style
#   datasets where DADA2's own paired-end merging performs adequately.
#   Requires --trunc-len-f and --trunc-len-r.
#
# --mode single: for datasets processed through fastp merge mode first
#   (11_fastp_merge_array.sh) then imported via 12_import_merged_reads.sh.
#   Reads are already merged single sequences, not R1/R2 pairs.
#   Requires --trunc-len (one value, not F/R).
#
# Usage (paired, unchanged from before):
#   sbatch --job-name=dada2_cirstea2022 \
#     09_dada2_denoise.sh \
#       --mode paired \
#       --demux-qza /rds/.../cirstea2022/cirstea2022_demux.qza \
#       --output-dir /rds/.../cirstea2022 \
#       --trunc-len-f 220 --trunc-len-r 190
#
# Usage (single, for fastp-merged datasets):
#   sbatch --job-name=dada2_tran2019 \
#     09_dada2_denoise.sh \
#       --mode single \
#       --demux-qza /rds/.../tran2019/merged_demux.qza \
#       --output-dir /rds/.../tran2019 \
#       --trunc-len 0
#
# Optional --max-ee-f / --max-ee-r (paired mode) or --max-ee (single mode);
# DADA2 default is 2.0 if not supplied. Use --trunc-len-f 0 --trunc-len-r 0
# (paired) or --trunc-len 0 (single) to disable hard truncation entirely and
# rely on expected-errors filtering instead — useful when truncation-length
# tuning hits an overlap ceiling.
#
# Truncation lengths are dataset-specific — always derive from that dataset's
# own demux.qzv quality plot and expected amplicon length, never reuse
# another dataset's values.
# =============================================================================

set -euo pipefail

# ---- defaults -----------------------------------------------------------
THREADS=8
OUTPUT_PREFIX=""   # optional prefix for output filenames, e.g. "tran2019_"
MODE="paired"      # paired (default, preserves prior behavior) or single
MAX_EE_F=2.0
MAX_EE_R=2.0
MAX_EE=2.0

# ---- parse args -----------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --demux-qza) DEMUX_QZA="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --trunc-len-f) TRUNC_LEN_F="$2"; shift 2 ;;
    --trunc-len-r) TRUNC_LEN_R="$2"; shift 2 ;;
    --trunc-len) TRUNC_LEN="$2"; shift 2 ;;
    --max-ee-f) MAX_EE_F="$2"; shift 2 ;;
    --max-ee-r) MAX_EE_R="$2"; shift 2 ;;
    --max-ee) MAX_EE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --output-prefix) OUTPUT_PREFIX="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

if [[ "${MODE}" != "paired" && "${MODE}" != "single" ]]; then
  echo "ERROR: --mode must be 'paired' or 'single', got: ${MODE}" >&2
  exit 1
fi

# ---- validate required args ------------------------------------------------
: "${DEMUX_QZA:?--demux-qza is required}"
: "${OUTPUT_DIR:?--output-dir is required}"

if [[ "${MODE}" == "paired" ]]; then
  : "${TRUNC_LEN_F:?--trunc-len-f is required for --mode paired}"
  : "${TRUNC_LEN_R:?--trunc-len-r is required for --mode paired}"
else
  : "${TRUNC_LEN:?--trunc-len is required for --mode single}"
fi

if [[ ! -f "$DEMUX_QZA" ]]; then
  echo "ERROR: demux artifact not found: $DEMUX_QZA" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

TABLE_QZA="${OUTPUT_DIR}/${OUTPUT_PREFIX}table.qza"
REPSEQS_QZA="${OUTPUT_DIR}/${OUTPUT_PREFIX}repseqs.qza"
STATS_QZA="${OUTPUT_DIR}/${OUTPUT_PREFIX}dada2stats.qza"
STATS_QZV="${OUTPUT_DIR}/${OUTPUT_PREFIX}dada2stats.qzv"

echo "=== DADA2 denoise (mode: ${MODE}) ==="
echo "Demux input:    $DEMUX_QZA"
if [[ "${MODE}" == "paired" ]]; then
  echo "Trunc len F:    $TRUNC_LEN_F"
  echo "Trunc len R:    $TRUNC_LEN_R"
  echo "Combined:       $((TRUNC_LEN_F + TRUNC_LEN_R)) bp"
  echo "Max EE F:       $MAX_EE_F"
  echo "Max EE R:       $MAX_EE_R"
else
  echo "Trunc len:      $TRUNC_LEN"
  echo "Max EE:         $MAX_EE"
fi
echo "Threads:        $THREADS"
echo "Output dir:     $OUTPUT_DIR"
echo "  table:        $TABLE_QZA"
echo "  rep-seqs:     $REPSEQS_QZA"
echo "  stats:        $STATS_QZA"
echo ""

# ---- module setup -----------------------------------------------------------
echo "Loading QIIME2..."
module purge
module load bear-apps/2023a
module load QIIME2/2025.4

export TMPDIR="/rds/projects/e/elhamsak-pd-thesis/tmp_qiime/${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"
echo "TMPDIR set to: ${TMPDIR}"

# ---- run DADA2 ---------------------------------------------------------------
if [[ "${MODE}" == "paired" ]]; then
  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "$DEMUX_QZA" \
    --p-trunc-len-f "$TRUNC_LEN_F" \
    --p-trunc-len-r "$TRUNC_LEN_R" \
    --p-max-ee-f "$MAX_EE_F" \
    --p-max-ee-r "$MAX_EE_R" \
    --p-n-threads "$THREADS" \
    --o-table "$TABLE_QZA" \
    --o-representative-sequences "$REPSEQS_QZA" \
    --o-denoising-stats "$STATS_QZA"
else
  qiime dada2 denoise-single \
    --i-demultiplexed-seqs "$DEMUX_QZA" \
    --p-trunc-len "$TRUNC_LEN" \
    --p-max-ee "$MAX_EE" \
    --p-n-threads "$THREADS" \
    --o-table "$TABLE_QZA" \
    --o-representative-sequences "$REPSEQS_QZA" \
    --o-denoising-stats "$STATS_QZA"
fi

echo ""
echo "=== DADA2 complete ==="

if [[ ! -s "${TABLE_QZA}" ]]; then
  echo "ERROR: DADA2 produced no output or empty feature table at ${TABLE_QZA}" >&2
  exit 1
fi

# ---- auto-generate stats visualisation for immediate inspection -------------
echo "Generating denoising stats visualisation..."
qiime metadata tabulate \
  --m-input-file "$STATS_QZA" \
  --o-visualization "$STATS_QZV"

echo "Stats visualisation: $STATS_QZV"
echo ""
echo "NOTE: view ${STATS_QZV} at https://view.qiime2.org to check"
echo "'percentage of input merged' (paired mode) or 'percentage of input"
echo "passed filter' (single mode) per sample. Expected healthy range:"
echo ">70% retained, <5% chimeric (per project benchmarks). If retention is"
echo "low, truncation length likely does not suit this dataset — revisit"
echo "--trunc-len(-f/-r), do not assume the run itself failed."

rm -rf "${TMPDIR:?}"/*