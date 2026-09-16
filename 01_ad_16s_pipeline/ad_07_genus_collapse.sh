#!/bin/bash
#SBATCH --job-name=genus_collapse
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/taxonomy/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/taxonomy/%x_%j.err
# =============================================================================
# /rds/projects/e/elhamsak-pd-thesis/scripts/ad_07_genus_collapse.sh — collapse a QIIME2 feature table to genus level (L6)
#
# Takes a DADA2 feature table + classify-sklearn taxonomy output and collapses
# to genus level (qiime taxa collapse --p-level 6), per the schedule doc's
# established cross-platform harmonisation approach. Exports a flat TSV
# (genera x samples) ready for cross-study merging.
#
# Usage:
#   sbatch --job-name=genus_collapse_cirstea 14_genus_collapse.sh \
#     --feature-table /rds/projects/.../cirstea2022/dada2_table.qza \
#     --taxonomy /rds/projects/.../cirstea2022/taxonomy/taxonomy.qza \
#     --output-dir /rds/projects/.../cirstea2022/genus_level
#
# Inputs expected:
#   --feature-table   FeatureTable[Frequency] .qza (from DADA2 denoise-paired/single)
#   --taxonomy        FeatureData[Taxonomy] .qza (from 13_taxonomy_assignment.sh)
#
# Output:
#   <output-dir>/genus_table.qza         (FeatureTable[Frequency], collapsed to L6)
#   <output-dir>/genus_table_export/feature-table.tsv  (flat TSV, genera x samples)
#
# NOTE: this script does NOT apply rare-taxon filtering or CLR transformation --
# those are separate, deliberate steps applied later during cross-study merging
# (per schedule doc Week 5 plan), so each dataset's raw genus-level counts are
# preserved here for inspection before any filtering decisions are made.
# =============================================================================

set -euo pipefail

# ---- parse args -----------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --feature-table) FEATURE_TABLE="$2"; shift 2 ;;
    --taxonomy) TAXONOMY="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --level) COLLAPSE_LEVEL="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

COLLAPSE_LEVEL="${COLLAPSE_LEVEL:-6}"  # 6 = genus, per established convention

# ---- validate required args ------------------------------------------------
: "${FEATURE_TABLE:?--feature-table is required}"
: "${TAXONOMY:?--taxonomy is required}"
: "${OUTPUT_DIR:?--output-dir is required}"

if [[ ! -f "${FEATURE_TABLE}" ]]; then
  echo "ERROR: feature table artifact not found: ${FEATURE_TABLE}" >&2
  exit 1
fi
if [[ ! -f "${TAXONOMY}" ]]; then
  echo "ERROR: taxonomy artifact not found: ${TAXONOMY}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

echo "=== Genus-level collapse (L${COLLAPSE_LEVEL}) ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "Feature table: ${FEATURE_TABLE}"
echo "Taxonomy: ${TAXONOMY}"
echo "Output dir: ${OUTPUT_DIR}"

# ---- module setup -----------------------------------------------------------
module purge
module load bear-apps/2023a
module load QIIME2/2025.4

export TMPDIR="/rds/projects/e/elhamsak-pd-thesis/tmp_qiime/${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"
echo "TMPDIR set to: ${TMPDIR}"

GENUS_TABLE="${OUTPUT_DIR}/genus_table.qza"
EXPORT_DIR="${OUTPUT_DIR}/genus_table_export"

# ---- Step 1: collapse to genus level ----------------------------------------
echo ""
echo "--- Step 1: qiime taxa collapse ---"
echo "Start: $(date)"

qiime taxa collapse \
  --i-table "${FEATURE_TABLE}" \
  --i-taxonomy "${TAXONOMY}" \
  --p-level "${COLLAPSE_LEVEL}" \
  --o-collapsed-table "${GENUS_TABLE}" \
  --verbose

echo "Step 1 complete: $(date)"

if [[ ! -s "${GENUS_TABLE}" ]]; then
  echo "ERROR: taxa collapse produced no output or empty file at ${GENUS_TABLE}" >&2
  exit 1
fi

# ---- Step 2: export to BIOM, then convert to flat TSV ----------------------
echo ""
echo "--- Step 2: export collapsed table ---"

qiime tools export \
  --input-path "${GENUS_TABLE}" \
  --output-path "${EXPORT_DIR}"

if [[ ! -f "${EXPORT_DIR}/feature-table.biom" ]]; then
  echo "ERROR: expected exported feature-table.biom not found at ${EXPORT_DIR}" >&2
  exit 1
fi

# Convert BIOM -> TSV (biom-format is available within the QIIME2 env)
biom convert \
  -i "${EXPORT_DIR}/feature-table.biom" \
  -o "${EXPORT_DIR}/feature-table.tsv" \
  --to-tsv

if [[ ! -f "${EXPORT_DIR}/feature-table.tsv" ]]; then
  echo "ERROR: BIOM-to-TSV conversion failed, expected file not found" >&2
  exit 1
fi

# ---- Sanity check: how many genera retained, how many samples ------------
N_GENERA=$(tail -n +3 "${EXPORT_DIR}/feature-table.tsv" | wc -l)
N_SAMPLES=$(head -n 2 "${EXPORT_DIR}/feature-table.tsv" | tail -n 1 | awk -F'\t' '{print NF-1}')

echo ""
echo "=== Done ==="
echo "Genera retained (pre-filtering): ${N_GENERA}"
echo "Samples in table: ${N_SAMPLES}"
echo "Output: ${GENUS_TABLE}"
echo "Exported TSV: ${EXPORT_DIR}/feature-table.tsv"
echo "Finished: $(date)"

echo ""
echo "NOTE: this table is UNFILTERED (no rare-taxon removal) and RAW COUNTS"
echo "(no CLR transformation). Apply rare-taxon filtering (<10% prevalence OR"
echo "mean rel. abundance <0.01%) and CLR transformation (pseudocount 0.5)"
echo "as a separate step during cross-study merging, per schedule doc Week 5 plan."

rm -rf "${TMPDIR:?}"