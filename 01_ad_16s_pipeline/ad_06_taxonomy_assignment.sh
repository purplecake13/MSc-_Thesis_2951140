#!/bin/bash
#SBATCH --job-name=taxonomy_assign
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/taxonomy/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/taxonomy/%x_%j.err
# =============================================================================
# /rds/projects/e/elhamsak-pd-thesis/scripts/ad_06_taxonomy_assignment.sh — QIIME2 classify-sklearn taxonomy assignment
#
# Assigns taxonomy to a DADA2-denoised feature table + representative sequences
# using a region-matched SILVA 138 classifier.
#
# Usage:
#   sbatch --job-name=taxonomy_cirstea 13_taxonomy_assignment.sh \
#     --rep-seqs /rds/projects/.../cirstea2022/repseqs.qza \
#     --classifier /rds/projects/e/elhamsak-ad-thesis/reference_genomes/silva138/silva-138-v3v4-classifier.qza \
#     --output-dir /rds/projects/.../cirstea2022/taxonomy
#
# For V1-V2 datasets (e.g. Yamashiro 2024), point --classifier at the custom-trained
# classifier instead:
#   --classifier /rds/projects/e/elhamsak-ad-thesis/reference_genomes/silva138/v1v2_classifier/silva-138-v1v2-classifier.qza
#
# Inputs expected (from upstream DADA2 step):
#   --rep-seqs    FeatureData[Sequence] .qza (representative ASV sequences)
#   --classifier  TaxonomicClassifier .qza (pre-trained or custom-trained, region-matched)
#
# Output:
#   <output-dir>/taxonomy.qza   (FeatureData[Taxonomy])
#   <output-dir>/taxonomy.qzv   (visualization summary)
#   <output-dir>/taxonomy_export/taxonomy.tsv  (flat TSV for downstream use)
# =============================================================================

set -euo pipefail

# ---- parse args -----------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --rep-seqs) REP_SEQS="$2"; shift 2 ;;
    --classifier) CLASSIFIER="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

THREADS="${THREADS:-8}"

# ---- validate required args ------------------------------------------------
: "${REP_SEQS:?--rep-seqs is required}"
: "${CLASSIFIER:?--classifier is required}"
: "${OUTPUT_DIR:?--output-dir is required}"

if [[ ! -f "${REP_SEQS}" ]]; then
  echo "ERROR: representative sequences artifact not found: ${REP_SEQS}" >&2
  exit 1
fi
if [[ ! -f "${CLASSIFIER}" ]]; then
  echo "ERROR: classifier artifact not found: ${CLASSIFIER}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

echo "=== Taxonomy assignment ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "Rep seqs: ${REP_SEQS}"
echo "Classifier: ${CLASSIFIER}"
echo "Output dir: ${OUTPUT_DIR}"

# ---- module setup -----------------------------------------------------------
module purge
module load bear-apps/2023a
module load QIIME2/2025.4

# ---- TMPDIR redirect (avoid node-local /tmp exhaustion, per known issue) ---
export TMPDIR="/rds/projects/e/elhamsak-pd-thesis/tmp_qiime/${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"
echo "TMPDIR set to: ${TMPDIR}"

TAXONOMY_QZA="${OUTPUT_DIR}/taxonomy.qza"
TAXONOMY_QZV="${OUTPUT_DIR}/taxonomy.qzv"
EXPORT_DIR="${OUTPUT_DIR}/taxonomy_export"

# ---- Step 1: classify-sklearn ------------------------------------------------
echo ""
echo "--- Step 1: qiime feature-classifier classify-sklearn ---"
echo "Start: $(date)"

qiime feature-classifier classify-sklearn \
  --i-classifier "${CLASSIFIER}" \
  --i-reads "${REP_SEQS}" \
  --p-n-jobs "${THREADS}" \
  --o-classification "${TAXONOMY_QZA}" \
  --verbose

echo "Step 1 complete: $(date)"

if [[ ! -s "${TAXONOMY_QZA}" ]]; then
  echo "ERROR: classify-sklearn produced no output or empty file at ${TAXONOMY_QZA}" >&2
  exit 1
fi

# ---- Step 2: generate visualization summary ---------------------------------
echo ""
echo "--- Step 2: qiime metadata tabulate (visualization) ---"

qiime metadata tabulate \
  --m-input-file "${TAXONOMY_QZA}" \
  --o-visualization "${TAXONOMY_QZV}"

# ---- Step 3: export to flat TSV for downstream genus collapse --------------
echo ""
echo "--- Step 3: export taxonomy to TSV ---"

qiime tools export \
  --input-path "${TAXONOMY_QZA}" \
  --output-path "${EXPORT_DIR}"

if [[ ! -f "${EXPORT_DIR}/taxonomy.tsv" ]]; then
  echo "ERROR: expected exported taxonomy.tsv not found at ${EXPORT_DIR}" >&2
  exit 1
fi

# ---- Sanity check: how many ASVs got a confident genus-level call? ---------
N_TOTAL=$(tail -n +2 "${EXPORT_DIR}/taxonomy.tsv" | wc -l)
N_GENUS=$(tail -n +2 "${EXPORT_DIR}/taxonomy.tsv" | awk -F'\t' '$2 ~ /g__[A-Za-z]/' | wc -l)

echo ""
echo "=== Done ==="
echo "Total ASVs classified: ${N_TOTAL}"
echo "ASVs with a genus-level call: ${N_GENUS} ($(awk "BEGIN {printf \"%.1f\", ${N_GENUS}/${N_TOTAL}*100}")%)"
echo "Output: ${TAXONOMY_QZA}"
echo "Exported TSV: ${EXPORT_DIR}/taxonomy.tsv"
echo "Finished: $(date)"

if [[ "${N_GENUS}" -eq 0 ]]; then
  echo "WARNING: zero ASVs received a genus-level classification. This likely indicates" >&2
  echo "a classifier/region mismatch -- double check the classifier matches this dataset's" >&2
  echo "amplicon region before trusting any downstream genus-level results." >&2
fi

# Clean up TMPDIR
rm -rf "${TMPDIR:?}"/*