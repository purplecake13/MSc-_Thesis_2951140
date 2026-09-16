#!/bin/bash
#SBATCH --job-name=v1v2_classifier_train
#SBATCH --account=elhamsak-pd-thesis
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/09_taxonomy/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/09_taxonomy/%x_%j.err

set -euo pipefail

echo "=========================================="
echo "V1-V2 custom SILVA 138 classifier training"
echo "Primers: 27F-mod / 338R (degenerate, matches Yamashiro 2024)"
echo "  Forward: AGRGTTTGATYMTGGCTCAG"
echo "  Reverse: TGCTGCCTCCCGTAGGAGT"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "=========================================="

# --- Paths (edit if your SILVA reference files live elsewhere) ---
PROJECT_DIR="/rds/projects/e/elhamsak-pd-thesis"
REF_DIR="/rds/projects/e/elhamsak-ad-thesis/reference_genomes/silva138"
OUT_DIR="/rds/projects/e/elhamsak-ad-thesis/reference_genomes/silva138/v1v2_classifier"
LOG_DIR="${PROJECT_DIR}/logs/09_taxonomy"

mkdir -p "${OUT_DIR}"
mkdir -p "${LOG_DIR}"

# --- Redirect TMPDIR to RDS project space instead of the compute node's local /tmp ---
# fit-classifier-naive-bayes stages substantial intermediate data while packaging
# the final classifier artifact. Node-local /tmp is often a small disk (or tmpfs)
# scoped to a single compute node, separate from and much smaller than RDS project
# storage. Pointing TMPDIR here avoids "No space left on device" failures even
# though the RDS project itself has plenty of free space (confirmed: 3.0T free).
export TMPDIR="${PROJECT_DIR}/tmp_qiime"
mkdir -p "${TMPDIR}"
echo "TMPDIR set to: ${TMPDIR}"

SEQS_IN="${REF_DIR}/silva-138-99-seqs.qza"
TAX_IN="${REF_DIR}/silva-138-99-tax.qza"
SEQS_V1V2="${OUT_DIR}/silva-138-v1v2-seqs.qza"
CLASSIFIER_OUT="${OUT_DIR}/silva-138-v1v2-classifier.qza"

# --- Sanity check inputs exist before burning walltime ---
if [[ ! -f "${SEQS_IN}" ]]; then
    echo "ERROR: Reference sequences not found at ${SEQS_IN}" >&2
    exit 1
fi
if [[ ! -f "${TAX_IN}" ]]; then
    echo "ERROR: Reference taxonomy not found at ${TAX_IN}" >&2
    exit 1
fi

# --- Load QIIME2 module (BlueBEAR) ---
module purge
module load bear-apps/2023a
module load QIIME2/2025.4

echo "QIIME2 module loaded:"
qiime --version

# --- Step 1: Extract V1-V2 region from full SILVA reference using 27F/338R primers ---
echo ""
echo "--- Step 1: extract-reads (V1-V2 region, 27F/338R primers) ---"
echo "Start: $(date)"

qiime feature-classifier extract-reads \
    --i-sequences "${SEQS_IN}" \
    --p-f-primer AGRGTTTGATYMTGGCTCAG \
    --p-r-primer TGCTGCCTCCCGTAGGAGT \
    --o-reads "${SEQS_V1V2}" \
    --verbose

echo "Step 1 complete: $(date)"

# Confirm output was actually produced (don't trust exit code alone)
if [[ ! -s "${SEQS_V1V2}" ]]; then
    echo "ERROR: extract-reads produced no output or empty file at ${SEQS_V1V2}" >&2
    exit 1
fi
echo "Step 1 output verified: ${SEQS_V1V2} ($(du -h "${SEQS_V1V2}" | cut -f1))"

# --- Step 2: Train naive Bayes classifier on extracted V1-V2 region ---
echo ""
echo "--- Step 2: fit-classifier-naive-bayes ---"
echo "Start: $(date)"

qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "${SEQS_V1V2}" \
    --i-reference-taxonomy "${TAX_IN}" \
    --o-classifier "${CLASSIFIER_OUT}" \
    --verbose

echo "Step 2 complete: $(date)"

# Confirm classifier output was actually produced
if [[ ! -s "${CLASSIFIER_OUT}" ]]; then
    echo "ERROR: fit-classifier-naive-bayes produced no output or empty file at ${CLASSIFIER_OUT}" >&2
    exit 1
fi
echo "Step 2 output verified: ${CLASSIFIER_OUT} ($(du -h "${CLASSIFIER_OUT}" | cut -f1))"

echo ""
echo "=========================================="
echo "V1-V2 classifier training complete."
echo "Classifier ready at: ${CLASSIFIER_OUT}"
echo "Finished: $(date)"
echo "=========================================="

# --- Clean up TMPDIR contents now that the run succeeded ---
# Leaving these around across reruns would slowly eat into RDS project space.
echo ""
echo "Cleaning up TMPDIR contents (${TMPDIR})..."
rm -rf "${TMPDIR:?}"/*
echo "TMPDIR cleaned."