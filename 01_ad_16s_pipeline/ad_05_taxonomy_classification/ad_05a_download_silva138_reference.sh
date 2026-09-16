#!/bin/bash
#SBATCH --job-name=silva138_download
#SBATCH --account=elhamsak-pd-thesis
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/09_taxonomy/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/09_taxonomy/%x_%j.err

set -euo pipefail

echo "=========================================="
echo "SILVA 138 reference download"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "=========================================="

REF_DIR="/rds/projects/e/elhamsak-ad-thesis/reference_genomes/silva138"
LOG_DIR="/rds/projects/e/elhamsak-pd-thesis/logs/09_taxonomy"

mkdir -p "${REF_DIR}"
mkdir -p "${LOG_DIR}"
cd "${REF_DIR}"

SEQS_URL="https://data.qiime2.org/2024.10/common/silva-138-99-seqs.qza"
TAX_URL="https://data.qiime2.org/2024.10/common/silva-138-99-tax.qza"

# Expected MD5s as published on the QIIME2 2024.10 data-resources page
SEQS_MD5_EXPECTED="de8886bb2c059b1e8752255d271f3010"
TAX_MD5_EXPECTED="f12d5b78bf4b1519721fe52803581c3d"

echo ""
echo "--- Downloading full-length SILVA 138 sequences ---"
wget --progress=dot:giga "${SEQS_URL}"

echo ""
echo "--- Downloading full-length SILVA 138 taxonomy ---"
wget --progress=dot:giga "${TAX_URL}"

# --- Verify downloads actually completed and are not truncated/corrupted ---
echo ""
echo "--- Verifying downloads ---"

SEQS_FILE="${REF_DIR}/silva-138-99-seqs.qza"
TAX_FILE="${REF_DIR}/silva-138-99-tax.qza"

if [[ ! -s "${SEQS_FILE}" ]]; then
    echo "ERROR: ${SEQS_FILE} missing or empty after download" >&2
    exit 1
fi
if [[ ! -s "${TAX_FILE}" ]]; then
    echo "ERROR: ${TAX_FILE} missing or empty after download" >&2
    exit 1
fi

echo "File sizes:"
du -h "${SEQS_FILE}" "${TAX_FILE}"

echo ""
echo "Computing MD5 checksums (this can take a few minutes for large files)..."
SEQS_MD5_ACTUAL=$(md5sum "${SEQS_FILE}" | awk '{print $1}')
TAX_MD5_ACTUAL=$(md5sum "${TAX_FILE}" | awk '{print $1}')

echo "  silva-138-99-seqs.qza expected MD5: ${SEQS_MD5_EXPECTED}"
echo "  silva-138-99-seqs.qza actual MD5:   ${SEQS_MD5_ACTUAL}"
echo "  silva-138-99-tax.qza  expected MD5: ${TAX_MD5_EXPECTED}"
echo "  silva-138-99-tax.qza  actual MD5:   ${TAX_MD5_ACTUAL}"

MD5_FAIL=0
if [[ "${SEQS_MD5_ACTUAL}" != "${SEQS_MD5_EXPECTED}" ]]; then
    echo "ERROR: seqs.qza MD5 mismatch -- download is corrupt or incomplete" >&2
    MD5_FAIL=1
fi
if [[ "${TAX_MD5_ACTUAL}" != "${TAX_MD5_EXPECTED}" ]]; then
    echo "ERROR: tax.qza MD5 mismatch -- download is corrupt or incomplete" >&2
    MD5_FAIL=1
fi

if [[ "${MD5_FAIL}" -eq 1 ]]; then
    echo "One or more MD5 checks failed. Do NOT proceed to classifier training until resolved." >&2
    exit 1
fi

echo ""
echo "Both checksums verified -- downloads are confirmed intact."

# --- Confirm QIIME2 can actually read these as valid artifacts ---
echo ""
echo "--- Confirming artifact validity with qiime tools peek ---"
module purge
module load bear-apps/2023a
module load QIIME2/2025.4

qiime tools peek "${SEQS_FILE}"
qiime tools peek "${TAX_FILE}"

echo ""
echo "=========================================="
echo "SILVA 138 reference download complete and verified."
echo "Finished: $(date)"
echo "=========================================="