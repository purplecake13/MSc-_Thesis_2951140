#!/bin/bash
#SBATCH --job-name=silva138_2_download
#SBATCH --account=elhamsak-pd-thesis
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/09_taxonomy/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/09_taxonomy/%x_%j.err
set -euo pipefail

echo "=========================================="
echo "SILVA 138.2 reference download (via RESCRIPt)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "=========================================="
echo "NOTE: this replaces the original SILVA 138 (2020) reference used to train"
echo "the V1-V2 classifier, to match the SILVA 138.2 (2025) taxonomy already"
echo "used by the pretrained full-length classifier -- fixes a genus/phylum"
echo "nomenclature mismatch identified via cross-classifier comparison on"
echo "Yamashiro2024 (e.g. Firmicutes/Bacillota, Prevotella/Segatella)."
echo "=========================================="

REF_DIR="/rds/projects/e/elhamsak-ad-thesis/reference_genomes/silva138_2"
LOG_DIR="/rds/projects/e/elhamsak-pd-thesis/logs/09_taxonomy"
mkdir -p "${REF_DIR}"
mkdir -p "${LOG_DIR}"
cd "${REF_DIR}"

export TMPDIR="/rds/projects/e/elhamsak-pd-thesis/tmp_qiime/${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"
echo "TMPDIR set to: ${TMPDIR}"

module purge
module load bear-apps/2023a
module load QIIME2/2025.4
echo "QIIME2 module loaded:"
qiime --version

SEQS_OUT="${REF_DIR}/silva-138-2-ssu-nr99-seqs.qza"
TAX_OUT="${REF_DIR}/silva-138-2-ssu-nr99-tax.qza"

echo ""
echo "--- Downloading + parsing SILVA 138.2 SSURef_NR99 via RESCRIPt ---"
echo "Start: $(date)"

qiime rescript get-silva-data \
    --p-version '138.2' \
    --p-target 'SSURef_NR99' \
    --o-silva-sequences "${SEQS_OUT}" \
    --o-silva-taxonomy "${TAX_OUT}" \
    --verbose

echo "Step complete: $(date)"

if [[ ! -s "${SEQS_OUT}" ]]; then
    echo "ERROR: ${SEQS_OUT} missing or empty after download" >&2
    exit 1
fi
if [[ ! -s "${TAX_OUT}" ]]; then
    echo "ERROR: ${TAX_OUT} missing or empty after download" >&2
    exit 1
fi

echo "File sizes:"
du -h "${SEQS_OUT}" "${TAX_OUT}"

echo ""
echo "--- Reverse-transcribing RNA sequences to DNA (required for extract-reads) ---"
echo "Start: $(date)"

SEQS_DNA_OUT="${REF_DIR}/silva-138-2-ssu-nr99-seqs-dna.qza"

qiime rescript reverse-transcribe \
    --i-rna-sequences "${SEQS_OUT}" \
    --o-dna-sequences "${SEQS_DNA_OUT}" \
    --verbose

echo "Step complete: $(date)"

if [[ ! -s "${SEQS_DNA_OUT}" ]]; then
    echo "ERROR: reverse-transcribe produced no output or empty file at ${SEQS_DNA_OUT}" >&2
    exit 1
fi
echo "DNA sequences: ${SEQS_DNA_OUT} ($(du -h "${SEQS_DNA_OUT}" | cut -f1))"

echo ""
echo "--- Confirming artifact validity with qiime tools peek ---"
qiime tools peek "${SEQS_OUT}"
qiime tools peek "${TAX_OUT}"
qiime tools peek "${SEQS_DNA_OUT}"

echo ""
echo "=========================================="
echo "SILVA 138.2 reference download complete and verified."
echo "Sequences: ${SEQS_OUT}"
echo "Sequences (DNA): ${SEQS_DNA_OUT}"
echo "Taxonomy:  ${TAX_OUT}"
echo "Finished: $(date)"
echo "=========================================="

rm -rf "${TMPDIR:?}"/*