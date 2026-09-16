#!/bin/bash
#SBATCH --job-name=v1v2_classifier_verify
#SBATCH --account=elhamsak-pd-thesis
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/09_taxonomy/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/09_taxonomy/%x_%j.err

set -euo pipefail

echo "=========================================="
echo "V1-V2 classifier verification"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "=========================================="

PROJECT_DIR="/rds/projects/e/elhamsak-pd-thesis"
OUT_DIR="/rds/projects/e/elhamsak-ad-thesis/reference_genomes/silva138/v1v2_classifier"
LOG_DIR="${PROJECT_DIR}/logs/taxonomy"

SEQS_V1V2="${OUT_DIR}/silva-138-v1v2-seqs.qza"
CLASSIFIER_OUT="${OUT_DIR}/silva-138-v1v2-classifier.qza"

mkdir -p "${LOG_DIR}"

# Redirect TMPDIR to RDS project space (see note in 12b) -- qiime tools export
# below also stages files and can hit the same node-local /tmp space limit.
export TMPDIR="${PROJECT_DIR}/tmp_qiime"
mkdir -p "${TMPDIR}"
echo "TMPDIR set to: ${TMPDIR}"

module purge
module load bear-apps/2023a
module load QIIME2/2025.4

# --- Check 1: artifact integrity ---
# qiime tools peek will fail loudly if the .qza is corrupt or incomplete,
# which a bare file-existence check would not catch.
echo ""
echo "--- Check 1: artifact integrity (qiime tools peek) ---"

if [[ ! -f "${SEQS_V1V2}" ]]; then
    echo "ERROR: extracted V1-V2 reference sequences not found at ${SEQS_V1V2}" >&2
    exit 1
fi
if [[ ! -f "${CLASSIFIER_OUT}" ]]; then
    echo "ERROR: trained classifier not found at ${CLASSIFIER_OUT}" >&2
    exit 1
fi

echo ""
echo "Peeking extracted V1-V2 sequences artifact:"
qiime tools peek "${SEQS_V1V2}"

echo ""
echo "Peeking trained classifier artifact:"
qiime tools peek "${CLASSIFIER_OUT}"

# --- Check 2: extracted region length sanity check ---
# V1-V2 amplicons (27F-mod/338R) are documented around ~310-430bp depending on
# exact primer variant and organism. If extracted reads are wildly outside this
# range, the primer match likely went wrong (e.g. matched a different region,
# or matched too permissively/restrictively due to the degenerate bases).
echo ""
echo "--- Check 2: extracted sequence length distribution ---"

EXPORT_DIR="${OUT_DIR}/seq_length_check_tmp"
mkdir -p "${EXPORT_DIR}"

qiime tools export \
    --input-path "${SEQS_V1V2}" \
    --output-path "${EXPORT_DIR}"

if [[ -f "${EXPORT_DIR}/dna-sequences.fasta" ]]; then
    echo "Sequence length summary (extracted V1-V2 region):"
    awk '/^>/{next} {print length($0)}' "${EXPORT_DIR}/dna-sequences.fasta" | \
        sort -n | \
        awk '
        {
            a[NR]=$1; sum+=$1; count++
        }
        END {
            if (count == 0) {
                print "WARNING: no sequences found in exported FASTA"
                exit 1
            }
            print "  Count:  " count
            print "  Min:    " a[1]
            print "  Max:    " a[count]
            print "  Mean:   " sum/count
            print "  Median: " a[int((count+1)/2)]
        }'
    echo ""
    echo "Expected range for 27F-mod/338R V1-V2 amplicons: roughly 280-430bp."
    echo "If min/max/mean fall well outside this, re-check primer sequences before proceeding."
else
    echo "ERROR: expected dna-sequences.fasta not found after export" >&2
    exit 1
fi

# Clean up temp export
rm -rf "${EXPORT_DIR}"

echo ""
echo "--- Check 3: feature count sanity ---"
echo "A V1-V2 extraction from full-length SILVA 138 should retain a substantial"
echo "fraction of reference sequences (tens of thousands to ~100k+ records likely)."
echo "If the count is suspiciously low (e.g. only hundreds), the primers may have"
echo "matched too few reference sequences -- check degenerate base handling."

echo ""
echo "=========================================="
echo "Verification complete: $(date)"
echo "Review the length distribution and feature counts above manually before"
echo "using this classifier on real Yamashiro 2024 ASVs."
echo "=========================================="

# --- Clean up TMPDIR contents ---
rm -rf "${TMPDIR:?}"/*