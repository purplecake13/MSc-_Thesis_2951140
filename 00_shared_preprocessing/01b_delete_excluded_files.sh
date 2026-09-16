#!/bin/bash
#SBATCH --job-name=filter_delete_files
#SBATCH --account=elhamsak-pd-thesis
#SBATCH --qos=bbdefault
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/02b_filter_delete/02b_filter_delete_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/02b_filter_delete/02b_filter_delete_%j.err
#
# =============================================================================
# 03_delete_excluded_samples.sh — delete downloaded FASTQ files for samples
# that fall OUTSIDE a kept-accessions list (e.g. non-gut, non-human, or MCI
# samples that were downloaded before a filter was applied).
#
# Designed to be used directly after 02_filter_studies.sh, using the SAME
# accessions list it produces ("kept" set) to identify what's safe to keep,
# and deleting everything else found on disk.
#
# SAFETY: this script ALWAYS does a dry run first and requires an explicit
# --confirm flag to actually delete anything. Running without --confirm just
# shows you what WOULD be deleted, with total size, so you can review before
# committing.
#
# Usage (dry run, default — shows what would be deleted, deletes nothing):
#   ./03_delete_excluded_samples.sh \
#     --accessions-list /path/to/dataset_accessions_<suffix>.txt \
#     --data-dir /path/to/downloaded/files
#
# Usage (actually delete):
#   ./03_delete_excluded_samples.sh \
#     --accessions-list /path/to/dataset_accessions_<suffix>.txt \
#     --data-dir /path/to/downloaded/files \
#     --confirm
#
# --accessions-list should be the file produced by 02_filter_studies.sh's
# Step 2 (one accession per line) — this is the KEEP list. Any FASTQ files
# in --data-dir whose accession is NOT in this list will be flagged for
# deletion.
# =============================================================================

set -euo pipefail

# ---- defaults -----------------------------------------------------------
CONFIRM=0

# ---- parse args -----------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --accessions-list) ACCESSIONS_LIST="$2"; shift 2 ;;
    --data-dir) DATA_DIR="$2"; shift 2 ;;
    --confirm) CONFIRM=1; shift ;;
    *) echo "ERROR: Unknown argument: $1" >&2; exit 1 ;;
  esac
done

# ---- validate required args ------------------------------------------------
: "${ACCESSIONS_LIST:?--accessions-list is required}"
: "${DATA_DIR:?--data-dir is required}"

if [[ ! -f "${ACCESSIONS_LIST}" ]]; then
  echo "ERROR: accessions list not found: ${ACCESSIONS_LIST}" >&2
  exit 1
fi
if [[ ! -d "${DATA_DIR}" ]]; then
  echo "ERROR: data directory not found: ${DATA_DIR}" >&2
  exit 1
fi

N_KEEP=$(wc -l < "${ACCESSIONS_LIST}")

echo "========================================"
echo "Delete excluded samples"
echo "Started: $(date)"
echo "========================================"
echo "Accessions list (KEEP set): ${ACCESSIONS_LIST} (${N_KEEP} accessions)"
echo "Data dir: ${DATA_DIR}"
echo "Mode: $([[ $CONFIRM -eq 1 ]] && echo 'CONFIRMED DELETE' || echo 'DRY RUN (nothing will be deleted)')"
echo ""

# ---- find files on disk whose accession is NOT in the keep list -----------
TO_DELETE=$(mktemp)

ls "${DATA_DIR}"/*.fastq.gz 2>/dev/null \
  | while read -r f; do
      ACC=$(basename "$f" | sed -E 's/_[12]\.fastq\.gz$//; s/\.fastq\.gz$//')
      if ! grep -qxF "${ACC}" "${ACCESSIONS_LIST}"; then
        echo "$f"
      fi
    done > "${TO_DELETE}"

N_TO_DELETE=$(wc -l < "${TO_DELETE}")

if [[ "${N_TO_DELETE}" -eq 0 ]]; then
  echo "Nothing to delete — all files on disk are in the keep list."
  rm -f "${TO_DELETE}"
  exit 0
fi

echo "Files identified for deletion (${N_TO_DELETE} files):"
cat "${TO_DELETE}"
echo ""

TOTAL_SIZE=$(du -ch $(cat "${TO_DELETE}") 2>/dev/null | tail -1 | cut -f1)
echo "Total size of files to delete: ${TOTAL_SIZE}"
echo ""

if [[ "${CONFIRM}" -eq 0 ]]; then
  echo "DRY RUN — no files deleted."
  echo "Review the list above carefully, then rerun with --confirm to actually delete."
  rm -f "${TO_DELETE}"
  exit 0
fi

# ---- actually delete, only reached if --confirm was passed -----------------
echo "Deleting ${N_TO_DELETE} files..."
while read -r f; do
  rm -f "$f"
  echo "  deleted: $f"
done < "${TO_DELETE}"

echo ""
echo "Done. Deleted ${N_TO_DELETE} files, freed approximately ${TOTAL_SIZE}."
echo "Finished: $(date)"

rm -f "${TO_DELETE}"