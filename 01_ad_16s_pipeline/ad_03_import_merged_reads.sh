#!/bin/bash
#SBATCH --job-name=import_merged
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/08_ad_qiime2/import_merged/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/08_ad_qiime2/import_merged/%x_%j.err
# =============================================================================
# /rds/projects/e/elhamsak-pd-thesis/scripts/ad_03_import_merged_reads.sh — wrap fastp-merged single-end FASTQs into a
# QIIME2 SampleData[SequencesWithQuality] artifact, ready for dada2 denoise-single
#
# BRIDGES 11_fastp_merge_array.sh's output (.fastq.gz per sample) to
# 09_dada2_denoise.sh's input requirement (--demux-qza). fastp's merge mode
# produces single combined reads (no longer R1/R2 pairs), so this must be
# imported as SampleData[SequencesWithQuality] (single-end), and the
# downstream DADA2 step must use denoise-single, NOT denoise-paired.
#
# Builds a QIIME2 manifest file automatically from the merged FASTQs in
# --merged-dir, then imports.
#
# Usage:
#   sbatch --job-name=import_merged_tran2019 12_import_merged_reads.sh \
#     --merged-dir /rds/projects/e/elhamsak-ad-thesis/ad_qiime2/tran2019/fastp_merge \
#     --output-dir /rds/projects/e/elhamsak-ad-thesis/ad_qiime2/tran2019
#
# Expects merged files named: <sample>_merged.fastq.gz
# (i.e. the default output naming from 11_fastp_merge_array.sh)
#
# Output:
#   <output-dir>/manifest.tsv      (QIIME2 manifest, auto-generated)
#   <output-dir>/merged_demux.qza  (SampleData[SequencesWithQuality])
#   <output-dir>/merged_demux.qzv  (summary visualization)
# =============================================================================

set -euo pipefail

# ---- parse args -----------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --merged-dir) MERGED_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --suffix) MERGED_SUFFIX="$2"; shift 2 ;;
    --sample-list) SAMPLE_LIST="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

MERGED_SUFFIX="${MERGED_SUFFIX:-_merged.fastq.gz}"

# ---- validate required args ------------------------------------------------
: "${MERGED_DIR:?--merged-dir is required}"
: "${OUTPUT_DIR:?--output-dir is required}"

if [[ ! -d "${MERGED_DIR}" ]]; then
  echo "ERROR: merged-dir not found: ${MERGED_DIR}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

MANIFEST="${OUTPUT_DIR}/manifest.tsv"
DEMUX_QZA="${OUTPUT_DIR}/merged_demux.qza"
DEMUX_QZV="${OUTPUT_DIR}/merged_demux.qzv"

echo "=== Import merged reads ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "Merged dir: ${MERGED_DIR}"
echo "Output dir: ${OUTPUT_DIR}"

# ---- Step 1: build manifest from merged FASTQs ------------------------------
echo ""
echo "--- Step 1: building QIIME2 manifest ---"

echo -e "sample-id\tabsolute-filepath" > "${MANIFEST}"

# Build the manifest from a sample list if provided, otherwise glob the
# merged-dir directly (less safe -- prefer --sample-list when available,
# since it gives explicit per-sample validation matching 06_fastp.sh /
# 07_bowtie2_host_removal.sh / 08_make_qiime_manifest.sh convention).
MISSING=0
WRITTEN=0

if [[ -n "${SAMPLE_LIST:-}" ]]; then
  if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERROR: sample list not found: ${SAMPLE_LIST}" >&2
    exit 1
  fi
  while IFS= read -r SAMPLE || [[ -n "$SAMPLE" ]]; do
    [[ -z "$SAMPLE" ]] && continue
    F="${MERGED_DIR}/${SAMPLE}${MERGED_SUFFIX}"
    if [[ ! -f "$F" ]]; then
      echo "WARNING: missing merged file for sample '${SAMPLE}': $F" >&2
      MISSING=$((MISSING + 1))
      continue
    fi
    printf "%s\t%s\n" "$SAMPLE" "$F" >> "${MANIFEST}"
    WRITTEN=$((WRITTEN + 1))
  done < "${SAMPLE_LIST}"
else
  echo "NOTE: no --sample-list provided, falling back to globbing ${MERGED_DIR}" >&2
  echo "for *${MERGED_SUFFIX} -- prefer passing --sample-list for explicit validation." >&2
  for f in "${MERGED_DIR}"/*"${MERGED_SUFFIX}"; do
    [[ -f "$f" ]] || continue
    SAMPLE=$(basename "$f" "${MERGED_SUFFIX}")
    printf "%s\t%s\n" "${SAMPLE}" "${f}" >> "${MANIFEST}"
    WRITTEN=$((WRITTEN + 1))
  done
fi

N_FOUND="${WRITTEN}"

if [[ "${MISSING}" -gt 0 ]]; then
  echo "" >&2
  echo "WARNING: ${MISSING} sample(s) skipped due to missing merged files." >&2
  echo "A manifest with fewer samples than expected will silently produce" >&2
  echo "a smaller demux artifact with no error -- check warnings above." >&2
fi

if [[ "${N_FOUND}" -eq 0 ]]; then
  echo "ERROR: no files matching *${MERGED_SUFFIX} found in ${MERGED_DIR}" >&2
  exit 1
fi

echo "Found ${N_FOUND} merged sample files."
echo "Manifest written to: ${MANIFEST}"

# ---- module setup -----------------------------------------------------------
module purge
module load bear-apps/2023a
module load QIIME2/2025.4

export TMPDIR="/rds/projects/e/elhamsak-pd-thesis/tmp_qiime/job_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"

# ---- Step 2: import as single-end sequences ---------------------------------
echo ""
echo "--- Step 2: qiime tools import ---"
echo "Start: $(date)"

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "${MANIFEST}" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --output-path "${DEMUX_QZA}"

echo "Step 2 complete: $(date)"

if [[ ! -s "${DEMUX_QZA}" ]]; then
  echo "ERROR: import produced no output or empty file at ${DEMUX_QZA}" >&2
  exit 1
fi

# ---- Step 3: summarize for inspection ----------------------------------------
echo ""
echo "--- Step 3: qiime demux summarize ---"

qiime demux summarize \
  --i-data "${DEMUX_QZA}" \
  --o-visualization "${DEMUX_QZV}"

echo ""
echo "=== Done ==="
echo "Imported ${N_FOUND} samples."
echo "Demux artifact: ${DEMUX_QZA}"
echo "Summary viz: ${DEMUX_QZV}"
echo "Finished: $(date)"
echo ""
echo "NEXT STEP: feed ${DEMUX_QZA} into DADA2 using denoise-SINGLE (not"
echo "denoise-paired) -- merged reads are no longer R1/R2 pairs. View"
echo "${DEMUX_QZV} at https://view.qiime2.org to choose an appropriate"
echo "--p-trunc-len for denoise-single based on the quality plot."

rm -rf "${TMPDIR:?}"/*