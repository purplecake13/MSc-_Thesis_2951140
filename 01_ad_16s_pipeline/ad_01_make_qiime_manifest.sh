#!/bin/bash
# =============================================================================
# /rds/projects/e/elhamsak-pd-thesis/scripts/ad_01_make_qiime_manifest.sh — generate a QIIME2 PairedEndFastqManifestPhred33V2
# manifest TSV from a fastp-trimmed data directory + sample list.
#
# Usage:
#   ./08_make_qiime_manifest.sh \
#     --data-dir /rds/projects/e/elhamsak-ad-thesis/ad_trimmed/tran2019 \
#     --sample-list /rds/projects/e/elhamsak-ad-thesis/ad_raw_data/tran2019/tran2019_samples.txt \
#     --output-manifest /rds/projects/e/elhamsak-ad-thesis/ad_raw_data/tran2019/manifest_tran2019.tsv
#
# sample-list = one sample basename per line (same convention as 06_fastp.sh /
#               07_bowtie2_host_removal.sh — built from R1 files only)
#
# Expects input files named: <data-dir>/<sample>_1.trimmed.fastq.gz
#                             <data-dir>/<sample>_2.trimmed.fastq.gz
# (i.e. the default output naming from 06_fastp.sh)
# =============================================================================

set -euo pipefail

# ---- defaults -----------------------------------------------------------
INPUT_SUFFIX_R1="_1.trimmed.fastq.gz"
INPUT_SUFFIX_R2="_2.trimmed.fastq.gz"

# ---- parse args -----------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir) DATA_DIR="$2"; shift 2 ;;
    --sample-list) SAMPLE_LIST="$2"; shift 2 ;;
    --output-manifest) OUTPUT_MANIFEST="$2"; shift 2 ;;
    --r1-suffix) INPUT_SUFFIX_R1="$2"; shift 2 ;;
    --r2-suffix) INPUT_SUFFIX_R2="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

# ---- validate required args ------------------------------------------------
: "${DATA_DIR:?--data-dir is required}"
: "${SAMPLE_LIST:?--sample-list is required}"
: "${OUTPUT_MANIFEST:?--output-manifest is required}"

if [[ ! -d "$DATA_DIR" ]]; then
  echo "ERROR: data directory not found: $DATA_DIR" >&2
  exit 1
fi

if [[ ! -f "$SAMPLE_LIST" ]]; then
  echo "ERROR: sample list not found: $SAMPLE_LIST" >&2
  exit 1
fi

N_SAMPLES=$(grep -c . "$SAMPLE_LIST" || true)
if [[ "$N_SAMPLES" -eq 0 ]]; then
  echo "ERROR: sample list is empty: $SAMPLE_LIST" >&2
  exit 1
fi

echo "=== Building QIIME2 manifest ==="
echo "Data dir:      $DATA_DIR"
echo "Sample list:   $SAMPLE_LIST ($N_SAMPLES samples)"
echo "R1 suffix:     $INPUT_SUFFIX_R1"
echo "R2 suffix:     $INPUT_SUFFIX_R2"
echo "Output:        $OUTPUT_MANIFEST"
echo ""

# ---- write header (PairedEndFastqManifestPhred33V2 requires this exact header,
#      tab-separated) ---------------------------------------------------------
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$OUTPUT_MANIFEST"

# ---- build rows, validating every file exists before writing ---------------
MISSING=0
WRITTEN=0

while IFS= read -r SAMPLE || [[ -n "$SAMPLE" ]]; do
  # skip blank lines defensively
  [[ -z "$SAMPLE" ]] && continue

  R1="${DATA_DIR}/${SAMPLE}${INPUT_SUFFIX_R1}"
  R2="${DATA_DIR}/${SAMPLE}${INPUT_SUFFIX_R2}"

  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "WARNING: missing files for sample '${SAMPLE}':" >&2
    [[ ! -f "$R1" ]] && echo "  missing: $R1" >&2
    [[ ! -f "$R2" ]] && echo "  missing: $R2" >&2
    MISSING=$((MISSING + 1))
    continue
  fi

  printf "%s\t%s\t%s\n" "$SAMPLE" "$R1" "$R2" >> "$OUTPUT_MANIFEST"
  WRITTEN=$((WRITTEN + 1))
done < "$SAMPLE_LIST"

echo ""
echo "=== Done ==="
echo "Samples written to manifest: $WRITTEN"
echo "Samples missing files (skipped): $MISSING"
echo "Manifest: $OUTPUT_MANIFEST"

if [[ "$MISSING" -gt 0 ]]; then
  echo ""
  echo "WARNING: $MISSING sample(s) were skipped due to missing files." >&2
  echo "Check the warnings above before importing into QIIME2 — a manifest" >&2
  echo "with fewer samples than expected will silently produce a smaller" >&2
  echo "feature table with no error." >&2
fi

if [[ "$WRITTEN" -eq 0 ]]; then
  echo "ERROR: no samples written to manifest — nothing to import." >&2
  exit 1
fi