#!/bin/bash
# /rds/projects/e/elhamsak-pd-thesis/scripts/pd_01_verify_host_removal.sh
# Generalized aggregate QC check across Bowtie2 host-removal outputs.
# Reusable for Bedarf, Boktor RUMC, Boktor TBC, Clasen, and Wallen.
# Run on the login node after each array job completes.
#
# USAGE:
#   ./pd_01_verify_host_removal.sh <dataset_key>
#
# Recognized dataset_key values (edit the case block below to add more):
#   bedarf2017 | boktor_rumc | boktor_tbc | clasen2024 | wallen2022
#
# /rds/projects/e/elhamsak-pd-thesis/scripts/pd_01_verify_host_removal.sh
set -uo pipefail

DATASET="${1:-}"
if [[ -z "$DATASET" ]]; then
  echo "Usage: $0 <dataset_key>"
  echo "  dataset_key one of: bedarf2017 boktor_rumc boktor_tbc clasen2024 wallen2022"
  exit 1
fi

# --- Path table -------------------------------------------------------
# Wallen lives under elhamsak-ad-thesis (storage balancing) - not a typo.
# Adjust sample-list filenames here if your actual naming differs.
PD_ROOT="/rds/projects/e/elhamsak-pd-thesis"
AD_ROOT="/rds/projects/e/elhamsak-ad-thesis"

case "$DATASET" in
  bedarf2017)
    OUT_DIR="${PD_ROOT}/pd_host_removed/bedarf2017"
    SAMPLE_LIST="${PD_ROOT}/pd_raw_data/bedarf2017/bedarf2017_samples.txt"
    ;;
  boktor_rumc)
    OUT_DIR="${PD_ROOT}/pd_host_removed/boktor2023/rumc"
    SAMPLE_LIST="${PD_ROOT}/pd_raw_data/boktor2023/rumc/boktor2023rumc_samples.txt"
    ;;
  boktor_tbc)
    OUT_DIR="${PD_ROOT}/pd_host_removed/boktor2023/tbc"
    SAMPLE_LIST="${PD_ROOT}/pd_raw_data/boktor2023/tbc/boktor2023tbc_samples.txt"
    ;;
  clasen2024)
    OUT_DIR="${PD_ROOT}/pd_host_removed/clasen2024"
    SAMPLE_LIST="${PD_ROOT}/pd_raw_data/clasen2024/clasen2024_samples.txt"
    ;;
  wallen2022)
    OUT_DIR="${AD_ROOT}/pd_host_removed/wallen2022"
    SAMPLE_LIST="${AD_ROOT}/pd_raw_data/wallen2022/wallen2022_samples.txt"
    ;;
  mao2021)
    OUT_DIR=/rds/projects/e/elhamsak-pd-thesis/pd_host_removed/mao2021
    SAMPLE_LIST=/rds/projects/e/elhamsak-pd-thesis/pd_raw_data/mao2021/mao2021_samples.txt
    ;;
  *)
    echo "Unknown dataset_key: ${DATASET}"
    echo "  Recognized: bedarf2017 boktor_rumc boktor_tbc clasen2024 wallen2022 mao2021"
    exit 1
    ;;
esac

if [[ ! -f "$SAMPLE_LIST" ]]; then
  echo "ERROR: sample list not found at ${SAMPLE_LIST}"
  echo "  Edit the path table in this script if your naming differs, then re-run."
  exit 1
fi

TOTAL=$(wc -l < "$SAMPLE_LIST")
MISSING_OUTPUT=0
ZERO_SIZE=0
OUT_OF_RANGE=0
MISSING_SUMMARY=0
UNPARSEABLE=0
LOW_THRESHOLD=0.5    # below this % is suspicious (near-zero, "already depleted" pattern)
HIGH_THRESHOLD=15    # above this % is suspicious (contamination / mislabeled sample type)

echo "=== ${DATASET} host-removal QC summary ==="
echo "OUT_DIR:      ${OUT_DIR}"
echo "SAMPLE_LIST:  ${SAMPLE_LIST}"
echo "Expected samples: ${TOTAL}"
echo ""

# Optional: collect rates for a quick distribution summary at the end
RATES_TMP=$(mktemp)

while IFS= read -r SAMPLE; do
  [[ -z "$SAMPLE" ]] && continue
  R1="${OUT_DIR}/${SAMPLE}_host_removed_1.fastq.gz"
  R2="${OUT_DIR}/${SAMPLE}_host_removed_2.fastq.gz"
  SUMMARY="${OUT_DIR}/${SAMPLE}_bowtie2_host_removed_summary.txt"

  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "MISSING OUTPUT: ${SAMPLE}"
    MISSING_OUTPUT=$((MISSING_OUTPUT+1))
    continue
  fi

  # exit code 0 does not mean success - check actual file sizes
  SIZE_R1=$(stat -c%s "$R1" 2>/dev/null || echo 0)
  SIZE_R2=$(stat -c%s "$R2" 2>/dev/null || echo 0)
  if [[ "$SIZE_R1" -lt 1000 || "$SIZE_R2" -lt 1000 ]]; then
    echo "ZERO/NEAR-ZERO SIZE: ${SAMPLE} (R1=${SIZE_R1}B R2=${SIZE_R2}B)"
    ZERO_SIZE=$((ZERO_SIZE+1))
    continue
  fi

  if [[ ! -f "$SUMMARY" ]]; then
    echo "MISSING SUMMARY: ${SAMPLE}"
    MISSING_SUMMARY=$((MISSING_SUMMARY+1))
    continue
  fi

  RATE=$(grep "overall alignment rate" "$SUMMARY" | grep -oE '[0-9]+\.[0-9]+' | head -1)
  if [[ -z "$RATE" ]]; then
    echo "COULD NOT PARSE RATE: ${SAMPLE}"
    UNPARSEABLE=$((UNPARSEABLE+1))
    continue
  fi

  echo "$RATE" >> "$RATES_TMP"

  # compare using awk (bash can't do float comparison natively)
  BELOW_LOW=$(awk -v r="$RATE" -v t="$LOW_THRESHOLD" 'BEGIN{print (r<t)?1:0}')
  ABOVE_HIGH=$(awk -v r="$RATE" -v t="$HIGH_THRESHOLD" 'BEGIN{print (r>t)?1:0}')
  if [[ "$BELOW_LOW" -eq 1 || "$ABOVE_HIGH" -eq 1 ]]; then
    echo "OUT OF RANGE: ${SAMPLE} -- ${RATE}%"
    OUT_OF_RANGE=$((OUT_OF_RANGE+1))
  fi
done < "$SAMPLE_LIST"

echo ""
echo "=== Summary: ${DATASET} ==="
echo "Total expected:      ${TOTAL}"
echo "Missing output:      ${MISSING_OUTPUT}"
echo "Zero/near-zero size: ${ZERO_SIZE}"
echo "Missing summary:     ${MISSING_SUMMARY}"
echo "Unparseable rate:    ${UNPARSEABLE}"
echo "Out-of-range rate:   ${OUT_OF_RANGE} (outside ${LOW_THRESHOLD}%-${HIGH_THRESHOLD}%)"

if [[ -s "$RATES_TMP" ]]; then
  echo ""
  echo "--- Alignment rate distribution (parsed samples) ---"
  sort -n "$RATES_TMP" | awk '
    { a[NR]=$1; sum+=$1 }
    END {
      n=NR
      if (n==0) { print "no data"; exit }
      printf "n=%d  min=%.2f%%  max=%.2f%%  mean=%.2f%%  median=%.2f%%\n",
        n, a[1], a[n], sum/n, (n%2==1 ? a[(n+1)/2] : (a[n/2]+a[n/2+1])/2)
    }'
fi
rm -f "$RATES_TMP"

echo ""
echo "If missing/zero-size counts are nonzero: check sacct for those SLURM array"
echo "task IDs for OOM/timeout/failed states before re-running."
echo "If out-of-range count is high: inspect those samples' summary files directly"
echo "  -- could indicate contamination, mislabeling, or (if near-zero across the"
echo "  whole dataset) the same pre-depletion pattern seen in Bedarf/Boktor."
echo "  For Wallen specifically: the paper's own pipeline reports <3% typical"
echo "  host-read removal, so near-zero here is NOT automatically the depletion"
echo "  pattern - it needs to be compared against that expected baseline, not"
echo "  against the Bedarf/Boktor 'fully pre-depleted' pattern."