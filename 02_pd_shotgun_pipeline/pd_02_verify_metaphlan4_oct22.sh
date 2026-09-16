#!/bin/bash
# pd_02_verify_metaphlan4_oct22.sh
# Usage: ./pd_02_verify_metaphlan4_oct22.sh <dataset_key>
# dataset_key: bedarf2017 | boktor_rumc | boktor_tbc | clasen2024 | wallen2022 | mao2021
#
# Checks, per sample:
#   1. Profile file exists
#   2. Profile file is non-trivial size (catches stub/header-only files)
#   3. Profile header contains #mpa_vOct22_CHOCOPhlAnSGB_202403
#      (catches accidental inclusion of stale vJan25 profiles)
#   4. Taxon row count -- flags LOW_TAXON (1-10 rows) separately from ZERO_TAXON,
#      since these are scientifically unusable but pass a naive non-empty check
#   5. Scans matching .err log for the two known silent-failure signatures:
#        - "BowTie2 output file detected" (stale .bz2 blocking re-run)
#        - "no sample found at line" (array size / sample-list mismatch)
#   6. Cross-checks output count against expected sample-list line count
#
# CHANGES FROM ORIGINAL pd_02_verify_metaphlan4.sh:
#   - All OUT_DIR paths updated pd_metaphlan4/ -> pd_metaphlan4_oct22/
#   - boktor_rumc SAMPLE_LIST updated to the filtered 139-sample list
#     (boktor2023rumc_samples_filtered.txt), matching the actual Oct22 submission
#   - Added VERSION_MISMATCH check (#3 above) — original script had no check
#     that profiles are actually the Oct22 index, which matters now that
#     vJan25 and vOct22 outputs exist as sibling directories

set -uo pipefail

DATASET="${1:-}"
if [ -z "$DATASET" ]; then
  echo "Usage: $0 <dataset_key>"
  echo "  dataset_key: bedarf2017 | boktor_rumc | boktor_tbc | clasen2024 | wallen2022 | mao2021"
  exit 1
fi

PD_BASE=/rds/projects/e/elhamsak-pd-thesis
AD_BASE=/rds/projects/e/elhamsak-ad-thesis
EXPECTED_HEADER="#mpa_vOct22_CHOCOPhlAnSGB_202403"

case "$DATASET" in
  bedarf2017)
    OUT_DIR=$PD_BASE/pd_metaphlan4_oct22/bedarf2017
    SAMPLE_LIST=$PD_BASE/pd_raw_data/bedarf2017/bedarf2017_merged_samples.txt
    LOG_DIR_PATTERN="$PD_BASE/logs/08_metaphlan4/mpa4_bedarf_oct22_*"
    ;;
  boktor_rumc)
    OUT_DIR=$PD_BASE/pd_metaphlan4_oct22/boktor2023/rumc
    SAMPLE_LIST=$PD_BASE/pd_raw_data/boktor2023/rumc/boktor2023rumc_samples_filtered.txt
    LOG_DIR_PATTERN="$PD_BASE/logs/08_metaphlan4/mpa4_boktor_rumc_oct22_*"
    ;;
  boktor_tbc)
    OUT_DIR=$PD_BASE/pd_metaphlan4_oct22/boktor2023/tbc
    SAMPLE_LIST=$PD_BASE/pd_raw_data/boktor2023/tbc/boktor2023tbc_samples.txt
    LOG_DIR_PATTERN="$PD_BASE/logs/08_metaphlan4/mpa4_boktor_tbc_oct22_*"
    ;;
  clasen2024)
    OUT_DIR=$PD_BASE/pd_metaphlan4_oct22/clasen2024
    SAMPLE_LIST=$PD_BASE/pd_raw_data/clasen2024/clasen2024_samples.txt
    LOG_DIR_PATTERN="$PD_BASE/logs/08_metaphlan4/mpa4_clasen_oct22_*"
    ;;
  wallen2022)
    OUT_DIR=$AD_BASE/pd_metaphlan4_oct22/wallen2022
    SAMPLE_LIST=$AD_BASE/pd_raw_data/wallen2022/wallen2022_samples.txt
    LOG_DIR_PATTERN="$PD_BASE/logs/08_metaphlan4/mpa4_wallen_oct22_*"
    ;;
  mao2021)
    OUT_DIR=$PD_BASE/pd_metaphlan4_oct22/mao2021
    SAMPLE_LIST=$PD_BASE/pd_raw_data/mao2021/mao2021_samples.txt
    LOG_DIR_PATTERN="$PD_BASE/logs/08_metaphlan4/mpa4_mao_oct22_*"
    ;;
  *)
    echo "Unknown dataset_key: $DATASET"
    exit 1
    ;;
esac

if [ ! -f "$SAMPLE_LIST" ]; then
  echo "ERROR: sample list not found: $SAMPLE_LIST"
  exit 1
fi

EXPECTED=$(wc -l < "$SAMPLE_LIST")

echo "=== $DATASET MetaPhlAn4 (Oct22) verification ==="
echo "OUT_DIR:      $OUT_DIR"
echo "SAMPLE_LIST:  $SAMPLE_LIST"
echo "Expected samples: $EXPECTED"
echo ""

MISSING_OUTPUT=0
ZERO_TAXON=0
LOW_TAXON=0
VERSION_MISMATCH=0
BZ2_COLLISION=0
SAMPLE_LIST_MISMATCH=0
OK=0

TAXON_COUNTS_FILE=$(mktemp)
LOW_TAXON_SAMPLES_FILE=$(mktemp)
MISSING_SAMPLES_FILE=$(mktemp)
VERSION_MISMATCH_FILE=$(mktemp)

while IFS= read -r sample; do
  [ -z "$sample" ] && continue

  PROFILE="$OUT_DIR/${sample}_metaphlan_profile.txt"

  if [ ! -f "$PROFILE" ]; then
    MISSING_OUTPUT=$((MISSING_OUTPUT+1))
    echo "$sample" >> "$MISSING_SAMPLES_FILE"
    continue
  fi

  # Check 3: confirm this is genuinely the Oct22-indexed profile, not a
  # stale vJan25 file sitting in the wrong directory
  if ! head -1 "$PROFILE" | grep -qF "$EXPECTED_HEADER"; then
    VERSION_MISMATCH=$((VERSION_MISMATCH+1))
    echo "$sample ($(head -1 "$PROFILE"))" >> "$VERSION_MISMATCH_FILE"
  fi

  # Count actual taxon rows: lines that are not comments/headers
  TAXON_ROWS=$(grep -vc "^#" "$PROFILE" 2>/dev/null || echo 0)

  echo "$TAXON_ROWS" >> "$TAXON_COUNTS_FILE"

  if [ "$TAXON_ROWS" -eq 0 ]; then
    ZERO_TAXON=$((ZERO_TAXON+1))
  elif [ "$TAXON_ROWS" -le 10 ]; then
    LOW_TAXON=$((LOW_TAXON+1))
    echo "$sample ($TAXON_ROWS rows)" >> "$LOW_TAXON_SAMPLES_FILE"
  else
    OK=$((OK+1))
  fi
done < "$SAMPLE_LIST"

echo "--- Output presence ---"
echo "Missing output files:     $MISSING_OUTPUT"
if [ "$MISSING_OUTPUT" -gt 0 ]; then
  echo "  (see: $MISSING_SAMPLES_FILE)"
fi
echo ""

echo "--- Database version check ---"
echo "Version mismatch (not Oct22 header): $VERSION_MISMATCH"
if [ "$VERSION_MISMATCH" -gt 0 ]; then
  echo "  Samples with unexpected header (see: $VERSION_MISMATCH_FILE):"
  cat "$VERSION_MISMATCH_FILE" | sed 's/^/    /'
  echo "  -> these are likely stale vJan25 profiles or corrupted headers -- do not merge as-is"
fi
echo ""

echo "--- Taxon row counts ---"
echo "OK (>10 taxon rows):      $OK"
echo "LOW_TAXON (1-10 rows):    $LOW_TAXON"
echo "ZERO_TAXON (0 rows):      $ZERO_TAXON"
if [ "$LOW_TAXON" -gt 0 ]; then
  echo "  Low-taxon samples (see: $LOW_TAXON_SAMPLES_FILE):"
  cat "$LOW_TAXON_SAMPLES_FILE" | sed 's/^/    /'
fi
if [ -s "$TAXON_COUNTS_FILE" ]; then
  echo "  Distribution: min=$(sort -n "$TAXON_COUNTS_FILE" | head -1)  max=$(sort -n "$TAXON_COUNTS_FILE" | tail -1)  median=$(sort -n "$TAXON_COUNTS_FILE" | awk '{a[NR]=$1} END{print a[int(NR/2)+1]}')"
fi
echo ""

echo "--- Known silent-failure signatures in logs ---"
if compgen -G "$LOG_DIR_PATTERN" > /dev/null; then
  for logdir in $LOG_DIR_PATTERN; do
    if [ -d "$logdir" ]; then
      BZ2_HITS=$(grep -l "BowTie2 output file detected" "$logdir"/*.err 2>/dev/null | wc -l)
      MISMATCH_HITS=$(grep -l "no sample found at line" "$logdir"/*.err 2>/dev/null | wc -l)
      BZ2_COLLISION=$((BZ2_COLLISION+BZ2_HITS))
      SAMPLE_LIST_MISMATCH=$((SAMPLE_LIST_MISMATCH+MISMATCH_HITS))
    fi
  done
fi
echo "Stale .bz2 collision errors:        $BZ2_COLLISION"
echo "  -> fix: find <output_dir> -name '*_bowtie2.bz2' -delete, then resubmit"
echo "Sample-list-mismatch errors:        $SAMPLE_LIST_MISMATCH"
echo "  -> fix: check --array range matches 'wc -l' of the sample list"
echo ""

echo "--- Cross-check ---"
ACTUAL_FILES=$(find "$OUT_DIR" -maxdepth 1 -name "*_metaphlan_profile.txt" 2>/dev/null | wc -l)
echo "Expected samples:          $EXPECTED"
echo "Profile files on disk:     $ACTUAL_FILES"
if [ "$ACTUAL_FILES" -ne "$EXPECTED" ]; then
  echo "  WARNING: file count does not match expected sample count"
  echo "  (if this is bedarf2017/mao2021/etc and count is HIGHER than expected,"
  echo "   check whether old vJan25 profiles are sitting in the same directory)"
fi
echo ""

echo "=== Verdict ==="
if [ "$MISSING_OUTPUT" -eq 0 ] && [ "$ZERO_TAXON" -eq 0 ] && [ "$LOW_TAXON" -eq 0 ] && [ "$VERSION_MISMATCH" -eq 0 ] && [ "$BZ2_COLLISION" -eq 0 ] && [ "$SAMPLE_LIST_MISMATCH" -eq 0 ]; then
  echo "CLEAN: all $EXPECTED samples have valid Oct22 profiles (>10 taxon rows), no known failure signatures detected."
  echo "-> safe to submit HUMAnN3 for this dataset."
else
  echo "ISSUES FOUND -- do not submit HUMAnN3 for this dataset until resolved:"
  [ "$MISSING_OUTPUT" -gt 0 ] && echo "  - $MISSING_OUTPUT samples missing output entirely"
  [ "$ZERO_TAXON" -gt 0 ] && echo "  - $ZERO_TAXON samples with zero taxon rows"
  [ "$LOW_TAXON" -gt 0 ] && echo "  - $LOW_TAXON samples with low (1-10) taxon rows -- needs manual review, not automatically excludable"
  [ "$VERSION_MISMATCH" -gt 0 ] && echo "  - $VERSION_MISMATCH samples not matching expected Oct22 header -- check for stale vJan25 files"
  [ "$BZ2_COLLISION" -gt 0 ] && echo "  - $BZ2_COLLISION stale-.bz2 collision failures -- clean up and resubmit"
  [ "$SAMPLE_LIST_MISMATCH" -gt 0 ] && echo "  - $SAMPLE_LIST_MISMATCH sample-list-mismatch failures -- fix array range and resubmit"
fi

rm -f "$TAXON_COUNTS_FILE" "$LOW_TAXON_SAMPLES_FILE" "$MISSING_SAMPLES_FILE" "$VERSION_MISMATCH_FILE"