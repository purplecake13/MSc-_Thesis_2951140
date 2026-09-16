#!/bin/bash
# pd_03_verify_humann3.sh
# Usage: ./pd_03_verify_humann3.sh <dataset_key> <job_name_prefix> [low_band_pct] [high_band_pct]
#
# dataset_key:       bedarf2017 | boktor_rumc | boktor_tbc | clasen2024 | wallen2022 | mao2021
# job_name_prefix:   the --job-name used at submission, e.g. humann3_bedarf_oct22
#                     (logs are flat: logs/09_humann3/{job_name_prefix}_{jobid}_{taskid}.out,
#                      so this must match exactly what you passed to --job-name)
#
# Checks, per sample (matched to its array task by LINE NUMBER in the sample
# list -- task N == line N, same convention already used by pd_03_humann3.sh
# and by the "no sample found at line" mismatch check):
#   1. genefamilies.tsv, pathabundance.tsv, pathcoverage.tsv all exist
#   2. Each file is non-trivial size (catches truncated/stub output)
#   3. genefamilies.tsv has a real number of rows (not just UNMAPPED)
#   4. Unaligned-read % pulled from the matching .out log (NOT computed by
#      dividing READS_UNMAPPED by a naive tsv sum -- READS_UNMAPPED is a raw
#      read count while every other row is CPM-normalized; dividing them
#      gives a meaningless number, as confirmed empirically this session)
#   5. Flags samples with unaligned % outside an expected band (default
#      15-45%) for manual review -- NOT an auto-exclusion, just a flag
#   6. Scans .err log for common humann3 failure signatures (missing
#      taxprofile version line, database path errors, OOM/timeout)
#   7. Cross-checks output count against expected sample-list line count
#
# LOG NAMING: logs/09_humann3/{job_name_prefix}_{SLURM_ARRAY_JOB_ID}_{task}.out
# Task N corresponds to LINE N of the sample list. If a dataset's array was
# submitted with --array=1-N against the same sample list used elsewhere in
# the pipeline, this mapping is correct; if you re-ran with a custom
# --array range or a reordered list, this will misattribute -- check the
# actual sbatch command used.
#
# NOTE ON UNITS: do not compute "%% unmapped" from genefamilies.tsv directly.
# READS_UNMAPPED there is a raw read count; every other row is Adjusted CPM.
# The correct figure is "Unaligned reads after translated alignment" from
# the job's .out log, which this script uses.

set -uo pipefail

DATASET="${1:-}"
JOB_NAME_PREFIX="${2:-}"
LOW_BAND="${3:-15}"
HIGH_BAND="${4:-45}"

if [ -z "$DATASET" ] || [ -z "$JOB_NAME_PREFIX" ]; then
  echo "Usage: $0 <dataset_key> <job_name_prefix> [low_band_pct] [high_band_pct]"
  echo "  dataset_key: bedarf2017 | boktor_rumc | boktor_tbc | clasen2024 | wallen2022 | mao2021"
  echo "  job_name_prefix: exact --job-name used at submission, e.g. humann3_bedarf_oct22"
  echo "  low_band_pct / high_band_pct: expected unaligned-read %% range (default 15-45)"
  exit 1
fi

PD_BASE=/rds/projects/e/elhamsak-pd-thesis
AD_BASE=/rds/projects/e/elhamsak-ad-thesis
LOG_DIR=$PD_BASE/logs/09_humann3

case "$DATASET" in
  bedarf2017)
    OUT_DIR=$PD_BASE/pd_humann3_oct22/bedarf2017
    SAMPLE_LIST=$PD_BASE/pd_raw_data/bedarf2017/bedarf2017_merged_samples.txt
    ;;
  boktor_rumc)
    OUT_DIR=$PD_BASE/pd_humann3_oct22/boktor2023/rumc
    SAMPLE_LIST=$PD_BASE/pd_raw_data/boktor2023/rumc/boktor2023rumc_samples_filtered.txt
    ;;
  boktor_tbc)
    OUT_DIR=$PD_BASE/pd_humann3_oct22/boktor2023/tbc
    SAMPLE_LIST=$PD_BASE/pd_raw_data/boktor2023/tbc/boktor2023tbc_samples.txt
    ;;
  clasen2024)
    OUT_DIR=$PD_BASE/pd_humann3_oct22/clasen2024
    SAMPLE_LIST=$PD_BASE/pd_raw_data/clasen2024/clasen2024_samples.txt
    ;;
  wallen2022)
    OUT_DIR=$AD_BASE/pd_humann3_oct22/wallen2022
    SAMPLE_LIST=$AD_BASE/pd_raw_data/wallen2022/wallen2022_samples.txt
    ;;
  mao2021)
    OUT_DIR=$PD_BASE/pd_humann3_oct22/mao2021
    SAMPLE_LIST=$PD_BASE/pd_raw_data/mao2021/mao2021_samples.txt
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

echo "=== $DATASET HUMAnN3 verification ==="
echo "OUT_DIR:          $OUT_DIR"
echo "SAMPLE_LIST:      $SAMPLE_LIST"
echo "LOG_DIR:          $LOG_DIR"
echo "JOB_NAME_PREFIX:  $JOB_NAME_PREFIX"
echo "Expected samples: $EXPECTED"
echo "Unaligned%% flag band: <${LOW_BAND}%% or >${HIGH_BAND}%% flagged for review"
echo ""

MISSING_FILES=0
EMPTY_GENEFAMILIES=0
UNALIGNED_NOT_FOUND=0
OUT_OF_BAND=0
LOG_NOT_FOUND=0
ERR_SIGNATURE=0
OK=0

MISSING_FILES_LOG=$(mktemp)
EMPTY_LOG=$(mktemp)
UNALIGNED_MISSING_LOG=$(mktemp)
OUT_OF_BAND_LOG=$(mktemp)
LOG_NOT_FOUND_LOG=$(mktemp)
ERR_SIGNATURE_LOG=$(mktemp)
UNALIGNED_VALUES_FILE=$(mktemp)

TASK_INDEX=0
while IFS= read -r sample; do
  [ -z "$sample" ] && continue
  TASK_INDEX=$((TASK_INDEX+1))

  GF="$OUT_DIR/${sample}_concat_2_genefamilies.tsv"
  RX="$OUT_DIR/${sample}_concat_3_reactions.tsv"
  PA="$OUT_DIR/${sample}_concat_4_pathabundance.tsv"
  # PC="$OUT_DIR/${sample}_pathcoverage.tsv"


  MISSING_THIS=0
  for f in "$GF" "$PA" "$RX"; do
    if [ ! -s "$f" ] || [ "$(stat -c%s "$f" 2>/dev/null || echo 0)" -lt 100 ]; then
      MISSING_THIS=1
    fi
  done
  if [ "$MISSING_THIS" -eq 1 ]; then
    MISSING_FILES=$((MISSING_FILES+1))
    echo "$sample (task $TASK_INDEX)" >> "$MISSING_FILES_LOG"
    continue
  fi

  ROW_COUNT=$(grep -vc "^#" "$GF" 2>/dev/null || echo 0)
  if [ "$ROW_COUNT" -le 1 ]; then
    EMPTY_GENEFAMILIES=$((EMPTY_GENEFAMILIES+1))
    echo "$sample (task $TASK_INDEX, $ROW_COUNT rows)" >> "$EMPTY_LOG"
    continue
  fi

  # OUT_LOG=$(find "$LOG_DIR" -maxdepth 1 -name "${JOB_NAME_PREFIX}_*_${TASK_INDEX}.out" 2>/dev/null | head -1)
  # ERR_LOG=$(find "$LOG_DIR" -maxdepth 1 -name "${JOB_NAME_PREFIX}_*_${TASK_INDEX}.err" 2>/dev/null | head -1)
  
  OUT_LOG=$(find "$LOG_DIR" -maxdepth 1 -name "${JOB_NAME_PREFIX}_*_${TASK_INDEX}.out" -printf '%T@ %p\n' 2>/dev/null | sort -rn | head -1 | cut -d' ' -f2-)
  ERR_LOG=$(find "$LOG_DIR" -maxdepth 1 -name "${JOB_NAME_PREFIX}_*_${TASK_INDEX}.err" -printf '%T@ %p\n' 2>/dev/null | sort -rn | head -1 | cut -d' ' -f2-)

  if [ -z "$OUT_LOG" ]; then
    LOG_NOT_FOUND=$((LOG_NOT_FOUND+1))
    echo "$sample (task $TASK_INDEX) -- no matching .out log for pattern ${JOB_NAME_PREFIX}_*_${TASK_INDEX}.out" >> "$LOG_NOT_FOUND_LOG"
  else
    UNALIGNED_PCT=$(grep -iE "Unaligned reads after translated alignment" "$OUT_LOG" 2>/dev/null | grep -oE '[0-9]+\.[0-9]+' | head -1)
    if [ -z "$UNALIGNED_PCT" ]; then
      UNALIGNED_NOT_FOUND=$((UNALIGNED_NOT_FOUND+1))
      echo "$sample (task $TASK_INDEX, log: $OUT_LOG)" >> "$UNALIGNED_MISSING_LOG"
    else
      echo "$sample $UNALIGNED_PCT" >> "$UNALIGNED_VALUES_FILE"
      BAND_CHECK=$(awk -v v="$UNALIGNED_PCT" -v lo="$LOW_BAND" -v hi="$HIGH_BAND" 'BEGIN{print (v<lo || v>hi) ? "1":"0"}')
      if [ "$BAND_CHECK" -eq 1 ]; then
        OUT_OF_BAND=$((OUT_OF_BAND+1))
        echo "$sample (task $TASK_INDEX, ${UNALIGNED_PCT}%)" >> "$OUT_OF_BAND_LOG"
      fi
    fi
  fi

  if [ -n "$ERR_LOG" ] && [ -s "$ERR_LOG" ]; then
    if grep -qE "does not contain the database version|CANCELLED|OUT_OF_MEMORY|DUE TO TIME LIMIT|No such file or directory" "$ERR_LOG" 2>/dev/null; then
      ERR_SIGNATURE=$((ERR_SIGNATURE+1))
      echo "$sample (task $TASK_INDEX, $ERR_LOG)" >> "$ERR_SIGNATURE_LOG"
    fi
  fi

  OK=$((OK+1))
done < "$SAMPLE_LIST"

echo "--- Output file presence ---"
echo "Missing/truncated outputs: $MISSING_FILES"
[ "$MISSING_FILES" -gt 0 ] && cat "$MISSING_FILES_LOG" | sed 's/^/    /'
echo ""

echo "--- Content check ---"
echo "Empty/near-empty genefamilies (<=1 row): $EMPTY_GENEFAMILIES"
[ "$EMPTY_GENEFAMILIES" -gt 0 ] && cat "$EMPTY_LOG" | sed 's/^/    /'
echo ""

echo "--- Log matching ---"
echo "Samples with no matching .out log found: $LOG_NOT_FOUND"
if [ "$LOG_NOT_FOUND" -gt 0 ]; then
  cat "$LOG_NOT_FOUND_LOG" | sed 's/^/    /'
  echo "  -> check JOB_NAME_PREFIX matches the actual --job-name used, and that"
  echo "     task-to-line mapping is correct (custom --array range or reordered list breaks this)"
fi
echo ""

echo "--- Unaligned-read %% (from .out log, correct source) ---"
echo "Samples with unaligned%% found: $((OK - UNALIGNED_NOT_FOUND - LOG_NOT_FOUND))"
echo "Samples with log found but no unaligned%% line: $UNALIGNED_NOT_FOUND"
[ "$UNALIGNED_NOT_FOUND" -gt 0 ] && cat "$UNALIGNED_MISSING_LOG" | sed 's/^/    /'
echo "Samples outside ${LOW_BAND}-${HIGH_BAND}%% band: $OUT_OF_BAND"
if [ "$OUT_OF_BAND" -gt 0 ]; then
  cat "$OUT_OF_BAND_LOG" | sed 's/^/    /'
fi
if [ -s "$UNALIGNED_VALUES_FILE" ]; then
  echo "  Distribution: min=$(sort -k2 -n "$UNALIGNED_VALUES_FILE" | head -1 | awk '{print $2}')%%  max=$(sort -k2 -n "$UNALIGNED_VALUES_FILE" | tail -1 | awk '{print $2}')%%"
fi
echo ""

echo "--- Known error signatures in .err logs ---"
echo "Error-log signature hits: $ERR_SIGNATURE"
if [ "$ERR_SIGNATURE" -gt 0 ]; then
  cat "$ERR_SIGNATURE_LOG" | sed 's/^/    /'
  echo "  -> 'does not contain the database version': taxprofile-dir points at a cleaned/rebuilt"
  echo "     file instead of raw MetaPhlAn4 output -- point --taxprofile-dir at the raw output dir"
  echo "  -> OOM/TIME LIMIT: increase --mem / --time in the sbatch header and resubmit that sample"
fi
echo ""

echo "--- Cross-check ---"
ACTUAL_FILES=$(find "$OUT_DIR" -maxdepth 1 -name "*_genefamilies.tsv" 2>/dev/null | wc -l)
echo "Expected samples:         $EXPECTED"
echo "genefamilies.tsv on disk: $ACTUAL_FILES"
[ "$ACTUAL_FILES" -ne "$EXPECTED" ] && echo "  WARNING: file count does not match expected sample count"
echo ""

echo "=== Verdict ==="
if [ "$MISSING_FILES" -eq 0 ] && [ "$EMPTY_GENEFAMILIES" -eq 0 ] && [ "$LOG_NOT_FOUND" -eq 0 ] && [ "$ERR_SIGNATURE" -eq 0 ] && [ "$ACTUAL_FILES" -eq "$EXPECTED" ]; then
  echo "CLEAN: all $EXPECTED samples have valid HUMAnN3 output, no known failure signatures detected."
  if [ "$OUT_OF_BAND" -gt 0 ]; then
    echo "NOTE: $OUT_OF_BAND sample(s) have unaligned%% outside the ${LOW_BAND}-${HIGH_BAND}%% band -- review before"
    echo "      treating as routine, but this does not block proceeding."
  fi
else
  echo "ISSUES FOUND -- review before treating this dataset's functional profiles as final:"
  [ "$MISSING_FILES" -gt 0 ] && echo "  - $MISSING_FILES samples missing/truncated output files"
  [ "$EMPTY_GENEFAMILIES" -gt 0 ] && echo "  - $EMPTY_GENEFAMILIES samples with empty/near-empty genefamilies output"
  [ "$LOG_NOT_FOUND" -gt 0 ] && echo "  - $LOG_NOT_FOUND samples with no matching .out log found (check JOB_NAME_PREFIX)"
  [ "$ERR_SIGNATURE" -gt 0 ] && echo "  - $ERR_SIGNATURE error-log signature hits -- see above"
  [ "$ACTUAL_FILES" -ne "$EXPECTED" ] && echo "  - output file count mismatch vs expected sample count"
fi

rm -f "$MISSING_FILES_LOG" "$EMPTY_LOG" "$UNALIGNED_MISSING_LOG" "$OUT_OF_BAND_LOG" "$LOG_NOT_FOUND_LOG" "$ERR_SIGNATURE_LOG" "$UNALIGNED_VALUES_FILE"
