#!/bin/bash
# verify_fastp_job.sh
# Reusable verification for fastp trimming or fastp-merge SLURM jobs.
# Checks: job states -> error logs -> file count/integrity -> JSON sanity.
#
# Usage:
#   ./verify_fastp_job.sh <dataset> <job_id> <expected_samples> <output_dir> <mode>
#
#   dataset          e.g. yildirim2022
#   job_id           SLURM job ID (parent array job ID, e.g. 48867688)
#   expected_samples number of samples expected (array size)
#   output_dir       directory containing *.trimmed.fastq.gz or merge output
#   mode             "paired" (R1/R2, e.g. trimming or merge output) or
#                     "single" (one file per sample, e.g. Chen2022)
#
# Example:
#   /rds/projects/e/elhamsak-pd-thesis/scripts/verify_fastp_job.sh yildirim2022 48867688 98 /rds/projects/e/elhamsak-ad-thesis/ad_trimmed/yildirim2022/with_sliding_window paired
#
# Adjust LOG_GLOB below if your SLURM log naming convention differs.

set -uo pipefail

DATASET="${1:?dataset name required}"
JOBID="${2:?SLURM job id required}"
EXPECTED="${3:?expected sample count required}"
OUTDIR="${4:?output directory required}"
MODE="${5:-paired}"   # paired | single

echo "=============================================="
echo " Verifying: $DATASET  (job $JOBID)"
echo "=============================================="

# --- Step 1: job states -------------------------------------------------
echo
echo "--- [1/4] SLURM job states (non-COMPLETED tasks) ---"
BAD_STATES=$(sacct -j "$JOBID" --format=JobID,State,ExitCode 2>/dev/null | grep -v COMPLETED | grep -v "^JobID" | grep -v "^----")
if [ -z "$BAD_STATES" ]; then
  echo "OK: all array tasks COMPLETED."
else
  echo "ISSUES FOUND:"
  echo "$BAD_STATES"
fi

# --- Step 2: error logs ---------------------------------------------------
echo
echo "--- [2/4] Error log scan ---"
# Adjust this glob if your logging directory structure differs.
LOG_GLOB="/rds/projects/e/elhamsak-pd-thesis/logs/*/*${JOBID}*.err"
MATCHES=$(grep -il -E "error|Error|ERROR|Traceback|No such file|cannot|not found" $LOG_GLOB 2>/dev/null)
if [ -z "$MATCHES" ]; then
  echo "OK: no error keywords found in matching logs (glob: $LOG_GLOB)"
  echo "NOTE: if this seems too easy, double-check the glob actually matched files:"
  ls $LOG_GLOB 2>/dev/null | head -3
else
  echo "ISSUES FOUND in:"
  echo "$MATCHES"
fi

# --- Step 3: file count + gzip integrity ----------------------------------
echo
echo "--- [3/4] File count and gzip integrity ---"
if [ "$MODE" == "single" ]; then
  FOUND=$(find "$OUTDIR" -maxdepth 1 -iname "*.trimmed.fastq.gz" 2>/dev/null | wc -l)
  echo "Expected: $EXPECTED   Found: $FOUND"
else
  FOUND_R1=$(find "$OUTDIR" -maxdepth 1 -iname "*_1.trimmed.fastq.gz" -o -iname "*_R1*.fastq.gz" -o -iname "*_1_merged.fastq.gz" -o -iname "*_merged.fastq.gz" 2>/dev/null | wc -l)
  FOUND_R2=$(find "$OUTDIR" -maxdepth 1 -iname "*_2.trimmed.fastq.gz" -o -iname "*_R2*.fastq.gz" 2>/dev/null | wc -l)
  echo "Expected: $EXPECTED   Found R1/merged-pattern: $FOUND_R1   Found R2: $FOUND_R2"
  echo "NOTE: file-naming patterns vary by step (trim vs merge) -- verify these counts look right for this stage."
fi

BAD_FILES=""
for f in "$OUTDIR"/*.fastq.gz; do
  [ -e "$f" ] || continue
  # use the gzip line when checking 16s, else use the line without gzip. 
  # if [ ! -s "$f" ] || ! gzip -t "$f" 2>/dev/null; then 
  if [ ! -s "$f" ] ; then # use when checking shotgun
    BAD_FILES="${BAD_FILES}${f}\n"
  fi
done

if [ -z "$BAD_FILES" ]; then
  echo "OK: no corrupt or empty files found."
else
  echo "CORRUPT/EMPTY FILES:"
  echo -e "$BAD_FILES"
  echo "$BAD_FILES" | sed -E 's#.*/##; s/_[12]?\.(trimmed\.)?fastq\.gz$//' | sort -u > "/tmp/${DATASET}_bad_samples.txt"
  echo "-> Affected sample IDs written to /tmp/${DATASET}_bad_samples.txt"
fi

# --- Step 4: JSON sanity check (near-zero reads) --------------------------
echo
echo "--- [4/4] fastp JSON sanity check (flagging near-zero read counts) ---"
JSON_COUNT=$(find "$OUTDIR"/reports -maxdepth 1 -iname "*.json" 2>/dev/null | wc -l)
echo "JSON reports found: $JSON_COUNT (expected ~$EXPECTED)"

if [ "$JSON_COUNT" -gt 0 ]; then
  FLAGGED=0
  for j in "$OUTDIR"/reports/*.json; do
    [ -e "$j" ] || continue
    python3 -c "
import json, sys
try:
    d = json.load(open('$j'))
    s = d.get('summary', {})
    before = s.get('before_filtering', {}).get('total_reads', 0)
    after = s.get('after_filtering', {}).get('total_reads', 0)
    name = '$j'.split('/')[-1]
    if before == 0:
        print(f'  FLAG (zero input): {name}  before={before} after={after}')
    elif after == 0:
        print(f'  FLAG (zero output): {name}  before={before} after={after}')
    elif after < before * 0.05:
        print(f'  FLAG (>95% loss): {name}  before={before} after={after}')
except Exception as e:
    print(f'  FLAG (unreadable JSON): $j  ({e})')
"
  done
  echo "(no output above these lines for a given file = looked fine)"
else
  echo "No JSON reports found -- cannot run sanity check. Check output directory path."
fi

echo
echo "=============================================="
echo " Done. Review any FLAG/ISSUE lines above before"
echo " treating $DATASET as verified for this stage."
echo "=============================================="