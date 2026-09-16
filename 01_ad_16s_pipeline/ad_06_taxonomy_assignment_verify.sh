#!/bin/bash
# ad_06_taxonomy_verify.sh
# 4-step verification for ad_06_taxonomy_assignment.sh jobs, plus
# taxonomy-quality diagnostics specific to this stage.
# Usage: ./ad_06_taxonomy_verify.sh <dataset_name> <job_id> <output_dir>

DATASET=$(echo "$1" | tr -d '[:space:]')
JOBID=$(echo "$2" | tr -d '[:space:]')
OUTDIR=$3
LOGDIR="/rds/projects/e/elhamsak-pd-thesis/logs/taxonomy"

if [[ -z "$DATASET" || -z "$JOBID" || -z "$OUTDIR" ]]; then
  echo "Usage: $0 <dataset_name> <job_id> <output_dir>"
  exit 1
fi

module purge
module load bear-apps/2023a
module load QIIME2/2025.4

echo "=================================================="
echo "Taxonomy verification: $DATASET (job $JOBID)"
echo "=================================================="

# --- Step 1: job state ---
echo -e "\n[1/4] Job state (sacct)"
sacct -j "$JOBID" --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS -P
STATE=$(sacct -j "$JOBID" --format=State -P --noheader | head -n1 | tr -d ' ')
[[ "$STATE" != "COMPLETED" ]] && echo "  -> WARNING: job state is '$STATE', not COMPLETED."

# --- Step 2: error log ---
echo -e "\n[2/4] Error log scan"
ERRLOG="${LOGDIR}/taxonomy_${DATASET}_${JOBID}.err"
if [[ -f "$ERRLOG" ]]; then
  echo "  Error log: $ERRLOG ($(stat -c%s "$ERRLOG") bytes)"
  grep -E "Traceback|Error|CRITICAL|Plugin error" "$ERRLOG" | head -n 20
  [[ $(grep -cE "Traceback|Error|CRITICAL|Plugin error" "$ERRLOG") -eq 0 ]] && echo "  -> none found"
  DEBUGLOG=$(grep -oE "/[^ ]*qiime2-q2cli-err-[^ ]*\.log" "$ERRLOG" | head -n1)
  if [[ -n "$DEBUGLOG" && -f "$DEBUGLOG" ]]; then
    echo -e "\n  --- Referenced debug log: $DEBUGLOG ---"
    cat "$DEBUGLOG"
  fi
else
  echo "  -> WARNING: expected error log not found at $ERRLOG"
fi

# --- Step 3: file integrity ---
echo -e "\n[3/4] Output artifact integrity"
for f in taxonomy.qza taxonomy.qzv taxonomy_export/taxonomy.tsv; do
  FPATH="${OUTDIR}/${f}"
  if [[ -f "$FPATH" ]]; then
    SIZE=$(stat -c%s "$FPATH")
    if [[ "$f" == *.qza || "$f" == *.qzv ]]; then
      if unzip -tq "$FPATH" > /dev/null 2>&1; then
        echo "  OK    $f  (${SIZE} bytes, zip-valid)"
      else
        echo "  FAIL  $f  (${SIZE} bytes) -- zip integrity check failed"
      fi
    else
      echo "  OK    $f  (${SIZE} bytes)"
    fi
  else
    echo "  MISSING  $f"
  fi
done

# --- Step 4: taxonomy quality diagnostics ---
echo -e "\n[4/4] Taxonomy quality diagnostics"
TSV="${OUTDIR}/taxonomy_export/taxonomy.tsv"

if [[ -f "$TSV" ]]; then
  N_TOTAL=$(tail -n +2 "$TSV" | wc -l)
  N_UNASSIGNED=$(tail -n +2 "$TSV" | awk -F'\t' '$2 ~ /^Unassigned/' | wc -l)
  N_DOMAIN_ONLY=$(tail -n +2 "$TSV" | awk -F'\t' '$2 ~ /^d__/ && $2 !~ /p__/' | wc -l)
  N_GENUS=$(tail -n +2 "$TSV" | awk -F'\t' '$2 ~ /g__[A-Za-z]/' | wc -l)

  echo "  Total ASVs classified:        $N_TOTAL"
  printf "  Unassigned:                   %d (%.1f%%)\n" "$N_UNASSIGNED" "$(awk "BEGIN{print ($N_UNASSIGNED/$N_TOTAL)*100}")"
  printf "  Domain-level only (weak):     %d (%.1f%%)\n" "$N_DOMAIN_ONLY" "$(awk "BEGIN{print ($N_DOMAIN_ONLY/$N_TOTAL)*100}")"
  printf "  Genus-level call achieved:    %d (%.1f%%)\n" "$N_GENUS" "$(awk "BEGIN{print ($N_GENUS/$N_TOTAL)*100}")"

  # Confidence score distribution (3rd column, if present)
  if awk -F'\t' 'NR==1{print $3}' "$TSV" | grep -qi confidence; then
    echo ""
    echo "  Confidence score distribution:"
    tail -n +2 "$TSV" | awk -F'\t' '{print $3}' | sort -n | awk '
      { a[NR]=$1; sum+=$1 }
      END {
        if (NR==0) { print "    no confidence values found"; exit }
        printf "    min=%.3f  median=%.3f  mean=%.3f  max=%.3f\n", a[1], a[int(NR/2)], sum/NR, a[NR]
      }'
    N_LOWCONF=$(tail -n +2 "$TSV" | awk -F'\t' '$3+0 < 0.7' | wc -l)
    printf "    ASVs below 0.7 confidence: %d (%.1f%%)\n" "$N_LOWCONF" "$(awk "BEGIN{print ($N_LOWCONF/$N_TOTAL)*100}")"
  fi

  echo ""
  if [[ "$N_GENUS" -eq 0 ]]; then
    echo "  *** WARNING: zero genus-level calls. Classifier/region mismatch likely. ***"
  elif (( $(awk "BEGIN{print ($N_UNASSIGNED/$N_TOTAL)*100 > 20}") )); then
    echo "  *** WARNING: >20% Unassigned. Investigate before trusting this dataset's taxonomy. ***"
  else
    echo "  No red flags detected by automated thresholds -- still worth a manual skim of taxonomy.qzv."
  fi
else
  echo "  -> WARNING: taxonomy.tsv not found, cannot run diagnostics"
fi

echo -e "\n=================================================="
echo "Verification complete for $DATASET"
echo "=================================================="