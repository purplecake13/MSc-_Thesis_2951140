#!/bin/bash
# ad_04_dada2_denoise_verify.sh
# 4-step verification for ad_04_dada2_denoise.sh jobs:
# job state -> error logs -> file integrity -> retention stats
# Usage: ./ad_04_dada2_denoise_verify.sh <dataset_name> <job_id> <output_dir>

module purge
module load bear-apps/2023a
module load QIIME2/2025.4

DATASET=$1
JOBID=$2
OUTDIR=$3
LOGDIR="/rds/projects/e/elhamsak-pd-thesis/logs/08_ad_qiime2/dada2_denoise"

if [[ -z "$DATASET" || -z "$JOBID" || -z "$OUTDIR" ]]; then
  echo "Usage: $0 <dataset_name> <job_id> <output_dir>"
  exit 1
fi

echo "=================================================="
echo "DADA2 verification: $DATASET (job $JOBID)"
echo "=================================================="

# --- Step 1: job state ---
echo -e "\n[1/4] Job state (sacct)"
sacct -j "$JOBID" --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS -P

STATE=$(sacct -j "$JOBID" --format=State -P --noheader | head -n1 | tr -d ' ')
if [[ "$STATE" != "COMPLETED" ]]; then
  echo "  -> WARNING: job state is '$STATE', not COMPLETED. Check error log below before trusting anything downstream."
fi

# --- Step 2: error log ---
echo -e "\n[2/4] Error log scan"
ERRLOG="${LOGDIR}/dada2_${DATASET}_${JOBID}.err"
if [[ -f "$ERRLOG" ]]; then
  ERRSIZE=$(stat -c%s "$ERRLOG")
  echo "  Error log: $ERRLOG (${ERRSIZE} bytes)"
  echo "  Non-boilerplate matches (Traceback, Error, CRITICAL, Plugin error):"
  grep -E "Traceback|Error|CRITICAL|Plugin error|Debug info" "$ERRLOG" | head -n 20
  [[ $(grep -cE "Traceback|Error|CRITICAL|Plugin error" "$ERRLOG") -eq 0 ]] && echo "  -> none found"

  # Surface the referenced QIIME2 debug log, if any, and its content
  DEBUGLOG=$(grep -oE "/[^ ]*qiime2-q2cli-err-[^ ]*\.log" "$ERRLOG" | head -n1)
  if [[ -n "$DEBUGLOG" && -f "$DEBUGLOG" ]]; then
    echo -e "\n  --- Referenced QIIME2 debug log: $DEBUGLOG ---"
    cat "$DEBUGLOG"
    echo "  --- end debug log ---"
  elif [[ -n "$DEBUGLOG" ]]; then
    echo "  -> Referenced debug log $DEBUGLOG not found on disk (may have been cleaned up)"
  fi
else
  echo "  -> WARNING: expected error log not found at $ERRLOG"
fi

# --- Step 3: file integrity ---
# NOTE: actual script output filenames confirmed from disk (Aug 2026 session):
#   table.qza, repseqs.qza (no hyphen), dada2stats.qza (+ .qzv companion)
echo -e "\n[3/4] Output artifact integrity"
for f in table.qza repseqs.qza dada2stats.qza; do
  FPATH="${OUTDIR}/${f}"
  if [[ -f "$FPATH" ]]; then
    SIZE=$(stat -c%s "$FPATH")
    if unzip -tq "$FPATH" > /dev/null 2>&1; then
      echo "  OK    $f  (${SIZE} bytes, zip-valid)"
    else
      echo "  FAIL  $f  (${SIZE} bytes) -- zip integrity check failed, likely truncated/corrupt"
    fi
  else
    echo "  MISSING  $f"
  fi
done

# --- Step 4: retention / read-count sanity ---
echo -e "\n[4/4] Denoising retention stats"
STATS_QZA="${OUTDIR}/dada2stats.qza"
STATS_DIR="${OUTDIR}/stats_export"

if [[ -f "$STATS_QZA" ]]; then
  mkdir -p "$STATS_DIR"
  qiime tools export --input-path "$STATS_QZA" --output-path "$STATS_DIR" 2>/dev/null

  if [[ -f "${STATS_DIR}/stats.tsv" ]]; then
    echo "  stats.tsv exported to ${STATS_DIR}/stats.tsv"
    echo ""
    column -t -s $'\t' "${STATS_DIR}/stats.tsv" | grep -v "^#q2:types" 2>/dev/null || cat "${STATS_DIR}/stats.tsv"
    echo ""
    echo "  Inspect the 'percentage of input non-chimeric' column above manually;"
    echo "  flag any sample under ~70% for review before trusting downstream merges."
  else
    echo "  -> WARNING: export produced no stats.tsv, inspect ${STATS_DIR} manually"
  fi
else
  echo "  -> WARNING: dada2stats.qza not found, cannot check retention"
fi

echo -e "\n=================================================="
echo "Verification complete for $DATASET"
echo "=================================================="