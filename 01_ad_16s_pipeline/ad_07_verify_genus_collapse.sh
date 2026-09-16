#!/bin/bash
# ad_07_verify_genus_collapse.sh
# Verifies qiime taxa collapse (genus-level) output across all AD 16S datasets,
# using the actual output structure produced by ad_07_genus_collapse.sh:
#   <genus_level_dir>/genus_table.qza
#   <genus_level_dir>/genus_table_export/feature-table.tsv
# and the genera/sample counts the collapse script itself prints to its .out log.
#
# Checks per dataset:
#   1. SLURM job state (COMPLETED)
#   2. genus_table.qza and feature-table.tsv both exist and are non-empty
#   3. Sample count in the genus-level TSV matches the pre-collapse table
#   4. Genera count is sensible: reduced vs pre-collapse feature count, and
#      not implausibly low (<5)
#   5. Feature IDs in the TSV actually look genus-level (g__ or D_5__ present)
#
# USAGE: ./ad_07_verify_genus_collapse.sh
# Requires: qiime2 env active (for exporting the PRE-collapse table only;
# the post-collapse side is read directly from the TSV, no qiime/biom needed).
set -uo pipefail

module purge
module load bear-apps/2023a
module load QIIME2/2025.4

WORK_TMP="/tmp/ad_genus_verify_$$"
mkdir -p "$WORK_TMP"
trap 'rm -rf "$WORK_TMP"' EXIT

AD_ROOT="/rds/projects/e/elhamsak-ad-thesis/ad_qiime2"
LOG_ROOT="/rds/projects/e/elhamsak-pd-thesis/logs/taxonomy"

# name | pre-collapse table.qza (relative to AD_ROOT) | genus_level dir (relative) | job_id
read -r -d '' DATASETS <<'EOF'
cirstea2022    cirstea2022/table.qza      cirstea2022/genus_level      51804822
yamashiro2024  yamashiro2024/table.qza    yamashiro2024/genus_level    51804828
tran2019       tran2019/table.qza         tran2019/genus_level         51804845
binyinli2019   binyinli2019/table.qza     binyinli2019/genus_level     51804859
ling2021       ling2021/table.qza         ling2021/genus_level         51804860
liu2019        liu2019/table.qza          liu2019/genus_level          51804868
ueda2021       ueda2021/table.qza         ueda2021/genus_level         51804871
yildirim2022   yildirim2022/table.qza     yildirim2022/genus_level     51804873
zhuang2018     zhuang2018/table.qza       zhuang2018/genus_level       51804894
EOF

printf "%-14s %-9s %-8s %-13s %-10s %-16s %s\n" \
  "DATASET" "JOBSTATE" "FILES" "SAMPLES(match)" "GENERA" "GENUS_LABELS" "OVERALL"
printf '%.0s-' {1..95}; echo

FAIL_COUNT=0

while read -r NAME TABLE_REL GENUS_DIR_REL JOBID; do
  [[ -z "$NAME" ]] && continue

  TABLE="${AD_ROOT}/${TABLE_REL}"
  GENUS_QZA="${AD_ROOT}/${GENUS_DIR_REL}/genus_table.qza"
  GENUS_TSV="${AD_ROOT}/${GENUS_DIR_REL}/genus_table_export/feature-table.tsv"

  # --- 1. SLURM state ---
  STATE=$(sacct -j "$JOBID" --format=State -X -n 2>/dev/null | tr -d '[:space:]')
  [[ -z "$STATE" ]] && STATE="NOTFOUND"

  # --- 2. file existence ---
  FILES_OK="OK"
  if [[ ! -s "$GENUS_QZA" || ! -s "$GENUS_TSV" ]]; then
    FILES_OK="MISSING"
    printf "%-14s %-9s %-8s %-13s %-10s %-16s %s\n" "$NAME" "$STATE" "$FILES_OK" "-" "-" "-" "FAIL"
    FAIL_COUNT=$((FAIL_COUNT+1))
    continue
  fi

  # --- genera count and sample count directly from the TSV ---
  # biom convert --to-tsv format: line1 = "# Constructed from biom file", line2 = header (#OTU ID <samples...>)
  N_SAMPLES_POST=$(sed -n '2p' "$GENUS_TSV" | awk -F'\t' '{print NF-1}')
  N_GENERA=$(tail -n +3 "$GENUS_TSV" | wc -l)

  # --- 3. pre-collapse sample count, via qiime export (only pre-side needed) ---
  PRE_DIR="${WORK_TMP}/${NAME}_pre"
  qiime tools export --input-path "$TABLE" --output-path "$PRE_DIR" >/dev/null 2>&1
  N_SAMPLES_PRE="?"
  N_FEATS_PRE="?"
  if [[ -f "${PRE_DIR}/feature-table.biom" ]]; then
    biom summarize-table -i "${PRE_DIR}/feature-table.biom" -o "${PRE_DIR}/summary.txt" 2>/dev/null
    N_SAMPLES_PRE=$(grep "Num samples" "${PRE_DIR}/summary.txt" | grep -oE '[0-9,]+' | tr -d ',')
    N_FEATS_PRE=$(grep "Num observations" "${PRE_DIR}/summary.txt" | grep -oE '[0-9,]+' | tr -d ',')
  fi

  SAMPLE_MATCH="?"
  if [[ "$N_SAMPLES_PRE" != "?" ]]; then
    if [[ "$N_SAMPLES_PRE" == "$N_SAMPLES_POST" ]]; then SAMPLE_MATCH="OK"; else SAMPLE_MATCH="MISMATCH"; fi
  fi

  # --- 4. genera count sanity vs pre-collapse feature count ---
  GENERA_FLAG="OK"
  if [[ "$N_FEATS_PRE" != "?" ]]; then
    if [[ "$N_GENERA" -ge "$N_FEATS_PRE" ]]; then GENERA_FLAG="NO_REDUCTION"; fi
  fi
  if [[ "$N_GENERA" -lt 5 ]]; then GENERA_FLAG="TOO_FEW"; fi

  # --- 5. genus-level label check: first column of TSV should contain g__ or D_5__ ---
  LABEL_CHECK="FAIL"
  FIRST_COL=$(tail -n +3 "$GENUS_TSV" | awk -F'\t' '{print $1}')
  if echo "$FIRST_COL" | grep -qiE 'g__|D_5__'; then
    LABEL_CHECK="OK"
  fi

  OVERALL="OK"
  if [[ "$STATE" != "COMPLETED" || "$SAMPLE_MATCH" == "MISMATCH" || "$GENERA_FLAG" != "OK" || "$LABEL_CHECK" != "OK" ]]; then
    OVERALL="CHECK"
    FAIL_COUNT=$((FAIL_COUNT+1))
  fi

  printf "%-14s %-9s %-8s %-13s %-10s %-16s %s\n" \
    "$NAME" "$STATE" "$FILES_OK" "${N_SAMPLES_PRE}->${N_SAMPLES_POST}(${SAMPLE_MATCH})" \
    "${N_GENERA}(${GENERA_FLAG})" "$LABEL_CHECK" "$OVERALL"

done <<< "$DATASETS"

echo ""
if [[ "$FAIL_COUNT" -eq 0 ]]; then
  echo "All datasets passed: job completed, files present, sample counts match pre-collapse, genera count reduced sensibly, genus-level labels confirmed."
else
  echo "${FAIL_COUNT} dataset(s) flagged CHECK/FAIL above - review before proceeding to cross-study merge."
  echo "Common causes: taxonomy.qza built against a different feature-table.qza than the one collapsed (feature ID mismatch),"
  echo "wrong --level passed to the collapse script (check the .out log's 'Genera retained' line), or for Yamashiro2024"
  echo "specifically, a high proportion of unclassified/family-only entries from its 81.2%-resolution V1-V2 classifier -"
  echo "inspect that dataset's actual TSV rows manually before assuming genuine failure:"
  echo "  cut -f1 <genus_level_dir>/genus_table_export/feature-table.tsv | tail -n +3 | less"
fi

echo ""
echo "Cross-check against each collapse job's own reported counts (it logs 'Genera retained' and 'Samples in table' directly):"
echo "  grep -E 'Genera retained|Samples in table' ${LOG_ROOT}/genus_collapse_*_<jobid>.out"
