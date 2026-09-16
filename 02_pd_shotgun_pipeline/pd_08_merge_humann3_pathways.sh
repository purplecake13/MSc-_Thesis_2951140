#!/usr/bin/env bash
# pd_09_merge_humann3_pathways.sh

module purge
module load bear-apps/2023a
module load humann/4.0.0a1-foss-2023a

EXCLUDE_LIST=/rds/projects/e/elhamsak-pd-thesis/pd_humann3_oct22/functional_analysis_exclusions.txt
ABX_FILE=/rds/projects/e/elhamsak-pd-thesis/cross_platform/excluded_samples_antibiotic.txt
MERGE_DIR=/rds/projects/e/elhamsak-pd-thesis/pd_humann3_merged

mkdir -p "$(dirname "$EXCLUDE_LIST")"
rm -rf "$MERGE_DIR/pathabundance_input"
mkdir -p "$MERGE_DIR/pathabundance_input"

# 1. Initialize exclusion file with the 8 oral contamination samples from clasen
cat << 'EOF' > "$EXCLUDE_LIST"
ERR15003744
ERR15003808
ERR15003814
ERR15003915
ERR15003916
ERR15003922
ERR15003954
ERR15003857
EOF

# 2. Append bare IDs from the antibiotic exclusion list
if [ -f "$ABX_FILE" ]; then
    # Extracts the final token after the last underscore (handles study_id_sample_id)
    sed -E 's/.*_//' "$ABX_FILE" >> "$EXCLUDE_LIST"
    sort -u "$EXCLUDE_LIST" -o "$EXCLUDE_LIST"
    echo "Combined exclusion list built: $(wc -l < "$EXCLUDE_LIST") unique sample IDs."
fi

# 3. Append Boktor severe-tier functional-analysis-only exclusions
#    (documented in sample_exclusions_log.tsv, reason=elevated_humann3_unaligned_pct_functional_only)
SAMPLE_EXCLUSIONS_LOG=/rds/projects/e/elhamsak-pd-thesis/sample_exclusions_log.tsv
if [ -f "$SAMPLE_EXCLUSIONS_LOG" ]; then
    awk -F'\t' '$4 == "elevated_humann3_unaligned_pct_functional_only" {print $1}' "$SAMPLE_EXCLUSIONS_LOG" >> "$EXCLUDE_LIST"
    sort -u "$EXCLUDE_LIST" -o "$EXCLUDE_LIST"
    echo "After Boktor functional-only exclusions: $(wc -l < "$EXCLUDE_LIST") unique sample IDs."
else
    echo "WARNING: $SAMPLE_EXCLUSIONS_LOG not found — Boktor severe-tier functional exclusions NOT applied."
fi

# 3. Process pathabundance files
# 3. Process pathabundance files
for f in /rds/projects/e/elhamsak-pd-thesis/pd_humann3_oct22/*/*_pathabundance.tsv \
         /rds/projects/e/elhamsak-pd-thesis/pd_humann3_oct22/*/*/*_pathabundance.tsv \
         /rds/projects/e/elhamsak-ad-thesis/pd_humann3_oct22/*/*_pathabundance.tsv; do
  sample=$(basename "$f" | sed -E 's/_concat_[0-9]+_pathabundance\.tsv//')
  if grep -qxF "$sample" "$EXCLUDE_LIST"; then
    echo "Excluding $sample from functional merge"
    continue
  fi
  ln -sf "$f" "$MERGE_DIR/pathabundance_input/${sample}_pathabundance.tsv"
done

# 4. Join tables via HUMAnN3
humann_join_tables -i "$MERGE_DIR/pathabundance_input" \
  -o "$MERGE_DIR/pd_pathabundance_joined.tsv" --file_name pathabundance

# 5. Normalise to CPM
humann_renorm_table -i "$MERGE_DIR/pd_pathabundance_joined.tsv" \
  -o "$MERGE_DIR/pd_pathabundance_cpm.tsv" --units cpm

# 6. Split stratified / unstratified
humann_split_stratified_table -i "$MERGE_DIR/pd_pathabundance_cpm.tsv" -o "$MERGE_DIR/"
