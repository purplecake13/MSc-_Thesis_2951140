#!/bin/bash
#SBATCH --job-name=filter_yamashiro
#SBATCH --account=elhamsak-pd-thesis
#SBATCH --qos=bbdefault
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/02c_filter_yamashiro_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/02c_filter_yamashiro_%j.err

# --- Paths ---
METADATA="/rds/projects/e/elhamsak-ad-thesis/ad_raw_data/yamashiro2024/yamashiro2024_metadata.txt"
ENA_REPORT="/rds/projects/e/elhamsak-ad-thesis/ad_raw_data/yamashiro2024/yamashiro2024_ena_report.tsv"
SCRIPT_DIR="/rds/projects/e/elhamsak-pd-thesis/scripts/03_download/ad"

# Output keep list path matching your standard naming conventions
KEEP_LIST="/rds/projects/e/elhamsak-ad-thesis/ad_raw_data/yamashiro2024/yamashiro2024_sample_list_no_MCI.txt"

echo "========================================"
echo "Generating Keep List for Yamashiro 2024"
echo "========================================"

# Temporary file to store S-number to Group mapping from metadata
MAPPING_TMP=$(mktemp)

# Step 1: Parse metadata file into an S-number -> Group mapping array layout
# Assumes columns: No. (which matches S-number suffix) and Group
awk 'NR>1 {print "S"$1"\t"$3}' "$METADATA" > "$MAPPING_TMP"

# Step 2: Read ENA Report, find Sample Title and Run Accession, cross-reference mapping, 
# and keep accessions where Group != "MCI"
awk -F'\t' -v mapping="$MAPPING_TMP" '
BEGIN {
    # Load metadata mapping into memory array: meta["S17"] = "MCI"
    while ((getline < mapping) > 0) {
        meta[$1] = $2
    }
    close(mapping)
}
NR==1 {
    # Dynamically find column indices from ENA report header
    for(i=1; i<=NF; i++) {
        if(tolower($i) ~ /sample[_-]title/) title_col=i
        if(tolower($i) ~ /run[_-]accession/) acc_col=i
    }
    next
}
{
    # Extract the S-number inside the parenthesis, e.g., "sample(S17)" -> "S17"
    match($title_col, /\(S[0-9]+\)/)
    if (RSTART > 0) {
        s_num = substr($title_col, RSTART+1, RLENGTH-2)
        group = meta[s_num]
        
        # Keep if the sample is Control or Dementia (exclude MCI)
        if (group != "MCI" && group != "") {
            print $acc_col
        }
    }
}' "$ENA_REPORT" | sort -u > "$KEEP_LIST"

echo "Keep list written to: $KEEP_LIST"
echo "Total samples kept (Control + Dementia): $(wc -l < "$KEEP_LIST")"

rm -f "$MAPPING_TMP"
