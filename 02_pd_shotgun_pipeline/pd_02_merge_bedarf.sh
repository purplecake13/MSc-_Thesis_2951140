#!/bin/bash
#SBATCH --job-name=merge_bedarf_post_qc
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/01_merge_bedarf/%x_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/01_merge_bedarf/%x_%j.err

# =============================================================================
# Merge Bowtie2 host-removed runs into biological samples for MetaPhlAn4
# =============================================================================

set -euo pipefail

# ---- Paths ------------------------------------------------------------------
RAW_DIR="/rds/projects/e/elhamsak-pd-thesis/pd_host_removed/bedarf2017"
OUT_DIR="/rds/projects/e/elhamsak-pd-thesis/pd_host_removed/bedarf2017_merged"
MAP_FILE="/rds/projects/e/elhamsak-pd-thesis/pd_raw_data/bedarf2017/bedarf2017_run_to_sample.tsv"
SAMPLE_LIST="/rds/projects/e/elhamsak-pd-thesis/pd_raw_data/bedarf2017/bedarf2017_merged_samples.txt"

mkdir -p "$OUT_DIR"
mkdir -p "$(dirname "$SAMPLE_LIST")"

echo "=== Starting Bedarf 2017 Post-QC Concatenation ==="

# Extract unique sample names from Column 2 (No header to skip)
SAMPLES=$(awk '{print $2}' "$MAP_FILE" | sort -u)

for SAMPLE in $SAMPLES; do
    echo "Processing sample: $SAMPLE"
    
    # Get all runs for this sample from Column 1
    RUNS=$(awk -v samp="$SAMPLE" '$2==samp {print $1}' "$MAP_FILE")
    
    R1_FILES=""
    R2_FILES=""
    
    # Build the list of files to concatenate using the Bowtie2 output suffix
    for RUN in $RUNS; do
        R1_FILES="$R1_FILES $RAW_DIR/${RUN}_host_removed_1.fastq.gz"
        R2_FILES="$R2_FILES $RAW_DIR/${RUN}_host_removed_2.fastq.gz"
    done
    
    # Concatenate directly into the output directory
    cat $R1_FILES > "${OUT_DIR}/${SAMPLE}_host_removed_1.fastq.gz"
    cat $R2_FILES > "${OUT_DIR}/${SAMPLE}_host_removed_2.fastq.gz"
done

echo "=== Concatenation Complete ==="

# ---- Generate the Sample List for MetaPhlAn4 -------------------------------
echo "Generating sample list for MetaPhlAn4..."
ls "${OUT_DIR}"/*_host_removed_1.fastq.gz | xargs -n 1 basename | sed 's/_host_removed_1.fastq.gz//' > "$SAMPLE_LIST"

echo "Sample list saved to: $SAMPLE_LIST"
echo "Total biological samples generated: $(wc -l < "$SAMPLE_LIST")"
