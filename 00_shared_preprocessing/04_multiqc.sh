#!/bin/bash
# =============================================================================

#SBATCH --job-name=multiqc
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
# Logs are redirected dynamically below (can't use shell vars in #SBATCH --output)

# =============================================================================
# multiqc.sh — Aggregate FastQC reports for one or more datasets
# =============================================================================
#
# USAGE:
#   # Run on a single dataset's QC folder:
#   sbatch multiqc.sh \
#          --qc-dir /rds/projects/e/elhamsak-pd-thesis/qc_reports/Tran_2019 \
#          --output-dir /rds/projects/e/elhamsak-pd-thesis/qc_reports/multiqc_Tran_2019 \
#          --report-name Tran_2019_multiqc
#
#   # Run across ALL AD datasets at once (recommended after all FastQC jobs finish):
#   sbatch multiqc.sh \
#          --qc-dir /rds/projects/e/elhamsak-pd-thesis/qc_reports \
#          --output-dir /rds/projects/e/elhamsak-pd-thesis/qc_reports/multiqc_AD_combined \
#          --report-name AD_combined_multiqc
#
#   # Run after FastQC array job completes (chain jobs with dependency):
#   FASTQC_JOB=$(sbatch --array=1-56%10 --parsable fastqc_array.sh \
#                       --data-dir .../Tran_2019 --output-dir .../qc_reports/Tran_2019)
#   sbatch --dependency=afterok:$FASTQC_JOB multiqc.sh \
#          --qc-dir .../qc_reports/Tran_2019 \
#          --output-dir .../qc_reports/multiqc_Tran_2019 \
#          --report-name Tran_2019_multiqc
#
# NOTES:
#   - MultiQC is loaded as a BlueBEAR module. Check the exact module name with:
#       module spider MultiQC
#     and update the MULTIQC_MODULE variable below if needed.
#   - The HTML report will be in --output-dir/<report-name>.html
#   - Logs go to a timestamped subdirectory so re-runs don't overwrite history.

# =============================================================================
# LOAD MULTIQC MODULE
# Check the exact name
# =============================================================================
module purge
module load bear-apps/2024a
module load MultiQC/1.32-foss-2024a

if ! command -v multiqc &> /dev/null; then
    echo "ERROR: multiqc not found after loading module '$MULTIQC_MODULE'."
    echo "Check available versions with: module spider MultiQC"
    exit 1
fi

echo "MultiQC version: $(multiqc --version)"

# =============================================================================
# PARSE ARGUMENTS
# =============================================================================
QC_DIR=""
OUTPUT_DIR=""
REPORT_NAME="multiqc_report"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --qc-dir)      QC_DIR="$2";      shift 2 ;;
        --output-dir)  OUTPUT_DIR="$2";  shift 2 ;;
        --report-name) REPORT_NAME="$2"; shift 2 ;;
        *)
            echo "ERROR: Unknown argument: $1"
            echo "Usage: sbatch multiqc.sh --qc-dir <path> --output-dir <path> [--report-name <name>]"
            exit 1
            ;;
    esac
done

# =============================================================================
# VALIDATE INPUTS
# =============================================================================
if [[ -z "$QC_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: Both --qc-dir and --output-dir are required."
    echo "Usage: sbatch multiqc.sh --qc-dir <path> --output-dir <path> [--report-name <name>]"
    exit 1
fi

if [[ ! -d "$QC_DIR" ]]; then
    echo "ERROR: QC directory does not exist: $QC_DIR"
    echo "Has the FastQC array job finished?"
    exit 1
fi

# =============================================================================
# SET UP TIMESTAMPED LOG DIRECTORY
# Each submission gets its own folder: logs/multiqc/YYYYMMDD_HHMM_<jobid>/
# =============================================================================
TIMESTAMP=$(date +"%Y%m%d_%H%M")
LOG_DIR="/rds/projects/e/elhamsak-pd-thesis/logs/05_multiqc/${TIMESTAMP}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
mkdir -p "$LOG_DIR"

exec > "${LOG_DIR}/multiqc.out" \
     2> "${LOG_DIR}/multiqc.err"

echo "====================================================="
echo "MultiQC Aggregation Job"
echo "Job ID       : ${SLURM_JOB_ID}"
echo "Node         : $(hostname)"
echo "Start time   : $(date)"
echo "QC input dir : $QC_DIR"
echo "Output dir   : $OUTPUT_DIR"
echo "Report name  : $REPORT_NAME"
echo "Log dir      : $LOG_DIR"
echo "====================================================="

# =============================================================================
# CHECK FOR FASTQC OUTPUT FILES IN THE INPUT DIRECTORY
# =============================================================================
FASTQC_FILE_COUNT=$(find "$QC_DIR" -name "*_fastqc.zip" | wc -l)
echo "FastQC zip files found: $FASTQC_FILE_COUNT"

if [[ $FASTQC_FILE_COUNT -eq 0 ]]; then
    echo "ERROR: No FastQC output files (*_fastqc.zip) found in $QC_DIR"
    echo "Check that the FastQC array job completed successfully before running MultiQC."
    exit 1
fi

# =============================================================================
# CREATE OUTPUT DIRECTORY AND RUN MULTIQC
# =============================================================================
mkdir -p "$OUTPUT_DIR"

multiqc \
    "$QC_DIR" \
    --outdir "$OUTPUT_DIR" \
    --filename "${REPORT_NAME}.html" \
    --force \
    --verbose

EXIT_CODE=$?

if [[ $EXIT_CODE -ne 0 ]]; then
    echo "ERROR: MultiQC failed with exit code $EXIT_CODE"
    exit $EXIT_CODE
fi

echo "-----------------------------------------------------"
echo "MultiQC complete."
echo "HTML report : ${OUTPUT_DIR}/${REPORT_NAME}.html"
echo "End time    : $(date)"
echo "====================================================="

# =============================================================================
# SUMMARY OF KEY FILES TO DOWNLOAD AND REVIEW
# =============================================================================
echo ""
echo "FILES TO DOWNLOAD FOR LOCAL REVIEW:"
echo "  scp <username>@bluebear.bham.ac.uk:${OUTPUT_DIR}/${REPORT_NAME}.html ."
echo ""
echo "WHAT TO CHECK IN THE REPORT:"
echo "  - Per Base Sequence Quality : want green (Q>28) across the full read length"
echo "  - Per Sequence Quality Score: want the peak at Q>30"
echo "  - Adapter Content           : if adapters detected, fastp will handle removal"
echo "  - Sequence Length Distribution: note if reads are all the same length (expected for Illumina)"
echo "  - % Duplicate reads         : 16S amplicon data is expected to have HIGH duplication — this is NORMAL"
echo "  - Per Base N Content        : should be near zero; high N content = sequencing failure"