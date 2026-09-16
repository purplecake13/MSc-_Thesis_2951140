#!/bin/bash
# =============================================================================
#SBATCH --job-name=fastqc_array
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
# Log paths are set dynamically below after parsing args — see LOG_DIR block.
# (--output and --error cannot reference shell variables, so we redirect manually.)
# =============================================================================
# fastqc_array.sh — Run FastQC on a directory of FASTQ files (SLURM array job)
# =============================================================================
#
# USAGE:
#   First, count your files:
#     ls /path/to/your/data/*.fastq.gz | wc -l
#
#   Then submit with the correct array size:
#     sbatch --array=1-<N>%10 fastqc_array.sh \
#            --data-dir /rds/projects/e/elhamsak-pd-thesis/raw_data/AD_16S/Tran_2019 \
#            --output-dir /rds/projects/e/elhamsak-pd-thesis/qc_reports/Tran_2019
#
#   Replace <N> with the number of .fastq.gz files in your data directory.
#   The %10 limits concurrent tasks to 10 — adjust if needed.
#
# EXAMPLES (one per dataset):
#   sbatch --array=1-56%10 fastqc_array.sh \
#          --data-dir .../AD_16S/Tran_2019 \
#          --output-dir .../qc_reports/Tran_2019
#
#   sbatch --array=1-59%10 fastqc_array.sh \
#          --data-dir .../PD_shotgun/Bedarf_2017 \
#          --output-dir .../qc_reports/Bedarf_2017
#
# NOTES:
#   - Matches both .fastq and .fastq.gz files automatically.
#   - Logs go to a timestamped subdirectory so re-runs don't overwrite history.
#   - Requires FastQC to be available as a BlueBEAR module (check: module spider fastqc).

# =============================================================================
# PARSE ARGUMENTS
# =============================================================================
DATA_DIR=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir)   DATA_DIR="$2";   shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        *)
            echo "ERROR: Unknown argument: $1"
            echo "Usage: sbatch --array=1-N%10 fastqc_array.sh --data-dir <path> --output-dir <path>"
            exit 1
            ;;
    esac
done

# =============================================================================
# VALIDATE INPUTS
# =============================================================================
if [[ -z "$DATA_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: Both --data-dir and --output-dir are required."
    echo "Usage: sbatch --array=1-N%10 fastqc_array.sh --data-dir <path> --output-dir <path>"
    exit 1
fi

if [[ ! -d "$DATA_DIR" ]]; then
    echo "ERROR: Data directory does not exist: $DATA_DIR"
    exit 1
fi

# =============================================================================
# SET UP TIMESTAMPED LOG DIRECTORY
# Each submission gets its own folder: logs/04_fastqc/YYYYMMDD_HHMMSS_<jobid>/
# =============================================================================
TIMESTAMP=$(date +"%Y%m%d_%H%M")
LOG_DIR="/rds/projects/e/elhamsak-pd-thesis/logs/04_fastqc/${TIMESTAMP}_${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}"
mkdir -p "$LOG_DIR"

# Redirect stdout and stderr to per-task log files within the timestamped folder
exec > "${LOG_DIR}/task_${SLURM_ARRAY_TASK_ID}.out" \
     2> "${LOG_DIR}/task_${SLURM_ARRAY_TASK_ID}.err"

echo "====================================================="
echo "FastQC Array Job"
echo "Job ID       : ${SLURM_ARRAY_JOB_ID}"
echo "Task ID      : ${SLURM_ARRAY_TASK_ID}"
echo "Node         : $(hostname)"
echo "Start time   : $(date)"
echo "Data dir     : $DATA_DIR"
echo "Output dir   : $OUTPUT_DIR"
echo "Log dir      : $LOG_DIR"
echo "====================================================="

# =============================================================================
# LOAD MODULE
# Check available version with: module spider fastqc
# =============================================================================
module purge
module load bear-apps/2023a
module load FastQC/0.11.9-Java-11 
echo "FastQC version: $(fastqc --version)"

# =============================================================================
# BUILD FILE LIST
# Matches both .fastq and .fastq.gz files (fasterq-dump can produce either)
# =============================================================================
mapfile -t FILES < <(find "$DATA_DIR" -maxdepth 1 \( -name "*.fastq" -o -name "*.fastq.gz" \) | sort)

TOTAL=${#FILES[@]}
echo "Total FASTQ files found: $TOTAL"

if [[ $TOTAL -eq 0 ]]; then
    echo "ERROR: No .fastq or .fastq.gz files found in $DATA_DIR"
    echo "Check the path and that fasterq-dump has been run."
    exit 1
fi

# Map SLURM task ID (1-indexed) to 0-indexed bash array
INDEX=$(( SLURM_ARRAY_TASK_ID - 1 ))

if [[ $INDEX -ge $TOTAL ]]; then
    echo "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} exceeds number of files ($TOTAL)."
    echo "Re-submit with --array=1-${TOTAL}%10"
    exit 1
fi

CURRENT_FILE="${FILES[$INDEX]}"
echo "Processing file: $CURRENT_FILE"

# =============================================================================
# CREATE OUTPUT DIRECTORY AND RUN FASTQC
# =============================================================================
mkdir -p "$OUTPUT_DIR"

fastqc \
    --threads 2 \
    --outdir "$OUTPUT_DIR" \
    "$CURRENT_FILE"

EXIT_CODE=$?

if [[ $EXIT_CODE -ne 0 ]]; then
    echo "ERROR: FastQC failed on $CURRENT_FILE with exit code $EXIT_CODE"
    exit $EXIT_CODE
fi

echo "-----------------------------------------------------"
echo "FastQC complete for: $CURRENT_FILE"
echo "Output written to  : $OUTPUT_DIR"
echo "End time           : $(date)"
echo "====================================================="