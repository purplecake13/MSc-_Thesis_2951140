#!/bin/bash
#SBATCH --job-name=ena_pd_download
#SBATCH --account=elhamsak-pd-thesis            # your BlueBEAR project account
#SBATCH --qos=bbdefault                         # standard queue
#SBATCH --time=48:00:00                         # 48 hours — generous for large datasets
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1                       # 4 parallel wget processes
#SBATCH --mem=8G                                # wget is not memory-intensive
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/03_download/03_download_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/03_download/03_download_%j.err

# ============================================================
# USAGE:
#   sbatch 03_download_pd_wave.sh \
#       --script /path/to/ena_download_script.sh \
#       --outdir /rds/projects/e/elhamsak-pd-thesis/raw_data/Bedarf_2017
#
# ARGUMENTS:
#   --script   Path to the ENA-provided shell script containing wget commands
#   --outdir   Directory where FASTQ files will be saved
#
# EXAMPLE:
#   sbatch 03_download.sh \
#       --script /rds/projects/e/elhamsak-pd-thesis/scripts/03_download/pd/bedarf_2017_ena_download.sh \
#       --outdir /rds/projects/e/elhamsak-pd-thesis/pd_raw_data/bedarf2017
#
# NOTE ON STORAGE CHECKS:
#   Storage before/after download is checked on whichever RDS project --outdir
#   actually points to (derived automatically from the path), not hardcoded to
#   elhamsak-pd-thesis. This matters because some datasets (e.g. Mao 2021,
#   Wallen 2022) are downloaded to elhamsak-ad-thesis instead, for storage
#   balancing reasons -- the script's own log should report on the correct
#   project either way.
# ============================================================

# ── Parse arguments ─────────────────────────────────────────
ENA_SCRIPT=""
OUTDIR=""

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --script) ENA_SCRIPT="$2"; shift ;;
        --outdir) OUTDIR="$2";     shift ;;
        *) echo "ERROR: Unknown argument: $1"; exit 1 ;;
    esac
    shift
done

# ── Validate inputs ─────────────────────────────────────────
if [[ -z "$ENA_SCRIPT" ]]; then
    echo "ERROR: --script argument is required"
    exit 1
fi

if [[ -z "$OUTDIR" ]]; then
    echo "ERROR: --outdir argument is required"
    exit 1
fi

if [[ ! -f "$ENA_SCRIPT" ]]; then
    echo "ERROR: ENA script not found: $ENA_SCRIPT"
    exit 1
fi

# ── Derive which RDS project --outdir lives under, for storage checks ───────
# Matches the pattern /rds/projects/<letter>/<project-name>/... and extracts
# just the /rds/projects/<letter>/<project-name> prefix. Falls back to the
# pd-thesis project (the script's own historical default) if the path doesn't
# match the expected RDS layout, so this never silently checks an empty path.
STORAGE_CHECK_DIR=$(echo "$OUTDIR" | grep -oE '^/rds/projects/[^/]+/[^/]+' || true)
if [[ -z "$STORAGE_CHECK_DIR" ]]; then
    echo "WARNING: could not derive RDS project root from --outdir ($OUTDIR)." >&2
    echo "Falling back to elhamsak-pd-thesis for storage checks." >&2
    STORAGE_CHECK_DIR="/rds/projects/e/elhamsak-pd-thesis"
fi

# ── Setup ───────────────────────────────────────────────────
module purge
module load bear-apps/2024a
module load wget/1.25.0-GCCcore-13.3.0

# Create output directory if it does not exist
mkdir -p "$OUTDIR"

# Log file for download progress — separate from SLURM log
PROGRESS_LOG="$OUTDIR/01_download_progress.log"

echo "========================================"    | tee -a "$PROGRESS_LOG"
echo "ENA download started: $(date)"              | tee -a "$PROGRESS_LOG"
echo "ENA script:    $ENA_SCRIPT"                 | tee -a "$PROGRESS_LOG"
echo "Output dir:    $OUTDIR"                     | tee -a "$PROGRESS_LOG"
echo "Storage check dir: $STORAGE_CHECK_DIR"      | tee -a "$PROGRESS_LOG"
echo "SLURM job ID:  $SLURM_JOB_ID"              | tee -a "$PROGRESS_LOG"
echo "========================================"    | tee -a "$PROGRESS_LOG"

# ── Count total files to download ───────────────────────────
# Extract only wget lines (skip blank lines and comments)
TOTAL=$(grep -c "^wget" "$ENA_SCRIPT" || true)
echo "Total files to download: $TOTAL"            | tee -a "$PROGRESS_LOG"
echo ""                                            | tee -a "$PROGRESS_LOG"

# ── Check storage before starting ───────────────────────────
echo "Storage check before download:"             | tee -a "$PROGRESS_LOG"
df -h "$STORAGE_CHECK_DIR"                         | tee -a "$PROGRESS_LOG"
echo ""                                            | tee -a "$PROGRESS_LOG"

# ── Extract URLs and download in parallel ───────────────────
# Extract just the URLs from the wget lines, then download
# using GNU parallel via xargs for controlled parallelism
# -nc = no-clobber (skip already downloaded files)
# -q  = quiet mode (less noise; progress tracked via our log instead)
# -P  = number of parallel downloads (matches --cpus-per-task)

grep "^wget" "$ENA_SCRIPT" \
    | awk '{print $NF}' \
    | xargs -P 1 -I{} bash -c '     # change P 4 to P 1
        URL="{}"
        FILENAME=$(basename "$URL")
        OUTFILE="'"$OUTDIR"'/$FILENAME"

        if [[ -f "$OUTFILE" ]]; then
            echo "[SKIP] Already exists: $FILENAME" | tee -a "'"$PROGRESS_LOG"'"
        else
            echo "[START] Downloading: $FILENAME" | tee -a "'"$PROGRESS_LOG"'"
            wget -nc -q \
                --tries=5 \
                --wait=30 \
                --retry-connrefused \
                --timeout=120 \
                -O "$OUTFILE" \
                "$URL"
            EXIT_CODE=$?
            if [[ $EXIT_CODE -eq 0 ]]; then
                SIZE=$(du -sh "$OUTFILE" | cut -f1)
                echo "[DONE] $FILENAME ($SIZE)" | tee -a "'"$PROGRESS_LOG"'"
            else
                echo "[FAILED] $FILENAME (exit code $EXIT_CODE)" | tee -a "'"$PROGRESS_LOG"'"
                rm -f "$OUTFILE"   # remove incomplete file so retry works
            fi
        fi
    '

# ── Post-download summary ────────────────────────────────────
echo ""                                                        | tee -a "$PROGRESS_LOG"
echo "========================================"                | tee -a "$PROGRESS_LOG"
echo "Download finished: $(date)"                             | tee -a "$PROGRESS_LOG"
echo ""                                                        | tee -a "$PROGRESS_LOG"

# Count outcomes
N_DONE=$(grep -c "^\[DONE\]"   "$PROGRESS_LOG" || true)
N_SKIP=$(grep -c "^\[SKIP\]"   "$PROGRESS_LOG" || true)
N_FAIL=$(grep -c "^\[FAILED\]" "$PROGRESS_LOG" || true)

echo "Successfully downloaded: $N_DONE / $TOTAL"             | tee -a "$PROGRESS_LOG"
echo "Skipped (already exist): $N_SKIP"                      | tee -a "$PROGRESS_LOG"
echo "Failed:                  $N_FAIL"                      | tee -a "$PROGRESS_LOG"
echo ""                                                        | tee -a "$PROGRESS_LOG"

# List any failed files so you know exactly what to retry
if [[ $N_FAIL -gt 0 ]]; then
    echo "FAILED FILES:"                                       | tee -a "$PROGRESS_LOG"
    grep "^\[FAILED\]" "$PROGRESS_LOG"                        | tee -a "$PROGRESS_LOG"
    echo ""                                                    | tee -a "$PROGRESS_LOG"
    echo "To retry failed files only, resubmit this job."    | tee -a "$PROGRESS_LOG"
    echo "wget -nc will skip files that completed."           | tee -a "$PROGRESS_LOG"
fi

# ── Verify paired-end completeness ──────────────────────────
# For paired-end samples, check that both _1 and _2 files exist
# Single-end files (no suffix) are expected to have no pair
echo ""                                                        | tee -a "$PROGRESS_LOG"
echo "Paired-end completeness check:"                         | tee -a "$PROGRESS_LOG"

PAIRED_INCOMPLETE=0
# Find all _1 files and check for matching _2
for F1 in "$OUTDIR"/*_1.fastq.gz; do
    [[ -f "$F1" ]] || continue   # skip if no _1 files exist
    BASENAME="${F1%_1.fastq.gz}"
    F2="${BASENAME}_2.fastq.gz"
    if [[ ! -f "$F2" ]]; then
        echo "[WARNING] Missing R2 for: $(basename $F1)"      | tee -a "$PROGRESS_LOG"
        PAIRED_INCOMPLETE=$((PAIRED_INCOMPLETE + 1))
    fi
done

if [[ $PAIRED_INCOMPLETE -eq 0 ]]; then
    echo "All paired-end files have matching R1 and R2."      | tee -a "$PROGRESS_LOG"
else
    echo "WARNING: $PAIRED_INCOMPLETE samples missing R2."    | tee -a "$PROGRESS_LOG"
fi

# ── Final storage check ──────────────────────────────────────
echo ""                                                        | tee -a "$PROGRESS_LOG"
echo "Storage check after download:"                          | tee -a "$PROGRESS_LOG"
df -h "$STORAGE_CHECK_DIR"                                      | tee -a "$PROGRESS_LOG"
echo ""                                                        | tee -a "$PROGRESS_LOG"
echo "Files in output directory:"                             | tee -a "$PROGRESS_LOG"
ls -lh "$OUTDIR"/*.fastq.gz 2>/dev/null | awk '{print $5, $9}' | tee -a "$PROGRESS_LOG"
echo ""                                                        | tee -a "$PROGRESS_LOG"
echo "Total size of downloaded files:"                        | tee -a "$PROGRESS_LOG"
du -sh "$OUTDIR"                                               | tee -a "$PROGRESS_LOG"

echo ""                                                        | tee -a "$PROGRESS_LOG"
echo "Progress log saved to: $PROGRESS_LOG"                   | tee -a "$PROGRESS_LOG"

# Exit with failure if any downloads failed, so SLURM marks the job accordingly
if [[ $N_FAIL -gt 0 ]]; then
    exit 1
fi

exit 0