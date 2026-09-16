#!/bin/bash
#SBATCH --job-name=filter_verify
#SBATCH --account=elhamsak-pd-thesis
#SBATCH --qos=bbdefault
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/rds/projects/e/elhamsak-pd-thesis/logs/02a_filter_studies/02_filter_verify_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-pd-thesis/logs/02a_filter_studies/02_filter_verify_%j.err

# ============================================================
# USAGE:
#   sbatch 02_filter_studies.sh \
#       --dataset cirstea2022 \
#       --ena-report /path/to/cirstea_2022_ena_report.tsv \
#       --ena-script /path/to/cirstea_2022_ena_download.sh \
#       --outdir /path/to/downloaded/files \
#       --sample-type fecal
#
# ARGUMENTS:
#   --dataset      Short name for the dataset (used in output filenames)
#   --ena-report   Path to the ENA metadata report TSV downloaded from ENA browser
#   --ena-script   Path to the original ENA wget download shell script
#   --outdir       Directory where FASTQ files were (or will be) downloaded
#   --sample-type  Keyword to filter on in the filter column (default: fecal)
#   --filter-col   Column name to filter on (default: "experiment alias")
#                  Use "scientific name" for Clasen_2024, which lacks an
#                  experiment alias column but encodes sample type there instead.
#                  The match is case-insensitive and partial (substring).
#   --exclude      If set, INVERTS the match: keep rows that DO NOT contain
#                  --sample-type in --filter-col, instead of keeping rows
#                  that DO contain it. Use this for exclusion filters like
#                  removing MCI samples (e.g. --sample-type "MCI" --exclude),
#                  as opposed to inclusion filters like gut-only or human-only
#                  (default behavior, no --exclude flag).
#
# WHAT THIS SCRIPT DOES:
#   Step 1 — Filter the ENA report to rows matching (or not matching, if
#            --exclude is set) --sample-type in --filter-col
#   Step 2 — Extract run accessions from the filtered report
#   Step 3 — Filter the original ENA wget script to matching accessions only
#   Step 4 — Compare expected accessions vs files actually on disk
#   Step 5 — Report what is missing, what is extra, and what is complete
#
# CHAINING MULTIPLE FILTERS (e.g. human-only AND not-MCI for the same dataset):
#   Each run starts from the --ena-report/--ena-script you point it at and
#   produces ONE filtered report + ONE filtered script as output. To apply a
#   SECOND filter on top of the first, point the second run's --ena-report
#   and --ena-script at the FIRST run's filtered OUTPUT files, not the
#   original unfiltered input. For example:
#
#     # Run 1: keep human samples only
#     sbatch 02_filter_studies.sh --dataset mystudy \
#       --ena-report mystudy_ena_report.tsv \
#       --ena-script mystudy_ena_download.sh \
#       --outdir ... --filter-col "host" --sample-type "Homo sapiens"
#     # -> produces mystudy_ena_report_Homo sapiens.tsv, mystudy_ena_download_Homo sapiens.sh
#
#     # Run 2: from run 1's output, also exclude MCI
#     sbatch 02_filter_studies.sh --dataset mystudy \
#       --ena-report "mystudy_ena_report_Homo sapiens.tsv" \
#       --ena-script "mystudy_ena_download_Homo sapiens.sh" \
#       --outdir ... --filter-col "diagnosis" --sample-type "MCI" --exclude
#
#   Running both filters against the ORIGINAL input independently (rather
#   than chaining) would each apply only ONE filter, not both combined.
#
# EXAMPLES:
#
#   Cirstea 2022 — uses "experiment alias" column (default), filter on "fecal":
#   sbatch 02_filter_studies.sh \
#       --dataset cirstea2022 \
#       --ena-report /rds/projects/e/elhamsak-ad-thesis/ad_raw_data/cirstea2022/cirstea_2022_ena_report.tsv \
#       --ena-script /rds/projects/e/elhamsak-pd-thesis/scripts/01_download/ad/cirstea_2022_ena_download.sh \
#       --outdir /rds/projects/e/elhamsak-ad-thesis/ad_raw_data/cirstea2022 \
#       --sample-type fecal \
#       --filter-col "experiment alias"
#
#   Clasen 2024 — uses "scientific name" column, filter on "fecal":
#   sbatch 02_filter_studies.sh \
#       --dataset clasen2024 \
#       --ena-report /rds/projects/e/elhamsak-pd-thesis/pd_raw_data/clasen_2024_ena_report.tsv \
#       --ena-script /rds/projects/e/elhamsak-pd-thesis/scripts/01_download/clasen_2024_ena_download.sh \
#       --outdir /rds/projects/e/elhamsak-pd-thesis/pd_raw_data/clasen2024 \
#       --sample-type fecal \
#       --filter-col "scientific name"
#
#   Human-only filter (study with mixed human/mouse samples):
#   sbatch 02_filter_studies.sh \
#       --dataset mystudy \
#       --ena-report ... --ena-script ... --outdir ... \
#       --filter-col "host" --sample-type "Homo sapiens"
#
#   MCI exclusion (Yamashiro 2024 and new AD studies):
#   sbatch 02_filter_studies.sh \
#       --dataset yamashiro2024 \
#       --ena-report ... --ena-script ... --outdir ... \
#       --filter-col "diagnosis" --sample-type "MCI" --exclude
# ============================================================


# ════════════════════════════════════════════════════════════
# ARGUMENT PARSING
# Read command-line arguments into named variables.
# All arguments are key-value pairs (--flag value), except
# --exclude which is a boolean switch (no value).
# ════════════════════════════════════════════════════════════

DATASET=""
ENA_REPORT=""
ENA_SCRIPT=""
OUTDIR=""
SAMPLE_TYPE="fecal"             # default filter keyword (matched as substring)
FILTER_COL="experiment alias"   # default column to filter on
                                 # override with --filter-col "scientific name" for Clasen_2024
EXCLUDE=0                       # 0 = keep rows matching SAMPLE_TYPE (default, original behavior)
                                 # 1 = keep rows NOT matching SAMPLE_TYPE (e.g. MCI exclusion)

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --dataset)     DATASET="$2";     shift ;;
        --ena-report)  ENA_REPORT="$2";  shift ;;
        --ena-script)  ENA_SCRIPT="$2";  shift ;;
        --outdir)      OUTDIR="$2";      shift ;;
        --sample-type) SAMPLE_TYPE="$2"; shift ;;
        --filter-col)  FILTER_COL="$2";  shift ;;
        --exclude)     EXCLUDE=1 ;;
        *) echo "ERROR: Unknown argument: $1"; exit 1 ;;
    esac
    shift
done


# ════════════════════════════════════════════════════════════
# INPUT VALIDATION
# Check all required arguments are present and that input
# files actually exist on disk before doing any work.
# ════════════════════════════════════════════════════════════

ERRORS=0

if [[ -z "$DATASET" ]]; then
    echo "ERROR: --dataset is required"; ERRORS=$((ERRORS+1))
fi
if [[ -z "$ENA_REPORT" ]]; then
    echo "ERROR: --ena-report is required"; ERRORS=$((ERRORS+1))
elif [[ ! -f "$ENA_REPORT" ]]; then
    echo "ERROR: ENA report not found: $ENA_REPORT"; ERRORS=$((ERRORS+1))
fi
if [[ -z "$ENA_SCRIPT" ]]; then
    echo "ERROR: --ena-script is required"; ERRORS=$((ERRORS+1))
elif [[ ! -f "$ENA_SCRIPT" ]]; then
    echo "ERROR: ENA script not found: $ENA_SCRIPT"; ERRORS=$((ERRORS+1))
fi
if [[ -z "$OUTDIR" ]]; then
    echo "ERROR: --outdir is required"; ERRORS=$((ERRORS+1))
fi

if [[ $ERRORS -gt 0 ]]; then
    echo "Exiting due to $ERRORS error(s)."
    exit 1
fi


# ════════════════════════════════════════════════════════════
# OUTPUT PATH SETUP
# Derive all output file paths from the dataset name and
# the directory containing the ENA script. Everything is
# co-located with the input script for easy reference.
#
# Filenames include a "_not" suffix when --exclude is set, so
# inclusion-filtered and exclusion-filtered outputs for the same
# --sample-type never collide (e.g. "..._MCI.tsv" for "keep MCI
# only" vs "..._not_MCI.tsv" for "exclude MCI").
# ════════════════════════════════════════════════════════════

SCRIPT_DIR=$(dirname "$ENA_SCRIPT")

if [[ "$EXCLUDE" -eq 1 ]]; then
    SUFFIX="not_${SAMPLE_TYPE}"
else
    SUFFIX="${SAMPLE_TYPE}"
fi

# Filtered ENA metadata report
FILTERED_REPORT="$SCRIPT_DIR/${DATASET}_ena_report_${SUFFIX}.tsv"

# List of run accessions (ERR/SRR IDs) to keep
ACCESSIONS_LIST="$SCRIPT_DIR/${DATASET}_accessions_${SUFFIX}.txt"

# Filtered download script (wget lines for matching samples only)
FILTERED_SCRIPT="$SCRIPT_DIR/${DATASET}_ena_download_${SUFFIX}.sh"

# Log file written alongside the downloaded files
VERIFY_LOG="$OUTDIR/${DATASET}_02_verify_${SUFFIX}.log"

mkdir -p "$OUTDIR"

# Print run header to both stdout and the log file
echo "========================================"   | tee "$VERIFY_LOG"
echo "Filter and verify: $DATASET"               | tee -a "$VERIFY_LOG"
echo "Started: $(date)"                          | tee -a "$VERIFY_LOG"
echo "SLURM job ID: $SLURM_JOB_ID"              | tee -a "$VERIFY_LOG"
echo ""                                           | tee -a "$VERIFY_LOG"
echo "ENA report:    $ENA_REPORT"                | tee -a "$VERIFY_LOG"
echo "ENA script:    $ENA_SCRIPT"                | tee -a "$VERIFY_LOG"
echo "Output dir:    $OUTDIR"                    | tee -a "$VERIFY_LOG"
echo "Filter column: $FILTER_COL"                | tee -a "$VERIFY_LOG"
echo "Sample type:   $SAMPLE_TYPE"               | tee -a "$VERIFY_LOG"
echo "Mode:          $([[ $EXCLUDE -eq 1 ]] && echo 'EXCLUDE (keep rows NOT matching)' || echo 'INCLUDE (keep rows matching)')" | tee -a "$VERIFY_LOG"
echo "========================================"   | tee -a "$VERIFY_LOG"
echo ""                                           | tee -a "$VERIFY_LOG"


# ════════════════════════════════════════════════════════════
# STEP 1 — FILTER ENA REPORT TO SAMPLE TYPE
#
# The ENA metadata TSV contains one row per run (FASTQ file).
# Different datasets use different columns to encode sample type,
# host species, or diagnosis group, e.g.:
#   - Cirstea 2022: "Experiment Alias" (e.g. "fecal_sample_PD_001")
#   - Clasen 2024:  "Scientific Name"  (e.g. "fecal metagenome")
#   - host species: a "host" or "scientific name" column (e.g. "Homo sapiens")
#   - diagnosis:    a "diagnosis" or similar column (e.g. "MCI")
#
# We find the requested column by matching its header
# case-insensitively, then keep (INCLUDE mode) or remove
# (EXCLUDE mode) rows where that column contains the
# --sample-type keyword (also case-insensitive, substring match).
# ════════════════════════════════════════════════════════════

if [[ "$EXCLUDE" -eq 1 ]]; then
    echo "── STEP 1: Filter ENA report, EXCLUDING rows matching '$SAMPLE_TYPE' in column '$FILTER_COL' ──" | tee -a "$VERIFY_LOG"
else
    echo "── STEP 1: Filter ENA report, KEEPING rows matching '$SAMPLE_TYPE' in column '$FILTER_COL' ──" | tee -a "$VERIFY_LOG"
fi

awk -F'\t' -v stype="$SAMPLE_TYPE" -v fcol="$FILTER_COL" -v exclude="$EXCLUDE" '
NR==1 {
    # Scan the header row for the requested filter column (case-insensitive exact match).
    # Storing the column index in `col` so data rows can use it.
    for(i=1; i<=NF; i++) {
        if(tolower($i) == tolower(fcol)) col=i
    }
    if (!col) {
        print "WARNING: Could not find column \"" fcol "\" — printing all rows"
    }
    print $0   # always print the header row
    next
}
{
    # If the column was not found, keep all rows (fail-open, matches prior behavior).
    if (!col) { print $0; next }

    is_match = tolower($col) ~ tolower(stype)

    # INCLUDE mode (exclude=0): keep rows that MATCH.
    # EXCLUDE mode (exclude=1): keep rows that DO NOT match.
    if ((exclude == 0 && is_match) || (exclude == 1 && !is_match)) {
        print $0
    }
}
' "$ENA_REPORT" > "$FILTERED_REPORT"

# Count samples before and after filtering (subtract 1 row for the header)
TOTAL_ROWS=$(wc -l < "$ENA_REPORT")
FILTERED_ROWS=$(wc -l < "$FILTERED_REPORT")
TOTAL_SAMPLES=$(( TOTAL_ROWS - 1 ))
FILTERED_SAMPLES=$(( FILTERED_ROWS - 1 ))
EXCLUDED_SAMPLES=$(( TOTAL_SAMPLES - FILTERED_SAMPLES ))

echo "Total samples in ENA report:         $TOTAL_SAMPLES"    | tee -a "$VERIFY_LOG"
echo "Samples kept after filtering:        $FILTERED_SAMPLES" | tee -a "$VERIFY_LOG"
echo "Samples excluded:                    $EXCLUDED_SAMPLES" | tee -a "$VERIFY_LOG"
echo "Filtered report saved to:            $FILTERED_REPORT"  | tee -a "$VERIFY_LOG"
echo ""                                                         | tee -a "$VERIFY_LOG"

# Hard stop if the filter matched nothing (INCLUDE mode) or removed everything
# (EXCLUDE mode) — likely a wrong column name or keyword either way.
if [[ $FILTERED_SAMPLES -eq 0 ]]; then
    echo "ERROR: No samples remained after filtering on '$SAMPLE_TYPE' in column '$FILTER_COL' (mode: $([[ $EXCLUDE -eq 1 ]] && echo exclude || echo include))." | tee -a "$VERIFY_LOG"
    echo ""                                                                    | tee -a "$VERIFY_LOG"
    echo "Troubleshooting — column names in your ENA report:"                 | tee -a "$VERIFY_LOG"
    # Print the header row with each column on its own line for easy reading
    head -1 "$ENA_REPORT" | tr '\t' '\n' | nl                                | tee -a "$VERIFY_LOG"
    echo ""                                                                    | tee -a "$VERIFY_LOG"
    echo "First 10 values found in column '$FILTER_COL':"                     | tee -a "$VERIFY_LOG"
    # Show sample values from the requested column to help diagnose wrong keywords
    awk -F'\t' -v fcol="$FILTER_COL" '
    NR==1 { for(i=1; i<=NF; i++) if(tolower($i)==tolower(fcol)) col=i; next }
    col { print $col }
    ' "$ENA_REPORT" | head -10                                                | tee -a "$VERIFY_LOG"
    exit 1
fi


# ════════════════════════════════════════════════════════════
# STEP 2 — EXTRACT RUN ACCESSIONS FROM FILTERED REPORT
#
# Each row in the ENA report corresponds to one sequencing run.
# We extract the run accession (ERR/SRR ID) column and save
# it as a plain text list — one accession per line.
# This list is used in Steps 3 and 4.
# ════════════════════════════════════════════════════════════

echo "── STEP 2: Extract run accessions ──" | tee -a "$VERIFY_LOG"

awk -F'\t' '
NR==1 {
    # Find the run accession column; its header may be "run accession" or "run_accession"
    for(i=1; i<=NF; i++) {
        if(tolower($i) ~ /run[_ ]accession/) col=i
    }
    if (!col) {
        print "ERROR: Could not find run accession column" > "/dev/stderr"
        exit 1
    }
    next   # skip the header row in the output
}
col { print $col }   # print one accession per line for data rows
' "$FILTERED_REPORT" > "$ACCESSIONS_LIST"

N_ACCESSIONS=$(wc -l < "$ACCESSIONS_LIST")
echo "Run accessions extracted:  $N_ACCESSIONS"               | tee -a "$VERIFY_LOG"
echo "Accession list saved to:   $ACCESSIONS_LIST"            | tee -a "$VERIFY_LOG"
echo ""                                                         | tee -a "$VERIFY_LOG"

# Show first 5 accessions as a quick sanity check
echo "First 5 accessions:"                                     | tee -a "$VERIFY_LOG"
head -5 "$ACCESSIONS_LIST"                                     | tee -a "$VERIFY_LOG"
echo ""                                                         | tee -a "$VERIFY_LOG"


# ════════════════════════════════════════════════════════════
# STEP 3 — FILTER ENA DOWNLOAD SCRIPT TO MATCHING ACCESSIONS
#
# The original ENA download script contains one wget line per
# FASTQ file (two lines per sample for paired-end data).
# We strip out all wget lines whose accession is NOT in our
# filtered list, producing a smaller download script that
# only fetches the samples we want to keep (per Step 1's mode).
#
# The shebang line (#!/bin/bash) is preserved if present.
# ════════════════════════════════════════════════════════════

echo "── STEP 3: Filter ENA download script ──" | tee -a "$VERIFY_LOG"

# Count how many wget lines are in the original (unfiltered) script
ORIGINAL_WGET=$(grep -c "^wget" "$ENA_SCRIPT" || true)

# Build the filtered script:
#   Line 1: preserve shebang (#!/bin/bash) if the original has one
#   Remaining lines: only wget lines whose URL contains an accession from our list
{
    head -1 "$ENA_SCRIPT" | grep "^#!" || true
    grep "^wget" "$ENA_SCRIPT" | grep -F -f "$ACCESSIONS_LIST"
} > "$FILTERED_SCRIPT"

chmod +x "$FILTERED_SCRIPT"

FILTERED_WGET=$(grep -c "^wget" "$FILTERED_SCRIPT" || true)
EXCLUDED_WGET=$(( ORIGINAL_WGET - FILTERED_WGET ))

echo "wget lines in original script:       $ORIGINAL_WGET"    | tee -a "$VERIFY_LOG"
echo "wget lines after filtering:          $FILTERED_WGET"    | tee -a "$VERIFY_LOG"
echo "wget lines excluded:                 $EXCLUDED_WGET"    | tee -a "$VERIFY_LOG"
echo "Filtered script saved to:            $FILTERED_SCRIPT"  | tee -a "$VERIFY_LOG"
echo ""                                                         | tee -a "$VERIFY_LOG"

# Sanity check: for paired-end data, expect 2 wget lines per accession.
# If fewer wget lines than accessions, some may have been missed.
if [[ $FILTERED_WGET -lt $N_ACCESSIONS ]]; then
    echo "WARNING: Fewer wget lines ($FILTERED_WGET) than accessions ($N_ACCESSIONS)." | tee -a "$VERIFY_LOG"
    echo "Some accessions may be missing from the download script."                    | tee -a "$VERIFY_LOG"
fi


# ════════════════════════════════════════════════════════════
# STEP 4 — VERIFY WHAT HAS ACTUALLY BEEN DOWNLOADED
#
# Compare the list of accessions we expect on disk (from Step 2)
# against what is actually there (by scanning filenames).
# Reports: missing files, unexpected extras, incomplete pairs,
# and suspiciously small files that may indicate failed downloads.
#
# NOTE: in EXCLUDE mode, "extra" files found on disk that are NOT
# in the expected (kept) list are exactly the samples you want to
# remove — see 03_delete_excluded_samples.sh for a script that
# deletes them directly using this same accession-list mechanism.
# ════════════════════════════════════════════════════════════

echo "── STEP 4: Verify downloaded files ──" | tee -a "$VERIFY_LOG"

# Use temp files for set comparisons between expected and found accessions
EXPECTED_FILE=$(mktemp)
FOUND_FILE=$(mktemp)

cp "$ACCESSIONS_LIST" "$EXPECTED_FILE"

# Extract accession IDs from FASTQ filenames already on disk.
# Strips the paired-end suffix and extension to get the bare accession:
#   ERR1234567_1.fastq.gz  →  ERR1234567
#   ERR1234567_2.fastq.gz  →  ERR1234567
#   ERR1234567.fastq.gz    →  ERR1234567
# sort -u ensures each accession appears only once even if both R1 and R2 are present.
ls "$OUTDIR"/*.fastq.gz 2>/dev/null \
    | xargs -I{} basename {} \
    | sed 's/_[12]\.fastq\.gz$//' \
    | sed 's/\.fastq\.gz$//' \
    | sort -u > "$FOUND_FILE"

N_FOUND=$(wc -l < "$FOUND_FILE")
echo "Accessions expected on disk:  $N_ACCESSIONS"            | tee -a "$VERIFY_LOG"
echo "Accessions found on disk:     $N_FOUND"                 | tee -a "$VERIFY_LOG"
echo ""                                                         | tee -a "$VERIFY_LOG"

# ── Missing accessions ───────────────────────────────────────
# comm -23: lines in expected (sorted) that are NOT in found (sorted)
MISSING=$(comm -23 <(sort "$EXPECTED_FILE") <(sort "$FOUND_FILE"))
N_MISSING=$(echo "$MISSING" | grep -c "." || true)

if [[ $N_MISSING -eq 0 ]]; then
    echo "✓ No missing accessions — all expected files are present." | tee -a "$VERIFY_LOG"
else
    echo "✗ MISSING ($N_MISSING accessions not yet downloaded):"     | tee -a "$VERIFY_LOG"
    echo "$MISSING"                                                    | tee -a "$VERIFY_LOG"
fi
echo ""                                                                | tee -a "$VERIFY_LOG"

# ── Extra (unexpected) files ─────────────────────────────────
# comm -13: lines in found that are NOT in expected
# These are samples on disk that don't belong to the kept set —
# e.g. non-gut, non-human, or MCI samples downloaded before this
# filter was applied. See 03_delete_excluded_samples.sh to remove them.
EXTRA=$(comm -13 <(sort "$EXPECTED_FILE") <(sort "$FOUND_FILE"))
N_EXTRA=$(echo "$EXTRA" | grep -c "." || true)

if [[ $N_EXTRA -eq 0 ]]; then
    echo "✓ No unexpected files — nothing extra on disk."        | tee -a "$VERIFY_LOG"
else
    echo "! EXTRA files on disk not in expected list ($N_EXTRA):" | tee -a "$VERIFY_LOG"
    echo "$EXTRA"                                                   | tee -a "$VERIFY_LOG"
    echo "These are samples outside the kept set for this filter." | tee -a "$VERIFY_LOG"
    echo "Use 03_delete_excluded_samples.sh to remove them and recover storage." | tee -a "$VERIFY_LOG"
fi
echo ""                                                             | tee -a "$VERIFY_LOG"

# ── Paired-end completeness check ───────────────────────────
# For every R1 file on disk, verify that the matching R2 exists.
# A missing R2 indicates a partially failed download — the sample
# cannot be processed by fastp or DADA2 without both reads.
echo "── Paired-end completeness check ──"                    | tee -a "$VERIFY_LOG"
PAIRED_INCOMPLETE=0
for F1 in "$OUTDIR"/*_1.fastq.gz; do
    [[ -f "$F1" ]] || continue
    BASENAME="${F1%_1.fastq.gz}"
    F2="${BASENAME}_2.fastq.gz"
    if [[ ! -f "$F2" ]]; then
        echo "  [WARNING] Missing R2 for: $(basename $F1)"    | tee -a "$VERIFY_LOG"
        PAIRED_INCOMPLETE=$((PAIRED_INCOMPLETE + 1))
    fi
done

if [[ $PAIRED_INCOMPLETE -eq 0 ]]; then
    echo "✓ All paired-end files have matching R1 and R2."    | tee -a "$VERIFY_LOG"
else
    echo "✗ $PAIRED_INCOMPLETE sample(s) missing R2 file."   | tee -a "$VERIFY_LOG"
fi
echo ""                                                         | tee -a "$VERIFY_LOG"

# ── File size sanity check ───────────────────────────────────
# FASTQ files from a normal sequencing run are usually hundreds
# of MB. A file under 10 MB almost certainly represents an
# incomplete or corrupted download that will fail QC.
echo "── File size check (flagging files < 10 MB) ──"        | tee -a "$VERIFY_LOG"
N_SMALL=0
for F in "$OUTDIR"/*.fastq.gz; do
    [[ -f "$F" ]] || continue
    SIZE=$(stat -c%s "$F" 2>/dev/null || echo 0)
    if [[ $SIZE -lt 10485760 ]]; then   # 10 MB = 10 * 1024 * 1024 bytes
        echo "  [WARNING] Suspiciously small: $(basename $F) ($(du -sh $F | cut -f1))" | tee -a "$VERIFY_LOG"
        N_SMALL=$((N_SMALL + 1))
    fi
done
if [[ $N_SMALL -eq 0 ]]; then
    echo "✓ All files are above 10 MB — no suspiciously small files." | tee -a "$VERIFY_LOG"
fi
echo ""                                                                 | tee -a "$VERIFY_LOG"

# ── Disk usage summary ───────────────────────────────────────
# Always print storage used and remaining — critical for knowing
# whether you can proceed to the next download wave.
echo "── Storage ──"                                          | tee -a "$VERIFY_LOG"
echo "Total size of $DATASET downloads:"                      | tee -a "$VERIFY_LOG"
du -sh "$OUTDIR"                                               | tee -a "$VERIFY_LOG"
echo ""                                                         | tee -a "$VERIFY_LOG"
echo "Project storage remaining:"                             | tee -a "$VERIFY_LOG"
df -h /rds/projects/e/elhamsak-pd-thesis/                     | tee -a "$VERIFY_LOG"

# Clean up temporary files used for set comparison
rm -f "$EXPECTED_FILE" "$FOUND_FILE"


# ════════════════════════════════════════════════════════════
# STEP 5 — FINAL SUMMARY AND ACTION INSTRUCTIONS
#
# Print a compact summary table and, if files are missing,
# give the exact command needed to download them.
# ════════════════════════════════════════════════════════════

echo ""                                                         | tee -a "$VERIFY_LOG"
echo "========================================"                 | tee -a "$VERIFY_LOG"
echo "SUMMARY: $DATASET"                                        | tee -a "$VERIFY_LOG"
echo "========================================"                 | tee -a "$VERIFY_LOG"
echo "Expected accessions (kept set): $N_ACCESSIONS"           | tee -a "$VERIFY_LOG"
echo "Found on disk:             $N_FOUND"                     | tee -a "$VERIFY_LOG"
echo "Missing:                   $N_MISSING"                   | tee -a "$VERIFY_LOG"
echo "Extra (excluded) on disk:  $N_EXTRA"                     | tee -a "$VERIFY_LOG"
echo "Incomplete pairs:          $PAIRED_INCOMPLETE"           | tee -a "$VERIFY_LOG"
echo "Small files flagged:       $N_SMALL"                     | tee -a "$VERIFY_LOG"
echo ""                                                         | tee -a "$VERIFY_LOG"

if [[ $N_MISSING -gt 0 ]]; then
    # Some files are not yet on disk — print the exact sbatch command to fetch them.
    # The filtered script already contains only the relevant accessions,
    # and wget -nc will skip any files that already exist, so it is safe to rerun.
    echo "ACTION REQUIRED: $N_MISSING files not yet downloaded." | tee -a "$VERIFY_LOG"
    echo "Use the filtered script to download only the missing accessions:" | tee -a "$VERIFY_LOG"
    echo ""                                                        | tee -a "$VERIFY_LOG"
    echo "  sbatch 03_download.sh \\"                             | tee -a "$VERIFY_LOG"
    echo "      --script $FILTERED_SCRIPT \\"                     | tee -a "$VERIFY_LOG"
    echo "      --outdir $OUTDIR"                                  | tee -a "$VERIFY_LOG"
    echo ""                                                        | tee -a "$VERIFY_LOG"
    echo "The -nc flag in the download script skips already-downloaded files." | tee -a "$VERIFY_LOG"
elif [[ $N_EXTRA -gt 0 ]]; then
    # All expected files are present, but there are extra files on disk
    # that fall outside the kept set (e.g. non-gut, non-human, or MCI samples
    # downloaded before this filter was applied).
    echo "ACTION OPTIONAL: $N_EXTRA excluded-set files on disk."  | tee -a "$VERIFY_LOG"
    echo "Run 03_delete_excluded_samples.sh with this same accession" | tee -a "$VERIFY_LOG"
    echo "list to remove them and recover storage."                | tee -a "$VERIFY_LOG"
else
    # All expected files present, nothing unexpected — ready for QC.
    echo "✓ ALL COMPLETE: All expected files are present,"        | tee -a "$VERIFY_LOG"
    echo "  no missing files, no unexpected extras."             | tee -a "$VERIFY_LOG"
    echo "  Ready to proceed to QC (FastQC + fastp)."           | tee -a "$VERIFY_LOG"
fi

echo ""                                                         | tee -a "$VERIFY_LOG"
echo "Finished: $(date)"                                       | tee -a "$VERIFY_LOG"
echo "Full log: $VERIFY_LOG"                                   | tee -a "$VERIFY_LOG"

exit 0