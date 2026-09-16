#!/bin/bash
# dada2_retention_summary.sh
# Summarizes per-dataset DADA2 retention (% non-chimeric) into buckets,
# and applies an absolute-read-count exclusion threshold (not just %),
# consistent with the project's existing 1,000-read minimum-depth convention.
#
# Usage: /rds/projects/e/elhamsak-pd-thesis/scripts/ad_04_dada2_denoise_retention_summary.sh <dataset1> <dataset2> ...
# Requires: each dataset's stats_export/stats.tsv to already exist
#           (produced by ad_04_dada2_denoise_verify.sh)

MIN_ABS_READS=1000   # same floor already used for Yildirim2022 sample exclusion

BASE_DIR="/rds/projects/e/elhamsak-ad-thesis/ad_qiime2"

if [[ $# -eq 0 ]]; then
  echo "Usage: $0 <dataset1> <dataset2> ..."
  exit 1
fi

printf "%-15s %6s %6s %6s %6s %6s %6s | %s\n" \
  "Dataset" "n" "<50%" "50-70" "70-90" ">=90%" "<${MIN_ABS_READS}reads" "Flagged samples (below abs. floor)"
echo "--------------------------------------------------------------------------------------------------"

for DS in "$@"; do
  STATS="${BASE_DIR}/${DS}/stats_export/stats.tsv"
  if [[ ! -f "$STATS" ]]; then
    echo "$DS: stats.tsv not found at $STATS -- skipping"
    continue
  fi

  # stats.tsv columns: sample-id input filtered pct_filter denoised non-chimeric pct_nonchim
  # skip header row(s); qiime metadata files sometimes have a #q2:types second line
  awk -F'\t' -v ds="$DS" -v floor="$MIN_ABS_READS" '
    NR==1 { next }
    $1 ~ /^#/ { next }
    {
      n++
      nonchim = $6
      pct = $7
      if (pct+0 < 50) lt50++
      else if (pct+0 < 70) lt70++
      else if (pct+0 < 90) lt90++
      else ge90++

      if (nonchim+0 < floor) {
        below_floor++
        flagged = flagged sep $1 "(" nonchim "reads," pct "%)"
        sep = ", "
      }
    }
    END {
      printf "%-15s %6d %6d %6d %6d %6d %6d | %s\n", \
        ds, n, lt50+0, lt70+0, lt90+0, ge90+0, below_floor+0, flagged
    }
  ' "$STATS"
done

echo ""
echo "Notes:"
echo "  - '<70%' and similar buckets are non-chimeric retention PERCENTAGE, for quick scanning."
echo "  - The actual exclusion recommendation uses ABSOLUTE non-chimeric read count"
echo "    (< ${MIN_ABS_READS} reads), consistent with the Yildirim2022 SRR14711477 precedent --"
echo "    a sample at 60% retention with 40,000 input reads still yields a healthy ~24,000 usable"
echo "    reads and should NOT be excluded on percentage alone. Only samples in the last column"
echo "    are exclusion candidates; treat the percentage columns as a review/QC signal, not a"
echo "    cutoff by themselves."
echo "  - Samples with low % but healthy absolute counts are fine to keep, but consider noting"
echo "    them in Methods/Limitations as a per-dataset quality observation (e.g. Ling2021, Cirstea2022"
echo "    show systematically lower retention than other datasets even where no sample is excluded)."