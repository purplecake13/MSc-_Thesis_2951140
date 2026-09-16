#!/usr/bin/env python3
"""
qc_01_antibiotic_filter.py

Applies the antibiotic-free-period exclusion rule agreed with Dr Sakhaee:
    - Prefer >=3 months (90 days) antibiotic-free period.
    - If the 3-month threshold would exclude >20% of samples, relax to
      >=1 month (30 days) instead.
    - Samples still failing the chosen threshold -- including fully
      undocumented/unknown antibiotic status -- are excluded. Nothing
      more lenient than 1 month is acceptable.

BEFORE RUNNING: confirm ANTIBIOTIC_COLUMN matches your actual
master_metadata.csv column. Check with:
    python -c "import pandas as pd; print(pd.read_csv('master_metadata.csv').columns.tolist())"
If your metadata stores this as a category rather than numeric days,
edit the `days` computation below accordingly.

Output:
    antibiotic_filter_summary.txt   -- drop % at both thresholds, decision, per-study breakdown
    excluded_samples_antibiotic.txt -- unique_sample_id list to exclude from ALL downstream steps
"""

import pandas as pd
import sys

METADATA_PATH = "/rds/projects/e/elhamsak-pd-thesis/cross_platform/master_metadata.csv"
OUT_SUMMARY = "/rds/projects/e/elhamsak-pd-thesis/cross_platform/antibiotic_filter_summary.txt"
OUT_EXCLUDED = "/rds/projects/e/elhamsak-pd-thesis/cross_platform/excluded_samples_antibiotic.txt"

# EDIT: column holding antibiotic-free period in weeks (numeric; NaN = unknown/undocumented)
ANTIBIOTIC_COLUMN = "antibiotic_free_period"

THRESH_3MO = 12
THRESH_1MO = 4
RELAX_TRIGGER_PCT = 20.0


def main():
    meta = pd.read_csv(METADATA_PATH)

    if "unique_sample_id" not in meta.columns:
        raise ValueError("master_metadata.csv has no unique_sample_id column.")

    if ANTIBIOTIC_COLUMN not in meta.columns:
        print(f"ERROR: '{ANTIBIOTIC_COLUMN}' not found in master_metadata.csv.")
        print(f"Available columns: {meta.columns.tolist()}")
        print("Edit ANTIBIOTIC_COLUMN at the top of this script, or if your data uses "
              "categories instead of numeric days, edit the filtering logic below.")
        sys.exit(1)

    n_total = len(meta)
    days = pd.to_numeric(meta[ANTIBIOTIC_COLUMN], errors="coerce")

    def pct_excluded(threshold_days):
        mask = days.isna() | (days < threshold_days)
        return mask, 100.0 * mask.sum() / n_total

    mask_3mo, pct_3mo = pct_excluded(THRESH_3MO)
    mask_1mo, pct_1mo = pct_excluded(THRESH_1MO)

    lines = []
    def log(msg):
        print(msg)
        lines.append(msg)

    log(f"Total samples in metadata: {n_total}")
    log(f"3-month (90d) threshold: would exclude {mask_3mo.sum()} samples ({pct_3mo:.1f}%)")
    log(f"1-month (30d) threshold: would exclude {mask_1mo.sum()} samples ({pct_1mo:.1f}%)")

    if pct_3mo <= RELAX_TRIGGER_PCT:
        chosen_threshold, chosen_mask = THRESH_3MO, mask_3mo
        log(f"\nDECISION: 3-month threshold drops {pct_3mo:.1f}% (<= {RELAX_TRIGGER_PCT}%) -- using it.")
    else:
        chosen_threshold, chosen_mask = THRESH_1MO, mask_1mo
        log(f"\nDECISION: 3-month threshold drops {pct_3mo:.1f}% (> {RELAX_TRIGGER_PCT}%) -- "
            f"relaxing to 1-month ({pct_1mo:.1f}% dropped instead).")

    excluded_ids = meta.loc[chosen_mask, "unique_sample_id"]
    log(f"\nFinal exclusion: {len(excluded_ids)} samples ({100*len(excluded_ids)/n_total:.1f}%) "
        f"at the {chosen_threshold}-day threshold.")

    if "study_id" in meta.columns:
        breakdown = meta.loc[chosen_mask, "study_id"].value_counts()
        log("\nExcluded samples by study (check no single dataset is driving this):")
        log(breakdown.to_string())

    excluded_ids.to_csv(OUT_EXCLUDED, index=False, header=False)
    with open(OUT_SUMMARY, "w") as f:
        f.write("\n".join(lines))

    log(f"\nWrote exclusion list to: {OUT_EXCLUDED}")
    log(f"Wrote summary to: {OUT_SUMMARY}")
    log("\nNEXT: apply this list via pd_07 / ad_04 scripts BEFORE cx_01's genus filtering "
        "(prevalence/abundance thresholds must be computed on the surviving sample set).")


if __name__ == "__main__":
    main()
    