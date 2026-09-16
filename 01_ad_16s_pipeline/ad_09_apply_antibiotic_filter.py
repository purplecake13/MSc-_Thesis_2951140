#!/usr/bin/env python3
"""
ad_09_apply_antibiotic_filter.py

AD's combined_ad_genus_table.tsv already contains real read counts -- no
units bug there. This script only removes antibiotic-excluded samples
before cx_01 filters/CLRs the AD side, so filtering thresholds are computed
on the correct final sample set.
"""

import pandas as pd

INPUT_PATH = "/rds/projects/e/elhamsak-ad-thesis/ad_genus_level/combined_ad_genus_table.tsv"
OUTPUT_PATH = INPUT_PATH  # overwrite in place -- cx_01's AD_INPUT path is unchanged
EXCLUDED_SAMPLES_PATH = "/rds/projects/e/elhamsak-pd-thesis/cross_platform/excluded_samples_antibiotic.txt"


def main():
    df = pd.read_csv(INPUT_PATH, sep="\t")
    sample_col, study_col = df.columns[0], df.columns[1]

    with open(EXCLUDED_SAMPLES_PATH) as f:
        excluded_ids = set(line.strip() for line in f if line.strip())

    df["_unique_id"] = df[study_col] + "_" + df[sample_col].astype(str)
    n_before = len(df)
    df_filtered = df[~df["_unique_id"].isin(excluded_ids)].drop(columns=["_unique_id"])
    n_after = len(df_filtered)

    print(f"AD table: {n_before} samples before filter, {n_after} after "
          f"({n_before - n_after} excluded).")

    df_filtered.to_csv(OUTPUT_PATH, sep="\t", index=False)
    print(f"Wrote filtered table to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
    