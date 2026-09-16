#!/usr/bin/env python3
"""
pd_07_build_combined_pd_table_final.py

Rebuilds combined_pd_genus_table.tsv (cx_01's PD_INPUT) using REAL per-sample
read counts (not relative abundance) with both:
  1. Antibiotic-free-period exclusion (qc_01)
  2. Oral contamination / saliva misclassification exclusions (clasen2024 7 samples)
applied BEFORE the table is built.

Output format matches cx_01's read_table() expectation exactly:
    sample_id \t study_id \t genus1 \t genus2 \t ...
"""

import os
import glob
import pandas as pd

# EDIT to match your actual directories
PD_INPUT_DIRS = {
    "bedarf2017":        "/rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/bedarf2017",
    "boktor2023_rumc":   "/rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/boktor2023/rumc",
    "boktor2023_tbc":    "/rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/boktor2023/tbc",
    "clasen2024":        "/rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/clasen2024",
    "wallen2022":        "/rds/projects/e/elhamsak-ad-thesis/pd_metaphlan4_oct22/wallen2022",
    "mao2021":           "/rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/mao2021",
}

EXCLUDED_SAMPLES_PATH = "/rds/projects/e/elhamsak-pd-thesis/cross_platform/excluded_samples_antibiotic.txt"
OUTPUT_PATH = "/rds/projects/e/elhamsak-pd-thesis/pd_genus_level/combined_pd_genus_table.tsv"

# Flagged oral-contamination / tissue-misclassification samples in clasen2024
ORAL_CONTAMINATED_SAMPLES = {
    "ERR15003744",
    "ERR15003808",
    "ERR15003814",
    "ERR15003915",
    "ERR15003916",
    "ERR15003922",
    "ERR15003954",
    "ERR15003857"
}


def parse_profile(filepath):
    genus_counts = {}
    sample_id = os.path.basename(filepath).replace("_metaphlan_profile.txt", "")
    with open(filepath) as f:
        header_seen = False
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#SampleID") or line.startswith("clade_name"):
                header_seen = True
                continue
            if line.startswith("#") or not header_seen:
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            clade_name = fields[0]
            try:
                read_count = float(fields[4])
            except (ValueError, IndexError):
                continue
            last_level = clade_name.split("|")[-1]
            if not last_level.startswith("g__") or len(last_level) <= 3:
                continue
            genus_counts[last_level] = genus_counts.get(last_level, 0.0) + read_count
    return sample_id, genus_counts


def main():
    excluded_ids = set()
    if os.path.exists(EXCLUDED_SAMPLES_PATH):
        with open(EXCLUDED_SAMPLES_PATH) as f:
            excluded_ids = set(line.strip() for line in f if line.strip())
        print(f"Loaded {len(excluded_ids)} excluded sample IDs from antibiotic filter.")
    else:
        print(f"WARNING: {EXCLUDED_SAMPLES_PATH} not found -- proceeding with NO "
              f"antibiotic exclusion applied.")

    records, sample_study = {}, {}
    n_excluded_abx = 0
    n_excluded_oral = 0

    for study_id, d in PD_INPUT_DIRS.items():
        if not os.path.isdir(d):
            print(f"WARNING: dir not found for {study_id}, skipping: {d}")
            continue
        files = glob.glob(os.path.join(d, "**", "*_metaphlan_profile.txt"), recursive=True)
        print(f"[{study_id}] found {len(files)} profile files")
        for fp in files:
            sample_id, genus_counts = parse_profile(fp)
            unique_id = f"{study_id}_{sample_id}"

            # Check antibiotic exclusions (both formatted and bare ID)
            if unique_id in excluded_ids or sample_id in excluded_ids:
                n_excluded_abx += 1
                continue

            # Check oral contamination exclusions
            if sample_id in ORAL_CONTAMINATED_SAMPLES:
                n_excluded_oral += 1
                print(f"  [EXCLUDED ORAL] {sample_id} ({study_id})")
                continue

            if not genus_counts:
                print(f"  WARNING: no genus rows parsed for {fp}")
                continue

            records[unique_id] = genus_counts
            sample_study[unique_id] = study_id

    print(f"\nExcluded {n_excluded_abx} samples via antibiotic filter.")
    print(f"Excluded {n_excluded_oral} samples via oral contamination filter.")
    print(f"Retained {len(records)} samples for the combined PD table.")

    df = pd.DataFrame.from_dict(records, orient="index").fillna(0.0)
    df.insert(0, "study_id", pd.Series(sample_study))
    df.index.name = "sample_id"
    
    assert df["study_id"].isna().sum() == 0, (
        f"NA study_id for {df['study_id'].isna().sum()} rows after build -- "
        "check that records/sample_study dicts are keyed identically."
    )

    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    df.to_csv(OUTPUT_PATH, sep="\t")
    print(f"\nWrote {df.shape[0]} samples x {df.shape[1]-1} genera to {OUTPUT_PATH}")
    print("\nNEXT: edit cx_01 (already_relative=True -> False for PD), then rerun "
          "cx_01 -> cx_02 -> cx_03.")


if __name__ == "__main__":
    main()
    