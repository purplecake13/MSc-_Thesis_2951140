#!/usr/bin/env python3
"""
pd_06_build_filtered_counts_for_dea.py

WHY THIS SCRIPT EXISTS (SCOPED, NOT A cx_01/cx_02/cx_03 REPLACEMENT)
------------------------------------------------------------------------
cx_01_filter_normalize.py treats PD's MetaPhlAn4 input as already-relative
abundance (already_relative=True), so the "filtered_raw" PD table it
writes -- which flows through cx_02 and cx_03 into
joined_genus_metadata.tsv -- has never contained real read counts.
ANCOM-BC2 needs genuine inter-sample depth variation to estimate its
per-sample bias-correction term; relative abundance destroys that
signal (confirmed: joined_genus_metadata.tsv's shotgun rows sum to
~90-100 regardless of sample, vs AD's 16S rows which vary 2,000-69,000
matching real depth).

Editing cx_01 directly and rerunning the full cx_01->cx_02->cx_03 chain
would also silently change inputs to diversity (cx_09) and ConQuR
(already run, feeds the RF/SHAP feature matrix) -- both already
completed and verified this session. Given the timeline, this script
instead builds a SEPARATE, scoped, count-based PD table -- same
filtering thresholds as cx_01 for methodological consistency, same
join logic as cx_03 -- that feeds ONLY cx_08's PD_vs_HC block. Nothing
upstream or downstream of cx_08 is touched.

STEPS
------------------------------------------------------------------------
1. Parse per-sample MetaPhlAn4 profiles (estimated_number_of_reads_from
   _the_clade column -- see rebuild_pd_genus_counts.py), tagging each
   sample with study_id from its input directory name.
2. Filter genera: keep only genus present in >=10% of samples AND with
   mean relative abundance >=0.01% -- SAME thresholds as cx_01's
   PREVALENCE_THRESHOLD / MEAN_ABUNDANCE_PCT_THRESHOLD, computed on
   relative abundance derived from these real counts (not the counts
   themselves) so the filter criterion is scale-independent -- but no
   CLR transform, since ANCOM-BC2 wants counts, not CLR values.
3. Join to master_metadata.csv: primary join on study_id + sample_id,
   fallback join on study_id + run_accession -- mirrors cx_03's join
   logic exactly (including the Boktor study_id naming fix), scoped to
   PD samples.

Output:
    /rds/projects/e/elhamsak-pd-thesis/cross_platform/pd_genus_counts_filtered_joined.tsv
    -- columns: unique_sample_id, study_id, disease_group, g__<genus>...
    -- ready to feed directly into cx_08's PD_vs_HC block.

Usage:
    python pd_06_build_filtered_counts_for_dea.py
"""

import os
import glob
import sys
import pandas as pd

# ---- EDIT THESE PATHS to match your actual per-study MetaPhlAn4 output dirs ----
# Directory basename is used as study_id -- must match master_metadata.csv's
# study_id spelling exactly (including the boktor2023_rumc/tbc convention).
PD_INPUT_DIRS = {
    "bedarf2017":        "/rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/bedarf2017",
    "boktor2023_rumc":   "/rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/boktor2023/rumc",
    "boktor2023_tbc":    "/rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/boktor2023/tbc",
    "clasen2024":        "/rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/clasen2024",
    "wallen2022":        "/rds/projects/e/elhamsak-ad-thesis/pd_metaphlan4_oct22/wallen2022",
    "mao2021":           "/rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/mao2021",
}

METADATA_PATH = "/rds/projects/e/elhamsak-pd-thesis/cross_platform/master_metadata.csv"
OUT_PATH = "/rds/projects/e/elhamsak-pd-thesis/cross_platform/pd_genus_counts_filtered_joined.tsv"
SUMMARY_PATH = "/rds/projects/e/elhamsak-pd-thesis/cross_platform/pd_genus_counts_filtered_joined_summary.txt"

PREVALENCE_THRESHOLD = 0.10          # matches cx_01
MEAN_ABUNDANCE_PCT_THRESHOLD = 0.01  # matches cx_01

log_lines = []
def log(msg):
    print(msg)
    log_lines.append(msg)


def parse_profile(filepath, study_id):
    """Same genus-level parse as rebuild_pd_genus_counts.py, but also
    returns study_id for downstream merging."""
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
    return sample_id, study_id, genus_counts


def build_count_matrix():
    records = {}
    sample_study = {}
    for study_id, d in PD_INPUT_DIRS.items():
        if not os.path.isdir(d):
            log(f"WARNING: dir not found for {study_id}, skipping: {d}")
            continue
        files = glob.glob(os.path.join(d, "**", "*_metaphlan_profile.txt"), recursive=True)
        log(f"[{study_id}] found {len(files)} profile files")
        for fp in files:
            sample_id, sid, genus_counts = parse_profile(fp, study_id)
            if not genus_counts:
                log(f"  WARNING: no genus rows parsed for {fp}")
                continue
            records[sample_id] = genus_counts
            sample_study[sample_id] = sid

    df = pd.DataFrame.from_dict(records, orient="index").fillna(0.0)
    df.index.name = "sample_id"
    study_series = pd.Series(sample_study, name="study_id")
    log(f"\nTotal samples parsed: {len(df)}, total genera (pre-filter): {df.shape[1]}")
    return df, study_series


def filter_genera(counts_df):
    """Compute per-sample relative abundance for filtering purposes only;
    return the subset of genus columns passing both thresholds. Mirrors
    cx_01's filter_genera logic exactly."""
    row_totals = counts_df.sum(axis=1)
    rel_abund = counts_df.div(row_totals, axis=0) * 100.0

    prevalence = (counts_df > 0).mean(axis=0)
    mean_abund = rel_abund.mean(axis=0)

    keep = (prevalence >= PREVALENCE_THRESHOLD) & (mean_abund >= MEAN_ABUNDANCE_PCT_THRESHOLD)
    kept_genera = keep[keep].index.tolist()
    log(f"Filtering kept {len(kept_genera)}/{counts_df.shape[1]} genera "
        f"(prevalence>={PREVALENCE_THRESHOLD*100:.0f}% AND mean_abund>="
        f"{MEAN_ABUNDANCE_PCT_THRESHOLD}% required)")
    return counts_df[kept_genera]


def join_metadata(counts_df, study_series):
    meta = pd.read_csv(METADATA_PATH)

    if "unique_sample_id" not in meta.columns:
        if {"study_id", "sample_id_local"}.issubset(meta.columns):
            meta["unique_sample_id"] = meta["study_id"] + "_" + meta["sample_id_local"].astype(str)
        else:
            raise ValueError("master_metadata.csv has neither unique_sample_id nor "
                              "study_id+sample_id_local -- cannot join.")

    counts_df = counts_df.copy()
    counts_df["study_id"] = study_series
    counts_df["unique_sample_id"] = counts_df["study_id"] + "_" + counts_df.index.astype(str)

    # ---- primary join ----
    primary = counts_df.merge(meta, on="unique_sample_id", how="inner", suffixes=("", "_meta"))
    log(f"Primary join (unique_sample_id): {len(primary)} samples matched.")

    matched_ids = set(primary["unique_sample_id"])
    unmatched = counts_df[~counts_df["unique_sample_id"].isin(matched_ids)]

    # ---- fallback join on study_id + run_accession (mirrors cx_03) ----
    fallback_result = pd.DataFrame()
    if "run_accession" in meta.columns and len(unmatched) > 0:
        meta_alt = meta[~meta["run_accession"].isna() & ~meta["run_accession"].astype(str).str.contains(";", na=False)].copy()
        meta_alt["alt_key"] = meta_alt["study_id"] + "_" + meta_alt["run_accession"].astype(str)
        meta_alt = meta_alt[~meta_alt["unique_sample_id"].isin(matched_ids)]

        fallback_result = unmatched.merge(
            meta_alt, left_on="unique_sample_id", right_on="alt_key",
            how="inner", suffixes=("", "_meta")
        )
        log(f"Fallback join (study_id + run_accession): {len(fallback_result)} additional samples matched.")

    joined = pd.concat([primary, fallback_result], ignore_index=True, sort=False)
    log(f"Total joined: {len(joined)} / {len(counts_df)} PD samples.")

    n_dropped = len(counts_df) - len(joined)
    if n_dropped > 0:
        pct = 100 * n_dropped / len(counts_df)
        log(f"Dropped {n_dropped} unmatched samples ({pct:.1f}%).")
        if pct > 10:
            log("*** WARNING: >10% unmatched -- check study_id spelling against master_metadata.csv.")

    if "disease_group" not in joined.columns:
        raise ValueError("Joined table has no disease_group column -- check master_metadata.csv.")
    na_disease = joined["disease_group"].isna().sum()
    log(f"disease_group NAs after join: {na_disease}")

    return joined


def main():
    counts_df, study_series = build_count_matrix()
    if counts_df.empty:
        log("ERROR: no samples parsed -- check PD_INPUT_DIRS paths.")
        sys.exit(1)

    filtered_df = filter_genera(counts_df)
    joined = join_metadata(filtered_df, study_series)

    genus_cols = [c for c in filtered_df.columns]
    keep_cols = ["unique_sample_id", "study_id", "disease_group"] + genus_cols
    keep_cols = [c for c in keep_cols if c in joined.columns]
    out_df = joined[keep_cols]

    out_df.to_csv(OUT_PATH, sep="\t", index=False)
    log(f"\nWrote {out_df.shape[0]} samples x {len(genus_cols)} genera to {OUT_PATH}")

    ct = out_df["disease_group"].value_counts()
    log(f"\ndisease_group counts:\n{ct.to_string()}")

    with open(SUMMARY_PATH, "w") as f:
        f.write("\n".join(log_lines))
    log(f"\nSummary written to: {SUMMARY_PATH}")


if __name__ == "__main__":
    main()
    