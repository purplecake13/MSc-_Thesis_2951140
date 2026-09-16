#!/usr/bin/env python3
"""
cx_02_combine_platforms.py

Combines the per-platform, per-dataset filtered + CLR-normalised genus
tables (PD shotgun + AD 16S) into a single samples x genera matrix, using
Option B from the schedule doc: keep the UNION of all genera across both
platforms, zero-filling any genus a given platform never detected, and add
a `platform` column so downstream batch correction (ConQuR/MMUPHin) and
models can use platform as an explicit covariate rather than conflating it
with disease signal.

Inputs (from cx_01_filter_normalize.py):
    /rds/projects/e/elhamsak-pd-thesis/filtered_normalized/pd_genus_filtered_clr.tsv
    /rds/projects/e/elhamsak-pd-thesis/filtered_normalized/ad_genus_filtered_clr.tsv

Both are expected to be samples (rows) x genera (columns), tab-separated,
with sample_id as the first column / index.

Output:
    /rds/projects/e/elhamsak-pd-thesis/cross_platform/combined_genus_clr_platform.tsv
        - all samples (PD + AD) as rows
        - union of all genera as columns (zero-filled where a platform
          never detected that genus)
        - a `platform` column ("shotgun" / "16S")
        - a `study_id`-equivalent `disease_arm` column ("PD" / "AD") for
          convenience, INFERRED from which input file the row came from
          -- NOT a substitute for the real study_id, which will come from
          the metadata master CSV join later.

Does NOT touch batch correction, metadata join, or ANCOM-BC2 -- this is
purely the matrix-alignment step. Batch correction is the next script.
"""

import pandas as pd
from pathlib import Path

# PD_CLR_PATH = Path("/rds/projects/e/elhamsak-pd-thesis/filtered_normalized/pd_genus_filtered_clr.tsv")
# AD_CLR_PATH = Path("/rds/projects/e/elhamsak-pd-thesis/filtered_normalized/ad_genus_filtered_clr.tsv")

# ConQuR needs pre-clr files. 
PD_CLR_PATH = Path("/rds/projects/e/elhamsak-pd-thesis/filtered_normalized/pd_genus_filtered_raw.tsv")
AD_CLR_PATH = Path("/rds/projects/e/elhamsak-pd-thesis/filtered_normalized/ad_genus_filtered_raw.tsv")

OUT_DIR = Path("/rds/projects/e/elhamsak-pd-thesis/cross_platform")
OUT_DIR.mkdir(parents=True, exist_ok=True)
# OUT_PATH = OUT_DIR / "combined_genus_clr_platform.tsv"
OUT_PATH = OUT_DIR / "combined_genus_raw_platform.tsv"
SUMMARY_PATH = OUT_DIR / "combine_platforms_summary.txt"


def load_clr_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    # Basic sanity: index should be sample IDs, columns genera, values numeric
    if df.shape[0] == 0 or df.shape[1] == 0:
        raise ValueError(f"{path} loaded as empty table -- check input file.")
        
    # Ensure study_id exists, then check that remaining taxonomic columns are numeric
    if "study_id" not in df.columns:
        raise ValueError(f"{path} missing required 'study_id' column.")
    taxa_cols = [c for c in df.columns if c != "study_id"]
    
    non_numeric = df[taxa_cols].select_dtypes(exclude="number").columns.tolist()
    if non_numeric:
        raise ValueError(
            f"{path} has non-numeric columns after index_col=0: {non_numeric}. "
            "Check that sample_id is really the first column."
        )
    return df


def main():
    pd_df = load_clr_table(PD_CLR_PATH)
    ad_df = load_clr_table(AD_CLR_PATH)

    # pd_genera = set(pd_df.columns)
    # ad_genera = set(ad_df.columns)
    # all_genera = sorted(pd_genera | ad_genera)
    
    # Extract genera columns (excluding study_id)
    pd_genera = set(pd_df.columns) - {"study_id"}
    ad_genera = set(ad_df.columns) - {"study_id"}
    
    # Define set comparisons for reporting
    shared = pd_genera & ad_genera
    pd_only = pd_genera - ad_genera
    ad_only = ad_genera - pd_genera
    
    all_genera = sorted(pd_genera | ad_genera)

    # Reindex both tables onto the full union of genera, zero-filling gaps.
    # Zero is the correct fill value here: these are CLR values, and a
    # genus a platform never detected is being treated as "not present /
    # not measurable on this platform" -- zero-filling a CLR matrix isn't
    # perfectly principled (true CLR zero != true abundance zero) but is
    # the standard pragmatic choice at this stage; ConQuR/MMUPHin will see
    # platform as a covariate and can absorb systematic platform-genus
    # detectability differences.
    
    # pd_aligned = pd_df.reindex(columns=all_genera, fill_value=0.0)
    # ad_aligned = ad_df.reindex(columns=all_genera, fill_value=0.0)
    pd_aligned = pd_df[list(pd_genera)].reindex(columns=all_genera, fill_value=0.0)
    ad_aligned = ad_df[list(ad_genera)].reindex(columns=all_genera, fill_value=0.0)
    
    # pd_aligned = pd_aligned.copy()
    #pd_aligned["platform"] = "shotgun"
    #pd_aligned["disease_arm"] = "PD"

    #ad_aligned = ad_aligned.copy()
    #ad_aligned["platform"] = "16S"
    #ad_aligned["disease_arm"] = "AD"
    
    pd_aligned["study_id"] = pd_df["study_id"].values
    pd_aligned["platform"] = "shotgun"
    pd_aligned["disease_arm"] = "PD"

    ad_aligned["study_id"] = ad_df["study_id"].values
    ad_aligned["platform"] = "16S"
    ad_aligned["disease_arm"] = "AD"

    combined = pd.concat([pd_aligned, ad_aligned], axis=0)
    combined.index.name = "sample_id"
    
    # combined = combined.rename_axis("sample_id_local")
    # combined["unique_sample_id"] = combined["study_id"].astype(str) + "_" + combined.index.astype(str)
    combined["unique_sample_id"] = combined.index.astype(str)

    # Check for sample_id collisions across platforms (should not happen,
    # but worth catching explicitly rather than silently overwriting rows)
    dup_ids = combined.index[combined.index.duplicated()].unique().tolist()

    combined.to_csv(OUT_PATH, sep="\t")

    summary_lines = [
        "=== Cross-platform combine summary ===",
        f"PD samples: {pd_df.shape[0]}, PD genera pre-union: {pd_df.shape[1]}",
        f"AD samples: {ad_df.shape[0]}, AD genera pre-union: {ad_df.shape[1]}",
        f"Union genera: {len(all_genera)}",
        f"  shared (both platforms): {len(shared)}",
        f"  PD-only (zero-filled for AD rows): {len(pd_only)}",
        f"  AD-only (zero-filled for PD rows): {len(ad_only)}",
        f"Combined table shape: {combined.shape[0]} samples x "
        f"{combined.shape[1] - 3} genera (+platform +disease_arm columns)",
        f"Output: {OUT_PATH}",
    ]
    if dup_ids:
        summary_lines.append(
            f"*** WARNING: {len(dup_ids)} sample_id collisions across "
            f"platforms -- inspect before trusting combined table: {dup_ids[:10]}"
            + (" ..." if len(dup_ids) > 10 else "")
        )
    else:
        summary_lines.append("No sample_id collisions between platforms -- OK.")

    summary_text = "\n".join(summary_lines)
    print(summary_text)
    SUMMARY_PATH.write_text(summary_text + "\n")


if __name__ == "__main__":
    main()
    