#!/usr/bin/env python3
"""
pd_05_merge_genus_level.py

Merges all 6 PD MetaPhlAn4 (Oct22) profiles into a single samples x genera
feature table at genus-level resolution, per the locked cross-platform
genus-collapse decision.

Usage:
    python pd_05_merge_genus_level.py

Reads the same dataset_key -> OUT_DIR / SAMPLE_LIST mapping used by
pd_02_verify_metaphlan4_oct22.sh, so this should only be run after that
verify script has returned CLEAN for a given dataset.

Genus-level row definition (per MetaPhlAn4 clade_name format):
    A row is genus-level if its pipe-delimited clade_name ends in "g__X"
    with no deeper "|s__..." or "|t__..." segment after it. This correctly
    includes "*_unclassified" genus rows and excludes species/SGB rows.

Respects sample_exclusions_log.tsv: any sample_id listed there is skipped
entirely, regardless of dataset, so exclusions logged during MetaPhlAn4/
HUMAnN3 QC automatically propagate into the merge without needing a
separate manual step.

Outputs (written to OUT_ROOT, created if missing):
    pd_genus_level/combined_pd_genus_table.tsv   samples x genera, filled 0
    pd_genus_level/per_dataset/<dataset>_genus_table.tsv   per-dataset tables
    pd_genus_level/sample_study_map.tsv          sample_id -> study_id
    pd_genus_level/merge_summary.txt             run summary / sanity stats

Does NOT filter rare taxa or apply CLR -- that's a separate downstream step
(rare-taxon filtering + normalisation), kept independent so this table can
be re-filtered without re-parsing all profiles.
"""

import os
import re
import sys
from collections import defaultdict

PD_BASE = "/rds/projects/e/elhamsak-pd-thesis"
AD_BASE = "/rds/projects/e/elhamsak-ad-thesis"
OUT_ROOT = f"{PD_BASE}/pd_genus_level"
EXCLUSIONS_LOG = f"{PD_BASE}/sample_exclusions_log.tsv"

# Same mapping as pd_02_verify_metaphlan4_oct22.sh -- keep in sync if that
# script's paths ever change.
DATASETS = {
    "bedarf2017": {
        "out_dir": f"{PD_BASE}/pd_metaphlan4_oct22/bedarf2017",
        "sample_list": f"{PD_BASE}/pd_raw_data/bedarf2017/bedarf2017_merged_samples.txt",
    },
    "boktor_rumc": {
        "out_dir": f"{PD_BASE}/pd_metaphlan4_oct22/boktor2023/rumc",
        "sample_list": f"{PD_BASE}/pd_raw_data/boktor2023/rumc/boktor2023rumc_samples_filtered.txt",
    },
    "boktor_tbc": {
        "out_dir": f"{PD_BASE}/pd_metaphlan4_oct22/boktor2023/tbc",
        "sample_list": f"{PD_BASE}/pd_raw_data/boktor2023/tbc/boktor2023tbc_samples.txt",
    },
    "clasen2024": {
        "out_dir": f"{PD_BASE}/pd_metaphlan4_oct22/clasen2024",
        "sample_list": f"{PD_BASE}/pd_raw_data/clasen2024/clasen2024_samples.txt",
    },
    "wallen2022": {
        "out_dir": f"{AD_BASE}/pd_metaphlan4_oct22/wallen2022",
        "sample_list": f"{AD_BASE}/pd_raw_data/wallen2022/wallen2022_samples.txt",
    },
    "mao2021": {
        "out_dir": f"{PD_BASE}/pd_metaphlan4_oct22/mao2021",
        "sample_list": f"{PD_BASE}/pd_raw_data/mao2021/mao2021_samples.txt",
    },
}

# Genus-level clade: ends in g__<name> with nothing deeper (no |s__ or |t__ after)
GENUS_ROW_RE = re.compile(r"g__[^|]+$")


def load_exclusions(path):
    excluded = set()
    if not os.path.isfile(path):
        print(f"NOTE: exclusions log not found at {path} -- proceeding with no exclusions applied.")
        return excluded
    with open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Expect sample_id as first tab-delimited field; tolerate a
            # header row by skipping anything that doesn't look like data.
            fields = line.split("\t")
            sample_id = fields[0].strip()
            if sample_id and sample_id.lower() not in ("sample_id", "sample"):
                excluded.add(sample_id)
    return excluded


def parse_genus_profile(profile_path):
    """Return {genus_name: relative_abundance} for one MetaPhlAn4 profile."""
    genus_abund = {}
    with open(profile_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            # vOct22 SGB-era output: clade_name, NCBI_tax_id, relative_abundance, [additional_species]
            if len(fields) < 3:
                continue
            clade_name = fields[0]
            try:
                rel_abund = float(fields[2])
            except ValueError:
                continue
            if GENUS_ROW_RE.search(clade_name):
                genus = clade_name.split("|")[-1]  # e.g. "g__Prevotella"
                # Guard against duplicate genus rows in a single profile
                # (shouldn't happen, but fail loudly rather than silently
                # overwrite if it does)
                if genus in genus_abund:
                    print(f"  WARNING: duplicate genus row '{genus}' in {profile_path} -- summing.")
                    genus_abund[genus] += rel_abund
                else:
                    genus_abund[genus] = rel_abund
    return genus_abund


def main():
    os.makedirs(OUT_ROOT, exist_ok=True)
    per_dataset_dir = os.path.join(OUT_ROOT, "per_dataset")
    os.makedirs(per_dataset_dir, exist_ok=True)

    excluded_samples = load_exclusions(EXCLUSIONS_LOG)
    print(f"Loaded {len(excluded_samples)} excluded sample_id(s) from {EXCLUSIONS_LOG}")

    combined = {}  # combined_sample_id -> {genus: abundance}
    sample_study_map = []  # (combined_sample_id, original_sample_id, dataset_key)
    all_genera = set()

    summary_lines = []
    summary_lines.append("=== PD genus-level merge summary ===\n")

    for dataset_key, paths in DATASETS.items():
        sample_list_path = paths["sample_list"]
        out_dir = paths["out_dir"]

        if not os.path.isfile(sample_list_path):
            print(f"ERROR: sample list not found for {dataset_key}: {sample_list_path}")
            sys.exit(1)

        with open(sample_list_path) as fh:
            samples = [s.strip() for s in fh if s.strip()]

        dataset_table = {}
        n_expected = len(samples)
        n_excluded = 0
        n_missing = 0
        n_ok = 0

        for sample in samples:
            if sample in excluded_samples:
                n_excluded += 1
                continue

            profile_path = os.path.join(out_dir, f"{sample}_metaphlan_profile.txt")
            if not os.path.isfile(profile_path):
                print(f"  WARNING [{dataset_key}]: expected profile missing, skipping: {profile_path}")
                n_missing += 1
                continue

            genus_abund = parse_genus_profile(profile_path)
            if not genus_abund:
                print(f"  WARNING [{dataset_key}]: no genus-level rows parsed for {sample} -- skipping.")
                n_missing += 1
                continue

            combined_sample_id = f"{dataset_key}_{sample}"
            dataset_table[sample] = genus_abund
            combined[combined_sample_id] = genus_abund
            sample_study_map.append((combined_sample_id, sample, dataset_key))
            all_genera.update(genus_abund.keys())
            n_ok += 1

        # Write per-dataset genus table for inspection/debugging
        dataset_genera = sorted({g for v in dataset_table.values() for g in v})
        per_dataset_path = os.path.join(per_dataset_dir, f"{dataset_key}_genus_table.tsv")
        with open(per_dataset_path, "w") as out:
            out.write("sample_id\t" + "\t".join(dataset_genera) + "\n")
            for sample, genus_abund in dataset_table.items():
                row = [f"{genus_abund.get(g, 0.0):.6f}" for g in dataset_genera]
                out.write(sample + "\t" + "\t".join(row) + "\n")

        line = (f"{dataset_key}: expected={n_expected} ok={n_ok} "
                f"excluded={n_excluded} missing/unparseable={n_missing} "
                f"genera_in_dataset={len(dataset_genera)}")
        print(line)
        summary_lines.append(line + "\n")

        if n_missing > 0:
            print(f"  -> {dataset_key} has {n_missing} sample(s) that should have passed "
                  f"pd_02_verify_metaphlan4_oct22.sh but failed to parse here. "
                  f"Investigate before trusting the combined table.")

    # Build combined wide table
    all_genera = sorted(all_genera)
    combined_path = os.path.join(OUT_ROOT, "combined_pd_genus_table.tsv")
    with open(combined_path, "w") as out:
        out.write("sample_id\tstudy_id\t" + "\t".join(all_genera) + "\n")
        for combined_sample_id, orig_sample, dataset_key in sample_study_map:
            genus_abund = combined[combined_sample_id]
            row = [f"{genus_abund.get(g, 0.0):.6f}" for g in all_genera]
            out.write(f"{combined_sample_id}\t{dataset_key}\t" + "\t".join(row) + "\n")

    # Sample -> study map, useful for joining metadata later
    map_path = os.path.join(OUT_ROOT, "sample_study_map.tsv")
    with open(map_path, "w") as out:
        out.write("combined_sample_id\toriginal_sample_id\tstudy_id\n")
        for combined_sample_id, orig_sample, dataset_key in sample_study_map:
            out.write(f"{combined_sample_id}\t{orig_sample}\t{dataset_key}\n")

    total_samples = len(sample_study_map)
    summary_lines.append(f"\nTOTAL combined samples: {total_samples}\n")
    summary_lines.append(f"TOTAL unique genera across all PD datasets: {len(all_genera)}\n")
    summary_lines.append(f"Combined table: {combined_path}\n")
    summary_lines.append(f"Sample/study map: {map_path}\n")

    summary_path = os.path.join(OUT_ROOT, "merge_summary.txt")
    with open(summary_path, "w") as out:
        out.writelines(summary_lines)

    print("\n" + "".join(summary_lines))
    print(f"Done. Combined table written to: {combined_path}")


if __name__ == "__main__":
    main()
    