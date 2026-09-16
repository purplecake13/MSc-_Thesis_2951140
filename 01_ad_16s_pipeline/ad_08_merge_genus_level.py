#!/usr/bin/env python3
"""
ad_08_merge_genus_level.py

Merges all 9 AD 16S genus-collapsed QIIME2 feature tables (produced by
ad_07_genus_collapse.sh) into a single samples x genera feature table,
the AD-side counterpart to pd_05_merge_genus_level.py.

Usage:
    python ad_08_merge_genus_level.py

Reads the same dataset -> genus_level directory mapping used by
ad_07_verify_genus_collapse.sh, so this should only be run after that
verify script shows all 9 datasets as OK.

INPUT FORMAT (biom convert --to-tsv output):
    line 1: "# Constructed from biom file" (comment, skipped)
    line 2: "#OTU ID<TAB>sample1<TAB>sample2..." (header)
    line 3+: "<taxonomy_string><TAB>count1<TAB>count2..."

Taxonomy strings are full lineages collapsed at genus level (L6), using
either SILVA-style "g__Genus" or Greengenes-style "D_5__Genus" prefixes
depending on dataset/classifier. This script extracts just the trailing
genus label from each, normalising both conventions to "g__Genus" so
genus keys line up with the PD-side MetaPhlAn4 table for the later
cross-platform combine.

VALUES ARE RAW COUNTS, not relative abundance (unlike the PD/MetaPhlAn4
table) -- this is expected and intentional. Per the locked design, CLR
normalisation is applied separately per platform, later, during the
cross-study merge step -- NOT here. Do not treat AD and PD values as
directly comparable until that normalisation has been applied.

Respects sample_exclusions_log.tsv (shared PD/AD log) the same way the
PD merge script does -- any sample_id listed there is skipped.

Outputs (written to OUT_ROOT):
    ad_genus_level/combined_ad_genus_table.tsv   samples x genera, raw counts, filled 0
    ad_genus_level/per_dataset/<dataset>_genus_table.tsv   per-dataset tables
    ad_genus_level/sample_study_map.tsv          sample_id -> study_id
    ad_genus_level/merge_summary.txt             run summary / sanity stats
"""

import os
import re
import sys

PD_BASE = "/rds/projects/e/elhamsak-pd-thesis"
AD_BASE = "/rds/projects/e/elhamsak-ad-thesis"
AD_ROOT = f"{AD_BASE}/ad_qiime2"
OUT_ROOT = f"{AD_BASE}/ad_genus_level"
EXCLUSIONS_LOG = f"{PD_BASE}/sample_exclusions_log.tsv"  # shared log, same as PD side

# Same 9 datasets / directory layout as ad_07_verify_genus_collapse.sh
DATASETS = {
    "cirstea2022": "cirstea2022/genus_level",
    "yamashiro2024": "yamashiro2024/genus_level",
    "tran2019": "tran2019/genus_level",
    "binyinli2019": "binyinli2019/genus_level",
    "ling2021": "ling2021/genus_level",
    "liu2019": "liu2019/genus_level",
    "ueda2021": "ueda2021/genus_level",
    "yildirim2022": "yildirim2022/genus_level",
    "zhuang2018": "zhuang2018/genus_level",
}

# Matches a trailing genus segment under either convention:
#   "...; g__Blautia"      (SILVA-style)
#   "...; D_5__Blautia"    (Greengenes-style naive-bayes classifier output)
# Also matches empty/unresolved genus calls, e.g. "...; g__" or "...; D_5__"
GENUS_SEGMENT_RE = re.compile(r"(?:g__|D_5__)([^;]*)$")


def load_exclusions(path):
    excluded = set()
    if not os.path.isfile(path):
        print(f"NOTE: exclusions log not found at {path} -- proceeding with no exclusions applied.")
        return excluded
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            sample_id = line.split("\t")[0].strip()
            if sample_id and sample_id.lower() not in ("sample_id", "sample"):
                excluded.add(sample_id)
    return excluded


def normalise_genus_label(taxonomy_string):
    """
    Extract the genus-level segment from a full QIIME2 taxonomy lineage
    string and normalise to 'g__Genus' regardless of source prefix
    convention. Returns None if no genus-level segment is found at all
    (shouldn't happen on a properly L6-collapsed table, but don't guess).
    """
    match = GENUS_SEGMENT_RE.search(taxonomy_string)
    if not match:
        return None
    genus_name = match.group(1).strip()
    return f"g__{genus_name}" if genus_name else "g__unclassified"


def parse_genus_tsv(tsv_path):
    """
    Parse a biom-converted feature-table.tsv into {sample_id: {genus: count}}.
    Returns (per_sample_dict, sample_order, warnings).
    """
    warnings = []
    with open(tsv_path) as fh:
        lines = fh.readlines()

    if len(lines) < 3:
        warnings.append(f"File has fewer than 3 lines (expected comment+header+data): {tsv_path}")
        return {}, [], warnings

    # line 0: "# Constructed from biom file" -- skip
    # line 1: header "#OTU ID\tsample1\tsample2..."
    header = lines[1].rstrip("\n").split("\t")
    sample_ids = header[1:]

    per_sample = {s: {} for s in sample_ids}

    for line_no, line in enumerate(lines[2:], start=3):
        line = line.rstrip("\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != len(sample_ids) + 1:
            warnings.append(f"{tsv_path}:{line_no} field count mismatch "
                             f"(expected {len(sample_ids)+1}, got {len(fields)}) -- skipping row.")
            continue
        taxonomy_string = fields[0]
        genus_label = normalise_genus_label(taxonomy_string)
        if genus_label is None:
            warnings.append(f"{tsv_path}:{line_no} no genus-level segment found in "
                             f"'{taxonomy_string}' -- skipping row.")
            continue
        for sample_id, raw_count in zip(sample_ids, fields[1:]):
            try:
                count = float(raw_count)
            except ValueError:
                count = 0.0
            if count == 0.0:
                continue
            if genus_label in per_sample[sample_id]:
                per_sample[sample_id][genus_label] += count
            else:
                per_sample[sample_id][genus_label] = count

    return per_sample, sample_ids, warnings


def main():
    os.makedirs(OUT_ROOT, exist_ok=True)
    per_dataset_dir = os.path.join(OUT_ROOT, "per_dataset")
    os.makedirs(per_dataset_dir, exist_ok=True)

    excluded_samples = load_exclusions(EXCLUSIONS_LOG)
    print(f"Loaded {len(excluded_samples)} excluded sample_id(s) from {EXCLUSIONS_LOG}")

    combined = {}  # combined_sample_id -> {genus: count}
    sample_study_map = []
    all_genera = set()
    summary_lines = ["=== AD genus-level merge summary ===\n"]

    for dataset_key, rel_dir in DATASETS.items():
        tsv_path = os.path.join(AD_ROOT, rel_dir, "genus_table_export", "feature-table.tsv")

        if not os.path.isfile(tsv_path):
            print(f"ERROR: genus table not found for {dataset_key}: {tsv_path}")
            print("  -> run/re-check ad_07_genus_collapse.sh and ad_07_verify_genus_collapse.sh for this dataset first.")
            sys.exit(1)

        per_sample, sample_ids, warnings = parse_genus_tsv(tsv_path)
        for w in warnings:
            print(f"  WARNING [{dataset_key}]: {w}")

        n_expected = len(sample_ids)
        n_excluded = 0
        n_ok = 0

        dataset_table = {}
        for sample_id, genus_counts in per_sample.items():
            if sample_id in excluded_samples:
                n_excluded += 1
                continue
            combined_sample_id = f"{dataset_key}_{sample_id}"
            dataset_table[sample_id] = genus_counts
            combined[combined_sample_id] = genus_counts
            sample_study_map.append((combined_sample_id, sample_id, dataset_key))
            all_genera.update(genus_counts.keys())
            n_ok += 1

        dataset_genera = sorted({g for v in dataset_table.values() for g in v})
        per_dataset_path = os.path.join(per_dataset_dir, f"{dataset_key}_genus_table.tsv")
        with open(per_dataset_path, "w") as out:
            out.write("sample_id\t" + "\t".join(dataset_genera) + "\n")
            for sample_id, genus_counts in dataset_table.items():
                row = [f"{genus_counts.get(g, 0.0):.1f}" for g in dataset_genera]
                out.write(sample_id + "\t" + "\t".join(row) + "\n")

        line = (f"{dataset_key}: expected={n_expected} ok={n_ok} "
                f"excluded={n_excluded} genera_in_dataset={len(dataset_genera)}")
        print(line)
        summary_lines.append(line + "\n")

    all_genera = sorted(all_genera)
    combined_path = os.path.join(OUT_ROOT, "combined_ad_genus_table.tsv")
    with open(combined_path, "w") as out:
        out.write("sample_id\tstudy_id\t" + "\t".join(all_genera) + "\n")
        for combined_sample_id, orig_sample, dataset_key in sample_study_map:
            genus_counts = combined[combined_sample_id]
            row = [f"{genus_counts.get(g, 0.0):.1f}" for g in all_genera]
            out.write(f"{combined_sample_id}\t{dataset_key}\t" + "\t".join(row) + "\n")

    map_path = os.path.join(OUT_ROOT, "sample_study_map.tsv")
    with open(map_path, "w") as out:
        out.write("combined_sample_id\toriginal_sample_id\tstudy_id\n")
        for combined_sample_id, orig_sample, dataset_key in sample_study_map:
            out.write(f"{combined_sample_id}\t{orig_sample}\t{dataset_key}\n")

    total_samples = len(sample_study_map)
    summary_lines.append(f"\nTOTAL combined samples: {total_samples}\n")
    summary_lines.append(f"TOTAL unique genera across all AD datasets: {len(all_genera)}\n")
    summary_lines.append("NOTE: values are RAW COUNTS, not relative abundance. "
                          "CLR normalisation happens later, per-platform, at the cross-study merge step.\n")
    summary_lines.append(f"Combined table: {combined_path}\n")
    summary_lines.append(f"Sample/study map: {map_path}\n")

    summary_path = os.path.join(OUT_ROOT, "merge_summary.txt")
    with open(summary_path, "w") as out:
        out.writelines(summary_lines)

    print("\n" + "".join(summary_lines))
    print(f"Done. Combined table written to: {combined_path}")


if __name__ == "__main__":
    main()
    