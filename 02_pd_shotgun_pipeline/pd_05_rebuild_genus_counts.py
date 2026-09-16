#!/usr/bin/env python3
"""
rebuild_pd_genus_counts.py

WHY THIS EXISTS
------------------------------------------------------------------------
joined_genus_metadata.tsv's shotgun rows were built from MetaPhlAn4's
relative_abundance column, pre-normalized to ~100 per sample regardless
of sequencing depth. ANCOM-BC2 needs real inter-sample depth variation
to estimate its per-sample bias-correction term; relative abundance
destroys that signal, which likely explains why PD_vs_HC ANCOM-BC2
found so few / possibly unreliable hits.

Per-sample MetaPhlAn4 profiles were run with -t rel_ab_w_read_stats,
so estimated_number_of_reads_from_the_clade is directly available --
no MetaPhlAn4 rerun or approximation needed. This script:

  1. Reads every per-sample *_metaphlan_profile.txt for all PD datasets
  2. Filters to genus-level rows (clade_name ending in g__<name>,
     i.e. no s__ species suffix)
  3. Pivots to a samples x genera matrix of estimated read counts
  4. Writes pd_genus_counts.tsv, ready to merge with metadata and feed
     into ANCOM-BC2 in place of the relative-abundance columns

Usage:
    python rebuild_pd_genus_counts.py \
        --input-dirs /rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/bedarf2017 \
                     /rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/boktor2023/rumc \
                     /rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/boktor2023/tbc \
                     /rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/clasen2024 \
                     /rds/projects/e/elhamsak-ad-thesis/pd_metaphlan4_oct22/wallen2022 \
                     /rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/mao2021 \
        --output /rds/projects/e/elhamsak-pd-thesis/cross_platform/pd_genus_counts.tsv

Adjust --input-dirs to match your actual per-study MetaPhlAn4 output
directories -- edit the DEFAULT_INPUT_DIRS list below or pass explicitly.
"""

import argparse
import glob
import os
import sys
import pandas as pd

DEFAULT_INPUT_DIRS = [
    "/rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/bedarf2017",
    "/rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/boktor2023_rumc",
    "/rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/boktor2023_tbc",
    "/rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/clasen2024",
    "/rds/projects/e/elhamsak-ad-thesis/pd_metaphlan4_oct22/wallen2022",
    "/rds/projects/e/elhamsak-pd-thesis/pd_metaphlan4_oct22/mao2021",
]


def parse_profile(filepath):
    """Parse one MetaPhlAn4 rel_ab_w_read_stats profile into a
    {genus_name: read_count} dict, summing any duplicate genus rows
    (shouldn't normally happen, but defensive)."""
    genus_counts = {}
    sample_id = os.path.basename(filepath).replace("_metaphlan_profile.txt", "")

    with open(filepath) as f:
        header_idx = None
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#SampleID") or line.startswith("clade_name"):
                header_idx = line.split("\t")
                continue
            if line.startswith("#"):
                continue
            if header_idx is None:
                continue

            fields = line.split("\t")
            if len(fields) < 5:
                continue

            clade_name = fields[0]
            try:
                read_count = float(fields[4])
            except (ValueError, IndexError):
                continue

            # Keep only genus-level rows: last taxonomic level present is
            # g__<name> with no s__ (species) suffix after it.
            levels = clade_name.split("|")
            last_level = levels[-1]
            if not last_level.startswith("g__"):
                continue
            if len(last_level) <= 3:
                # g__ with empty name -- skip
                continue

            genus_name = last_level  # e.g. "g__Prevotella"
            genus_counts[genus_name] = genus_counts.get(genus_name, 0.0) + read_count

    return sample_id, genus_counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dirs", nargs="+", default=DEFAULT_INPUT_DIRS,
                     help="Directories containing *_metaphlan_profile.txt files (searched recursively)")
    ap.add_argument("--output", required=True, help="Output TSV path (samples x genera, read counts)")
    ap.add_argument("--pattern", default="*_metaphlan_profile.txt",
                     help="Filename glob pattern for per-sample profiles")
    args = ap.parse_args()

    all_files = []
    for d in args.input_dirs:
        if not os.path.isdir(d):
            print(f"WARNING: input dir not found, skipping: {d}", file=sys.stderr)
            continue
        found = glob.glob(os.path.join(d, "**", args.pattern), recursive=True)
        print(f"[{os.path.basename(d.rstrip('/'))}] found {len(found)} profile files")
        all_files.extend(found)

    if not all_files:
        print("ERROR: no profile files found across any input dir. Check --input-dirs.", file=sys.stderr)
        sys.exit(1)

    print(f"\nTotal profile files to parse: {len(all_files)}")

    records = {}
    failed = []
    for fp in all_files:
        try:
            sample_id, genus_counts = parse_profile(fp)
            if not genus_counts:
                failed.append((fp, "no genus-level rows parsed"))
                continue
            records[sample_id] = genus_counts
        except Exception as e:
            failed.append((fp, str(e)))

    if failed:
        print(f"\nWARNING: {len(failed)} files failed to parse or yielded no genus rows:")
        for fp, reason in failed[:20]:
            print(f"  {fp}: {reason}")
        if len(failed) > 20:
            print(f"  ... and {len(failed) - 20} more")

    print(f"\nSuccessfully parsed {len(records)} samples")

    df = pd.DataFrame.from_dict(records, orient="index").fillna(0.0)
    df.index.name = "unique_sample_id"

    # Sanity check: row sums should now vary with sequencing depth,
    # not cluster near 100
    row_sums = df.sum(axis=1)
    print(f"\nRow sum range: {row_sums.min():.1f} to {row_sums.max():.1f}")
    print(f"Row sum mean: {row_sums.mean():.1f}, median: {row_sums.median():.1f}")
    if row_sums.max() < 200:
        print("WARNING: row sums still look like relative abundance, not counts -- "
              "check that field index 4 (estimated_number_of_reads_from_the_clade) "
              "was parsed correctly.")

    df.to_csv(args.output, sep="\t")
    print(f"\nWrote {df.shape[0]} samples x {df.shape[1]} genera to {args.output}")
    print("\nNEXT STEP: this table's sample IDs are the raw MetaPhlAn accession-based "
          "IDs (e.g. SRR19064316), not your unique_sample_id join key. Join this to "
          "master_metadata.tsv on run_accession/sample_id per your existing fallback-join "
          "convention, then restrict cx_08's PD_vs_HC block to source genus columns from "
          "this table instead of joined_genus_metadata.tsv's relative-abundance columns.")


if __name__ == "__main__":
    main()
    