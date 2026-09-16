#!/usr/bin/env python3
"""
cx_01_filter_normalize.py

Applies rare-taxon filtering and CLR normalisation SEPARATELY per platform
to the PD (MetaPhlAn4, relative abundance) and AD (QIIME2, raw counts)
combined genus-level tables, per the locked design:
    - filter: remove genus if prevalence <10% of samples OR mean relative
      abundance <0.01%
    - CLR: clr(x) = log(x / geometric_mean(x)), calculated pseudocount 
      (half the minimum non-zero value in each platform's filtered table) before CLR
    - filtering and CLR applied per-platform, independently, before any
      cross-platform merge

Usage:
    python cx_01_filter_normalize.py

Inputs:
    /rds/projects/e/elhamsak-pd-thesis/pd_genus_level/combined_pd_genus_table.tsv
    /rds/projects/e/elhamsak-ad-thesis/ad_genus_level/combined_ad_genus_table.tsv

Outputs (written to OUT_ROOT):
    filtered_normalized/pd_genus_filtered_clr.tsv
    filtered_normalized/ad_genus_filtered_clr.tsv
    filtered_normalized/pd_genus_filtered_raw.tsv   (post-filter, pre-CLR, for inspection)
    filtered_normalized/ad_genus_filtered_raw.tsv
    filtered_normalized/filter_normalize_summary.txt
"""

import os
import math

PD_BASE = "/rds/projects/e/elhamsak-pd-thesis"
AD_BASE = "/rds/projects/e/elhamsak-ad-thesis"

PD_INPUT = f"{PD_BASE}/pd_genus_level/combined_pd_genus_table.tsv"
AD_INPUT = f"{AD_BASE}/ad_genus_level/combined_ad_genus_table.tsv"
OUT_ROOT = f"{PD_BASE}/filtered_normalized"

# Excluded from MAIN analysis but preserved in ad_genus_level/per_dataset/
# for supplementary APOE4 work: tran2019's prospective design is incompatible
# with the symptomatic-AD case definition used elsewhere (same logic as the
# Ferreiro2025/Jia2025 exclusion from the AD shotgun arm).
MAIN_ANALYSIS_EXCLUDE_STUDIES = {"tran2019", "ueda2021"}

PREVALENCE_THRESHOLD = 0.10      # <10% of samples -> candidate for removal
MEAN_ABUNDANCE_PCT_THRESHOLD = 0.01  # <0.01% mean relative abundance -> candidate for removal


def read_table(path, exclude_studies=None):
    """Read a sample_id/study_id/genus... TSV into (samples, study_ids, genera, matrix[sample][genus])."""
    exclude_studies = exclude_studies or set()
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        genera = header[2:]
        samples = []
        study_ids = []
        matrix = {}
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            sample_id = fields[0]
            study_id = fields[1]
            if study_id in exclude_studies:
                continue
            values = [float(v) for v in fields[2:]]
            samples.append(sample_id)
            study_ids.append(study_id)
            matrix[sample_id] = dict(zip(genera, values))
    return samples, study_ids, genera, matrix


def compute_relative_abundance_pct(matrix, samples, genera, already_relative):
    """
    Return {sample: {genus: rel_abund_pct}}.
    If already_relative (PD/MetaPhlAn), values are used as-is (assumed
    already 0-100 percent scale). Otherwise (AD/counts), convert to percent
    of that sample's total.
    """
    rel = {}
    for s in samples:
        row = matrix[s]
        if already_relative:
            rel[s] = row
        else:
            total = sum(row.values())
            if total == 0:
                rel[s] = {g: 0.0 for g in genera}
            else:
                rel[s] = {g: (v / total) * 100.0 for g, v in row.items()}
    return rel


def filter_genera(samples, genera, rel_abund):
    """Return the list of genera passing BOTH prevalence and mean-abundance thresholds."""
    n = len(samples)
    kept = []
    for g in genera:
        present = sum(1 for s in samples if rel_abund[s].get(g, 0.0) > 0)
        prevalence = present / n if n > 0 else 0.0
        mean_abund = sum(rel_abund[s].get(g, 0.0) for s in samples) / n if n > 0 else 0.0
        # Remove if prevalence <10% OR mean abundance <0.01% -> keep only if BOTH pass
        if prevalence >= PREVALENCE_THRESHOLD and mean_abund >= MEAN_ABUNDANCE_PCT_THRESHOLD:
            kept.append(g)
    return kept


def clr_transform(matrix, samples, genera, pseudocount=None):
    """
    Per-sample CLR: clr(x_i) = log((x_i + pc) / geometric_mean(x + pc)).
    If pseudocount is None, uses half the minimum non-zero value observed
    across the filtered table -- a standard data-driven choice that scales
    correctly to each platform's native units (percent for PD, raw counts
    for AD), rather than imposing a fixed value that only suits one scale.
    """
    if pseudocount is None:
        nonzero_vals = [matrix[s].get(g, 0.0) for s in samples for g in genera
                         if matrix[s].get(g, 0.0) > 0]
        if not nonzero_vals:
            raise ValueError("No non-zero values found -- cannot compute data-driven pseudocount.")
        pseudocount = min(nonzero_vals) / 2.0
        print(f"  Data-driven pseudocount: {pseudocount:.8f} (half of min non-zero value)")

    clr = {}
    for s in samples:
        row = matrix[s]
        vals = [row.get(g, 0.0) + pseudocount for g in genera]
        log_vals = [math.log(v) for v in vals]
        gm_log = sum(log_vals) / len(log_vals)
        clr[s] = {g: (log_vals[i] - gm_log) for i, g in enumerate(genera)}
    return clr, pseudocount


def write_table(path, samples, study_ids, genera, matrix, fmt="{:.6f}"):
    with open(path, "w") as out:
        out.write("sample_id\tstudy_id\t" + "\t".join(genera) + "\n")
        for s, study in zip(samples, study_ids):
            row = [fmt.format(matrix[s].get(g, 0.0)) for g in genera]
            out.write(f"{s}\t{study}\t" + "\t".join(row) + "\n")


def summarize_clr(matrix, samples, genera, label):
    """Quick sanity stats on CLR output -- flat printout, not a file write."""
    all_vals = [matrix[s][g] for s in samples for g in genera]
    if not all_vals:
        return f"{label}: no values to summarize."
    mean_v = sum(all_vals) / len(all_vals)
    min_v = min(all_vals)
    max_v = max(all_vals)
    return (f"{label} CLR summary: n_values={len(all_vals)} mean={mean_v:.4f} "
            f"min={min_v:.4f} max={max_v:.4f}")


def process_platform(label, input_path, already_relative, exclude_studies=None):
    print(f"\n=== Processing {label} ===")
    samples, study_ids, genera, matrix = read_table(input_path, exclude_studies=exclude_studies)
    print(f"{label}: {len(samples)} samples, {len(genera)} genera (pre-filter)")

    rel_abund = compute_relative_abundance_pct(matrix, samples, genera, already_relative)
    kept_genera = filter_genera(samples, genera, rel_abund)
    n_dropped = len(genera) - len(kept_genera)
    print(f"{label}: filtering kept {len(kept_genera)}/{len(genera)} genera "
          f"(dropped {n_dropped}, prevalence>={PREVALENCE_THRESHOLD*100:.0f}% "
          f"AND mean_abund>={MEAN_ABUNDANCE_PCT_THRESHOLD}% required to keep)")

    # Filtered raw (post-filter, original units, pre-CLR) -- for inspection
    filtered_raw_path = os.path.join(OUT_ROOT, f"{label.lower()}_genus_filtered_raw.tsv")
    write_table(filtered_raw_path, samples, study_ids, kept_genera, matrix, fmt="{:.6f}")

    # CLR on filtered data -- data-driven pseudocount, computed separately per platform
    clr_matrix, used_pseudocount = clr_transform(matrix, samples, kept_genera)
    print(f"{label}: used pseudocount={used_pseudocount:.8f}")
    clr_path = os.path.join(OUT_ROOT, f"{label.lower()}_genus_filtered_clr.tsv")
    write_table(clr_path, samples, study_ids, kept_genera, clr_matrix, fmt="{:.6f}")

    summary = summarize_clr(clr_matrix, samples, kept_genera, label)
    print(summary)

    return {
        "label": label,
        "n_samples": len(samples),
        "n_genera_pre_filter": len(genera),
        "n_genera_post_filter": len(kept_genera),
        "n_dropped": n_dropped,
        "clr_summary": summary,
        "filtered_raw_path": filtered_raw_path,
        "clr_path": clr_path,
        "pseudocount": used_pseudocount,
    }


def main():
    os.makedirs(OUT_ROOT, exist_ok=True)

    pd_result = process_platform("PD", PD_INPUT, already_relative=False)
    ad_result = process_platform("AD", AD_INPUT, already_relative=False, exclude_studies=MAIN_ANALYSIS_EXCLUDE_STUDIES)

    summary_path = os.path.join(OUT_ROOT, "filter_normalize_summary.txt")
    with open(summary_path, "w") as out:
        out.write("=== Rare-taxon filtering + per-platform CLR summary ===\n\n")
        if MAIN_ANALYSIS_EXCLUDE_STUDIES:
            out.write(
                f"AD side excludes {sorted(MAIN_ANALYSIS_EXCLUDE_STUDIES)} from main "
                f"analysis (per MAIN_ANALYSIS_EXCLUDE_STUDIES) -> AD n={ad_result['n_samples']}.\n\n"
            )
        else:
            out.write(f"AD side uses all datasets, no exclusions applied -> AD n={ad_result['n_samples']}.\n\n")
        for result in (pd_result, ad_result):
            out.write(f"{result['label']}:\n")
            out.write(f"  samples: {result['n_samples']}\n")
            out.write(f"  genera pre-filter: {result['n_genera_pre_filter']}\n")
            out.write(f"  genera post-filter: {result['n_genera_post_filter']} "
                       f"(dropped {result['n_dropped']})\n")
            out.write(f"  {result['clr_summary']}\n")
            out.write(f"  pseudocount used: {result['pseudocount']:.8f} "
                       f"(data-driven: half of min non-zero value, computed separately per platform)\n")
            out.write(f"  filtered raw table: {result['filtered_raw_path']}\n")
            out.write(f"  CLR table: {result['clr_path']}\n\n")
    print(f"\nSummary written to: {summary_path}")


if __name__ == "__main__":
    main()
    