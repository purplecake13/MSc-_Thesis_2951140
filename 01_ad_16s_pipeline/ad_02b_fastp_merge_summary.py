#!/usr/bin/env python3
"""
/rds/projects/e/elhamsak-pd-thesis/scripts/ad_02b_fastp_merge_summary.py — aggregate fastp merge-mode JSON reports

Run after all 11_fastp_merge_array.sh array tasks complete. Produces a
per-sample summary table (total reads, merged reads, merge rate %) plus
an overall average, in the same format used to compare the DADA2 and
VSEARCH-merge baselines for Tran2019.

Usage:
    python3 11b_fastp_merge_summary.py --json-dir /path/to/fastp_merge/output --out summary.tsv
"""

import argparse
import json
import glob
import os
import csv


def main():
    parser = argparse.ArgumentParser(description="Summarise fastp merge-mode JSON reports")
    parser.add_argument("--json-dir", required=True, help="Directory containing *_fastp_merge.json files")
    parser.add_argument("--out", default="fastp_merge_summary.tsv", help="Output TSV path")
    args = parser.parse_args()

    json_files = sorted(glob.glob(os.path.join(args.json_dir, "*_fastp_merge.json")))
    if not json_files:
        print(f"ERROR: no *_fastp_merge.json files found in {args.json_dir}")
        return

    rows = []
    for jf in json_files:
        sample = os.path.basename(jf).replace("_fastp_merge.json", "")
        with open(jf) as f:
            data = json.load(f)

        # fastp 0.24.0 JSON structure (confirmed against real output, 23 June 2026):
        # summary.before_filtering.total_reads = R1 + R2 combined count (i.e. 2x read pairs)
        # summary.after_filtering.total_reads  = merged single-read output count, when -m
        #   merge mode is used (confirmed by after_filtering.read1_mean_length jumping to
        #   ~400+, consistent with a merged V3-V4 amplicon read, not a raw ~230bp read)
        summary = data.get("summary", {})
        total_reads_before = summary.get("before_filtering", {}).get("total_reads", None)
        merged_reads = summary.get("after_filtering", {}).get("total_reads", None)

        # total input READ PAIRS = total_reads_before / 2 (fastp counts R1+R2 separately)
        input_pairs = (total_reads_before / 2) if total_reads_before else None
        merge_rate = (100 * merged_reads / input_pairs) if (input_pairs and merged_reads is not None) else None

        rows.append({
            "sample": sample,
            "input_read_pairs": int(input_pairs) if input_pairs else "NA",
            "merged_reads": merged_reads if merged_reads is not None else "NA",
            "merge_rate_pct": round(merge_rate, 2) if merge_rate is not None else "NA",
        })

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample", "input_read_pairs", "merged_reads", "merge_rate_pct"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    valid_rates = [r["merge_rate_pct"] for r in rows if r["merge_rate_pct"] != "NA"]
    if valid_rates:
        avg_rate = sum(valid_rates) / len(valid_rates)
        print(f"Summary written to: {args.out}")
        print(f"Samples processed:  {len(rows)}")
        print(f"Mean merge rate:    {avg_rate:.2f}%")
        print(f"Min / Max:          {min(valid_rates):.2f}% / {max(valid_rates):.2f}%")
        print()
        print("VERIFICATION STEP (do this before trusting the numbers above):")
        print("Cross-check the first sample's JSON-derived merged-read count against an")
        print("actual line count of its merged FASTQ output, to confirm")
        print("'after_filtering.total_reads' really represents merged reads and not")
        print("merged+unmerged combined. Run:")
        first_sample = rows[0]["sample"]
        print(f"  zcat {args.json_dir}/{first_sample}_merged.fastq.gz | wc -l")
        print("  # divide by 4 (FASTQ = 4 lines per read) and compare to this script's")
        print(f"  # reported merged_reads value for {first_sample} ({rows[0]['merged_reads']})")
        print("  # — these should match exactly. If they don't, the field interpretation")
        print("  # is wrong and the merge rates above should not be trusted yet.")
    else:
        print("WARNING: no valid merge rates extracted. Inspect one JSON file directly:")
        print(f"  python3 -m json.tool {json_files[0]} | less")
        print("and check the key names under 'summary.before_filtering' and 'summary.after_filtering'.")


if __name__ == "__main__":
    main()