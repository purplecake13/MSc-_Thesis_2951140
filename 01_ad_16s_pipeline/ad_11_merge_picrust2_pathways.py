#!/usr/bin/env python
import pandas as pd
import os

studies = ["binyinli2019", "cirstea2022", "ling2021", "liu2019", "ueda2021",
           "yamashiro2024", "yildirim2022", "zhuang2018"]
base = "/rds/projects/e/elhamsak-ad-thesis/ad_picrust2"
out_dir = "/rds/projects/e/elhamsak-ad-thesis/ad_picrust2_merged"
os.makedirs(out_dir, exist_ok=True)

dfs = []
for s in studies:
    # ADJUST this path if PICRUSt2 output structure differs --
    # default PICRUSt2 pathway output location is pathways_out/path_abun_unstrat.tsv.gz
    path = f"{base}/{s}/pathways_out/path_abun_unstrat.tsv.gz"
    if not os.path.exists(path):
        print(f"WARNING: not found for {s}: {path}")
        continue
    df = pd.read_csv(path, sep="\t", index_col=0)
    print(f"{s}: {df.shape[0]} pathways x {df.shape[1]} samples")
    dfs.append(df)

if not dfs:
    raise SystemExit("No pathway tables found for any study -- check the path pattern above.")

merged = pd.concat(dfs, axis=1, join="outer").fillna(0)
print(f"\nMerged: {merged.shape[0]} pathways x {merged.shape[1]} samples")

# CPM normalise (matches your PD-side normalisation for consistency)
cpm = merged.div(merged.sum(axis=0), axis=1) * 1e6

merged.to_csv(f"{out_dir}/ad_pathabundance_joined.tsv", sep="\t")
cpm.to_csv(f"{out_dir}/ad_pathabundance_cpm.tsv", sep="\t")
print("Wrote joined + CPM tables.")
