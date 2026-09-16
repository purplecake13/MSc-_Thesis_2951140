#!/usr/bin/env bash
# ad_11_merge_picrust2_pathways.sh
# Merges path_abun_unstrat.tsv.gz across all 8 AD PICRUSt2 outputs into one
# pathway x sample table (outer join on pathway ID, zero-fill missing).

STUDIES=(binyinli2019 cirstea2022 ling2021 liu2019 ueda2021 yamashiro2024 yildirim2022 zhuang2018)
BASE=/rds/projects/e/elhamsak-ad-thesis/ad_picrust2
OUT=/rds/projects/e/elhamsak-ad-thesis/ad_picrust2_merged
mkdir -p "$OUT"

# Quick pre-flight: confirm expected file exists for at least one study before
# invoking Python, so a path-naming mismatch fails fast and obviously.
CHECK_PATH="$BASE/binyinli2019/pathways_out/path_abun_unstrat.tsv.gz"
if [ ! -f "$CHECK_PATH" ]; then
  echo "WARNING: expected file not found at $CHECK_PATH"
  echo "Run: find $BASE/binyinli2019 -iname '*path_abun*'  to locate the actual filename/structure."
fi

/rds/projects/e/elhamsak-pd-thesis/envs/py_analysis/bin/python /rds/projects/e/elhamsak-pd-thesis/scripts/ad_11_merge_picrust2_pathways.py
