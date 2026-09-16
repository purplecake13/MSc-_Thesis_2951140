#!/bin/bash
MAPFILE=/rds/bear-apps/2023a/EL8-ice/software/MetaPhlAn/4.1.1-foss-2023a/lib/python3.11/site-packages/metaphlan/utils/mpa_vOct22_CHOCOPhlAnSGB_202212_SGB2GTDB.tsv
MPAFILE=/rds/projects/e/elhamsak-pd-thesis/cross_platform/merged_pd_metaphlan.txt

# for ggb in GGB9350 GGB36267 GGB9345 GGB33512 GGB9699 GGB9619 GGB9695 GGB9694 GGB9781 GGB9770 GGB3746 GGB45432 GGB9775 GGB2980; do
# for ggb in GGB2980 GGB33512 GGB36267 GGB3746 GGB45432 GGB9345 GGB9350 GGB9619 GGB9694 GGB9695 GGB9699 GGB9770 GGB9775 GGB9781; do
for ggb in GGB9350 GGB36267 GGB9345 GGB33512 GGB9699 GGB9619 GGB9695 GGB9694 GGB9781 GGB9770 GGB3746 GGB45432 GGB9775 GGB2980 GGB9770; do
  echo "=== $ggb ==="
  sgbs=$(grep "g__${ggb}|" "$MPAFILE" | grep -oP 't__SGB\K[0-9]+' | sort -u)
  genera=""
  for sgb in $sgbs; do
    hit=$(grep -P "^SGB${sgb}\t" "$MAPFILE")
    if [ -z "$hit" ]; then
      echo "  SGB${sgb}: NOT FOUND in mapping file"
    else
      genus=$(echo "$hit" | grep -oP 'g__\K[^;]*')
      echo "  SGB${sgb} -> g__${genus:-<blank>}"
      genera="${genera}${genus}|"
    fi
  done
  uniq_genera=$(echo "$genera" | tr '|' '\n' | sort -u | grep -v '^$')
  n_uniq=$(echo "$uniq_genera" | grep -c .)
  if [ "$n_uniq" -eq 1 ]; then
    echo "  --> CONSISTENT: $uniq_genera"
  elif [ "$n_uniq" -eq 0 ]; then
    echo "  --> UNRESOLVED: no genus assignment found for any SGB"
  else
    echo "  --> AMBIGUOUS: multiple genera ($uniq_genera)"
  fi
  echo ""
done
