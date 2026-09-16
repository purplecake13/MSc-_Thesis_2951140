#!/usr/bin/env Rscript
# cx_12_pd_vs_ad_pathway_comparison.R
#
# Direct PD-vs-AD comparison at the FUNCTIONAL (pathway) level, mirroring
# cx_08e's genus-level approach.
#
# *** CAVEAT -- STRONGER than the genus-level comparison (cx_08d/cx_08e) ***
# This confounds TWO things simultaneously:
#   1. Platform: PD is shotgun, AD is 16S (same confound as genus-level test)
#   2. Data type: PD pathways are HUMAnN3-derived (measured directly from
#      shotgun reads), while AD pathways are PICRUSt2-derived (PREDICTED
#      from 16S taxonomy via reference genome inference, not measured).
# A difference here could reflect genuine biology, platform effects, OR
# systematic differences between measured and predicted functional profiles.
# Treat as strictly exploratory/hypothesis-generating; state both caveats
# explicitly wherever this table is used.
#
# Method: pathways are matched by MetaCyc ID (PD table embeds descriptions
# in the ID string, e.g. "P164-PWY: purine nucleobases degradation I
# (anaerobic)" -- stripped to bare ID before matching against AD's bare
# PICRUSt2 IDs). Both tables are already CPM-normalised (comparable within
# each platform's own scale), so no additional row-total step is needed
# here (unlike the genus-level test, which used raw counts).
# Wilcoxon rank-sum test (PD vs AD, two-sided) per shared pathway; BH
# correction across the shared-pathway set.
#
# Output:
#   .../differential_abundance/pd_vs_ad_pathway_test.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

PD_PATHWAY_PATH <- "/rds/projects/e/elhamsak-pd-thesis/pd_humann3_merged/pd_pathabundance_cpm_unstratified.tsv"
AD_PATHWAY_PATH <- "/rds/projects/e/elhamsak-ad-thesis/ad_picrust2_merged/ad_pathabundance_cpm.tsv"
OUT_PATH <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/differential_abundance/pd_vs_ad_pathway_test.csv"

# ---- Load PD pathways (HUMAnN3, measured) ----------------------------------
# pd_raw <- read_tsv(PD_PATHWAY_PATH, comment = "#", col_names = TRUE, show_col_types = FALSE)
pd_raw <- read_tsv(PD_PATHWAY_PATH, comment = "", col_names = TRUE, show_col_types = FALSE)
colnames(pd_raw)[1] <- "pathway_full"
colnames(pd_raw) <- gsub("_concat_Abundance$", "", colnames(pd_raw))

# Strip embedded description -> bare MetaCyc ID (e.g. "P164-PWY: purine..." -> "P164-PWY")
pd_raw <- pd_raw %>%
  mutate(pathway_id = trimws(sub(":.*$", "", pathway_full))) %>%
  filter(!pathway_id %in% c("UNMAPPED", "UNINTEGRATED")) %>%
  select(-pathway_full) %>%
  distinct(pathway_id, .keep_all = TRUE)  # guard against rare ID collisions after stripping

pd_mat <- pd_raw %>% select(-pathway_id) %>% as.matrix()
rownames(pd_mat) <- pd_raw$pathway_id

# ---- Load AD pathways (PICRUSt2, predicted) --------------------------------
ad_raw <- read_tsv(AD_PATHWAY_PATH, show_col_types = FALSE)
colnames(ad_raw)[1] <- "pathway_id"
ad_raw <- ad_raw %>%
  filter(!pathway_id %in% c("UNMAPPED", "UNINTEGRATED")) %>%
  distinct(pathway_id, .keep_all = TRUE)

ad_mat <- ad_raw %>% select(-pathway_id) %>% as.matrix()
rownames(ad_mat) <- ad_raw$pathway_id

# ---- Intersect on pathway ID ------------------------------------------------
shared_pathways <- intersect(rownames(pd_mat), rownames(ad_mat))
cat(sprintf("PD pathways: %d | AD pathways: %d | Shared: %d\n",
            nrow(pd_mat), nrow(ad_mat), length(shared_pathways)))

if (length(shared_pathways) < 10) {
  stop("Fewer than 10 shared pathways found -- check ID-stripping regex against actual PD pathway_full values (head(pd_raw$pathway_id) before the intersect) before proceeding.")
}

cat("\n*** CAVEAT: PD (HUMAnN3, measured, shotgun) vs AD (PICRUSt2, PREDICTED, 16S) ***\n")
cat("*** comparison confounds platform AND measured-vs-predicted data type ***\n\n")

# ---- Test each shared pathway ------------------------------------------------
results <- lapply(shared_pathways, function(p) {
  pd_vals <- as.numeric(pd_mat[p, ])
  ad_vals <- as.numeric(ad_mat[p, ])
  test <- tryCatch(
    wilcox.test(pd_vals, ad_vals),
    error = function(e) NULL
  )
  if (is.null(test)) return(NULL)
  data.frame(
    pathway_id = p,
    median_PD = median(pd_vals),
    median_AD = median(ad_vals),
    p_value = test$p.value
  )
}) %>% bind_rows()

results$q_value <- p.adjust(results$p_value, method = "BH")
results <- results %>% arrange(q_value)

write_csv(results, OUT_PATH)
cat(sprintf("\n%d/%d shared pathways significant at q<0.05 for direct PD-vs-AD comparison.\n",
            sum(results$q_value < 0.05), nrow(results)))
cat(sprintf("Results written to: %s\n", OUT_PATH))
cat("\n*** Report BOTH caveats above wherever this table is used. ***\n")
