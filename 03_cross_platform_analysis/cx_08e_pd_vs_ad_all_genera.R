#!/usr/bin/env Rscript
# cx_08e_pd_vs_ad_all_genera.R
#
# Extension of cx_08d: rather than restricting to the 7 genera pre-selected
# as independently significant in both PD_vs_HC and AD_vs_HC, this tests
# ALL union genera (283) directly between PD and AD samples. Answers a
# broader question: "which genera differentiate PD from AD directly,"
# rather than "do the genera we already flagged differ."
#
# *** CAVEAT (same as cx_08d, flag prominently in Results/Discussion) ***
# No dataset contains both PD and AD patients -- inherently confounded
# with platform: 100% of PD samples are shotgun, 100% of AD samples are
# 16S. Cannot fully separate a genuine PD-vs-AD disease effect from a
# technology/protocol effect.
#
# Method: relative abundance (% of each sample's FULL genus-level total,
# not just the tested subset -- see cx_08d's row_totals bug/fix for why
# this matters); Wilcoxon rank-sum test (PD vs AD, two-sided) per genus;
# BH correction across all tested genera.
#
# Output:
#   .../differential_abundance/pd_vs_ad_all_genera_test.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

JOINED_PATH <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/joined_genus_metadata.tsv"
OUT_PATH <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/differential_abundance/pd_vs_ad_all_genera_test.csv"

joined <- read_tsv(JOINED_PATH, show_col_types = FALSE)
disease_data <- joined %>% filter(disease_group %in% c("PD", "AD"))

genus_cols <- grep("^g__", colnames(joined), value = TRUE)
cat(sprintf("Testing %d union genera (PD vs AD, direct comparison)\n", length(genus_cols)))
cat("\n*** CAVEAT: PD (shotgun) vs AD (16S) comparison is confounded with platform ***\n\n")

full_mat <- as.matrix(disease_data[, genus_cols])
row_totals <- rowSums(full_mat, na.rm = TRUE)

n_zero_total <- sum(row_totals == 0)
if (n_zero_total > 0) {
  cat(sprintf("WARNING: %d samples have zero total genus abundance -- excluding from this test.\n",
              n_zero_total))
}
keep <- row_totals > 0
disease_data <- disease_data[keep, ]
full_mat <- full_mat[keep, , drop = FALSE]
row_totals <- row_totals[keep]

rel_abund <- sweep(full_mat, 1, row_totals, "/") * 100

# Drop genera with near-zero prevalence in this PD/AD subset -- testing an
# all-zero or near-all-zero column wastes a BH correction "slot" on a
# guaranteed-null result and dilutes power on genuinely testable genera.
prevalence <- colMeans(rel_abund > 0)
testable <- names(prevalence[prevalence >= 0.05])  # present in >=5% of PD+AD samples
cat(sprintf("Dropped %d genera present in <5%% of samples (untestable); %d remain.\n",
            length(genus_cols) - length(testable), length(testable)))

results <- lapply(testable, function(g) {
  vals <- rel_abund[, g]
  test <- tryCatch(
    wilcox.test(vals ~ disease_data$disease_group),
    error = function(e) NULL
  )
  if (is.null(test)) return(NULL)
  data.frame(
    taxon = g,
    median_PD = median(vals[disease_data$disease_group == "PD"]),
    median_AD = median(vals[disease_data$disease_group == "AD"]),
    p_value = test$p.value
  )
}) %>% bind_rows()

results$q_value <- p.adjust(results$p_value, method = "BH")
results <- results %>% arrange(q_value)

write_csv(results, OUT_PATH)
cat(sprintf("\n%d/%d union genera significant at q<0.05 for direct PD-vs-AD comparison.\n",
            sum(results$q_value < 0.05), nrow(results)))
cat(sprintf("Results written to: %s\n", OUT_PATH))
cat("\n*** Report the platform-confound caveat above wherever this table is used. ***\n")
