#!/usr/bin/env Rscript
# cx_08d_pd_vs_ad_shared_significance.R
#
# Per Sakhaee (supervisor meeting): test whether abundance of the shared
# significant genera (significant in BOTH PD_vs_HC and AD_vs_HC, same
# direction -- cx_08's shared_genera.csv) differs directly between PD and
# AD samples.
#
# *** CAVEAT (flag prominently in Results/Discussion, raise with Sakhaee
# if not already discussed) ***
# No dataset in either arm contains both PD and AD patients, so this
# comparison is inherently confounded with platform: 100% of PD samples
# are shotgun, 100% of AD samples are 16S. This test cannot fully separate
# a genuine PD-vs-AD disease effect from a technology/protocol effect.
# Proceeding as instructed, with this limitation stated explicitly in
# every output.
#
# Method: relative abundance (not raw counts -- not comparable across
# platforms) for each shared genus; Wilcoxon rank-sum test (PD vs AD,
# two-sided); BH correction across the shared-genera set.
#
# Output:
#   .../differential_abundance/pd_vs_ad_shared_taxa_test.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

JOINED_PATH <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/joined_genus_metadata.tsv"
SHARED_PATH <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/differential_abundance/shared_genera.csv"
OUT_PATH <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/differential_abundance/pd_vs_ad_shared_taxa_test.csv"

shared <- read_csv(SHARED_PATH, show_col_types = FALSE)
shared_taxa <- shared$taxon
cat(sprintf("Testing %d shared genera: %s\n", length(shared_taxa), paste(shared_taxa, collapse=", ")))
cat("\n*** CAVEAT: PD (shotgun) vs AD (16S) comparison is confounded with platform ***\n\n")

joined <- read_tsv(JOINED_PATH, show_col_types = FALSE)
disease_data <- joined %>% filter(disease_group %in% c("PD", "AD"))

# mat <- as.matrix(disease_data[, shared_taxa])
# row_totals <- rowSums(mat)
# rel_abund <- sweep(mat, 1, row_totals, "/") * 100  # relative abundance %, comparable across platforms

mat <- as.matrix(disease_data[, shared_taxa])

# row_totals must reflect each sample's FULL genus-level abundance, not just
# the 7 shared genera being tested here -- otherwise "relative abundance" is
# silently redefined as "% of these 7 genera" instead of "% of the sample",
# and any sample with all-zero counts across just these 7 genera causes a
# 0/0 division, which propagates as NaN -> NA through median() (this is what
# caused the median_AD = NA / non-NA p-value mismatch).
genus_cols <- setdiff(colnames(joined), c("unique_sample_id", "sample_id",
                                           "study_id", "disease_group",
                                           "platform", "disease_arm"))
# Find which "genus_cols" are actually not numeric -- these are leftover
# metadata columns that need to be added to the exclusion list above.
col_classes <- sapply(disease_data[, genus_cols], class)
non_numeric <- genus_cols[col_classes != "numeric" & col_classes != "integer"]
cat("Non-numeric columns caught in genus_cols (need to be excluded):\n")
print(non_numeric)
genus_cols <- grep("^g__", colnames(joined), value = TRUE)
                                           
full_mat <- as.matrix(disease_data[, genus_cols])
row_totals <- rowSums(full_mat, na.rm = TRUE)

n_zero_total <- sum(row_totals == 0)
if (n_zero_total > 0) {
  cat(sprintf("WARNING: %d samples have zero total genus abundance across all genera -- excluding from this test.\n",
              n_zero_total))
}
keep <- row_totals > 0
disease_data <- disease_data[keep, ]
mat <- mat[keep, , drop = FALSE]
row_totals <- row_totals[keep]

rel_abund <- sweep(mat, 1, row_totals, "/") * 100  # relative abundance %, comparable across platforms


results <- lapply(shared_taxa, function(g) {
  vals <- rel_abund[, g]
  test <- wilcox.test(vals ~ disease_data$disease_group)
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
cat("Results:\n")
print(results)
cat(sprintf("\n%d/%d shared genera significant at q<0.05 for direct PD-vs-AD comparison.\n",
            sum(results$q_value < 0.05), nrow(results)))
cat("\n*** Report the platform-confound caveat above wherever this table is used. ***\n")

