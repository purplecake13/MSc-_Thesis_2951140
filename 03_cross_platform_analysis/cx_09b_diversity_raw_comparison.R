#!/usr/bin/env Rscript
# cx_09b_diversity_pd_raw_comparison.R
#
# WHY THIS SCRIPT EXISTS
# ------------------------------------------------------------------------
# cx_08's original PD_vs_HC ANCOM-BC2 run (ConQuR-corrected table +
# study_id covariate) found only 1/74 significant genera, vs 50/233 on
# raw data + study_id covariate -- caused by double-removing study-
# associated variance (ConQuR already adjusts for it; adding study_id
# again in the same model removes it twice).
#
# cx_09's diversity analysis (PERMANOVA via adonis2) may have applied
# the same pattern for PD_vs_HC: ConQuR-corrected distance matrix, with
# `study` included as a PERMANOVA term. This script re-runs PD_vs_HC
# alpha + beta diversity on the RAW (uncorrected) shotgun-only table,
# with study as a covariate, and compares R^2/p-values against the
# original ConQuR-based result -- same comparison logic as cx_08b,
# applied to diversity instead of DA.
#
# Does NOT touch AD_vs_HC diversity (already run on raw platform-matched
# data, consistent with this script's approach already).
#
# Output:
#   .../cross_platform/diversity/
#       pd_vs_hc_raw_alpha_diversity.csv
#       pd_vs_hc_raw_permanova_summary.txt
#       pd_vs_hc_raw_pcoa.png

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(vegan)
  library(ggplot2)
})

JOINED_PATH <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/joined_genus_metadata.tsv"
OUT_DIR <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/diversity"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
SUMMARY_PATH <- file.path(OUT_DIR, "pd_vs_hc_raw_permanova_summary.txt")

log_lines <- c("=== cx_09b PD_vs_HC raw-vs-ConQuR diversity comparison ===")
add_log <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

# ---- load raw shotgun-only PD subset (relative abundance -- fine for
# Bray-Curtis; note this is NOT the counts-based table from pd_06, since
# diversity metrics want relative abundance, not raw counts) ----
joined <- read_tsv(JOINED_PATH, show_col_types = FALSE)
genus_cols <- grep("^g__", names(joined), value = TRUE)

pd_data <- joined %>% filter(disease_group %in% c("PD", "HC"), platform == "shotgun")
add_log("[PD_vs_HC RAW] n=%d (PD=%d, HC=%d)",
        nrow(pd_data), sum(pd_data$disease_group == "PD"), sum(pd_data$disease_group == "HC"))

pd_mat <- as.matrix(pd_data[, genus_cols])
rownames(pd_mat) <- pd_data$unique_sample_id

# ============================================================================
# Alpha diversity: Shannon, richness
# ============================================================================
add_log("\n########## Alpha diversity (raw) ##########")

shannon <- diversity(pd_mat, index = "shannon")
richness <- specnumber(pd_mat)

alpha_df <- data.frame(
  unique_sample_id = pd_data$unique_sample_id,
  disease_group = pd_data$disease_group,
  study_id = pd_data$study_id,
  shannon = shannon,
  richness = richness
)
write_csv(alpha_df, file.path(OUT_DIR, "pd_vs_hc_raw_alpha_diversity.csv"))

kw_shannon <- kruskal.test(shannon ~ disease_group, data = alpha_df)
kw_richness <- kruskal.test(richness ~ disease_group, data = alpha_df)
add_log("Shannon Kruskal-Wallis: chi-sq=%.3f, p=%.4g", kw_shannon$statistic, kw_shannon$p.value)
add_log("Richness Kruskal-Wallis: chi-sq=%.3f, p=%.4g", kw_richness$statistic, kw_richness$p.value)

# ============================================================================
# Beta diversity: Bray-Curtis PCoA + PERMANOVA + betadisper
# ============================================================================
add_log("\n########## Beta diversity (raw) ##########")

bray_dist <- vegdist(pd_mat, method = "bray")

pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2)
pcoa_df <- data.frame(
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2],
  disease_group = pd_data$disease_group,
  study_id = pd_data$study_id
)
p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = disease_group)) +
  geom_point(alpha = 0.6) +
  labs(title = "PD_vs_HC Bray-Curtis PCoA (raw, uncorrected)",
       x = "PCoA1", y = "PCoA2") +
  theme_minimal()
ggsave(file.path(OUT_DIR, "pd_vs_hc_raw_pcoa.png"), p, width = 7, height = 5, dpi = 300)

# PERMANOVA: disease_group + study_id, on RAW (uncorrected) data.
# Compare this R^2 against the original ConQuR-based cx_09 result --
# if raw shows a MUCH higher disease_group R^2, that's the same
# double-correction pattern found in cx_08's ANCOM-BC2 fix.
# permanova <- adonis2(bray_dist ~ disease_group + study_id, data = pd_data, permutations = 999)
permanova <- adonis2(bray_dist ~ disease_group + study_id, data = pd_data, permutations = 999, by = "margin")
add_log("\nPERMANOVA (raw + study_id covariate):")
permanova_str <- capture.output(print(permanova))
log_lines <<- c(log_lines, permanova_str)
cat(paste(permanova_str, collapse = "\n"), "\n")

bd <- betadisper(bray_dist, pd_data$disease_group)
bd_test <- permutest(bd, permutations = 999)
add_log("\nbetadisper (disease_group): F=%.3f, p=%.4g",
        bd_test$tab$F[1], bd_test$tab$`Pr(>F)`[1])

writeLines(log_lines, SUMMARY_PATH)
add_log("\nWrote outputs to: %s", OUT_DIR)
add_log("\nCOMPARISON STEP: manually compare disease_group R^2 and p-value above")
add_log("against cx_09's original ConQuR-based PD_vs_HC PERMANOVA result")
add_log("(R^2=0.018, betadisper p=3.5e-24 per session notes). A substantially")
add_log("higher raw R^2 would indicate the same double-correction pattern found")
add_log("in cx_08's ANCOM-BC2 fix; a similar R^2 would suggest ConQuR wasn't")
add_log("meaningfully suppressing signal for this particular metric.")
