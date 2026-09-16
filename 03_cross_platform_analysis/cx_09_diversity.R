#!/usr/bin/env Rscript
# cx_09_diversity.R
#
# Alpha and beta diversity, split the same way as cx_08 (ANCOM-BC2) and
# for the same reason: PD_vs_HC has no platform confound (100% shotgun
# both sides), so the ConQuR-corrected table is valid. AD_vs_HC is
# confounded with platform in the pooled table (100% of AD is 16S; the
# ConQuR batch reference has zero AD samples -- see cx_07), so this uses
# the RAW, platform-matched 16S-only subset instead.
#
# ALPHA DIVERSITY
#   Shannon index + observed richness, computed on RELATIVE ABUNDANCE
#   (not CLR -- CLR is for beta diversity/Aitchison below; diversity
#   indices should be computed on the un-transformed compositional data).
#   Kruskal-Wallis disease_group comparison per subset (2 groups per
#   subset here, so KW reduces to a rank-sum test -- fine, kept for
#   consistency with the original 3-class plan and reusable if you later
#   add a pooled alpha-diversity summary across all three groups).
#
# BETA DIVERSITY
#   Aitchison distance = Euclidean distance on CLR-transformed data
#   (pseudocount 0.5 before CLR, matching your existing normalisation
#   convention). PCoA via cmdscale. PERMANOVA (adonis2) partitioning
#   disease_group and study_id. This PERMANOVA is a REPORTED RESULT here
#   (Results: Beta Diversity), not a correction-diagnostic like the ones
#   in cx_05/cx_06/cx_07.
#
# NOTE ON NEGATIVE VALUES: ConQuR's regression-based correction can
# produce small negative abundances (same issue as cx_08). These are
# floored to 0 before diversity indices/CLR, and the count logged.
#
# Output:
#   .../cross_platform/diversity_analysis/
#       alpha_diversity_<subset>.csv        (per-sample Shannon + richness)
#       alpha_diversity_<subset>_boxplot.png
#       alpha_diversity_<subset>_kruskal.txt
#       pcoa_aitchison_<subset>.png
#       permanova_aitchison_<subset>.txt
#       diversity_summary.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(vegan)
  library(ggplot2)
})

options(bitmapType = "cairo")  # headless-safe, same reason as cx_05/06/07

JOINED_PATH  <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/joined_genus_metadata.tsv"
CONQUR_PATH  <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/conqur_corrected_table.tsv"
OUT_DIR      <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/diversity_analysis"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
SUMMARY_PATH <- file.path(OUT_DIR, "diversity_summary.txt")

log_lines <- c("=== cx_09_diversity.R summary ===")
add_log <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

floor_negatives <- function(mat, label) {
  n_neg <- sum(mat < 0)
  if (n_neg > 0) {
    add_log("[%s] %d negative cells found (out of %d) -- flooring to 0.",
            label, n_neg, length(mat))
    mat[mat < 0] <- 0
  }
  mat
}

clr_transform <- function(mat, pseudocount = 0.5) {
  m <- mat + pseudocount
  log_m <- log(m)
  gm <- rowMeans(log_m)
  sweep(log_m, 1, gm, FUN = "-")
}

# ============================================================================
# Per-subset pipeline: alpha diversity, Aitchison beta diversity, PCoA, PERMANOVA
# ============================================================================
run_diversity <- function(tax_mat, meta, subset_label, out_dir) {
  tax_mat <- floor_negatives(tax_mat, subset_label)

  # ---- alpha diversity (on relative abundance, un-transformed) --------------
  shannon   <- vegan::diversity(tax_mat, index = "shannon")
  richness  <- vegan::specnumber(tax_mat)

  alpha_df <- tibble(
    unique_sample_id = rownames(tax_mat),
    shannon           = shannon,
    richness          = richness
  ) %>%
    left_join(meta, by = "unique_sample_id")

  write_csv(alpha_df, file.path(out_dir, sprintf("alpha_diversity_%s.csv", subset_label)))

  kw_shannon  <- kruskal.test(shannon  ~ disease_group, data = alpha_df)
  kw_richness <- kruskal.test(richness ~ disease_group, data = alpha_df)

  kw_lines <- c(
    sprintf("=== Kruskal-Wallis: %s ===", subset_label),
    sprintf("Shannon:  chi-sq=%.3f, df=%d, p=%.4g", kw_shannon$statistic, kw_shannon$parameter, kw_shannon$p.value),
    sprintf("Richness: chi-sq=%.3f, df=%d, p=%.4g", kw_richness$statistic, kw_richness$parameter, kw_richness$p.value)
  )
  writeLines(kw_lines, file.path(out_dir, sprintf("alpha_diversity_%s_kruskal.txt", subset_label)))
  add_log("[%s] Shannon KW p=%.4g | Richness KW p=%.4g",
          subset_label, kw_shannon$p.value, kw_richness$p.value)

  group_medians <- alpha_df %>%
    group_by(disease_group) %>%
    summarise(median_shannon = median(shannon), median_richness = median(richness), n = n(), .groups = "drop")
  add_log("[%s] group medians:\n%s", subset_label,
          paste(capture.output(print(group_medians)), collapse = "\n"))

  alpha_long <- alpha_df %>%
    select(disease_group, shannon, richness) %>%
    tidyr::pivot_longer(cols = c(shannon, richness), names_to = "metric", values_to = "value")

  p_alpha <- ggplot(alpha_long, aes(x = disease_group, y = value, fill = disease_group)) +
    geom_boxplot(outlier.alpha = 0.4) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(title = sprintf("Alpha diversity -- %s", gsub("_", " ", subset_label)),
         x = NULL, y = NULL) +
    theme_minimal() +
    theme(legend.position = "none")
  ggsave(file.path(out_dir, sprintf("alpha_diversity_%s_boxplot.png", subset_label)),
         p_alpha, width = 7, height = 4.5)

  # ---- beta diversity: Aitchison (Euclidean on CLR) --------------------------
  clr_mat <- clr_transform(as.matrix(tax_mat))
  aitch_dist <- dist(clr_mat, method = "euclidean")

  pcoa <- cmdscale(aitch_dist, k = 2, eig = TRUE)
  var_explained <- round(100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]), 1)
  pcoa_df <- tibble(
    unique_sample_id = rownames(clr_mat),
    PCo1 = pcoa$points[, 1],
    PCo2 = pcoa$points[, 2]
  ) %>%
    left_join(meta, by = "unique_sample_id")

  p_pcoa <- ggplot(pcoa_df, aes(PCo1, PCo2, color = disease_group)) +
    geom_point(size = 1.8, alpha = 0.75) +
    labs(title = sprintf("Aitchison PCoA -- %s", gsub("_", " ", subset_label)),
         x = sprintf("PCo1 (%.1f%%)", var_explained[1]),
         y = sprintf("PCo2 (%.1f%%)", var_explained[2])) +
    theme_minimal()
  ggsave(file.path(out_dir, sprintf("pcoa_aitchison_%s.png", subset_label)),
         p_pcoa, width = 6.5, height = 5.5)

  # ---- PERMANOVA (reported result, not a correction-diagnostic here) --------
  perm_meta <- meta %>% filter(unique_sample_id %in% rownames(clr_mat))
  perm_meta <- perm_meta[match(rownames(clr_mat), perm_meta$unique_sample_id), ]
  stopifnot(identical(perm_meta$unique_sample_id, rownames(clr_mat)))

  perm_full  <- adonis2(aitch_dist ~ disease_group + study_id, data = perm_meta)
  perm_disease_only <- adonis2(aitch_dist ~ disease_group, data = perm_meta)
  bd <- betadisper(aitch_dist, perm_meta$disease_group)
  bd_test <- anova(bd)

  perm_lines <- c(
    sprintf("=== PERMANOVA (Aitchison): %s ===", subset_label),
    sprintf("disease_group alone: R2=%.3f, p=%.4g",
            perm_disease_only$R2[1], perm_disease_only$`Pr(>F)`[1]),
    sprintf("disease_group (adjusted for study_id): R2=%.3f, p=%.4g",
            perm_full$R2[1], perm_full$`Pr(>F)`[1]),
    sprintf("study_id (residual after disease_group): R2=%.3f, p=%.4g",
            perm_full$R2[2], perm_full$`Pr(>F)`[2]),
    sprintf("betadisper(disease_group) p=%.4g -- want NOT significant (else adonis2",
            bd_test$`Pr(>F)`[1]),
    "  result may partly reflect dispersion rather than centroid shift)"
  )
  writeLines(perm_lines, file.path(out_dir, sprintf("permanova_aitchison_%s.txt", subset_label)))
  add_log("[%s] PERMANOVA disease_group R2=%.3f (p=%.4g), study_id residual R2=%.3f, betadisper p=%.4g",
          subset_label, perm_disease_only$R2[1], perm_disease_only$`Pr(>F)`[1],
          perm_full$R2[2], bd_test$`Pr(>F)`[1])

  invisible(list(alpha = alpha_df, pcoa = pcoa_df))
}

# ============================================================================
# PD_vs_HC on ConQuR-corrected table
# ============================================================================
add_log("\n########## PD_vs_HC (ConQuR-corrected) ##########")
conqur <- read_tsv(CONQUR_PATH, show_col_types = FALSE)
genus_cols <- grep("^g__", names(conqur), value = TRUE)

pd_meta <- conqur %>%
  filter(disease_group %in% c("PD", "HC")) %>%
  transmute(unique_sample_id, study_id = factor(study_id),
            disease_group = factor(disease_group, levels = c("HC", "PD")))
pd_tax <- conqur %>%
  filter(disease_group %in% c("PD", "HC")) %>%
  select(all_of(genus_cols)) %>% as.data.frame()
rownames(pd_tax) <- pd_meta$unique_sample_id

run_diversity(pd_tax, pd_meta, "PD_vs_HC", OUT_DIR)

# ============================================================================
# AD_vs_HC on RAW, platform-matched (16S-only) table
# ============================================================================
add_log("\n########## AD_vs_HC (raw, platform-matched 16S-only) ##########")
joined <- read_tsv(JOINED_PATH, show_col_types = FALSE)
stopifnot(all(genus_cols %in% names(joined)))

ad_meta <- joined %>%
  filter(disease_group %in% c("AD", "HC"), platform == "16S") %>%
  transmute(unique_sample_id, study_id = factor(study_id),
            disease_group = factor(disease_group, levels = c("HC", "AD")))
ad_tax <- joined %>%
  filter(disease_group %in% c("AD", "HC"), platform == "16S") %>%
  select(all_of(genus_cols)) %>% as.data.frame()
rownames(ad_tax) <- ad_meta$unique_sample_id

run_diversity(ad_tax, ad_meta, "AD_vs_HC", OUT_DIR)

writeLines(log_lines, SUMMARY_PATH)
add_log("\nWrote all diversity outputs to: %s", OUT_DIR)
add_log("\nCHECKPOINT: compare PD_vs_HC and AD_vs_HC Aitchison PERMANOVA disease_group R2")
add_log("against the genus-level Bray-Curtis R2 from cx_06/cx_07 -- they needn't match")
add_log("exactly (different distance metrics) but should tell a broadly consistent story.")
