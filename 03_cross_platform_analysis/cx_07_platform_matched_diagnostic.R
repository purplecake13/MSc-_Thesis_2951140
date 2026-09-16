#!/usr/bin/env Rscript
# cx_07_platform_matched_diagnostic.R
#
# WHY THIS SCRIPT EXISTS
# ------------------------------------------------------------------------
# The AD_vs_HC subset in cx_06 pools ALL HC samples (from both AD 16S
# studies and PD shotgun studies) against ALL AD samples (100% 16S).
# That means disease_group and platform are almost perfectly aliased for
# the AD side: every AD sample is 16S, but only ~39% of the pooled HC
# samples are. A clean AD vs HC split in that subset is equally well
# explained by "16S vs shotgun" as by "AD vs HC" -- you cannot tell which
# from the pooled PCoA/PERMANOVA alone.
#
# This script isolates the honest comparison: AD vs HC restricted to
# platform == "16S" only (i.e. AD studies' own paired controls, dropping
# every borrowed shotgun-HC sample). If disease separation survives when
# platform is held constant, that's real signal. If disease R2 collapses
# toward the platform-confounded pooled number, the "clean split" seen in
# cx_06 was substantially a platform artifact.
#
# Runs the same four-axis diagnostic as cx_06 (disease R2, residual
# study R2, betadisper, PCoA) for THREE versions of AD_vs_HC, side by side:
#   1. raw_before        -- no batch correction, all platforms pooled
#   2. raw_before_16Sonly -- no batch correction, 16S-only (platform-matched)
#   3. <method>_pooled    -- corrected table, all platforms pooled (= cx_06's number)
#   4. <method>_16Sonly   -- corrected table, 16S-only (platform-matched)
#
# for each of MMUPHin and ConQuR (whichever corrected tables exist).
#
# Output:
#   .../cross_platform/platform_matched_diagnostic/platform_matched_summary.txt
#   .../cross_platform/platform_matched_diagnostic/pcoa_<version>.png

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(vegan)
  library(ggplot2)
})

options(bitmapType = "cairo")  # headless-safe, same reason as cx_05/cx_06

JOINED_PATH <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/joined_genus_metadata.tsv"

CORRECTED_TABLES <- list(
  MMUPHin = "/rds/projects/e/elhamsak-pd-thesis/cross_platform/batch_corrected_table.tsv",
  ConQuR  = "/rds/projects/e/elhamsak-pd-thesis/cross_platform/conqur_corrected_table.tsv"
)

OUT_DIR <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/platform_matched_diagnostic"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
SUMMARY_PATH <- file.path(OUT_DIR, "platform_matched_summary.txt")

log_lines <- c("=== cx_07_platform_matched_diagnostic.R summary ===")
add_log <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

plot_pcoa <- function(mat, group, title) {
  bray <- vegdist(mat, method = "bray")
  pcoa <- cmdscale(bray, k = 2, eig = TRUE)
  df <- data.frame(PCo1 = pcoa$points[, 1], PCo2 = pcoa$points[, 2], group = group)
  var_explained <- round(100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]), 1)
  ggplot(df, aes(PCo1, PCo2, color = group)) +
    geom_point(size = 1.8, alpha = 0.75) +
    labs(title = title,
         x = sprintf("PCo1 (%.1f%%)", var_explained[1]),
         y = sprintf("PCo2 (%.1f%%)", var_explained[2])) +
    theme_minimal()
}

run_diagnostics <- function(tax_mat, meta, version_label, out_dir) {
  bray <- vegdist(tax_mat, method = "bray")

  perm_disease <- adonis2(bray ~ disease_group, data = meta)
  perm_full    <- adonis2(bray ~ disease_group + study_id, data = meta)
  bd           <- betadisper(bray, meta$disease_group)
  bd_test      <- anova(bd)

  n_by_group <- table(meta$disease_group)
  add_log("[%s] n=%d (%s)", version_label, nrow(meta),
          paste(names(n_by_group), n_by_group, sep = "=", collapse = ", "))
  add_log("[%s] disease_group R2 = %.3f (p=%.3f)",
          version_label, perm_disease$R2[1], perm_disease$`Pr(>F)`[1])
  add_log("[%s] residual study_id R2 (after disease_group) = %.3f",
          version_label, perm_full$R2[2])
  add_log("[%s] betadisper(disease_group) p = %.3f -- want NOT significant",
          version_label, bd_test$`Pr(>F)`[1])

  ggsave(file.path(out_dir, sprintf("pcoa_%s.png", version_label)),
         plot_pcoa(tax_mat, meta$disease_group,
                    sprintf("AD vs HC -- %s", gsub("_", " ", version_label))),
         width = 6, height = 5)

  invisible(list(disease_R2 = perm_disease$R2[1], disease_p = perm_disease$`Pr(>F)`[1]))
}

subset_and_run <- function(tax_all, meta_all, platform_restrict, version_label, out_dir) {
  meta_sub <- meta_all %>% filter(disease_group %in% c("AD", "HC"))
  if (platform_restrict) {
    meta_sub <- meta_sub %>% filter(platform == "16S")
  }
  meta_sub <- meta_sub %>%
    mutate(disease_group = droplevels(factor(disease_group)),
           study_id      = droplevels(factor(study_id)))
  tax_sub <- tax_all[meta_sub$unique_sample_id, , drop = FALSE]
  run_diagnostics(tax_sub, meta_sub, version_label, out_dir)
}

# ---- raw (pre-correction) ----------------------------------------------------
joined <- read_tsv(JOINED_PATH, show_col_types = FALSE)
genus_cols <- grep("^g__", names(joined), value = TRUE)
add_log("Raw joined table: %d samples, %d genus columns.", nrow(joined), length(genus_cols))

raw_tax <- joined %>% select(all_of(genus_cols)) %>% as.data.frame()
rownames(raw_tax) <- joined$unique_sample_id
raw_meta <- joined %>%
  transmute(unique_sample_id, study_id, disease_group, platform)

add_log("\nPlatform composition of pooled HC (raw, all studies):")
add_log(paste(capture.output(print(table(
  (raw_meta %>% filter(disease_group == "HC"))$platform
))), collapse = "\n"))

add_log("\n########## RAW (no batch correction) ##########")
subset_and_run(raw_tax, raw_meta, platform_restrict = FALSE, "raw_pooled",  OUT_DIR)
subset_and_run(raw_tax, raw_meta, platform_restrict = TRUE,  "raw_16Sonly", OUT_DIR)

# ---- corrected tables ---------------------------------------------------------
for (method_label in names(CORRECTED_TABLES)) {
  path <- CORRECTED_TABLES[[method_label]]
  if (!file.exists(path)) {
    add_log("\n[%s] corrected table not found at %s -- skipping.", method_label, path)
    next
  }
  add_log("\n########## %s ##########", method_label)
  corrected <- read_tsv(path, show_col_types = FALSE)
  corr_genus_cols <- intersect(genus_cols, names(corrected))
  if (length(corr_genus_cols) < length(genus_cols)) {
    add_log("*** NOTE: %d/%d genus columns matched for %s -- using the intersection.",
            length(corr_genus_cols), length(genus_cols), method_label)
  }
  corr_tax <- corrected %>% select(all_of(corr_genus_cols)) %>% as.data.frame()
  rownames(corr_tax) <- corrected$unique_sample_id
  corr_meta <- corrected %>%
    transmute(unique_sample_id, study_id, disease_group, platform)

  subset_and_run(corr_tax, corr_meta, platform_restrict = FALSE,
                 sprintf("%s_pooled", method_label), OUT_DIR)
  subset_and_run(corr_tax, corr_meta, platform_restrict = TRUE,
                 sprintf("%s_16Sonly", method_label), OUT_DIR)
}

writeLines(log_lines, SUMMARY_PATH)
add_log("\nWrote summary + PCoA plots to: %s", OUT_DIR)
add_log("\nCHECKPOINT -- for each method (raw, MMUPHin, ConQuR), compare pooled vs 16Sonly:")
add_log("  - If disease_group R2 HOLDS (similar or only modestly lower) 16Sonly vs pooled:")
add_log("    the pooled signal is largely real disease separation. Good.")
add_log("  - If disease_group R2 COLLAPSES toward ~0 in 16Sonly:")
add_log("    the pooled 'clean split' was substantially a platform artifact, not disease")
add_log("    biology -- do not report the pooled AD_vs_HC number as a clean disease result;")
add_log("    use the 16Sonly (platform-matched) number instead, and note the pooled result")
add_log("    was confounded in Methods/Limitations.")
add_log("  - n for 16Sonly comparisons is smaller (AD=323 vs HC=~294, its own paired")
add_log("    controls only) -- expect wider PCoA scatter and check p-values, not just R2.")
