#!/usr/bin/env Rscript
# cx_06_subset_permanova.R
#
# The pooled 3-group (PD/AD/HC) PERMANOVA on the full corrected table is
# misleading: no study contains both AD and PD patients, so disease_group
# is partially confounded with study_id for the AD-vs-PD contrast
# specifically, and that confound alone will drag the pooled disease R2
# down regardless of whether correction worked. PD-vs-HC and AD-vs-HC are
# NOT confounded this way (HC appears across both cohort types), so this
# script isolates those two honest, identifiable contrasts and reports
# before/after R2 for each -- separately -- rather than one polluted number.
#
# Runs against whichever corrected table(s) already exist (MMUPHin and/or
# ConQuR) -- safe to run now with just MMUPHin done, and rerun once ConQuR
# finishes to add it to the same comparison log.
#
# For each subset (PD-vs-HC, AD-vs-HC) x each available method:
#   - adonis2(dist ~ disease_group)              -- headline R2, before/after
#   - adonis2(dist ~ disease_group + study_id)    -- partitions how much
#     study_id still explains WITHIN the subset after correction (want this
#     small/non-significant post-correction; a large residual study_id term
#     means correction didn't fully homogenize even the honest subsets)
#   - betadisper(disease_group)                    -- dispersion check, since
#     a significant adonis2 result can partly reflect unequal dispersion
#     rather than a true centroid (mean composition) shift
#   - PCoA plot, before and after, coloured by disease_group
#
# Output:
#   .../cross_platform/subset_permanova/subset_permanova_summary.txt
#   .../cross_platform/subset_permanova/pcoa_<subset>_<method>_<before|after>.png

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(vegan)
  library(ggplot2)
})

options(bitmapType = "cairo")  # headless-safe, same reason as the correction scripts

JOINED_PATH <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/joined_genus_metadata.tsv"

# path => display label. Add/remove rows here as more methods finish.
CORRECTED_TABLES <- list(
  MMUPHin = "/rds/projects/e/elhamsak-pd-thesis/cross_platform/batch_corrected_table.tsv",
  ConQuR  = "/rds/projects/e/elhamsak-pd-thesis/cross_platform/conqur_corrected_table.tsv"
)

OUT_DIR <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/subset_permanova"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
SUMMARY_PATH <- file.path(OUT_DIR, "subset_permanova_summary.txt")

log_lines <- c("=== cx_06_subset_permanova.R summary ===")
add_log <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

# ---- load raw (pre-correction) table ----------------------------------------
joined <- read_tsv(JOINED_PATH, show_col_types = FALSE)
genus_cols <- grep("^g__", names(joined), value = TRUE)
add_log("Raw joined table: %d samples, %d genus columns.", nrow(joined), length(genus_cols))

raw_tax <- joined %>% select(all_of(genus_cols)) %>% as.data.frame()
rownames(raw_tax) <- joined$unique_sample_id
raw_meta <- joined %>%
  transmute(unique_sample_id, study_id = factor(study_id), disease_group = factor(disease_group))

# ---- subset definitions ------------------------------------------------------
subsets <- list(
  "PD_vs_HC" = c("PD", "HC"),
  "AD_vs_HC" = c("AD", "HC")
)

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

run_subset_diagnostics <- function(tax_mat, meta, subset_name, method_label, out_dir) {
  bray <- vegdist(tax_mat, method = "bray")

  perm_disease <- adonis2(bray ~ disease_group, data = meta)
  perm_full    <- adonis2(bray ~ disease_group + study_id, data = meta)
  bd           <- betadisper(bray, meta$disease_group)
  bd_test      <- anova(bd)

  add_log("[%s | %s] n=%d (%s)", subset_name, method_label, nrow(meta),
          paste(names(table(meta$disease_group)), table(meta$disease_group),
                sep = "=", collapse = ", "))
  add_log("[%s | %s] disease_group R2 = %.3f (p=%.3f)",
          subset_name, method_label, perm_disease$R2[1], perm_disease$`Pr(>F)`[1])
  add_log("[%s | %s] residual study_id R2 (after disease_group) = %.3f (p=%.3f)",
          subset_name, method_label, perm_full$R2[2], perm_full$`Pr(>F)`[2])
  add_log("[%s | %s] betadisper(disease_group) p = %.3f -- want NOT significant",
          subset_name, method_label, bd_test$`Pr(>F)`[1])

  ggsave(file.path(out_dir, sprintf("pcoa_%s_%s.png", subset_name, method_label)),
         plot_pcoa(tax_mat, meta$disease_group,
                    sprintf("%s -- %s", gsub("_", " ", subset_name), method_label)),
         width = 6, height = 5)

  invisible(list(disease_R2 = perm_disease$R2[1], study_R2 = perm_full$R2[2]))
}

# ---- BEFORE correction: run once per subset (same raw data regardless of method) ----
add_log("\n########## BEFORE CORRECTION ##########")
for (subset_name in names(subsets)) {
  groups <- subsets[[subset_name]]
  meta_sub <- raw_meta %>% filter(disease_group %in% groups) %>%
    mutate(disease_group = droplevels(disease_group), study_id = droplevels(study_id))
  tax_sub <- raw_tax[meta_sub$unique_sample_id, , drop = FALSE]
  run_subset_diagnostics(tax_sub, meta_sub, subset_name, "raw_before", OUT_DIR)
}

# ---- AFTER correction: once per available method --------------------------------
for (method_label in names(CORRECTED_TABLES)) {
  path <- CORRECTED_TABLES[[method_label]]
  if (!file.exists(path)) {
    add_log("\n[%s] corrected table not found yet at %s -- skipping (rerun this script once it exists).",
            method_label, path)
    next
  }
  add_log("\n########## AFTER CORRECTION: %s ##########", method_label)
  corrected <- read_tsv(path, show_col_types = FALSE)
  corr_genus_cols <- intersect(genus_cols, names(corrected))
  if (length(corr_genus_cols) < length(genus_cols)) {
    add_log("*** NOTE: %d/%d genus columns matched between raw and %s corrected table -- using the intersection.",
            length(corr_genus_cols), length(genus_cols), method_label)
  }
  corr_tax <- corrected %>% select(all_of(corr_genus_cols)) %>% as.data.frame()
  rownames(corr_tax) <- corrected$unique_sample_id
  corr_meta <- corrected %>%
    transmute(unique_sample_id, study_id = factor(study_id), disease_group = factor(disease_group))

  for (subset_name in names(subsets)) {
    groups <- subsets[[subset_name]]
    meta_sub <- corr_meta %>% filter(disease_group %in% groups) %>%
      mutate(disease_group = droplevels(disease_group), study_id = droplevels(study_id))
    tax_sub <- corr_tax[meta_sub$unique_sample_id, , drop = FALSE]
    run_subset_diagnostics(tax_sub, meta_sub, subset_name, method_label, OUT_DIR)
  }
}

writeLines(log_lines, SUMMARY_PATH)
add_log("\nWrote summary + PCoA plots to: %s", OUT_DIR)
add_log("\nCHECKPOINT: for each subset, compare raw_before -> MMUPHin/ConQuR:")
add_log("  - disease_group R2 should hold (not collapse) after correction")
add_log("  - residual study_id R2 should drop toward ~0 after correction")
add_log("  - betadisper p should not newly become significant after correction")
