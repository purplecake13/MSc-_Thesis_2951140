#!/usr/bin/env Rscript
# cx_05_conqur_full.R
#
# Full-scale batch correction using ConQuR -- run in PARALLEL with
# cx_05_mmuphin_full.R purely for comparison at full scale. The 2-study
# pilot favoured MMUPHin (much cleaner study R2 drop, no new betadisper
# artifact, no zeroed samples) -- this run checks whether that holds, or
# reverses, once ConQuR has all 14 studies to work with instead of 2.
#
# Runs on the full joined table (~1839 samples, 14 studies, both platforms).
#   batch      = study_id
#   covariates = disease_group + platform   (Arm A / coarse decision, same
#                as the MMUPHin full run -- see cx_05_mmuphin_full.R header)
#   batch_ref  = largest study by sample count (chosen programmatically
#                below -- with 14 studies there's no single obvious pilot-
#                style choice, so "largest" is the least arbitrary default;
#                note in your log if you want to override this)
#
# Output:
#   .../cross_platform/conqur_corrected_table.tsv
#   .../cross_platform/full_conqur/  (before/after PCoA x2, PERMANOVA +
#                                     betadisper diagnostics, summary txt)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ConQuR)
  library(foreach)   # ConQuR dependency for %do% -- missing this broke cx_04 pilot originally
  library(vegan)
  library(ggplot2)
})

# SLURM nodes are headless -- no X11 display, so the default png/plot
# device errors out with "unable to open connection to X11 display" the
# first time ggsave() tries to render. Cairo doesn't need X11.
options(bitmapType = "cairo")

JOINED_PATH  <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/joined_genus_metadata.tsv"
OUT_DIR      <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform"
FULL_DIR     <- file.path(OUT_DIR, "full_conqur")
dir.create(FULL_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_TABLE    <- file.path(OUT_DIR, "conqur_corrected_table.tsv")
SUMMARY_PATH <- file.path(FULL_DIR, "conqur_full_summary.txt")

log_lines <- c("=== cx_05_conqur_full.R summary ===")
add_log <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

joined <- read_tsv(JOINED_PATH, show_col_types = FALSE)
add_log("Loaded joined table: %d samples.", nrow(joined))

stopifnot(!any(is.na(joined$disease_group)))
stopifnot(!any(is.na(joined$platform)))
stopifnot(!any(is.na(joined$study_id)))

n_studies <- length(unique(joined$study_id))
add_log("Studies: %d", n_studies)
add_log(paste(capture.output(print(table(joined$study_id, joined$disease_group))), collapse = "\n"))

# genus columns: select by g__ prefix, NOT a hardcoded exclusion list --
# same fix as cx_04/cx_05_mmuphin_full.
genus_cols <- grep("^g__", names(joined), value = TRUE)
add_log("Genus columns detected: %d", length(genus_cols))
if (length(genus_cols) < 50) {
  stop("Only ", length(genus_cols), " g__ columns found -- check column naming ",
       "before trusting this run (expected ~292, per the pilot-validated cx_02 output).")
}

tax_tab <- joined %>% select(all_of(genus_cols)) %>% as.data.frame()
rownames(tax_tab) <- joined$unique_sample_id
# ConQuR takes tax_tab as samples x genera directly -- NOT transposed and
# NOT renormalized (unlike the MMUPHin script; ConQuR works on the raw/
# filtered counts as-is).
stopifnot(!any(is.na(tax_tab)))

meta_df <- joined %>%
  transmute(
    study_id      = factor(study_id),
    disease_group = factor(disease_group),
    platform      = factor(platform)  # kept for PCoA/PERMANOVA plots below, NOT used as a covariate
  ) %>%
  as.data.frame()
rownames(meta_df) <- joined$unique_sample_id

batchid <- meta_df$study_id
# platform dropped as a covariate: it's structurally confounded with
# study_id (every study uses exactly one platform), so batch correction on
# study_id already absorbs any platform effect. ConQuR's regularized fit
# won't hard-error on this the way MMUPHin's adjust_batch does, but the
# platform coefficient would be unidentifiable/arbitrary rather than
# meaningful -- same underlying problem, just silent instead of caught.
covar <- meta_df %>% select(disease_group)
stopifnot(!any(is.na(covar)))

batch_sizes <- table(batchid)
batch_ref <- names(batch_sizes)[which.max(batch_sizes)]
add_log("batch_ref chosen as largest study: %s (n=%d)", batch_ref, max(batch_sizes))
add_log(paste(capture.output(print(batch_sizes)), collapse = "\n"))

add_log("Running ConQuR: batch=study_id, covariates=disease_group, batch_ref=%s ...", batch_ref)
corrected_mat <- ConQuR(
  tax_tab    = tax_tab,
  batchid    = batchid,
  covariates = covar,
  batch_ref  = batch_ref
)

saveRDS(corrected_mat, file.path(FULL_DIR, "full_corrected.rds"))

# ---- empty-row check --------------------------------------------------------
empty_before <- sum(rowSums(tax_tab) == 0)
empty_after  <- sum(rowSums(corrected_mat) == 0)
add_log("Empty rows in tax_tab (pre-correction): %d", empty_before)
add_log("Empty rows in corrected table (post-ConQuR): %d", empty_after)
if (empty_after > empty_before) {
  zeroed_ids <- rownames(corrected_mat)[rowSums(corrected_mat) == 0 & rowSums(tax_tab) > 0]
  add_log("*** WARNING: %d previously-nonzero sample(s) zeroed out post-correction: %s",
          length(zeroed_ids), paste(head(zeroed_ids, 10), collapse = ", "))
}

# ---- PCoA before/after, by study and by disease_group -----------------------
plot_pcoa <- function(mat, group, title) {
  bray <- vegdist(mat, method = "bray")
  pcoa <- cmdscale(bray, k = 2, eig = TRUE)
  df <- data.frame(PCo1 = pcoa$points[, 1], PCo2 = pcoa$points[, 2], group = group)
  var_explained <- round(100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]), 1)
  ggplot(df, aes(PCo1, PCo2, color = group)) +
    geom_point(size = 1.6, alpha = 0.7) +
    labs(title = title,
         x = sprintf("PCo1 (%.1f%%)", var_explained[1]),
         y = sprintf("PCo2 (%.1f%%)", var_explained[2])) +
    theme_minimal()
}

ggsave(file.path(FULL_DIR, "full_pcoa_before_study.png"),
       plot_pcoa(tax_tab, meta_df$study_id, "Full ConQuR run: before correction (by study)"),
       width = 6.5, height = 5.5)
ggsave(file.path(FULL_DIR, "full_pcoa_after_study.png"),
       plot_pcoa(corrected_mat, meta_df$study_id, "Full ConQuR run: after ConQuR (by study)"),
       width = 6.5, height = 5.5)
ggsave(file.path(FULL_DIR, "full_pcoa_before_disease.png"),
       plot_pcoa(tax_tab, meta_df$disease_group, "Full ConQuR run: before correction (by disease_group)"),
       width = 6.5, height = 5.5)
ggsave(file.path(FULL_DIR, "full_pcoa_after_disease.png"),
       plot_pcoa(corrected_mat, meta_df$disease_group, "Full ConQuR run: after ConQuR (by disease_group)"),
       width = 6.5, height = 5.5)

# ---- PERMANOVA ----------------------------------------------------------------
perm_study_before   <- adonis2(vegdist(tax_tab, method = "bray") ~ study_id, data = meta_df)
perm_study_after    <- adonis2(vegdist(corrected_mat, method = "bray") ~ study_id, data = meta_df)
perm_disease_before <- adonis2(vegdist(tax_tab, method = "bray") ~ disease_group, data = meta_df)
perm_disease_after  <- adonis2(vegdist(corrected_mat, method = "bray") ~ disease_group, data = meta_df)

add_log("\n=== PERMANOVA: study_id R2 before=%.3f after=%.3f ===",
        perm_study_before$R2[1], perm_study_after$R2[1])
add_log("=== PERMANOVA: disease_group R2 before=%.3f after=%.3f ===",
        perm_disease_before$R2[1], perm_disease_after$R2[1])

# ---- betadisper (study_id dispersion) -----------------------------------------
bd_before <- betadisper(vegdist(tax_tab, method = "bray"), meta_df$study_id)
bd_after  <- betadisper(vegdist(corrected_mat, method = "bray"), meta_df$study_id)
bd_before_test <- anova(bd_before)
bd_after_test  <- anova(bd_after)
add_log("\nbetadisper (study_id) p-value: before=%.3f after=%.3f -- want NOT newly significant.",
        bd_before_test$`Pr(>F)`[1], bd_after_test$`Pr(>F)`[1])

# ---- write outputs -------------------------------------------------------------
corrected_df <- as.data.frame(corrected_mat)
corrected_df$unique_sample_id <- rownames(corrected_df)
corrected_df <- corrected_df %>%
  left_join(joined %>% select(unique_sample_id, study_id, disease_group, platform),
            by = "unique_sample_id")

write_tsv(corrected_df, OUT_TABLE)
writeLines(log_lines, SUMMARY_PATH)

add_log("\nWrote ConQuR-corrected table to: %s", OUT_TABLE)
add_log("Wrote PCoA plots + summary to: %s", FULL_DIR)
add_log("\nCHECKPOINT: compare mmuphin_full_summary.txt vs conqur_full_summary.txt --")
add_log("same four axes as the pilot: study R2 drop, disease R2 held,")
add_log("no new betadisper artifact, no zeroed samples.")
