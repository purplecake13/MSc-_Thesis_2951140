#!/usr/bin/env Rscript
# cx_08b_pd_sensitivity.R
#
# The primary PD_vs_HC ANCOM-BC2 model (cx_08, fix_formula = disease_group
# + study_id, on ConQuR-corrected data) returned only 1/74 significant
# genera despite n=1516 -- backwards from what sample size alone predicts
# vs AD_vs_HC's n=617. cx_06 shows why this is plausible, not obviously a
# bug: residual study_id R2 stays at 0.847 even AFTER ConQuR correction for
# PD_vs_HC (vs 0.670 raw) -- correction did not meaningfully reduce
# study's RELATIVE dominance over disease here. Wallen2022 alone
# contributes 490/756 (65%) of all PD samples and 234/760 (31%) of all HC
# samples -- an extreme single-study dominance for a 14-study table.
#
# This script runs 3 variants on the SAME PD_vs_HC ConQuR table to
# separate "genuine near-null disease signal" from "modeling artifact of
# double-adjusting for batch on top of incomplete correction, dominated by
# one oversized study":
#
#   A) baseline       : fix_formula = disease_group + study_id, ALL studies
#                       (= cx_08's primary result, rerun here for a clean
#                       side-by-side log)
#   B) no_study_covar : fix_formula = disease_group only, ALL studies.
#                       Tests whether including study_id on top of an
#                       already-corrected table is what's crushing power.
#                       *** EXPLORATORY ONLY *** -- relies entirely on
#                       ConQuR having removed batch effect, which cx_06
#                       shows it did NOT fully do here. A hit appearing
#                       ONLY in this variant may be residual-batch leaking
#                       through as a false "disease" effect, not a
#                       genuine PD-specific genus.
#   C) drop_wallen    : fix_formula = disease_group + study_id, EXCLUDING
#                       wallen2022. Tests whether one oversized, high-
#                       leverage study is driving the residual study_id
#                       dominance seen in the baseline.
#
# Output: .../differential_abundance/pd_sensitivity/
#   ancombc2_PD_vs_HC_<variant>_sig.csv  (one per variant)
#   pd_sensitivity_summary.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(phyloseq)
  library(ANCOMBC)
})

CONQUR_PATH <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/conqur_corrected_table.tsv"
OUT_DIR     <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/differential_abundance/pd_sensitivity"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
SUMMARY_PATH <- file.path(OUT_DIR, "pd_sensitivity_summary.txt")

log_lines <- c("=== cx_08b_pd_sensitivity.R summary ===")
add_log <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

floor_negatives <- function(mat, label) {
  n_neg <- sum(mat < 0)
  if (n_neg > 0) {
    add_log("[%s] %d negative cells found -- flooring to 0 before ANCOM-BC2.", label, n_neg)
    mat[mat < 0] <- 0
  }
  mat
}

run_ancombc2_variant <- function(tax_mat, meta, fix_formula, label, out_dir) {
  otu <- otu_table(t(as.matrix(tax_mat)), taxa_are_rows = TRUE)
  sdata <- sample_data(meta %>% column_to_rownames("unique_sample_id"))
  ps <- phyloseq(otu, sdata)

  add_log("[%s] n=%d (HC=%d, PD=%d), fix_formula='%s'",
          label, nsamples(ps), sum(meta$disease_group == "HC"),
          sum(meta$disease_group == "PD"), fix_formula)

  out <- ancombc2(
    data         = ps,
    fix_formula  = fix_formula,
    p_adj_method = "BH",
    group        = "disease_group",
    struc_zero   = TRUE,
    neg_lb       = TRUE,
    alpha        = 0.05,
    n_cl         = 1
  )

  res <- out$res
  lfc_col  <- grep("^lfc_disease_group",  names(res), value = TRUE)
  q_col    <- grep("^q_disease_group",    names(res), value = TRUE)
  diff_col <- grep("^diff_disease_group", names(res), value = TRUE)

  res_cols <- res %>%
    transmute(taxon = taxon, lfc = .data[[lfc_col]], q = .data[[q_col]],
              diff_sig = .data[[diff_col]])

  n_q_only <- sum(res_cols$q < 0.05, na.rm = TRUE)
  n_robust <- sum(res_cols$diff_sig, na.rm = TRUE)
  add_log("[%s] %d genera at q<0.05 (raw); %d pass robust diff_ call (of %d tested)",
          label, n_q_only, n_robust, nrow(res))

  sig <- res_cols %>% filter(diff_sig) %>% arrange(q)
  write_csv(sig, file.path(out_dir, sprintf("ancombc2_PD_vs_HC_%s_sig.csv", label)))
  sig
}

conqur <- read_tsv(CONQUR_PATH, show_col_types = FALSE)
genus_cols <- grep("^g__", names(conqur), value = TRUE)

pd_full <- conqur %>% filter(disease_group %in% c("PD", "HC"))
add_log("Wallen2022 share of PD_vs_HC subset: PD=%d/%d (%.0f%%), HC=%d/%d (%.0f%%)",
        sum(pd_full$study_id == "wallen2022" & pd_full$disease_group == "PD"),
        sum(pd_full$disease_group == "PD"),
        100 * sum(pd_full$study_id == "wallen2022" & pd_full$disease_group == "PD") / sum(pd_full$disease_group == "PD"),
        sum(pd_full$study_id == "wallen2022" & pd_full$disease_group == "HC"),
        sum(pd_full$disease_group == "HC"),
        100 * sum(pd_full$study_id == "wallen2022" & pd_full$disease_group == "HC") / sum(pd_full$disease_group == "HC"))

make_meta_tax <- function(df) {
  meta <- df %>%
    transmute(unique_sample_id, study_id = factor(study_id),
              disease_group = factor(disease_group, levels = c("HC", "PD")))
  tax <- df %>% select(all_of(genus_cols)) %>% as.data.frame()
  rownames(tax) <- meta$unique_sample_id
  tax <- floor_negatives(tax, "PD_vs_HC")
  list(meta = meta, tax = tax)
}

# ---- A) baseline: disease_group + study_id, all studies ---------------------
add_log("\n########## A) baseline (disease_group + study_id, all studies) ##########")
d_a <- make_meta_tax(pd_full)
sig_a <- run_ancombc2_variant(d_a$tax, d_a$meta, "disease_group + study_id", "baseline", OUT_DIR)

# ---- B) no study_id covariate, all studies -----------------------------------
add_log("\n########## B) no_study_covar (disease_group only, all studies) ##########")
add_log("*** CAUTION: results here are EXPLORATORY ONLY. Dropping study_id relies")
add_log("entirely on ConQuR having removed batch effect, which cx_06 shows it did")
add_log("NOT fully do for PD_vs_HC (residual study R2 = 0.847). A hit appearing")
add_log("ONLY in this variant (not in A or C) may be a residual-batch false")
add_log("positive, not a genuine disease effect -- do not report as PD-specific")
add_log("without corroboration from ALDEx2 or the other variants.")
sig_b <- run_ancombc2_variant(d_a$tax, d_a$meta, "disease_group", "no_study_covar", OUT_DIR)

# ---- C) drop wallen2022, disease_group + study_id -----------------------------
add_log("\n########## C) drop_wallen (disease_group + study_id, wallen2022 excluded) ##########")
pd_nowallen <- pd_full %>% filter(study_id != "wallen2022")
d_c <- make_meta_tax(pd_nowallen %>% mutate(study_id = droplevels(factor(study_id))))
sig_c <- run_ancombc2_variant(d_c$tax, d_c$meta, "disease_group + study_id", "drop_wallen", OUT_DIR)

# ---- comparison ---------------------------------------------------------------
add_log("\n########## Comparison across variants ##########")
add_log("A) baseline (all studies, +study_id):      %d significant", nrow(sig_a))
add_log("B) no_study_covar (all studies, no covar):  %d significant", nrow(sig_b))
add_log("C) drop_wallen (13 studies, +study_id):     %d significant", nrow(sig_c))

overlap_ac <- intersect(sig_a$taxon, sig_c$taxon)
add_log("Overlap between A and C (baseline vs drop_wallen): %d genera: %s",
        length(overlap_ac), paste(overlap_ac, collapse = ", "))

writeLines(log_lines, SUMMARY_PATH)
add_log("\nWrote all outputs to: %s", OUT_DIR)
add_log("\nCHECKPOINT: if C >> A, Wallen's dominance is likely suppressing power in the")
add_log("full model -- consider whether drop_wallen should become the PRIMARY")
add_log("PD_vs_HC result (with the full-data model reported as a robustness check),")
add_log("not the other way around. If B >> A with little overlap to C, treat B's")
add_log("extra hits as residual-batch artifacts, not genuine PD-specific genera.")
