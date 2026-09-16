#!/usr/bin/env Rscript
# cx_08b_ancombc2_pd_raw_comparison.R
#
# WHY THIS SCRIPT EXISTS
# ------------------------------------------------------------------------
# cx_08 ran PD_vs_HC on the ConQuR-corrected table with fix_formula =
# disease_group + study_id, and found only 1/74 significant genera --
# suspiciously low given the well-replicated PD gut signature in the
# literature (down Prevotella/Faecalibacterium, up Akkermansia) and a
# sample size of n=1516.
#
# Hypothesis: PD_vs_HC has NO platform confound (100% shotgun both sides,
# unlike AD_vs_HC), so ConQuR correction isn't structurally necessary here
# the way it might be thought to be for AD. Running ANCOM-BC2 with
# fix_formula = disease_group + study_id on TOP OF an already
# batch-corrected table may double-remove study-associated variance,
# taking real disease signal with it (disease status covaries with study
# -- e.g. Wallen2022 alone is ~48% of the PD_vs_HC subset).
#
# This script re-runs PD_vs_HC on the RAW, uncorrected, shotgun-only
# subset of joined_genus_metadata.tsv, with study_id as a covariate
# (same logic already applied to AD_vs_HC in cx_08: correct for platform
# confound only where a platform confound exists; otherwise use raw data
# + study covariate). Compare hit counts and direction against cx_08's
# ConQuR-based result before deciding which is the primary PD_vs_HC
# model for the thesis.
#
# Output:
#   .../cross_platform/differential_abundance/
#       ancombc2_PD_vs_HC_RAW_sig.csv
#       ancombc2_PD_vs_HC_comparison_summary.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(phyloseq)
  library(ANCOMBC)
})

JOINED_PATH  <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/joined_genus_metadata.tsv"
CONQUR_SIG_PATH <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/differential_abundance/ancombc2_PD_vs_HC_sig.csv"
OUT_DIR      <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/differential_abundance"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
SUMMARY_PATH <- file.path(OUT_DIR, "ancombc2_PD_vs_HC_comparison_summary.txt")

log_lines <- c("=== cx_08b PD_vs_HC raw-vs-ConQuR comparison ===")
add_log <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

floor_negatives <- function(mat, label) {
  n_neg <- sum(mat < 0)
  if (n_neg > 0) {
    add_log("[%s] %d negative cells found (out of %d) -- flooring to 0 before ANCOM-BC2.",
            label, n_neg, length(mat))
    mat[mat < 0] <- 0
  }
  mat
}

# Same run_ancombc2 as cx_08, with the diff_sig + passed_ss fix already applied.
run_ancombc2 <- function(tax_mat, meta, label, out_dir) {
  otu <- otu_table(t(as.matrix(tax_mat)), taxa_are_rows = TRUE)
  sdata <- sample_data(meta %>% column_to_rownames("unique_sample_id"))
  ps <- phyloseq(otu, sdata)

  add_log("[%s] running ancombc2: n=%d, fix_formula='disease_group + study_id', group='disease_group'",
          label, nsamples(ps))

  out <- ancombc2(
    data          = ps,
    fix_formula   = "disease_group + study_id",
    p_adj_method  = "BH",
    group         = "disease_group",
    struc_zero    = TRUE,
    neg_lb        = TRUE,
    alpha         = 0.05,
    n_cl          = 1
  )

  res <- out$res
  lfc_col  <- grep("^lfc_disease_group",  names(res), value = TRUE)
  q_col    <- grep("^q_disease_group",    names(res), value = TRUE)
  diff_col <- grep("^diff_disease_group", names(res), value = TRUE)
  ss_col   <- grep("^passed_ss_disease_group", names(res), value = TRUE)
  if (length(lfc_col) != 1 || length(q_col) != 1 || length(diff_col) != 1) {
    stop("[", label, "] expected exactly one lfc/q/diff disease_group column, found: ",
         paste(lfc_col, collapse=","), " / ", paste(q_col, collapse=","),
         " / ", paste(diff_col, collapse=","))
  }

  res_cols <- res %>%
    transmute(
      taxon      = taxon,
      lfc        = .data[[lfc_col]],
      q          = .data[[q_col]],
      diff_sig   = .data[[diff_col]],
      passed_ss  = if (length(ss_col) == 1) .data[[ss_col]] else NA
    ) %>%
    mutate(robust_sig = diff_sig & !is.na(passed_ss) & passed_ss)

  n_q_only <- sum(res_cols$q < 0.05, na.rm = TRUE)
  n_robust <- sum(res_cols$robust_sig, na.rm = TRUE)
  add_log("[%s] %d genera at q<0.05 (raw); %d pass diff_ AND passed_ss (of %d tested)",
          label, n_q_only, n_robust, nrow(res))

  sig <- res_cols %>%
    filter(robust_sig) %>%
    arrange(q)

  write_csv(sig, file.path(out_dir, sprintf("ancombc2_%s_sig.csv", label)))
  sig
}

# ============================================================================
# PD_vs_HC on RAW, shotgun-only table (no ConQuR correction)
# ============================================================================
add_log("\n########## PD_vs_HC (RAW, shotgun-only, uncorrected) ##########")
joined <- read_tsv(JOINED_PATH, show_col_types = FALSE)
genus_cols <- grep("^g__", names(joined), value = TRUE)

# CONFIRM this platform label matches your data before trusting results --
# run: print(unique(joined$platform)) separately if unsure.
PD_PLATFORM_LABEL <- "shotgun"

pd_meta_raw <- joined %>%
  filter(disease_group %in% c("PD", "HC"), platform == PD_PLATFORM_LABEL) %>%
  transmute(unique_sample_id, study_id = factor(study_id),
            disease_group = factor(disease_group, levels = c("HC", "PD")))
pd_tax_raw <- joined %>%
  filter(disease_group %in% c("PD", "HC"), platform == PD_PLATFORM_LABEL) %>%
  select(all_of(genus_cols)) %>% as.data.frame()
rownames(pd_tax_raw) <- pd_meta_raw$unique_sample_id
pd_tax_raw <- floor_negatives(pd_tax_raw, "PD_vs_HC_RAW")

add_log("[PD_vs_HC_RAW] n=%d (PD=%d, HC=%d)",
        nrow(pd_meta_raw), sum(pd_meta_raw$disease_group == "PD"),
        sum(pd_meta_raw$disease_group == "HC"))

pd_sig_raw <- run_ancombc2(pd_tax_raw, pd_meta_raw, "PD_vs_HC_RAW", OUT_DIR)

# ============================================================================
# Compare against cx_08's ConQuR-based result
# ============================================================================
add_log("\n########## Comparison: RAW vs ConQuR-corrected PD_vs_HC ##########")

if (file.exists(CONQUR_SIG_PATH)) {
  pd_sig_conqur <- read_csv(CONQUR_SIG_PATH, show_col_types = FALSE)
  add_log("ConQuR-corrected model: %d significant genera", nrow(pd_sig_conqur))
  add_log("Raw (uncorrected) model: %d significant genera", nrow(pd_sig_raw))

  overlap <- intersect(pd_sig_conqur$taxon, pd_sig_raw$taxon)
  only_conqur <- setdiff(pd_sig_conqur$taxon, pd_sig_raw$taxon)
  only_raw <- setdiff(pd_sig_raw$taxon, pd_sig_conqur$taxon)

  add_log("Overlap (significant in both): %d -- %s", length(overlap), paste(overlap, collapse=", "))
  add_log("Significant only in ConQuR model: %d", length(only_conqur))
  add_log("Significant only in raw model: %d -- %s", length(only_raw), paste(only_raw, collapse=", "))

  # Literature check: PD expected down Prevotella/Faecalibacterium/Segatella,
  # up Akkermansia
  lit_genera <- c("g__Prevotella", "g__Faecalibacterium", "g__Segatella", "g__Akkermansia")
  add_log("\nLiterature check genera present in raw-model hits:")
  for (g in lit_genera) {
    hit <- pd_sig_raw %>% filter(taxon == g)
    if (nrow(hit) > 0) {
      add_log("  %s: lfc=%.3f, q=%.2e (%s)", g, hit$lfc, hit$q,
              ifelse(hit$lfc < 0, "down in PD -- matches expected direction",
                     "up in PD -- check against expected direction"))
    } else {
      add_log("  %s: not significant in raw model", g)
    }
  }
} else {
  add_log("WARNING: could not find cx_08 ConQuR sig file at %s -- comparison skipped.", CONQUR_SIG_PATH)
}

writeLines(log_lines, SUMMARY_PATH)
add_log("\nWrote comparison summary to: %s", SUMMARY_PATH)
add_log("\nDECISION POINT: if the raw model recovers substantially more of the")
add_log("expected PD literature signature than the ConQuR model, consider using")
add_log("RAW + study_id covariate as the primary PD_vs_HC model in the thesis,")
add_log("consistent with the AD_vs_HC approach (batch-correct only where a")
add_log("platform confound exists; otherwise raw + study covariate). Flag this")
add_log("decision explicitly to Sakhaee before finalising Methods language.")
