#!/usr/bin/env Rscript
# cx_08_ancombc2.R
#
# WHY TWO SEPARATE MODELS, BOTH ON RAW DATA + study_id COVARIATE
# ------------------------------------------------------------------------
# UPDATED (this version): both PD_vs_HC and AD_vs_HC now run on RAW,
# uncorrected genus tables with study_id as a covariate. This replaces
# the earlier version, which sourced PD_vs_HC from the ConQuR-corrected
# table.
#
# Why the change: a direct comparison (cx_08b) showed the ConQuR-based
# PD_vs_HC model found only 1/74 significant genera, vs 50/233 on the
# raw table with the same fix_formula. Literature-expected PD signals
# (down Prevotella, down Faecalibacterium, up Akkermansia) were present
# and correctly directioned in the raw model but absent from the ConQuR
# model. Root cause: PD_vs_HC has no platform confound (100% shotgun
# both sides), so ConQuR correction was not structurally necessary here.
# Stacking ConQuR correction with a study_id covariate in the same
# ANCOM-BC2 model double-removed study-associated variance -- and
# disease status covaries with study (Wallen2022 alone is ~48% of the
# PD_vs_HC subset) -- taking real disease signal with it.
#
# This mirrors the logic already applied to AD_vs_HC (cx_07 finding):
# batch-correct only where a genuine platform confound exists within
# the subset being tested; otherwise use raw data with study_id as a
# covariate. Neither disease-arm DA model uses ConQuR any longer.
# ConQuR's remaining role in the thesis is the pooled multiclass
# PD/AD/HC feature matrix for Random Forest/SHAP, where platform IS
# confounded with disease across the full three-group comparison.
#
#   PD_vs_HC : RAW joined table, restricted to platform == "shotgun"
#   AD_vs_HC : RAW joined table, restricted to platform == "16S"
#
# fix_formula = disease_group + study_id (age/sex excluded: 6/9 AD
# datasets have zero demographic coverage; see decisions log).
#
# SIGNIFICANCE FILTER FIX (this version): the previous version filtered
# on diff_disease_group* alone. ANCOM-BC2's own guidance is that a
# robust call requires diff_* TRUE AND passed_ss_* TRUE (the pseudo-
# count sensitivity check) -- diff_* alone reproduced raw q<0.05 exactly
# and overcounted hits (e.g. AD_vs_HC dropped from 80 to 19 genera once
# passed_ss was correctly required). Fixed here.
#
# GGB/unresolved-genus note: many raw-model PD hits are MetaPhlAn4 GGB-
# coded genome bins without a resolved genus name. Run
# sgb_to_gtdb_profile.py on the source MetaPhlAn4 profiles (see
# decisions log) to resolve these before finalising the PD-specific
# gene list for the thesis -- this script does not do that remapping.
#
# Output:
#   .../cross_platform/differential_abundance/
#       ancombc2_PD_vs_HC_sig.csv
#       ancombc2_AD_vs_HC_sig.csv
#       pd_specific_genera.csv   (sig in PD, not in AD)
#       ad_specific_genera.csv   (sig in AD, not in PD)
#       shared_genera.csv        (sig in both, same direction)
#       ancombc2_summary.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(phyloseq)
  library(ANCOMBC)
})

JOINED_PATH   <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/joined_genus_metadata.tsv"
OUT_DIR       <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/differential_abundance"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
SUMMARY_PATH  <- file.path(OUT_DIR, "ancombc2_summary.txt")

log_lines <- c("=== cx_08_ancombc2.R summary ===")
add_log <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

# Neither raw table should contain negatives, but keep this defensively --
# cheap to check, and protects against any unexpected upstream artifact.
floor_negatives <- function(mat, label) {
  n_neg <- sum(mat < 0)
  if (n_neg > 0) {
    add_log("[%s] %d negative cells found (out of %d) -- flooring to 0 before ANCOM-BC2.",
            label, n_neg, length(mat))
    mat[mat < 0] <- 0
  }
  mat
}

run_ancombc2 <- function(tax_mat, meta, label, out_dir) {
  # phyloseq wants taxa x samples
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
  n_diff_only <- sum(res_cols$diff_sig, na.rm = TRUE)
  n_robust <- sum(res_cols$robust_sig, na.rm = TRUE)
  add_log("[%s] %d genera at q<0.05 (raw); %d pass diff_ alone; %d pass diff_ AND passed_ss (robust, of %d tested)",
          label, n_q_only, n_diff_only, n_robust, nrow(res))
  if (n_diff_only > n_robust) {
    add_log("[%s] *** %d genera passed diff_ but FAILED the sensitivity check -- excluded from 'sig'.",
            label, n_diff_only - n_robust)
  }

  # 'sig' uses the robust diff_ AND passed_ss call, not diff_ alone.
  sig <- res_cols %>%
    filter(robust_sig) %>%
    arrange(q)

  write_csv(sig, file.path(out_dir, sprintf("ancombc2_%s_sig.csv", label)))
  sig
}

joined <- read_tsv(JOINED_PATH, show_col_types = FALSE)
genus_cols <- grep("^g__", names(joined), value = TRUE)

# # ============================================================================
# # PD_vs_HC on RAW, shotgun-only table (no ConQuR correction -- see header)
# # ============================================================================
# add_log("\n########## PD_vs_HC (RAW, shotgun-only, uncorrected) ##########")

# pd_meta <- joined %>%
#   filter(disease_group %in% c("PD", "HC"), platform == "shotgun") %>%
#   transmute(unique_sample_id, study_id = factor(study_id),
#             disease_group = factor(disease_group, levels = c("HC", "PD")))
# pd_tax <- joined %>%
#   filter(disease_group %in% c("PD", "HC"), platform == "shotgun") %>%
#   select(all_of(genus_cols)) %>% as.data.frame()
# rownames(pd_tax) <- pd_meta$unique_sample_id
# pd_tax <- floor_negatives(pd_tax, "PD_vs_HC")

# add_log("[PD_vs_HC] n=%d (PD=%d, HC=%d)",
#         nrow(pd_meta), sum(pd_meta$disease_group == "PD"), sum(pd_meta$disease_group == "HC"))

# pd_sig <- run_ancombc2(pd_tax, pd_meta, "PD_vs_HC", OUT_DIR)

# ============================================================================
# PD_vs_HC on RAW, shotgun-only, FILTERED COUNT table (see pd_06 script --
# built separately from cx_01/cx_02/cx_03 to avoid touching the shared
# upstream pipeline that diversity/ConQuR already depend on)
# ============================================================================
add_log("\n########## PD_vs_HC (RAW COUNTS, shotgun-only, filtered) ##########")

PD_COUNTS_PATH <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/pd_genus_counts_filtered_joined.tsv"
pd_counts_raw <- read_tsv(PD_COUNTS_PATH, show_col_types = FALSE)
pd_genus_cols <- grep("^g__", names(pd_counts_raw), value = TRUE)

pd_meta <- pd_counts_raw %>%
  filter(disease_group %in% c("PD", "HC")) %>%
  transmute(unique_sample_id, study_id = factor(study_id),
            disease_group = factor(disease_group, levels = c("HC", "PD")))
pd_tax <- pd_counts_raw %>%
  filter(disease_group %in% c("PD", "HC")) %>%
  select(all_of(pd_genus_cols)) %>% as.data.frame()
rownames(pd_tax) <- pd_meta$unique_sample_id
pd_tax <- floor_negatives(pd_tax, "PD_vs_HC")

add_log("[PD_vs_HC] n=%d (PD=%d, HC=%d)",
        nrow(pd_meta), sum(pd_meta$disease_group == "PD"), sum(pd_meta$disease_group == "HC"))

pd_sig <- run_ancombc2(pd_tax, pd_meta, "PD_vs_HC", OUT_DIR)

# ============================================================================
# AD_vs_HC on RAW, platform-matched (16S-only) table
# ============================================================================
add_log("\n########## AD_vs_HC (raw, platform-matched 16S-only) ##########")

ad_meta <- joined %>%
  filter(disease_group %in% c("AD", "HC"), platform == "16S") %>%
  transmute(unique_sample_id, study_id = factor(study_id),
            disease_group = factor(disease_group, levels = c("HC", "AD")))
ad_tax <- joined %>%
  filter(disease_group %in% c("AD", "HC"), platform == "16S") %>%
  select(all_of(genus_cols)) %>% as.data.frame()
rownames(ad_tax) <- ad_meta$unique_sample_id
ad_tax <- floor_negatives(ad_tax, "AD_vs_HC")

add_log("[AD_vs_HC] n=%d (AD=%d, HC=%d) -- platform-matched, no shotgun-HC included",
        nrow(ad_meta), sum(ad_meta$disease_group == "AD"), sum(ad_meta$disease_group == "HC"))

ad_sig <- run_ancombc2(ad_tax, ad_meta, "AD_vs_HC", OUT_DIR)

# ============================================================================
# PD-specific / AD-specific / shared signatures -- built by INTERSECTION,
# never from a pooled 3-group model (see cx_08 header history).
# ============================================================================
add_log("\n########## Signature intersection ##########")

pd_genera <- pd_sig$taxon
ad_genera <- ad_sig$taxon

pd_specific <- pd_sig %>% filter(!taxon %in% ad_genera)
ad_specific <- ad_sig %>% filter(!taxon %in% pd_genera)

shared_taxa <- intersect(pd_genera, ad_genera)
shared <- pd_sig %>%
  filter(taxon %in% shared_taxa) %>%
  rename(lfc_pd = lfc, q_pd = q) %>%
  inner_join(ad_sig %>% rename(lfc_ad = lfc, q_ad = q), by = "taxon") %>%
  mutate(same_direction = sign(lfc_pd) == sign(lfc_ad))

n_pd_ggb <- sum(grepl("^g__GGB", pd_specific$taxon))
n_ad_ggb <- sum(grepl("^g__GGB", ad_specific$taxon))

add_log("PD-specific genera (sig in PD, not in AD): %d (%d unresolved GGB-coded)", nrow(pd_specific), n_pd_ggb)
add_log("AD-specific genera (sig in AD, not in PD): %d (%d unresolved GGB-coded)", nrow(ad_specific), n_ad_ggb)
add_log("Genera significant in both: %d (%d same direction, %d opposite direction)",
        nrow(shared), sum(shared$same_direction), sum(!shared$same_direction))
add_log("NOTE: 'shared' here means independently significant in two separate")
add_log("2-group models run on raw, platform-restricted subsets (shotgun-only for")
add_log("PD, 16S-only for AD), each with study_id as a covariate -- report this")
add_log("construction explicitly in Methods; it is not the same as a joint 3-group")
add_log("model's shared-effect term.")
if (n_pd_ggb > 0) {
  add_log("ACTION: %d PD-specific hits are unresolved GGB-coded genome bins.", n_pd_ggb)
  add_log("Run sgb_to_gtdb_profile.py on source MetaPhlAn4 profiles to resolve")
  add_log("genus names before finalising the PD-specific gene list for the thesis.")
}

write_csv(pd_specific, file.path(OUT_DIR, "pd_specific_genera.csv"))
write_csv(ad_specific, file.path(OUT_DIR, "ad_specific_genera.csv"))
write_csv(shared,      file.path(OUT_DIR, "shared_genera.csv"))

writeLines(log_lines, SUMMARY_PATH)
add_log("\nWrote all outputs to: %s", OUT_DIR)
add_log("\nCHECKPOINT: sanity-check pd_specific/ad_specific against published direction")
add_log("(PD: down Prevotella/Faecalibacterium, up Akkermansia; AD: down butyrate")
add_log("producers) before treating results as final.")
