#!/usr/bin/env Rscript
# cx_10_da_bar_charts.R
#
# Figures for differential abundance (taxa) AND differential pathway results.
# For each arm/feature-type: bar chart of top |LFC| hits, a volcano plot,
# and a cross-disease (PD vs AD) comparison heatmap for features significant
# in at least one arm.
#
# GGB placeholder genus names are resolved to their GTDB/SGB-consistent
# names via a lookup table (see GGB_LOOKUP below) before plotting.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
})

# ---------------------------------------------------------------------------
# CONFIG — adjust paths to match your actual output locations
# ---------------------------------------------------------------------------

TAXA_DIR <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/differential_abundance"

PD_PATHWAY_DIR <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/maaslin2_pd_pathways"
AD_PATHWAY_DIR <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/maaslin2_ad_pathways"

FIG_DIR <- file.path(TAXA_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# GGB -> resolved name lookup (from cx_08c_lookup_ggb.sh output)
# Extend this table any time cx_08c turns up new GGB codes.
# Keys should match however GGB appears in your taxon strings, WITHOUT
# a g__ prefix (that's stripped before lookup).
# ---------------------------------------------------------------------------

GGB_LOOKUP <- c(
  "GGB2980"   = "Copromorpha",
  "GGB33512"  = "Scatomorpha",
  "GGB36267"  = "HGM13222",
  "GGB3746"   = "Lachnospira",
  "GGB45432"  = "CAJFUR01",
  "GGB9345"   = "SFEL01",
  "GGB9350"   = "UBA1685",
  "GGB9619"   = "Dysosmobacter",
  "GGB9694"   = "Scatomorpha",
  "GGB9695"   = "Scatomorpha",
  "GGB9699"   = "Faecousia",
  "GGB9770"   = "Limiplasma",
  "GGB9775"   = "Spyradocola",
  "GGB9781"   = "Pullichristensenella"
)

resolve_ggb <- function(taxon_clean) {
  # taxon_clean has already had "g__" stripped. Replace any exact GGB code
  # match using the lookup table; leave everything else untouched.
  ifelse(taxon_clean %in% names(GGB_LOOKUP),
         GGB_LOOKUP[taxon_clean],
         taxon_clean)
}

# MaAsLin2 all_results.tsv / significant_results.tsv have one row PER
# METADATA VARIABLE per feature (disease_group, age, sex, study_id all show
# up as separate rows for the same pathway). We only ever want the disease
# effect row, or every plot silently mixes in age/sex/study_id coefficients.
DISEASE_METADATA_VAR <- "disease_group"

# ---------------------------------------------------------------------------
# Generic reader: normalises ANCOM-BC2 (taxon/lfc/q_val, CSV) and MaAsLin2
# (feature/coef/qval, TSV, multi-row-per-feature) output into a common
# {feature, lfc, qval} shape. File type is inferred from the extension.
# ---------------------------------------------------------------------------

read_sig <- function(path, feature_type = c("taxon", "pathway")) {
  feature_type <- match.arg(feature_type)

  is_tsv <- grepl("\\.tsv$", path, ignore.case = TRUE)
  df <- if (is_tsv) read_tsv(path, show_col_types = FALSE) else read_csv(path, show_col_types = FALSE)

  if (all(c("taxon", "lfc") %in% names(df))) {
    df <- df %>% rename(feature = taxon, lfc = lfc)
    qcol <- intersect(c("q_val", "qval", "q.value", "q_value", "padj", "p_adj", "q"), names(df))[1]
    if (is.na(qcol)) {
      stop(sprintf(
        "Couldn't find a q-value column in %s. Actual columns: %s\nAdd the real column name to the qcol candidate list in read_sig().",
        path, paste(names(df), collapse = ", ")
      ))
    }
    df <- df %>% rename(qval = all_of(qcol))
  } else if (all(c("feature", "coef") %in% names(df))) {
    # MaAsLin2 format: filter down to the disease effect row per feature
    if ("metadata" %in% names(df)) {
      n_before <- n_distinct(df$feature)
      df <- df %>% filter(metadata == DISEASE_METADATA_VAR)
      if (nrow(df) == 0) {
        stop(sprintf(
          "No rows with metadata == '%s' in %s. Available metadata values: %s",
          DISEASE_METADATA_VAR, path, paste(unique(df$metadata), collapse = ", ")
        ))
      }
      cat(sprintf("  [%s] kept %d disease_group rows (from %d features x all metadata vars)\n",
                  basename(path), nrow(df), n_before))
    }
    df <- df %>% rename(lfc = coef)
    qcol <- intersect(c("qval", "q_val", "q.value", "q_value", "padj", "p_adj", "q"), names(df))[1]
    if (is.na(qcol)) {
      stop(sprintf(
        "Couldn't find a q-value column in %s. Actual columns: %s\nAdd the real column name to the qcol candidate list in read_sig().",
        path, paste(names(df), collapse = ", ")
      ))
    }
    df <- df %>% rename(qval = all_of(qcol))
  } else {
    stop(sprintf("Unrecognised column format in %s: %s", path, paste(names(df), collapse = ", ")))
  }

  if (feature_type == "pathway") {
    # MaAsLin2 sanitizes "ID: description" into one string, joining the ID
    # and description with ".." and turning hyphens within the ID into ".",
    # e.g. "UDPNAGSYN.PWY..UDP.N.acetyl.D.glucosamine.biosynthesis.I".
    # Keep only the ID portion and restore the hyphen: "UDPNAGSYN-PWY".
    df <- df %>%
      mutate(
        feature_clean = str_remove(feature, "^pwy_|^PWY_"),
        feature_clean = str_split_fixed(feature_clean, "\\.\\.", 2)[, 1],
        feature_clean = str_replace_all(feature_clean, "\\.", "-")
      )
  } else {
    df <- df %>%
      mutate(feature_clean = str_remove(feature, "^g__")) %>%
      mutate(feature_clean = resolve_ggb(feature_clean))
  }

  df
}

# ---------------------------------------------------------------------------
# 1. Bar chart — top 15 by |LFC|
# ---------------------------------------------------------------------------

make_bar_chart <- function(sig_path, label, disease_label, out_path,
                            feature_type = "taxon") {
  sig <- read_sig(sig_path, feature_type)

  top_n <- sig %>%
    mutate(abs_lfc = abs(lfc)) %>%
    arrange(desc(abs_lfc)) %>%
    slice_head(n = 15) %>%
    mutate(
      direction = if_else(lfc > 0, sprintf("Up in %s", disease_label), "Down"),
      feature_clean = factor(feature_clean, levels = rev(feature_clean))
    )

  p <- ggplot(top_n, aes(x = feature_clean, y = lfc, fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = setNames(c("#1b9e77", "#d95f02"),
                                         c(sprintf("Up in %s", disease_label), "Down"))) +
    labs(title = sprintf("%s: Top significant %ss by log-fold-change", label,
                          feature_type),
         x = NULL, y = "Log-fold-change (vs HC)", fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(out_path, p, width = 8, height = 6, dpi = 300)
  cat(sprintf("Wrote %s (%d %ss shown of %d significant)\n",
              out_path, nrow(top_n), feature_type, nrow(sig)))
}

# ---------------------------------------------------------------------------
# 2. Volcano plot — all tested features, not just the significant-only subset.
#    NOTE: this needs the FULL results table (all features tested), not the
#    sig-only CSV, so effect sizes near zero / non-significant points show up.
#    If you only have the sig-only CSVs on disk, point full_path at the same
#    file for now and flag in the caption that only significant features are
#    shown (comment below marks where to change this).
# ---------------------------------------------------------------------------

make_volcano <- function(full_path, label, disease_label, out_path,
                          feature_type = "taxon", q_thresh = 0.05,
                          label_top_n = 10) {
  # label_top_n: how many of the most-significant points get a text label
  # via ggrepel. Set to 0 for a fully unlabeled plot (relies on colour/
  # legend only) if a labeled version looks too cluttered for a given panel.
  df <- read_sig(full_path, feature_type)

  df <- df %>%
    mutate(
      neg_log10_q = -log10(pmax(qval, 1e-300)),  # avoid Inf on qval == 0
      sig = qval < q_thresh,
      direction = case_when(
        sig & lfc > 0 ~ sprintf("Up in %s (q<%.2g)", disease_label, q_thresh),
        sig & lfc < 0 ~ sprintf("Down (q<%.2g)", q_thresh),
        TRUE ~ "Not significant"
      )
    )

  label_df <- df %>%
    filter(sig) %>%
    arrange(qval) %>%
    slice_head(n = label_top_n)

  cols <- setNames(
    c("#1b9e77", "#d95f02", "grey70"),
    c(sprintf("Up in %s (q<%.2g)", disease_label, q_thresh),
      sprintf("Down (q<%.2g)", q_thresh),
      "Not significant")
  )

  p <- ggplot(df, aes(x = lfc, y = neg_log10_q, color = direction)) +
    geom_point(alpha = 0.7, size = 1.8) +
    geom_hline(yintercept = -log10(q_thresh), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
    ggrepel::geom_text_repel(
      data = label_df, aes(label = feature_clean),
      size = 3, max.overlaps = 20, show.legend = FALSE
    ) +
    scale_color_manual(values = cols) +
    labs(title = sprintf("%s: Volcano plot (%ss)", label, feature_type),
         x = "Log-fold-change (vs HC)", y = expression(-log[10](q)), color = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(out_path, p, width = 8, height = 7, dpi = 300)
  cat(sprintf("Wrote %s\n", out_path))
}

# ---------------------------------------------------------------------------
# 3. Cross-disease heatmap — union of features significant in PD and/or AD,
#    one column per arm, cell = LFC (NA/grey where not significant in that
#    arm). Useful for spotting mirror patterns like the O-antigen one.
# ---------------------------------------------------------------------------

make_cross_disease_heatmap <- function(pd_path, ad_path, out_path,
                                        feature_type = "taxon", top_n = 51, 
                                        exclude = character(0)) {
  pd_sig <- read_sig(pd_path, feature_type) %>% filter(!feature_clean %in% exclude) %>% select(feature_clean, lfc) %>% rename(PD = lfc)
  ad_sig <- read_sig(ad_path, feature_type) %>% filter(!feature_clean %in% exclude) %>% select(feature_clean, lfc) %>% rename(AD = lfc)

  combined <- full_join(pd_sig, ad_sig, by = "feature_clean")

  # rank by whichever arm has the larger |LFC| for that feature, then trim
  combined <- combined %>%
    mutate(rank_val = pmax(abs(PD), abs(AD), na.rm = TRUE)) %>%
    arrange(desc(rank_val)) %>%
    slice_head(n = top_n)

  mat <- as.matrix(combined %>% select(PD, AD))
  rownames(mat) <- combined$feature_clean

  col_fun <- colorRamp2(c(-max(abs(mat), na.rm = TRUE), 0, max(abs(mat), na.rm = TRUE)),
                         c("#d95f02", "white", "#1b9e77"))

  png(out_path, width = 6, height = max(4, 0.25 * nrow(mat) + 1.5), units = "in", res = 300,
      type = "cairo")
  ht <- Heatmap(
    mat, name = "LFC vs HC",
    col = col_fun,
    na_col = "grey90",
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_gp = grid::gpar(fontsize = 9),
    column_title = sprintf("PD vs AD: %s log-fold-change comparison", feature_type),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.na(mat[i, j])) {
        grid::grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = grid::gpar(fontsize = 7))
      }
    }
  )
  draw(ht)
  dev.off()
  cat(sprintf("Wrote %s (%d %ss shown)\n", out_path, nrow(mat), feature_type))
}

# ---------------------------------------------------------------------------
# RUN — Taxa (genera)
# ---------------------------------------------------------------------------

pd_taxa_sig <- file.path(TAXA_DIR, "ancombc2_PD_vs_HC_sig.csv")
ad_taxa_sig <- file.path(TAXA_DIR, "ancombc2_AD_vs_HC_sig.csv")

make_bar_chart(pd_taxa_sig, "PD_vs_HC", "PD",
               file.path(FIG_DIR, "pd_vs_hc_top_genera_lfc.png"), "taxon")
make_bar_chart(ad_taxa_sig, "AD_vs_HC", "AD",
               file.path(FIG_DIR, "ad_vs_hc_top_genera_lfc.png"), "taxon")

make_volcano(pd_taxa_sig, "PD_vs_HC", "PD",
             file.path(FIG_DIR, "pd_vs_hc_genera_volcano.png"), "taxon")
make_volcano(ad_taxa_sig, "AD_vs_HC", "AD",
             file.path(FIG_DIR, "ad_vs_hc_genera_volcano.png"), "taxon")

make_cross_disease_heatmap(pd_taxa_sig, ad_taxa_sig, exclude = "Granulicatella",
                            file.path(FIG_DIR, "pd_ad_genera_heatmap.png"), "taxon")

# ---------------------------------------------------------------------------
# RUN — Pathways
# Bar chart + heatmap use significant_results.tsv (already q<0.25 or q<0.05
# filtered by MaAsLin2 itself, per whatever threshold you ran it with).
# Volcano uses all_results.tsv so non-significant points show up too.
# ---------------------------------------------------------------------------

pd_pathway_sig <- file.path(PD_PATHWAY_DIR, "significant_results.tsv")
ad_pathway_sig <- file.path(AD_PATHWAY_DIR, "significant_results.tsv")
pd_pathway_all <- file.path(PD_PATHWAY_DIR, "all_results.tsv")
ad_pathway_all <- file.path(AD_PATHWAY_DIR, "all_results.tsv")

if (all(file.exists(pd_pathway_sig, ad_pathway_sig, pd_pathway_all, ad_pathway_all))) {
  make_bar_chart(pd_pathway_sig, "PD_vs_HC", "PD",
                 file.path(FIG_DIR, "pd_vs_hc_top_pathways_lfc.png"), "pathway")
  make_bar_chart(ad_pathway_sig, "AD_vs_HC", "AD",
                 file.path(FIG_DIR, "ad_vs_hc_top_pathways_lfc.png"), "pathway")

  # NOTE: for AD, MaAsLin2's strongest q was 0.061 (zero FDR-significant at
  # q<0.05). If AD significant_results.tsv is empty under the threshold you
  # ran MaAsLin2 with, the AD pathway bar chart / heatmap calls below will
  # error on an empty table — expected until the q<0.05 vs q<0.25 reporting
  # decision comes back from your supervisor. Comment out the AD pathway
  # calls if you hit that and only need the PD pathway figures for now.

  make_volcano(pd_pathway_all, "PD_vs_HC", "PD",
               file.path(FIG_DIR, "pd_vs_hc_pathways_volcano.png"), "pathway")
  make_volcano(ad_pathway_all, "AD_vs_HC", "AD",
               file.path(FIG_DIR, "ad_vs_hc_pathways_volcano.png"), "pathway")

  make_cross_disease_heatmap(pd_pathway_sig, ad_pathway_sig,
                              file.path(FIG_DIR, "pd_ad_pathways_heatmap.png"), "pathway")
} else {
  cat("\n[SKIPPED] Pathway figures — one or more files not found:\n")
  cat(" ", pd_pathway_sig, "\n")
  cat(" ", ad_pathway_sig, "\n")
  cat(" ", pd_pathway_all, "\n")
  cat(" ", ad_pathway_all, "\n")
}

cat("\nDone. Figures in:", FIG_DIR, "\n")
