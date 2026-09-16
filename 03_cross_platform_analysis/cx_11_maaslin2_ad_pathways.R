#!/usr/bin/env Rscript
# cx_11_maaslin2_ad_pathways.R  
suppressPackageStartupMessages({
  library(Maaslin2)
  library(dplyr)
  library(readr)
})

options(bitmapType = "cairo")

PATHWAY_PATH <- "/rds/projects/e/elhamsak-ad-thesis/ad_picrust2_merged/ad_pathabundance_cpm.tsv"
METADATA_PATH <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/joined_genus_metadata.tsv"
OUT_DIR <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/maaslin2_ad_pathways"

pathways <- read_tsv(PATHWAY_PATH, show_col_types = FALSE)
pathway_ids <- pathways[[1]]
mat <- as.data.frame(t(pathways[, -1]))
colnames(mat) <- pathway_ids
rownames(mat) <- colnames(pathways)[-1]

meta <- read_tsv(METADATA_PATH, show_col_types = FALSE) %>%
  filter(disease_group %in% c("AD", "HC"), platform == "16S") %>%
  # unique_sample_id is "study_id_SRR..." / "study_id_ERR..." -- strip the
  # study-prefix so it matches the bare run accession columns in the
  # PICRUSt2-merged pathway table (confirmed via direct inspection: e.g.
  # "cirstea2022_SRR16303365" -> "SRR16303365")
  mutate(bare_run_id = sub("^[^_]+_", "", unique_sample_id))

meta <- as.data.frame(meta)
rownames(meta) <- meta$bare_run_id

common <- intersect(rownames(mat), rownames(meta))
cat(sprintf("Pathway table samples: %d | Metadata samples: %d | Overlap: %d\n",
            nrow(mat), nrow(meta), length(common)))

if (length(common) < 0.8 * min(nrow(mat), nrow(meta))) {
  stop(sprintf("Overlap (%d) is far below expected -- check bare_run_id construction before proceeding.", length(common)))
}

mat <- mat[common, , drop = FALSE]
meta <- meta[common, , drop = FALSE]

cat(sprintf("Rows going into Maaslin2: mat=%d, meta=%d\n", nrow(mat), nrow(meta)))

# Delete any stale output directory before rerunning, so Maaslin2 can't
# possibly reuse old files:
unlink(OUT_DIR, recursive = TRUE)

fit <- Maaslin2(
  input_data = mat,
  input_metadata = meta,
  output = OUT_DIR,
  fixed_effects = c("disease_group", "age", "sex", "study_id"),
  reference = c("disease_group,HC", "study_id,binyinli2019"),
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  standardize = FALSE
)

cat("MaAsLin2 (AD pathways, DRAFT for supervisor preview) complete. Results in:", OUT_DIR, "\n")
