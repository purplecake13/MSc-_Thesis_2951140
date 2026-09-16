#!/usr/bin/env Rscript
# cx_10_maaslin2_pd_pathways.R
suppressPackageStartupMessages({
  library(Maaslin2)
  library(dplyr)
  library(readr)
})

PATHWAY_PATH <- "/rds/projects/e/elhamsak-pd-thesis/pd_humann3_merged/pd_pathabundance_cpm_unstratified.tsv"
METADATA_PATH <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/joined_genus_metadata.tsv"
OUT_DIR <- "/rds/projects/e/elhamsak-pd-thesis/cross_platform/maaslin2_pd_pathways"

header_line <- readLines(PATHWAY_PATH, n = 1)
header_line <- sub("^#\\s*", "", header_line)   # strip leading '# ' marker only, keep the rest
col_names_vec <- strsplit(header_line, "\t")[[1]]

pathways <- read_tsv(PATHWAY_PATH, skip = 1, col_names = col_names_vec, show_col_types = FALSE)
# pathways <- read_tsv(PATHWAY_PATH, comment = "#", col_names = TRUE, show_col_types = FALSE)
colnames(pathways)[1] <- "pathway"
# Pathway table sample columns have NO study prefix, just "<bare_id>_concat_Abundance"
colnames(pathways) <- gsub("_concat_Abundance$", "", colnames(pathways))

pathway_ids <- pathways[[1]]
mat <- as.data.frame(t(pathways[, -1]))
colnames(mat) <- pathway_ids
rownames(mat) <- colnames(pathways)[-1]

meta <- read_tsv(METADATA_PATH, show_col_types = FALSE) %>%
  filter(disease_group %in% c("PD", "HC"), platform == "shotgun") %>%
  # unique_sample_id = "<study_id>_<bare_id>" -- strip by known study_id length
  # (safe even for study names containing underscores, e.g. boktor2023_rumc)
  mutate(bare_run_id = substring(unique_sample_id, nchar(study_id) + 2))

meta <- as.data.frame(meta)
rownames(meta) <- meta$bare_run_id

common <- intersect(rownames(mat), rownames(meta))
cat(sprintf("Pathway table samples: %d | Metadata samples: %d | Overlap: %d\n",
            nrow(mat), nrow(meta), length(common)))

if (length(common) < 0.8 * min(nrow(mat), nrow(meta))) {
  stop(sprintf("Overlap (%d) far below expected -- check bare_run_id / pathway column stripping.", length(common)))
}

mat <- mat[common, , drop = FALSE]
meta <- meta[common, , drop = FALSE]

fit <- Maaslin2(
  input_data = mat,
  input_metadata = meta,
  output = OUT_DIR,
  fixed_effects = c("disease_group", "age", "sex", "study_id"),
  reference = c("disease_group,HC", "study_id,mao2021"),  # multi-level study_id needs a reference
  normalization = "NONE",
  transform = "LOG",
  min_prevalence = 0.1,
  standardize = FALSE,
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

cat("MaAsLin2 (PD pathways) complete. Results in:", OUT_DIR, "\n")
