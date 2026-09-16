#!/usr/bin/env Rscript
# install_deps.R
#
# One-off setup: installs everything cx_04/cx_05 (ConQuR + MMUPHin) need
# into your personal R library. Run ONCE on a LOGIN node (needs internet;
# compute nodes on BlueBEAR typically don't have it) -- do not put this
# inside an sbatch job.
#
# Usage:
#   module load bear-apps/2025a          # confirm exact name via `module avail`
#   module load R/4.5.0-gfbf-2025a       # confirm exact name via `module avail R/4.5`
#   Rscript install_deps.R
#
# NOTE: the CRAN package list below is my best-known approximation of
# ConQuR's dependency tree, not a guaranteed-complete list pulled from its
# current DESCRIPTION file. The remotes::install_github(dependencies=TRUE)
# call at the end will pull in anything this list misses automatically --
# treat the manual list as a head start, not the authority.

# ---- personal library setup -------------------------------------------------
# Ensures installs land in a writable location and match where R will look
# for packages automatically in future sessions/jobs (no need to set
# R_LIBS_USER by hand if this lines up with R's own default path).
user_lib <- Sys.getenv("R_LIBS_USER")
if (identical(user_lib, "")) {
  user_lib <- file.path(Sys.getenv("HOME"),
                         "R", paste0(R.version$platform, "-library"),
                         paste(R.version$major,
                               sub("\\..*$", "", R.version$minor), sep = "."))
}
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
cat("R version:", R.version.string, "\n")
cat("Installing into:", user_lib, "\n\n")

options(repos = c(CRAN = "https://cran.ma.imperial.ac.uk"))

# ---- CRAN dependencies -------------------------------------------------------
cran_pkgs <- c(
  "dplyr", "readr", "ggplot2", "vegan", "tibble", "tidyr", "purrr",
  "forcats", "magrittr", "ade4", "doParallel", "foreach", "GUniFrac",
  "gam", "glmnet", "MASS", "matrixStats", "cluster", "gridExtra"
)
missing_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) {
  cat("Installing missing CRAN packages:", paste(missing_cran, collapse = ", "), "\n")
  install.packages(missing_cran, dependencies = TRUE)
} else {
  cat("All listed CRAN dependencies already present.\n")
}

# ---- BiocManager + MMUPHin (Bioconductor) -----------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("MMUPHin", quietly = TRUE)) {
  cat("Installing MMUPHin via BiocManager...\n")
  BiocManager::install("MMUPHin", update = FALSE, ask = FALSE)
} else {
  cat("MMUPHin already present.\n")
}

# ---- ConQuR (GitHub only -- not on CRAN or Bioconductor) --------------------
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("ConQuR", quietly = TRUE)) {
  cat("Installing ConQuR from GitHub (wdl2459/ConQuR)...\n")
  # force = TRUE: remotes caches the checked-out git SHA separately from
  # whether the library install actually succeeded. If a prior attempt
  # failed partway (e.g. due to a lock file or a missing dependency),
  # remotes will otherwise see "SHA unchanged" and skip reinstalling even
  # though the package isn't actually present -- force bypasses that.
  remotes::install_github("wdl2459/ConQuR", dependencies = TRUE,
                           upgrade = "never", force = TRUE)
} else {
  cat("ConQuR already present.\n")
}

# ---- final verification ------------------------------------------------------
cat("\n=== Verification ===\n")
check_pkgs <- c("ConQuR", "MMUPHin", "foreach", "vegan", "dplyr", "readr", "ggplot2")
results <- sapply(check_pkgs, requireNamespace, quietly = TRUE)
for (pkg in check_pkgs) {
  cat(sprintf("%-10s %s\n", pkg, if (results[[pkg]]) "OK" else "STILL MISSING"))
}
if (!all(results)) {
  cat("\nSome packages are still missing -- re-run this script; it's safe to ",
      "run repeatedly (skips anything already installed). If a specific ",
      "package keeps failing, the error above it will usually name the real ",
      "missing system/CRAN dependency to add to cran_pkgs.\n", sep = "")
}