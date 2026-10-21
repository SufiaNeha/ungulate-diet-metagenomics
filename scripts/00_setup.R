# ------------------------------------------------------------
# 1. General settings
# ------------------------------------------------------------

options(stringsAsFactors = FALSE)

set.seed(345)

# ------------------------------------------------------------
# 2. Packages required for manuscript analyses
# ------------------------------------------------------------

# CRAN packages
cran_packages <- c(
  "vegan",          # Diversity, Bray-Curtis, PERMANOVA, NMDS
  "iNEXT",          # Sample coverage
  "spaa",           # Pianka's dietary overlap
  "ggplot2",        # Figures
  "dplyr",          # Data manipulation
  "tidyr",          # Data reshaping
  "readr",          # CSV import/export
  "readxl",         # Metadata import
  "tibble",         # Data frames
  "car",            # Levene's test
  "patchwork",      # Combine figure panels
  "ggrepel",        # NMDS genus-vector labels
  "scales",         # Figure axis formatting
  "VennDiagram",    # Dietary genus overlap diagrams
  "gridExtra",      # Arrange Venn panels
  "RColorBrewer"    # Color palettes
)

# Bioconductor packages
bioc_packages <- c(
  "phyloseq",       # Taxonomic abundance data
  "ANCOMBC"         # ANCOM-BC2
)

# Additional package used for pairwise PERMANOVA
github_packages <- c(
  "pairwiseAdonis"
)

all_packages <- c(
  cran_packages,
  bioc_packages,
  github_packages
)

# ------------------------------------------------------------
# 3. Check required packages
# ------------------------------------------------------------

missing_packages <- all_packages[
  !vapply(
    all_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {

  message(
    "Missing packages:\n",
    paste(missing_packages, collapse = ", ")
  )

  message(
    "\nInstall CRAN packages using install.packages().",
    "\nInstall phyloseq and ANCOMBC using BiocManager.",
    "\nInstall pairwiseAdonis from its GitHub repository."
  )

}

# Stop if core packages needed by all scripts are missing

core_packages <- c(
  "phyloseq",
  "vegan",
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "readxl",
  "tibble",
  "patchwork"
)

core_missing <- core_packages[
  !vapply(
    core_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(core_missing) > 0) {
  stop(
    "Install missing core packages: ",
    paste(core_missing, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 4. Load core packages
# ------------------------------------------------------------

suppressPackageStartupMessages({

  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(tibble)
  library(patchwork)

})

# Other packages will be called using package::function()
# in the relevant analysis scripts.

# ------------------------------------------------------------
# 5. Create output directories
# ------------------------------------------------------------

output_dirs <- c(
  "results/tables",
  "results/figures",
  "results/objects"
)

for (directory in output_dirs) {

  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

}

# ------------------------------------------------------------
# 6. Study factors
# ------------------------------------------------------------

host_levels <- c(
  "Aoudad",
  "Desert bighorn sheep",
  "Mule deer"
)

season_levels <- c(
  "Spring",
  "Winter"
)

# Variables included in manuscript statistical models

model_variables <- c(
  "Host",
  "Season",
  "Age",
  "Sex",
  "Body_condition"
)

# ------------------------------------------------------------
# 7. Shared helper functions
# ------------------------------------------------------------

# Export tables

write_table <- function(x, filename) {

  readr::write_csv(
    as.data.frame(x),
    file.path("results/tables", filename)
  )

}

# Extract samples x taxa matrix from phyloseq

get_mat <- function(ps) {

  mat <- as(
    phyloseq::otu_table(ps),
    "matrix"
  )

  if (phyloseq::taxa_are_rows(ps)) {
    mat <- t(mat)
  }

  mat
}

# Extract metadata

get_meta <- function(ps) {

  metadata <- as(
    phyloseq::sample_data(ps),
    "data.frame"
  )

  metadata$SampleID <- rownames(metadata)

  metadata
}

# Convert counts to proportions

make_prop <- function(mat) {

  totals <- rowSums(mat)

  if (any(!is.finite(totals)) ||
      any(totals <= 0)) {
    stop("All samples must have positive total counts")
  }

  sweep(mat, 1, totals, "/")
}

