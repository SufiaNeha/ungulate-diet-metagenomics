options(stringsAsFactors = FALSE)

packages <- c(
  "phyloseq", "vegan", "ggplot2",
  "dplyr", "tidyr", "readxl",
  "readr", "tibble", "iNEXT",
  "car", "ANCOMBC", "ggrepel",
  "patchwork", "scales"
)

missing <- packages[
  !vapply(
    packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing)) {
  stop(
    "Install packages: ",
    paste(missing, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(patchwork)
})

set.seed(345)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

write_table <- function(x, name) {
  readr::write_csv(
    as.data.frame(x),
    file.path("results/tables", name)
  )
}

get_mat <- function(ps) {
  x <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) x <- t(x)
  x
}

get_meta <- function(ps) {
  x <- as(phyloseq::sample_data(ps), "data.frame")
  x$SampleID <- rownames(x)
  x
}

make_prop <- function(x) {
  sweep(x, 1, rowSums(x), "/")
}
