
source("scripts/00_setup.R")

# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

otu_file <- "data/otu.tab"
tax_file <- "data/tax_out.tab"
metadata_file <- "data/metadata.xlsx"

input_files <- c(otu_file, tax_file, metadata_file)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing input files: ",
    paste(missing_files, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 2. Import abundance data
# ------------------------------------------------------------

otu_long <- readr::read_tsv(
  otu_file,
  show_col_types = FALSE,
  col_types = readr::cols(
    .default = readr::col_character()
  )
)

required_columns <- c("sample", "Var1", "Freq")

if (!all(required_columns %in% names(otu_long))) {
  stop("otu.tab must contain sample, Var1, and Freq")
}

otu_long <- otu_long %>%
  dplyr::transmute(
    SampleID = trimws(sample),
    TaxonID = trimws(Var1),
    Count = suppressWarnings(as.numeric(Freq))
  )

if (
  anyNA(otu_long$SampleID) ||
  anyNA(otu_long$TaxonID) ||
  any(otu_long$SampleID == "") ||
  any(otu_long$TaxonID == "")
) {
  stop("Missing sample or taxon identifiers")
}

if (
  anyNA(otu_long$Count) ||
  any(!is.finite(otu_long$Count)) ||
  any(otu_long$Count < 0)
) {
  stop("Invalid abundance values")
}

if (any(abs(otu_long$Count - round(otu_long$Count)) > 1e-8)) {
  stop("Abundance values must be whole-number counts")
}

otu_long <- otu_long %>%
  dplyr::group_by(SampleID, TaxonID) %>%
  dplyr::summarise(
    Count = sum(Count),
    .groups = "drop"
  )

otu_wide <- otu_long %>%
  tidyr::pivot_wider(
    names_from = TaxonID,
    values_from = Count,
    values_fill = 0
  )

otu_df <- as.data.frame(otu_wide)

rownames(otu_df) <- otu_df$SampleID
otu_df$SampleID <- NULL

otu_mat <- as.matrix(otu_df)
storage.mode(otu_mat) <- "numeric"

if (nrow(otu_mat) == 0 || ncol(otu_mat) == 0) {
  stop("The abundance matrix is empty")
}

# ------------------------------------------------------------
# 3. Import taxonomy
# ------------------------------------------------------------

tax_raw <- read.delim(
  tax_file,
  header = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = "character"
)

if (ncol(tax_raw) != 8) {
  stop("tax_out.tab must contain eight columns")
}

colnames(tax_raw) <- c(
  "TaxonID",
  "Species",
  "Genus",
  "Family",
  "Order",
  "Class",
  "Phylum",
  "Kingdom"
)

tax_raw$TaxonID <- trimws(tax_raw$TaxonID)

if (
  anyNA(tax_raw$TaxonID) ||
  any(tax_raw$TaxonID == "")
) {
  stop("Taxonomy contains missing taxon IDs")
}

# Verify consistency of duplicate taxonomy records

rank_columns <- c(
  "Species",
  "Genus",
  "Family",
  "Order",
  "Class",
  "Phylum",
  "Kingdom"
)

tax_signatures <- apply(
  tax_raw[, rank_columns, drop = FALSE],
  1,
  function(x) paste(
    ifelse(is.na(x), "<NA>", x),
    collapse = "\r"
  )
)

tax_check <- data.frame(
  TaxonID = tax_raw$TaxonID,
  Signature = tax_signatures
)

conflicting_ids <- tax_check %>%
  dplyr::group_by(TaxonID) %>%
  dplyr::summarise(
    N = dplyr::n_distinct(Signature),
    .groups = "drop"
  ) %>%
  dplyr::filter(N > 1)

if (nrow(conflicting_ids) > 0) {
  stop(
    "Conflicting taxonomy for IDs: ",
    paste(conflicting_ids$TaxonID, collapse = ", ")
  )
}

tax_raw <- tax_raw %>%
  dplyr::distinct(TaxonID, .keep_all = TRUE)

tax_df <- as.data.frame(tax_raw)

rownames(tax_df) <- tax_df$TaxonID
tax_df$TaxonID <- NULL

missing_taxa <- setdiff(
  colnames(otu_mat),
  rownames(tax_df)
)

if (length(missing_taxa) > 0) {
  stop(
    "Taxonomy missing for: ",
    paste(missing_taxa, collapse = ", ")
  )
}

tax_df <- tax_df[
  colnames(otu_mat),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 4. Filter plant assignments
# ------------------------------------------------------------

plant_keep <- !is.na(tax_df$Kingdom) &
  trimws(tax_df$Kingdom) == "Viridiplantae"

algae <- !is.na(tax_df$Phylum) &
  grepl(
    "Chlorophyta",
    tax_df$Phylum,
    ignore.case = TRUE
  )

keep_taxa <- plant_keep & !algae

exclusion_log <- data.frame(
  TaxonID = rownames(tax_df)[!keep_taxa],
  Reason = ifelse(
    algae[!keep_taxa],
    "Chlorophyta",
    "Non-Viridiplantae"
  ),
  stringsAsFactors = FALSE
)

otu_mat <- otu_mat[
  ,
  keep_taxa,
  drop = FALSE
]

tax_df <- tax_df[
  keep_taxa,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 5. Apply verified taxonomic exclusions
# ------------------------------------------------------------

exclusion_file <- "config/taxonomy_exclusions.csv"

if (file.exists(exclusion_file)) {

  exclusions <- readr::read_csv(
    exclusion_file,
    col_types = readr::cols(
      TaxonID = readr::col_character(),
      Reason = readr::col_character()
    ),
    show_col_types = FALSE
  )

  if (!all(c("TaxonID", "Reason") %in% names(exclusions))) {
    stop("Exclusion file requires TaxonID and Reason")
  }

  exclusions$TaxonID <- trimws(exclusions$TaxonID)

  if (
    anyNA(exclusions$TaxonID) ||
    any(exclusions$TaxonID == "") ||
    anyNA(exclusions$Reason) ||
    any(trimws(exclusions$Reason) == "")
  ) {
    stop("Exclusion records require TaxonID and Reason")
  }

  exclude_ids <- intersect(
    colnames(otu_mat),
    exclusions$TaxonID
  )

  exclusion_log <- dplyr::bind_rows(
    exclusion_log,
    exclusions %>%
      dplyr::filter(TaxonID %in% exclude_ids)
  )

  otu_mat <- otu_mat[
    ,
    !colnames(otu_mat) %in% exclude_ids,
    drop = FALSE
  ]

  tax_df <- tax_df[
    colnames(otu_mat),
    ,
    drop = FALSE
  ]

} else {

  warning(
    "config/taxonomy_exclusions.csv not found. ",
    "Regional taxonomic exclusions have not been applied."
  )

}

write_table(
  exclusion_log,
  "taxonomic_exclusion_log.csv"
)

# ------------------------------------------------------------
# 6. Remove unobserved taxa
# ------------------------------------------------------------

otu_mat <- otu_mat[
  ,
  colSums(otu_mat) > 0,
  drop = FALSE
]

tax_df <- tax_df[
  colnames(otu_mat),
  ,
  drop = FALSE
]

if (ncol(otu_mat) == 0) {
  stop("No plant taxa remain after filtering")
}

zero_samples <- rownames(otu_mat)[
  rowSums(otu_mat) == 0
]

if (length(zero_samples) > 0) {
  stop(
    "Samples with zero plant reads: ",
    paste(zero_samples, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 7. Import metadata
# ------------------------------------------------------------

metadata <- readxl::read_excel(
  metadata_file
) %>%
  as.data.frame()

colnames(metadata)[1] <- "SampleID"

metadata$SampleID <- trimws(
  as.character(metadata$SampleID)
)

if (
  anyNA(metadata$SampleID) ||
  any(metadata$SampleID == "") ||
  anyDuplicated(metadata$SampleID)
) {
  stop("Invalid or duplicated metadata sample IDs")
}

rownames(metadata) <- metadata$SampleID

# ------------------------------------------------------------
# 8. Validate analysis variables
# ------------------------------------------------------------

required_metadata <- c(
  "Host",
  "Season",
  "Age",
  "Sex",
  "Body_condition"
)

missing_metadata <- setdiff(
  required_metadata,
  names(metadata)
)

if (length(missing_metadata) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 9. Standardize host species
# ------------------------------------------------------------

host_raw <- tolower(
  trimws(as.character(metadata$Host))
)

metadata$Host <- dplyr::case_when(
  host_raw %in% c(
    "ammotragus lervia",
    "aoudad"
  ) ~ "Aoudad",

  host_raw %in% c(
    "ovis canadensis",
    "ovis canadensis mexicana",
    "desert bighorn sheep"
  ) ~ "Desert bighorn sheep",

  host_raw %in% c(
    "odocoileus hemionus",
    "mule deer"
  ) ~ "Mule deer",

  TRUE ~ NA_character_
)

if (anyNA(metadata$Host)) {
  stop("Unrecognized host species in metadata")
}

metadata$Host <- factor(
  metadata$Host,
  levels = host_levels
)

# ------------------------------------------------------------
# 10. Standardize seasons
# ------------------------------------------------------------

season_raw <- tolower(
  trimws(as.character(metadata$Season))
)

metadata$Season <- dplyr::case_when(
  season_raw == "spring" ~ "Spring",
  season_raw == "winter" ~ "Winter",
  TRUE ~ NA_character_
)

if (anyNA(metadata$Season)) {
  stop("Unrecognized season in metadata")
}

metadata$Season <- factor(
  metadata$Season,
  levels = season_levels
)

# ------------------------------------------------------------
# 11. Preserve demographic variables
# ------------------------------------------------------------

metadata$Sex <- factor(
  trimws(as.character(metadata$Sex))
)

# Preserve Age and Body_condition in their original
# measurement scales.

metadata$Age <- metadata$Age

metadata$Body_condition <- metadata$Body_condition

metadata$Group <- interaction(
  metadata$Host,
  metadata$Season,
  sep = " - ",
  drop = TRUE
)

# ------------------------------------------------------------
# 12. Match metadata to abundance data
# ------------------------------------------------------------

missing_metadata_samples <- setdiff(
  rownames(otu_mat),
  rownames(metadata)
)

if (length(missing_metadata_samples) > 0) {
  stop(
    "Samples missing from metadata: ",
    paste(missing_metadata_samples, collapse = ", ")
  )
}

extra_metadata_samples <- setdiff(
  rownames(metadata),
  rownames(otu_mat)
)

if (length(extra_metadata_samples) > 0) {
  warning(
    "Metadata-only samples excluded: ",
    paste(extra_metadata_samples, collapse = ", ")
  )
}

metadata <- metadata[
  rownames(otu_mat),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(otu_mat),
    rownames(metadata)
  )
)

# ------------------------------------------------------------
# 13. Construct phyloseq objects
# ------------------------------------------------------------

ps_full <- phyloseq::phyloseq(
  phyloseq::otu_table(
    otu_mat,
    taxa_are_rows = FALSE
  ),
  phyloseq::tax_table(
    as.matrix(tax_df)
  ),
  phyloseq::sample_data(metadata)
)

# Normalize missing rank labels without inventing
# taxonomic identities.

tax_matrix <- as(
  phyloseq::tax_table(ps_full),
  "matrix"
)

for (rank in c("Genus", "Family")) {

  values <- trimws(
    as.character(tax_matrix[, rank])
  )

  values[
    is.na(values) | values == ""
  ] <- NA_character_

  tax_matrix[, rank] <- values
}

phyloseq::tax_table(ps_full) <- phyloseq::tax_table(
  tax_matrix
)

ps_genus <- phyloseq::tax_glom(
  ps_full,
  taxrank = "Genus",
  NArm = TRUE
)

ps_family <- phyloseq::tax_glom(
  ps_full,
  taxrank = "Family",
  NArm = TRUE
)

# ------------------------------------------------------------
# 14. Validate taxonomic aggregation
# ------------------------------------------------------------

if (phyloseq::ntaxa(ps_genus) == 0) {
  stop("No classified genera remain")
}

if (phyloseq::ntaxa(ps_family) == 0) {
  stop("No classified families remain")
}

# ------------------------------------------------------------
# 15. Export metadata and sample summaries
# ------------------------------------------------------------

sample_summary <- metadata %>%
  dplyr::count(
    Host,
    Season,
    name = "Samples",
    .drop = FALSE
  )

write_table(
  sample_summary,
  "sample_counts_host_season.csv"
)

depth_table <- data.frame(
  SampleID = rownames(otu_mat),
  PlantReads = rowSums(otu_mat)
)

write_table(
  depth_table,
  "plant_read_depths.csv"
)

write_table(
  metadata,
  "curated_metadata.csv"
)

# ------------------------------------------------------------
# 16. Save objects for downstream analyses
# ------------------------------------------------------------

saveRDS(
  list(
    full = ps_full,
    genus = ps_genus,
    family = ps_family
  ),
  "results/phyloseq_objects.rds"
)

saveRDS(
  list(
    raw = get_mat(ps_genus),
    family_raw = get_mat(ps_family),
    metadata = get_meta(ps_genus)
  ),
  "results/analysis_objects.rds"
)

saveRDS(
  ps_full,
  "results/objects/phyloseq_full.rds"
)

saveRDS(
  ps_genus,
  "results/objects/phyloseq_genus.rds"
)

saveRDS(
  ps_family,
  "results/objects/phyloseq_family.rds"
)

# ------------------------------------------------------------
# 17. Data validation report
# ------------------------------------------------------------

validation <- data.frame(
  Metric = c(
    "Samples",
    "Taxonomic IDs",
    "Plant genera",
    "Plant families",
    "Minimum plant reads",
    "Maximum plant reads"
  ),
  Value = c(
    phyloseq::nsamples(ps_full),
    phyloseq::ntaxa(ps_full),
    phyloseq::ntaxa(ps_genus),
    phyloseq::ntaxa(ps_family),
    min(phyloseq::sample_sums(ps_full)),
    max(phyloseq::sample_sums(ps_full))
  )
)

write_table(
  validation,
  "data_validation_summary.csv"
)

print(sample_summary)
print(validation)

# ------------------------------------------------------------
# 18. Completion
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "Data import and curation completed\n",
  "Samples: ", phyloseq::nsamples(ps_full), "\n",
  "Genera: ", phyloseq::ntaxa(ps_genus), "\n",
  "Families: ", phyloseq::ntaxa(ps_family), "\n",
  "========================================\n",
  sep = ""
)
