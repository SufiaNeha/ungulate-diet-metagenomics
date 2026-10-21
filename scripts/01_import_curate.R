source("scripts/00_setup.R")

# --------------------------------------------------
# Read original OTU table
# --------------------------------------------------

long <- read.delim(
  "data/otu.tab",
  check.names = FALSE
)

stopifnot(all(
  c("sample", "Var1", "Freq") %in% names(long)
))

long <- long %>%
  mutate(
    sample = as.character(sample),
    Var1 = as.character(Var1),
    Freq = as.numeric(Freq)
  )

if (anyNA(long$Freq) || any(long$Freq < 0)) {
  stop("Invalid abundance counts")
}

wide <- long %>%
  group_by(sample, Var1) %>%
  summarise(Freq = sum(Freq), .groups = "drop") %>%
  pivot_wider(
    names_from = Var1,
    values_from = Freq,
    values_fill = 0
  )

otu <- as.matrix(wide[, -1, drop = FALSE])
rownames(otu) <- wide$sample
storage.mode(otu) <- "numeric"

# --------------------------------------------------
# Read taxonomy
# --------------------------------------------------

tax <- read.delim(
  "data/tax_out.tab",
  header = FALSE,
  check.names = FALSE
)

stopifnot(ncol(tax) == 8)

names(tax) <- c(
  "TaxID", "Species", "Genus", "Family",
  "Order", "Class", "Phylum", "Kingdom"
)

tax$TaxID <- as.character(tax$TaxID)

if (anyDuplicated(tax$TaxID)) {
  stop("Duplicate taxonomy IDs")
}

colnames(otu) <- paste0("tax.", colnames(otu))
rownames(tax) <- paste0("tax.", tax$TaxID)

if (!all(colnames(otu) %in% rownames(tax))) {
  stop("Taxonomy missing for some OTU IDs")
}

tax <- tax[colnames(otu), , drop = FALSE]

# --------------------------------------------------
# Optional abundance filter
# --------------------------------------------------

# Use MIN_CELL_COUNT=19 only if this matches
# the final manuscript analysis.

min_count <- as.numeric(
  Sys.getenv("MIN_CELL_COUNT", "0")
)

if (min_count > 0) {
  otu[otu < min_count] <- 0
}

# --------------------------------------------------
# Taxonomy curation
# --------------------------------------------------

cur <- read_csv(
  "config/taxonomy_curation.csv",
  show_col_types = FALSE,
  col_types = cols(.default = col_character())
)

stopifnot(all(
  c(
    "TaxID", "Include", "Reason",
    "VerifiedGenus", "VerifiedFamily"
  ) %in% names(cur)
))

cur$TaxID <- ifelse(
  grepl("^tax\\.", cur$TaxID),
  cur$TaxID,
  paste0("tax.", cur$TaxID)
)

if (anyDuplicated(cur$TaxID)) {
  stop("Duplicate curated taxids")
}

if (!setequal(cur$TaxID, colnames(otu))) {
  stop("Curation file must cover every OTU taxid")
}

cur <- cur[
  match(colnames(otu), cur$TaxID),
  ,
  drop = FALSE
]

if (anyNA(cur$Include) ||
    any(!cur$Include %in% c("TRUE", "FALSE"))) {
  stop("Use TRUE or FALSE for every taxid")
}

fix_genus <- !is.na(cur$VerifiedGenus) &
  nzchar(cur$VerifiedGenus)

fix_family <- !is.na(cur$VerifiedFamily) &
  nzchar(cur$VerifiedFamily)

tax$Genus[fix_genus] <- cur$VerifiedGenus[fix_genus]
tax$Family[fix_family] <- cur$VerifiedFamily[fix_family]

# Retain curated terrestrial plant taxa
keep <- cur$Include == "TRUE" &
  tax$Kingdom == "Viridiplantae" &
  tax$Phylum != "Chlorophyta" &
  !is.na(tax$Genus) &
  nzchar(tax$Genus) &
  !is.na(tax$Family) &
  nzchar(tax$Family)

if (anyNA(keep)) {
  stop("Incomplete taxonomy assignments")
}

otu <- otu[, keep, drop = FALSE]

tax <- tax[
  keep,
  c(
    "Kingdom", "Phylum", "Class",
    "Order", "Family", "Genus", "Species"
  ),
  drop = FALSE
]

if (any(rowSums(otu) == 0)) {
  stop("Some samples have zero curated plant reads")
}

# --------------------------------------------------
# Read metadata
# --------------------------------------------------

meta <- as.data.frame(
  readxl::read_excel("data/metadata.xlsx")
)

if (!"SampleID" %in% names(meta)) {
  names(meta)[1] <- "SampleID"
}

stopifnot(all(
  c("SampleID", "Host", "Season") %in% names(meta)
))

meta$SampleID <- as.character(meta$SampleID)

if (anyDuplicated(meta$SampleID)) {
  stop("Duplicate sample IDs in metadata")
}

if (!setequal(rownames(otu), meta$SampleID)) {
  stop("OTU sample IDs do not match metadata")
}

meta <- meta[
  match(rownames(otu), meta$SampleID),
  ,
  drop = FALSE
]

rownames(meta) <- meta$SampleID

# --------------------------------------------------
# Standardize host names
# --------------------------------------------------

host_map <- c(
  "Ammotragus lervia" = "Aoudad",
  "Ovis canadensis" = "Desert bighorn sheep",
  "Odocoileus hemionus" = "Mule deer"
)

idx <- match(meta$Host, names(host_map))

meta$Host[!is.na(idx)] <- unname(
  host_map[idx[!is.na(idx)]]
)

meta$Season <- tools::toTitleCase(
  tolower(trimws(as.character(meta$Season)))
)

if (!all(meta$Season %in% c("Spring", "Winter"))) {
  stop("Unexpected season label")
}

meta$Host <- factor(
  meta$Host,
  levels = c(
    "Aoudad",
    "Desert bighorn sheep",
    "Mule deer"
  )
)

if (anyNA(meta$Host)) {
  stop("Unexpected host label")
}

meta$Season <- factor(
  meta$Season,
  levels = c("Spring", "Winter")
)

meta$Group <- interaction(
  meta$Host,
  meta$Season,
  sep = " - ",
  drop = TRUE
)

# --------------------------------------------------
# Create phyloseq objects
# --------------------------------------------------

ps <- phyloseq(
  otu_table(otu, taxa_are_rows = FALSE),
  tax_table(as.matrix(tax)),
  sample_data(meta)
)

ps_genus <- tax_glom(
  ps,
  taxrank = "Genus",
  NArm = TRUE
)

ps_family <- tax_glom(
  ps,
  taxrank = "Family",
  NArm = TRUE
)

saveRDS(
  list(
    ps = ps,
    genus = ps_genus,
    family = ps_family
  ),
  "results/phyloseq_objects.rds"
)

write_table(
  data.frame(
    SampleID = sample_names(ps),
    CuratedPlantReads = sample_sums(ps)
  ),
  "curated_read_counts.csv"
)

cat(
  "Samples:", nsamples(ps),
  "\nGenera:", ntaxa(ps_genus),
  "\nFamilies:", ntaxa(ps_family),
  "\n"
)

if (nsamples(ps) != 46) {
  warning("Manuscript reports 46 samples")
}
