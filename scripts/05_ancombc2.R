source("scripts/00_setup.R")

set.seed(345)

alpha_threshold <- 0.05
top_n_genera <- 10

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

dir.create(
  "results/objects",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Load genus-level count data
# ------------------------------------------------------------

objects <- readRDS(
  "results/phyloseq_objects.rds"
)

ps_genus <- objects$genus

if (!inherits(ps_genus, "phyloseq")) {
  stop("Genus object must be a phyloseq object")
}

metadata <- get_meta(ps_genus)

required_metadata <- c(
  "SampleID",
  "Host",
  "Season"
)

if (!all(required_metadata %in% names(metadata))) {
  stop("Required metadata columns are missing")
}

if (
  anyNA(metadata$SampleID) ||
  anyDuplicated(metadata$SampleID)
) {
  stop("Invalid or duplicate sample IDs")
}

metadata$Host <- factor(
  metadata$Host,
  levels = c(
    "Aoudad",
    "Desert bighorn sheep",
    "Mule deer"
  )
)

metadata$Season <- factor(
  metadata$Season,
  levels = c("Spring", "Winter")
)

if (
  anyNA(metadata$Host) ||
  anyNA(metadata$Season)
) {
  stop("Invalid Host or Season values")
}

metadata$Group <- paste(
  as.character(metadata$Host),
  as.character(metadata$Season),
  sep = " - "
)

counts <- get_mat(ps_genus)

if (!setequal(
  rownames(counts),
  metadata$SampleID
)) {
  stop("Count matrix and metadata do not match")
}

if (
  anyNA(counts) ||
  any(counts < 0) ||
  any(abs(counts - round(counts)) > 1e-8)
) {
  stop("ANCOM-BC2 requires nonnegative integer counts")
}

# ------------------------------------------------------------
# 2. Load pairwise PERMANOVA results
# ------------------------------------------------------------

pairwise_file <- paste0(
  "results/tables/",
  "PERMANOVA_pairwise_all.csv"
)

if (!file.exists(pairwise_file)) {
  stop("Pairwise PERMANOVA results not found")
}

pairwise_results <- readr::read_csv(
  pairwise_file,
  show_col_types = FALSE
)

required_pairwise <- c(
  "Group1",
  "Group2",
  "P_adjusted"
)

if (!all(
  required_pairwise %in% names(pairwise_results)
)) {
  stop("Pairwise PERMANOVA columns are missing")
}

significant_pairs <- pairwise_results %>%
  dplyr::filter(
    !is.na(P_adjusted),
    P_adjusted < alpha_threshold
  ) %>%
  dplyr::arrange(
    P_adjusted,
    Group1,
    Group2
  ) %>%
  dplyr::distinct(
    Group1,
    Group2,
    .keep_all = TRUE
  )

if (nrow(significant_pairs) == 0) {
  stop("No significant pairwise PERMANOVA comparisons")
}

observed_groups <- unique(metadata$Group)

unknown_groups <- setdiff(
  unique(c(
    significant_pairs$Group1,
    significant_pairs$Group2
  )),
  observed_groups
)

if (length(unknown_groups) > 0) {
  stop(
    "Unknown comparison groups: ",
    paste(unknown_groups, collapse = ", ")
  )
}

significant_pairs$ContrastID <- sprintf(
  "Contrast_%02d",
  seq_len(nrow(significant_pairs))
)

significant_pairs$Contrast <- paste(
  significant_pairs$Group1,
  "vs",
  significant_pairs$Group2
)

write_table(
  significant_pairs,
  "ANCOMBC2_significant_contrasts.csv"
)

# ------------------------------------------------------------
# 3. Prepare genus taxonomy
# ------------------------------------------------------------

taxonomy <- as(
  phyloseq::tax_table(ps_genus),
  "matrix"
)

if (!"Genus" %in% colnames(taxonomy)) {
  stop("Genus rank missing from taxonomy")
}

genus_lookup <- data.frame(
  TaxonID = rownames(taxonomy),
  Genus = trimws(
    as.character(taxonomy[, "Genus"])
  ),
  stringsAsFactors = FALSE
)

if (
  anyNA(genus_lookup$Genus) ||
  any(genus_lookup$Genus == "")
) {
  stop("Missing genus labels")
}

# ------------------------------------------------------------
# 4. Run ANCOM-BC2 for one pairwise contrast
# ------------------------------------------------------------

run_contrast <- function(
    group1,
    group2,
    contrast_id
) {

  selected_ids <- metadata$SampleID[
    metadata$Group %in% c(group1, group2)
  ]

  comparison_metadata <- metadata[
    match(selected_ids, metadata$SampleID),
    ,
    drop = FALSE
  ]

  comparison_metadata$ContrastGroup <- factor(
    comparison_metadata$Group,
    levels = c(group1, group2)
  )

  if (any(table(
    comparison_metadata$ContrastGroup
  ) < 2)) {
    warning(
      "Skipping ", contrast_id,
      ": insufficient samples"
    )
    return(NULL)
  }

  ps_subset <- phyloseq::prune_samples(
    selected_ids,
    ps_genus
  )

  ps_subset <- phyloseq::prune_taxa(
    phyloseq::taxa_sums(ps_subset) > 0,
    ps_subset
  )

  if (phyloseq::ntaxa(ps_subset) < 2) {
    warning(
      "Skipping ", contrast_id,
      ": insufficient taxa"
    )
    return(NULL)
  }

  comparison_metadata <- comparison_metadata[
    match(
      phyloseq::sample_names(ps_subset),
      comparison_metadata$SampleID
    ),
    ,
    drop = FALSE
  ]

  rownames(comparison_metadata) <-
    comparison_metadata$SampleID

  phyloseq::sample_data(ps_subset) <-
    phyloseq::sample_data(comparison_metadata)

  # The first group is the reference.
  # Positive LFC indicates higher abundance
  # in the second group.

  set.seed(345)

  fit <- ANCOMBC::ancombc2(
    data = ps_subset,
    tax_level = NULL,
    fix_formula = "ContrastGroup",
    rand_formula = NULL,
    p_adj_method = "BH",
    pseudo_sens = TRUE,
    prv_cut = 0,
    lib_cut = 0,
    s0_perc = 0.05,
    group = "ContrastGroup",
    struc_zero = FALSE,
    neg_lb = FALSE,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE,
    n_cl = 1,
    verbose = FALSE
  )

  result <- as.data.frame(fit$res)

  if (nrow(result) == 0) {
    warning("Empty ANCOM-BC2 result: ", contrast_id)
    return(NULL)
  }

  lfc_column <- grep(
    "^lfc_ContrastGroup",
    names(result),
    value = TRUE
  )

  q_column <- grep(
    "^q_ContrastGroup",
    names(result),
    value = TRUE
  )

  diff_column <- grep(
    "^diff_ContrastGroup",
    names(result),
    value = TRUE
  )

  if (
    length(lfc_column) != 1 ||
    length(q_column) != 1 ||
    length(diff_column) > 1
  ) {
    stop(
      "Unexpected ANCOM-BC2 output structure: ",
      contrast_id
    )
  }

  if (!"taxon" %in% names(result)) {
    result$taxon <- rownames(result)
  }

  result$TaxonID <- as.character(result$taxon)

  result <- result %>%
    dplyr::left_join(
      genus_lookup,
      by = "TaxonID"
    )

  if (anyNA(result$Genus)) {
    stop("Could not match genus labels")
  }

  result$LFC <- as.numeric(
    result[[lfc_column]]
  )

  result$Q <- as.numeric(
    result[[q_column]]
  )

  if (length(diff_column) == 1) {
    result$Differential <- as.logical(
      result[[diff_column]]
    )
  } else {
    result$Differential <- (
      !is.na(result$Q) &
      result$Q < alpha_threshold
    )
  }

  result$ContrastID <- contrast_id
  result$Group1 <- group1
  result$Group2 <- group2

  result$Contrast <- paste(
    group1,
    "vs",
    group2
  )

  result$EnrichedIn <- ifelse(
    result$LFC > 0,
    group2,
    group1
  )

  write_table(
    result,
    paste0(
      "ANCOMBC2_",
      contrast_id,
      "_all.csv"
    )
  )

  saveRDS(
    fit,
    file.path(
      "results/objects",
      paste0(
        "ANCOMBC2_",
        contrast_id,
        ".rds"
      )
    )
  )

  result %>%
    dplyr::select(
      ContrastID,
      Contrast,
      Group1,
      Group2,
      TaxonID,
      Genus,
      LFC,
      Q,
      Differential,
      EnrichedIn
    )
}

# ------------------------------------------------------------
# 5. Run all significant comparisons
# ------------------------------------------------------------

contrast_results <- vector(
  "list",
  nrow(significant_pairs)
)

for (i in seq_len(nrow(significant_pairs))) {

  comparison <- significant_pairs[
    i,
    ,
    drop = FALSE
  ]

  message(
    "Running ANCOM-BC2: ",
    comparison$Contrast
  )

  contrast_results[[i]] <- run_contrast(
    group1 = comparison$Group1,
    group2 = comparison$Group2,
    contrast_id = comparison$ContrastID
  )
}

all_results <- dplyr::bind_rows(
  contrast_results
)

if (nrow(all_results) == 0) {
  stop("No ANCOM-BC2 results were generated")
}

write_table(
  all_results,
  "ANCOMBC2_all_contrasts.csv"
)

# ------------------------------------------------------------
# 6. Identify significant genera
# ------------------------------------------------------------

significant_results <- all_results %>%
  dplyr::filter(
    !is.na(Q),
    !is.na(LFC),
    !is.na(Differential),
    Differential,
    Q < alpha_threshold
  ) %>%
  dplyr::arrange(
    ContrastID,
    Q,
    dplyr::desc(abs(LFC))
  )

write_table(
  significant_results,
  "ANCOMBC2_significant_genera.csv"
)

# ------------------------------------------------------------
# 7. Select genera for Figure 8
# ------------------------------------------------------------

plot_data <- significant_results %>%
  dplyr::group_by(ContrastID) %>%
  dplyr::arrange(
    Q,
    dplyr::desc(abs(LFC)),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(
    n = top_n_genera
  ) %>%
  dplyr::ungroup()

write_table(
  plot_data,
  "Figure_8_plot_data.csv"
)

# ------------------------------------------------------------
# 8. Plotting function
# ------------------------------------------------------------

make_ancombc2_panel <- function(contrast_id) {

  comparison <- significant_pairs %>%
    dplyr::filter(
      ContrastID == contrast_id
    )

  panel_data <- plot_data %>%
    dplyr::filter(
      ContrastID == contrast_id
    )

  if (nrow(panel_data) == 0) {

    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0,
          y = 0,
          label = "No significant genera",
          size = 5
        ) +
        ggplot2::xlim(-1, 1) +
        ggplot2::ylim(-1, 1) +
        ggplot2::labs(
          title = comparison$Contrast
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            face = "bold",
            hjust = 0.5
          ),
          plot.background = ggplot2::element_rect(
            fill = "white",
            color = NA
          )
        )
    )
  }

  panel_data <- panel_data %>%
    dplyr::arrange(LFC) %>%
    dplyr::mutate(
      PlotID = factor(
        TaxonID,
        levels = TaxonID
      )
    )

  genus_labels <- stats::setNames(
    panel_data$Genus,
    panel_data$TaxonID
  )

  group_colors <- stats::setNames(
    c("#0072B2", "#D55E00"),
    c(
      comparison$Group1,
      comparison$Group2
    )
  )

  ggplot2::ggplot(
    panel_data,
    ggplot2::aes(
      x = PlotID,
      y = LFC,
      fill = EnrichedIn
    )
  ) +
    ggplot2::geom_col(
      width = 0.75
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.5
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_x_discrete(
      labels = genus_labels
    ) +
    ggplot2::scale_fill_manual(
      values = group_colors,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = comparison$Contrast,
      x = NULL,
      y = "Log fold change",
      fill = "Higher abundance in"
    ) +
    ggplot2::theme_classic(
      base_size = 12
    ) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.y = ggplot2::element_text(
        face = "italic"
      ),
      legend.title = ggplot2::element_text(
        face = "bold"
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = 11
      ),
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      )
    )
}

# ------------------------------------------------------------
# 9. Create Figure 8 panels
# ------------------------------------------------------------

panel_ids <- significant_pairs$ContrastID

panel_letters <- LETTERS[
  seq_along(panel_ids)
]

figure_panels <- lapply(
  panel_ids,
  make_ancombc2_panel
)

names(figure_panels) <- panel_letters

# ------------------------------------------------------------
# 10. Export individual panels
# ------------------------------------------------------------

for (i in seq_along(figure_panels)) {

  save_figure(
    figure_panels[[i]],
    paste0(
      "Figure_8",
      panel_letters[i],
      "_ANCOMBC2.png"
    ),
    width = 9,
    height = 7
  )
}

# ------------------------------------------------------------
# 11. Combine Figure 8
# ------------------------------------------------------------

figure8 <- patchwork::wrap_plots(
  figure_panels,
  ncol = 2
) +
  patchwork::plot_annotation(
    tag_levels = "A"
  )

figure_height <- max(
  8,
  7 * ceiling(length(figure_panels) / 2)
)

save_figure(
  figure8,
  "Figure_8_ANCOMBC2_differential_abundance.png",
  width = 17,
  height = figure_height
)

# ------------------------------------------------------------
# 12. Export panel mapping
# ------------------------------------------------------------

panel_mapping <- significant_pairs %>%
  dplyr::mutate(
    Panel = panel_letters
  ) %>%
  dplyr::select(
    Panel,
    ContrastID,
    Group1,
    Group2,
    Contrast,
    P_adjusted
  )

write_table(
  panel_mapping,
  "Figure_8_panel_mapping.csv"
)

# ------------------------------------------------------------
# 13. Save analysis objects
# ------------------------------------------------------------

saveRDS(
  list(
    significant_pairs = significant_pairs,
    all_results = all_results,
    significant_results = significant_results,
    plot_data = plot_data,
    panel_mapping = panel_mapping,
    figure8 = figure8
  ),
  "results/objects/ancombc2_objects.rds"
)

# ------------------------------------------------------------
# 14. Completion
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "ANCOM-BC2 analysis completed\n",
  "Contrasts tested: ",
  nrow(significant_pairs), "\n",
  "Significant genus results: ",
  nrow(significant_results), "\n",
  "Figure 8 panels: ",
  length(figure_panels), "\n",
  "No predetermined plant genera\n",
  "Figures: PNG, 600 DPI, white background\n",
  "========================================\n",
  sep = ""
)
