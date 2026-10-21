source("scripts/00_setup.R")

# ------------------------------------------------------------
# 1. Load curated data
# ------------------------------------------------------------

objects <- readRDS(
  "results/phyloseq_objects.rds"
)

ps_genus <- objects$genus
ps_family <- objects$family

metadata <- get_meta(ps_genus)

host_order <- c(
  "Aoudad",
  "Desert bighorn sheep",
  "Mule deer"
)

season_order <- c(
  "Spring",
  "Winter"
)

metadata$Host <- factor(
  metadata$Host,
  levels = host_order
)

metadata$Season <- factor(
  metadata$Season,
  levels = season_order
)

if (
  anyNA(metadata$Host) ||
  anyNA(metadata$Season)
) {
  stop("Missing or unrecognized Host/Season values")
}

if (
  anyNA(metadata$SampleID) ||
  anyDuplicated(metadata$SampleID)
) {
  stop("Missing or duplicate SampleID values")
}

# ------------------------------------------------------------
# 2. Aggregate counts by taxonomic name
# ------------------------------------------------------------

aggregate_taxa <- function(ps, rank) {

  mat <- get_mat(ps)

  taxonomy <- as(
    phyloseq::tax_table(ps),
    "matrix"
  )

  if (!rank %in% colnames(taxonomy)) {
    stop("Taxonomic rank not found: ", rank)
  }

  labels <- trimws(
    as.character(
      taxonomy[colnames(mat), rank]
    )
  )

  if (
    anyNA(labels) ||
    any(labels == "")
  ) {
    stop(
      "Missing taxonomic labels at rank: ",
      rank
    )
  }

  aggregated <- rowsum(
    t(mat),
    group = labels,
    reorder = TRUE
  )

  t(aggregated)
}

genus_counts <- aggregate_taxa(
  ps_genus,
  "Genus"
)

family_counts <- aggregate_taxa(
  ps_family,
  "Family"
)

# ------------------------------------------------------------
# 3. Match sample order
# ------------------------------------------------------------

if (!setequal(
  rownames(genus_counts),
  metadata$SampleID
)) {
  stop("Genus counts and metadata do not match")
}

if (!setequal(
  rownames(family_counts),
  metadata$SampleID
)) {
  stop("Family counts and metadata do not match")
}

genus_counts <- genus_counts[
  metadata$SampleID,
  ,
  drop = FALSE
]

family_counts <- family_counts[
  metadata$SampleID,
  ,
  drop = FALSE
]

if (
  any(rowSums(genus_counts) <= 0) ||
  any(rowSums(family_counts) <= 0)
) {
  stop("Zero-count samples detected")
}

# ------------------------------------------------------------
# 4. Sample-level relative abundances
# ------------------------------------------------------------

genus_prop <- make_prop(genus_counts)
family_prop <- make_prop(family_counts)

if (
  any(!is.finite(genus_prop)) ||
  any(!is.finite(family_prop))
) {
  stop("Non-finite relative abundances detected")
}

# ------------------------------------------------------------
# 5. Shared plotting theme
# ------------------------------------------------------------

composition_theme <- function() {

  ggplot2::theme_classic(
    base_size = 14
  ) +
    ggplot2::theme(

      axis.title = ggplot2::element_text(
        face = "bold",
        color = "black"
      ),

      axis.text = ggplot2::element_text(
        face = "bold",
        color = "#555555"
      ),

      legend.title = ggplot2::element_text(
        face = "bold",
        color = "black"
      ),

      legend.text = ggplot2::element_text(
        color = "black"
      ),

      plot.title = ggplot2::element_text(
        face = "bold",
        color = "black"
      ),

      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 17
      ),

      panel.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),

      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      )
    )
}

# ------------------------------------------------------------
# 6. Helper: Aggregate proportions by category
# ------------------------------------------------------------

aggregate_proportions <- function(
    prop_matrix,
    categories
) {

  if (
    length(categories) != ncol(prop_matrix) ||
    anyNA(categories)
  ) {
    stop("Invalid category assignments")
  }

  result <- t(
    rowsum(
      t(prop_matrix),
      group = categories,
      reorder = TRUE
    )
  )

  if (
    any(abs(rowSums(result) - 1) > 1e-8)
  ) {
    stop(
      "Aggregated proportions do not sum to one"
    )
  }

  result
}

# ------------------------------------------------------------
# 7. Helper: Mean relative abundance by host
# ------------------------------------------------------------

mean_by_host <- function(prop_matrix) {

  result <- dplyr::bind_rows(
    lapply(
      host_order,
      function(host_name) {

        ids <- metadata$SampleID[
          metadata$Host == host_name
        ]

        if (length(ids) == 0) {
          stop("No samples for ", host_name)
        }

        means <- colMeans(
          prop_matrix[
            ids,
            ,
            drop = FALSE
          ]
        )

        data.frame(
          Host = host_name,
          Category = names(means),
          MeanAbundance = unname(means),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  result$Host <- factor(
    result$Host,
    levels = host_order
  )

  result
}

# ------------------------------------------------------------
# 8. Functional group assignments
# ------------------------------------------------------------

# Required input:
# config/plant_functional_groups.csv
#
# Columns:
# Genus,FunctionalGroup
#
# Functional groups:
# Browse, Forb, Grass, Succulent

functional_file <- "config/plant_functional_groups.csv"

if (!file.exists(functional_file)) {
  stop(
    "Missing functional group mapping: ",
    functional_file
  )
}

functional_map <- readr::read_csv(
  functional_file,
  show_col_types = FALSE
)

if (!all(
  c("Genus", "FunctionalGroup") %in%
    names(functional_map)
)) {
  stop(
    "Mapping requires Genus and ",
    "FunctionalGroup columns"
  )
}

functional_map <- functional_map %>%
  dplyr::mutate(
    Genus = trimws(as.character(Genus)),
    FunctionalGroup = trimws(
      as.character(FunctionalGroup)
    )
  )

functional_map$FunctionalGroup <- dplyr::recode(
  functional_map$FunctionalGroup,
  "Browses" = "Browse",
  "Forbs" = "Forb",
  "Grasses" = "Grass",
  "Succulents" = "Succulent"
)

functional_levels <- c(
  "Browse",
  "Forb",
  "Grass",
  "Succulent"
)

if (
  anyNA(functional_map$Genus) ||
  any(functional_map$Genus == "") ||
  anyDuplicated(functional_map$Genus)
) {
  stop("Invalid or duplicate genus assignments")
}

if (
  anyNA(functional_map$FunctionalGroup) ||
  any(!functional_map$FunctionalGroup %in%
        functional_levels)
) {
  stop("Invalid functional group assignments")
}

missing_genera <- setdiff(
  colnames(genus_prop),
  functional_map$Genus
)

if (length(missing_genera) > 0) {
  stop(
    "Missing functional group assignments for: ",
    paste(missing_genera, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 9. Calculate functional group composition
# ------------------------------------------------------------

functional_categories <- functional_map$FunctionalGroup[
  match(
    colnames(genus_prop),
    functional_map$Genus
  )
]

functional_prop <- aggregate_proportions(
  genus_prop,
  functional_categories
)

functional_summary <- mean_by_host(
  functional_prop
) %>%
  dplyr::rename(
    FunctionalGroup = Category
  )

functional_summary$FunctionalGroup <- factor(
  functional_summary$FunctionalGroup,
  levels = functional_levels
)

write_table(
  functional_summary,
  "Figure_4A_functional_groups.csv"
)

# ------------------------------------------------------------
# 10. Figure 4A: Functional groups
# ------------------------------------------------------------

functional_colors <- c(
  "Browse" = "#BE70EF",
  "Forb" = "#00BFC4",
  "Grass" = "#7CAE00",
  "Succulent" = "#F8766D"
)

p_functional <- ggplot2::ggplot(
  functional_summary,
  ggplot2::aes(
    x = Host,
    y = MeanAbundance,
    fill = FunctionalGroup
  )
) +
  ggplot2::geom_col(
    width = 0.8
  ) +
  ggplot2::scale_fill_manual(
    values = functional_colors,
    breaks = c(
      "Succulent",
      "Grass",
      "Forb",
      "Browse"
    ),
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Functional group"
  ) +
  composition_theme()

# ------------------------------------------------------------
# 11. Select plant families
# ------------------------------------------------------------

# Families with overall mean relative abundance
# >= 2% are displayed individually.
# Remaining families are combined as Other.

family_threshold <- 0.02

family_mean <- sort(
  colMeans(family_prop),
  decreasing = TRUE
)

selected_families <- names(
  family_mean[
    family_mean >= family_threshold
  ]
)

if (length(selected_families) == 0) {
  stop("No families meet the 2% threshold")
}

family_selection_table <- data.frame(
  Family = names(family_mean),
  OverallMeanAbundance = unname(family_mean),
  Selected = names(family_mean) %in%
    selected_families,
  stringsAsFactors = FALSE
)

write_table(
  family_selection_table,
  "Figure_4B_family_selection.csv"
)

# ------------------------------------------------------------
# 12. Combine remaining families as Other
# ------------------------------------------------------------

family_categories <- ifelse(
  colnames(family_prop) %in%
    selected_families,
  colnames(family_prop),
  "Other"
)

family_plot_prop <- aggregate_proportions(
  family_prop,
  family_categories
)

family_summary <- mean_by_host(
  family_plot_prop
) %>%
  dplyr::rename(
    Family = Category
  )

family_order <- c(
  selected_families,
  "Other"
)

family_summary$Family <- factor(
  family_summary$Family,
  levels = family_order
)

family_totals <- family_summary %>%
  dplyr::group_by(Host) %>%
  dplyr::summarise(
    Total = sum(MeanAbundance),
    .groups = "drop"
  )

if (
  any(abs(family_totals$Total - 1) > 1e-8)
) {
  stop("Family proportions do not sum to one")
}

write_table(
  family_summary,
  "Figure_4B_family_composition.csv"
)

# ------------------------------------------------------------
# 13. Generate family colors
# ------------------------------------------------------------

family_colors <- stats::setNames(
  grDevices::hcl.colors(
    length(selected_families),
    palette = "Dynamic"
  ),
  selected_families
)

family_colors <- c(
  family_colors,
  "Other" = "#A6CEE3"
)

# ------------------------------------------------------------
# 14. Figure 4B: Plant families
# ------------------------------------------------------------

p_family <- ggplot2::ggplot(
  family_summary,
  ggplot2::aes(
    x = Host,
    y = MeanAbundance,
    fill = Family
  )
) +
  ggplot2::geom_col(
    width = 0.8
  ) +
  ggplot2::scale_fill_manual(
    values = family_colors,
    breaks = rev(family_order),
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Family"
  ) +
  composition_theme()

# ------------------------------------------------------------
# 15. Combine and export Figure 4
# ------------------------------------------------------------

figure4 <- (
  p_functional +
    p_family
) +
  patchwork::plot_layout(
    ncol = 2,
    widths = c(1, 1.2)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  )

save_figure(
  figure4,
  "Figure_4_functional_groups_families.png",
  width = 18,
  height = 8
)

save_figure(
  p_functional,
  "Figure_4A_functional_groups.png",
  width = 9,
  height = 7
)

save_figure(
  p_family,
  "Figure_4B_plant_families.png",
  width = 10,
  height = 7
)

# ------------------------------------------------------------
# 16. Select top 10 genera
# ------------------------------------------------------------

get_top_genera <- function(
    host_name,
    season_name,
    n = 10
) {

  ids <- metadata$SampleID[
    metadata$Host == host_name &
      metadata$Season == season_name
  ]

  if (length(ids) == 0) {
    stop(
      "No samples for ",
      host_name,
      " - ",
      season_name
    )
  }

  mean_abundance <- colMeans(
    genus_prop[
      ids,
      ,
      drop = FALSE
    ]
  )

  mean_abundance <- sort(
    mean_abundance[
      mean_abundance > 0
    ],
    decreasing = TRUE
  )

  names(
    head(mean_abundance, n)
  )
}

# ------------------------------------------------------------
# 17. Calculate host-season genus composition
# ------------------------------------------------------------

calculate_genus_composition <- function(
    host_name,
    season_name
) {

  ids <- metadata$SampleID[
    metadata$Host == host_name &
      metadata$Season == season_name
  ]

  mean_abundance <- colMeans(
    genus_prop[
      ids,
      ,
      drop = FALSE
    ]
  )

  top_genera <- get_top_genera(
    host_name,
    season_name,
    n = 10
  )

  selected_abundance <- mean_abundance[
    top_genera
  ]

  other_abundance <- sum(
    mean_abundance[
      !names(mean_abundance) %in%
        top_genera
    ]
  )

  data.frame(
    Host = host_name,
    Season = season_name,
    Genus = c(
      top_genera,
      "Other"
    ),
    RelativeAbundance = c(
      unname(selected_abundance),
      other_abundance
    ),
    Rank = c(
      seq_along(top_genera),
      length(top_genera) + 1L
    ),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# 18. Generate Figure 5 composition data
# ------------------------------------------------------------

figure5_data <- dplyr::bind_rows(
  lapply(
    host_order,
    function(host_name) {

      dplyr::bind_rows(
        lapply(
          season_order,
          function(season_name) {

            calculate_genus_composition(
              host_name,
              season_name
            )
          }
        )
      )
    }
  )
)

figure5_data$Host <- factor(
  figure5_data$Host,
  levels = host_order
)

figure5_data$Season <- factor(
  figure5_data$Season,
  levels = season_order
)

write_table(
  figure5_data,
  "Figure_5_genus_composition.csv"
)

top_genera_table <- figure5_data %>%
  dplyr::filter(
    Genus != "Other"
  ) %>%
  dplyr::arrange(
    Host,
    Season,
    Rank
  )

write_table(
  top_genera_table,
  "Figure_5_top_10_genera.csv"
)

# ------------------------------------------------------------
# 19. Validate genus proportions
# ------------------------------------------------------------

genus_totals <- figure5_data %>%
  dplyr::group_by(
    Host,
    Season
  ) %>%
  dplyr::summarise(
    Total = sum(RelativeAbundance),
    .groups = "drop"
  )

if (
  any(abs(genus_totals$Total - 1) > 1e-8)
) {
  stop("Genus proportions do not sum to one")
}

# ------------------------------------------------------------
# 20. Automatic genus color assignment
# ------------------------------------------------------------

# The same genus receives the same color
# across spring and winter within a host.

make_host_palette <- function(host_name) {

  host_data <- figure5_data %>%
    dplyr::filter(
      Host == host_name,
      Genus != "Other"
    )

  genera <- unique(
    as.character(host_data$Genus)
  )

  colors <- grDevices::hcl.colors(
    length(genera),
    palette = "Dynamic"
  )

  c(
    stats::setNames(
      colors,
      genera
    ),
    "Other" = "#A6CEE3"
  )
}

# ------------------------------------------------------------
# 21. Plot seasonal genus composition by host
# ------------------------------------------------------------

make_host_figure <- function(host_name) {

  plot_data <- figure5_data %>%
    dplyr::filter(
      Host == host_name
    )

  palette <- make_host_palette(
    host_name
  )

  genus_levels <- c(
    setdiff(
      names(palette),
      "Other"
    ),
    "Other"
  )

  plot_data$Genus <- factor(
    plot_data$Genus,
    levels = genus_levels
  )

  plot_data$Season <- factor(
    plot_data$Season,
    levels = season_order
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = Season,
      y = RelativeAbundance,
      fill = Genus
    )
  ) +
    ggplot2::geom_col(
      width = 0.65
    ) +
    ggplot2::scale_fill_manual(
      values = palette,
      breaks = rev(genus_levels),
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(
        accuracy = 1
      ),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Relative abundance",
      fill = "Genus"
    ) +
    composition_theme() +
    ggplot2::theme(
      legend.text = ggplot2::element_text(
        face = "italic",
        size = 10
      ),
      axis.text.x = ggplot2::element_text(
        face = "bold",
        color = "black"
      )
    )
}

# ------------------------------------------------------------
# 22. Figure 5A: Aoudad
# ------------------------------------------------------------

p_aoudad <- make_host_figure(
  "Aoudad"
)

# ------------------------------------------------------------
# 23. Figure 5B: Desert bighorn sheep
# ------------------------------------------------------------

p_bighorn <- make_host_figure(
  "Desert bighorn sheep"
)

# ------------------------------------------------------------
# 24. Figure 5C: Mule deer
# ------------------------------------------------------------

p_mule <- make_host_figure(
  "Mule deer"
)

# ------------------------------------------------------------
# 25. Export individual Figure 5 panels
# ------------------------------------------------------------

save_figure(
  p_aoudad,
  "Figure_5A_aoudad.png",
  width = 12,
  height = 6
)

save_figure(
  p_bighorn,
  "Figure_5B_desert_bighorn_sheep.png",
  width = 12,
  height = 6
)

save_figure(
  p_mule,
  "Figure_5C_mule_deer.png",
  width = 12,
  height = 6
)

# ------------------------------------------------------------
# 26. Combine and export Figure 5
# ------------------------------------------------------------

figure5 <- (
  p_aoudad /
    p_bighorn /
    p_mule
) +
  patchwork::plot_annotation(
    tag_levels = "A"
  )

save_figure(
  figure5,
  "Figure_5_seasonal_genus_composition.png",
  width = 15,
  height = 19
)

# ------------------------------------------------------------
# 27. Save analysis objects
# ------------------------------------------------------------

saveRDS(
  list(
    functional_summary = functional_summary,
    family_summary = family_summary,
    family_selection = family_selection_table,
    genus_summary = figure5_data,
    top_genera = top_genera_table,
    figure4 = figure4,
    figure5 = figure5
  ),
  "results/objects/composition_objects.rds"
)

# ------------------------------------------------------------
# 28. Completion
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "Dietary composition analysis completed\n",
  "Figure 4: Functional groups and families\n",
  "Figure 5: Seasonal genus composition\n",
  "Family selection: Overall mean >= 2%\n",
  "Genus selection: Top 10 per host-season\n",
  "No predetermined plant genera or families\n",
  "Figures: PNG, 600 DPI, white background\n",
  "========================================\n",
  sep = ""
)

