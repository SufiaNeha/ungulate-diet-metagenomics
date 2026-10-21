
source("scripts/00_setup.R")

# ------------------------------------------------------------
# 1. Load curated data
# ------------------------------------------------------------

objects <- readRDS("results/phyloseq_objects.rds")

ps_genus <- objects$genus
ps_family <- objects$family

genus_counts <- get_mat(ps_genus)
family_counts <- get_mat(ps_family)

metadata <- get_meta(ps_genus)

metadata <- metadata[
  rownames(genus_counts),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(genus_counts),
    rownames(metadata)
  )
)

family_counts <- family_counts[
  rownames(genus_counts),
  ,
  drop = FALSE
]

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

if (anyNA(metadata$Host) || anyNA(metadata$Season)) {
  stop("Invalid Host or Season metadata")
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

  labels <- trimws(
    as.character(
      taxonomy[colnames(mat), rank]
    )
  )

  invalid <- is.na(labels) | labels == ""

  if (any(invalid)) {
    stop(
      "Missing ", rank, " labels for ",
      sum(invalid), " taxa"
    )
  }

  t(
    rowsum(
      t(mat),
      group = labels,
      reorder = TRUE
    )
  )
}

genus_counts <- aggregate_taxa(
  ps_genus,
  "Genus"
)

family_counts <- aggregate_taxa(
  ps_family,
  "Family"
)

genus_prop <- make_prop(genus_counts)
family_prop <- make_prop(family_counts)

# ------------------------------------------------------------
# 3. Shared colors and plot theme
# ------------------------------------------------------------

host_order <- c(
  "Aoudad",
  "Desert bighorn sheep",
  "Mule deer"
)

season_order <- c("Spring", "Winter")

theme_composition <- function() {

  theme_classic(base_size = 15) +
    theme(
      axis.title = element_text(
        face = "bold",
        color = "black"
      ),
      axis.text = element_text(
        face = "bold",
        color = "#555555"
      ),
      legend.title = element_text(
        face = "bold",
        color = "black"
      ),
      legend.text = element_text(
        face = "bold",
        color = "black"
      ),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      panel.background = element_rect(
        fill = "white",
        color = NA
      ),
      plot.background = element_rect(
        fill = "white",
        color = NA
      )
    )
}

# ------------------------------------------------------------
# 4. Functional group classification
# ------------------------------------------------------------

# Mapping file columns:
# Genus,FunctionalGroup
#
# Functional groups:
# Browse, Forb, Grass, Succulent

functional_file <- "config/plant_functional_groups.csv"

if (!file.exists(functional_file)) {
  stop(
    "Missing file: ",
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
    "Functional group mapping requires ",
    "Genus and FunctionalGroup"
  )
}

functional_map <- functional_map %>%
  dplyr::mutate(
    Genus = trimws(as.character(Genus)),
    FunctionalGroup = trimws(
      as.character(FunctionalGroup)
    )
  )

# Normalize singular and plural labels.

functional_map$FunctionalGroup <- dplyr::recode(
  functional_map$FunctionalGroup,
  "Browse" = "Browse",
  "Browses" = "Browse",
  "Forb" = "Forb",
  "Forbs" = "Forb",
  "Grass" = "Grass",
  "Grasses" = "Grass",
  "Succulent" = "Succulent",
  "Succulents" = "Succulent"
)

functional_levels <- c(
  "Browse",
  "Forb",
  "Grass",
  "Succulent"
)

if (anyDuplicated(functional_map$Genus)) {
  stop("Duplicate genus assignments in mapping")
}

if (any(
  !functional_map$FunctionalGroup %in%
    functional_levels
)) {
  stop("Invalid functional group classification")
}

missing_genera <- setdiff(
  colnames(genus_prop),
  functional_map$Genus
)

if (length(missing_genera) > 0) {
  stop(
    "Functional group assignments missing for: ",
    paste(missing_genera, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 5. Functional group relative abundance
# ------------------------------------------------------------

genus_long <- as.data.frame(
  genus_prop,
  check.names = FALSE
) %>%
  tibble::rownames_to_column("SampleID") %>%
  tidyr::pivot_longer(
    cols = -SampleID,
    names_to = "Genus",
    values_to = "RelativeAbundance"
  ) %>%
  dplyr::left_join(
    metadata %>%
      dplyr::select(SampleID, Host, Season),
    by = "SampleID"
  )

functional_long <- genus_long %>%
  dplyr::left_join(
    functional_map,
    by = "Genus"
  ) %>%
  dplyr::group_by(
    SampleID,
    Host,
    FunctionalGroup
  ) %>%
  dplyr::summarise(
    RelativeAbundance = sum(RelativeAbundance),
    .groups = "drop"
  )

functional_summary <- functional_long %>%
  dplyr::group_by(
    Host,
    FunctionalGroup
  ) %>%
  dplyr::summarise(
    MeanAbundance = mean(RelativeAbundance),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    Host,
    FunctionalGroup = functional_levels,
    fill = list(MeanAbundance = 0)
  ) %>%
  dplyr::mutate(
    FunctionalGroup = factor(
      FunctionalGroup,
      levels = functional_levels
    )
  )

write_table(
  functional_summary,
  "Figure_4A_functional_group_data.csv"
)

# ------------------------------------------------------------
# 6. Figure 4A — Functional groups
# ------------------------------------------------------------

functional_colors <- c(
  "Browse" = "#BE70EF",
  "Forb" = "#00BFC4",
  "Grass" = "#7CAE00",
  "Succulent" = "#F8766D"
)

p_functional <- ggplot(
  functional_summary,
  aes(
    x = Host,
    y = MeanAbundance,
    fill = FunctionalGroup
  )
) +
  geom_col(
    width = 0.8,
    color = NA
  ) +
  scale_fill_manual(
    values = functional_colors,
    breaks = c(
      "Succulent",
      "Grass",
      "Forb",
      "Browse"
    )
  ) +
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Functional group"
  ) +
  theme_composition()

# ------------------------------------------------------------
# 7. Family relative abundance
# ------------------------------------------------------------

family_long <- as.data.frame(
  family_prop,
  check.names = FALSE
) %>%
  tibble::rownames_to_column("SampleID") %>%
  tidyr::pivot_longer(
    cols = -SampleID,
    names_to = "Family",
    values_to = "RelativeAbundance"
  ) %>%
  dplyr::left_join(
    metadata %>%
      dplyr::select(SampleID, Host, Season),
    by = "SampleID"
  )

# Families displayed in Figure 4B

figure4_families <- c(
  "Fabaceae",
  "Poaceae",
  "Amaranthaceae",
  "Asteraceae",
  "Cupressaceae",
  "Malvaceae",
  "Brassicaceae",
  "Fagaceae",
  "Cactaceae",
  "Acanthaceae"
)

family_long <- family_long %>%
  dplyr::mutate(
    Family = ifelse(
      Family %in% figure4_families,
      Family,
      "Other"
    )
  ) %>%
  dplyr::group_by(
    SampleID,
    Host,
    Family
  ) %>%
  dplyr::summarise(
    RelativeAbundance = sum(RelativeAbundance),
    .groups = "drop"
  )

family_summary <- family_long %>%
  dplyr::group_by(
    Host,
    Family
  ) %>%
  dplyr::summarise(
    MeanAbundance = mean(RelativeAbundance),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    Host,
    Family = c(
      figure4_families,
      "Other"
    ),
    fill = list(MeanAbundance = 0)
  )

family_summary$Family <- factor(
  family_summary$Family,
  levels = c(
    figure4_families,
    "Other"
  )
)

write_table(
  family_summary,
  "Figure_4B_family_data.csv"
)

# ------------------------------------------------------------
# 8. Figure 4B — Plant families
# ------------------------------------------------------------

family_colors <- c(
  "Fabaceae" = "#FFFF99",
  "Poaceae" = "#6A3D9A",
  "Amaranthaceae" = "#CAB2D6",
  "Asteraceae" = "#FF7F00",
  "Cupressaceae" = "#FDBF6F",
  "Malvaceae" = "#E31A1C",
  "Brassicaceae" = "#FB9A99",
  "Fagaceae" = "#33A02C",
  "Cactaceae" = "#B2DF8A",
  "Acanthaceae" = "#1F78B4",
  "Other" = "#A6CEE3"
)

p_family <- ggplot(
  family_summary,
  aes(
    x = Host,
    y = MeanAbundance,
    fill = Family
  )
) +
  geom_col(
    width = 0.8,
    color = NA
  ) +
  scale_fill_manual(
    values = family_colors,
    breaks = rev(
      c(
        figure4_families,
        "Other"
      )
    )
  ) +
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Family"
  ) +
  theme_composition()

# ------------------------------------------------------------
# 9. Combine and export Figure 4
# ------------------------------------------------------------

figure4 <- (
  p_functional +
    p_family
) +
  patchwork::plot_layout(
    ncol = 2,
    widths = c(1, 1.15)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  )

figure4 <- figure4 &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 18
    )
  )

save_figure(
  figure4,
  "Figure_4_functional_groups_families.png",
  width = 18,
  height = 8
)

# ------------------------------------------------------------
# 10. Figure 5 — Genus composition
# ------------------------------------------------------------

# Genus lists follow the six panels in Figure 5.

figure5_genera <- list(

  Aoudad = list(

    Spring = c(
      "Aristida",
      "Astragalus",
      "Hordeum",
      "Prosopis",
      "Juniperus",
      "Bidens",
      "Lepidium",
      "Quercus",
      "Sida",
      "Lycium"
    ),

    Winter = c(
      "Aristida",
      "Astragalus",
      "Hordeum",
      "Prosopis",
      "Atriplex",
      "Bidens",
      "Lepidium",
      "Descurainia",
      "Dyschoriste",
      "Sida"
    )
  ),

  "Desert bighorn sheep" = list(

    Spring = c(
      "Muhlenbergia",
      "Astragalus",
      "Sphaeralcea",
      "Prosopis",
      "Atriplex",
      "Brickellia",
      "Lepidium",
      "Quercus",
      "Rhus",
      "Chenopodium"
    ),

    Winter = c(
      "Muhlenbergia",
      "Astragalus",
      "Sphaeralcea",
      "Prosopis",
      "Atriplex",
      "Brickellia",
      "Lepidium",
      "Quercus",
      "Triticum",
      "Chenopodium"
    )
  ),

  "Mule deer" = list(

    Spring = c(
      "Acaciella",
      "Prosopis",
      "Atriplex",
      "Juniperus",
      "Gossypium",
      "Chrysothamnus",
      "Digitaria",
      "Artemisia",
      "Triticum",
      "Euphorbia"
    ),

    Winter = c(
      "Acaciella",
      "Prosopis",
      "Atriplex",
      "Juniperus",
      "Gossypium",
      "Chrysothamnus",
      "Opuntia",
      "Digitaria",
      "Yucca",
      "Brassica"
    )
  )
)

# ------------------------------------------------------------
# 11. Figure 5 colors
# ------------------------------------------------------------

# Consistent colors for the displayed genus ranks.

genus_rank_colors <- c(
  "#B15928",
  "#6A3D9A",
  "#CAB2D6",
  "#E9965B",
  "#FF7F00",
  "#FDBF6F",
  "#E31A1C",
  "#FB9A99",
  "#1F78B4",
  "#A6CEE3"
)

# ------------------------------------------------------------
# 12. Calculate seasonal mean composition
# ------------------------------------------------------------

get_season_composition <- function(
    host_name,
    season_name
) {

  selected_samples <- metadata$SampleID[
    metadata$Host == host_name &
      metadata$Season == season_name
  ]

  if (length(selected_samples) == 0) {
    stop(
      "No samples for ",
      host_name,
      " - ",
      season_name
    )
  }

  selected_matrix <- genus_prop[
    selected_samples,
    ,
    drop = FALSE
  ]

  mean_abundance <- colMeans(
    selected_matrix
  )

  selected_genera <- figure5_genera[
    [host_name]
  ][[season_name]]

  absent_genera <- setdiff(
    selected_genera,
    names(mean_abundance)
  )

  if (length(absent_genera) > 0) {
    warning(
      "Genera absent from curated data for ",
      host_name, " - ", season_name, ": ",
      paste(absent_genera, collapse = ", ")
    )
  }

  values <- stats::setNames(
    numeric(length(selected_genera)),
    selected_genera
  )

  present <- intersect(
    selected_genera,
    names(mean_abundance)
  )

  values[present] <- mean_abundance[present]

  other <- max(
    0,
    1 - sum(values)
  )

  data.frame(
    Host = host_name,
    Season = season_name,
    Genus = c(
      selected_genera,
      "Other"
    ),
    RelativeAbundance = c(
      unname(values),
      other
    ),
    Rank = c(
      seq_along(selected_genera),
      11L
    ),
    stringsAsFactors = FALSE
  )
}

figure5_data <- dplyr::bind_rows(
  lapply(
    host_order,
    function(host_name) {

      dplyr::bind_rows(
        lapply(
          season_order,
          function(season_name) {

            get_season_composition(
              host_name,
              season_name
            )
          }
        )
      )
    }
  )
)

write_table(
  figure5_data,
  "Figure_5_genus_composition_data.csv"
)

# ------------------------------------------------------------
# 13. Plot each Figure 5 panel
# ------------------------------------------------------------

make_figure5_panel <- function(
    host_name,
    season_name
) {

  plot_data <- figure5_data %>%
    dplyr::filter(
      Host == host_name,
      Season == season_name
    )

  genus_order <- c(
    plot_data$Genus[
      plot_data$Genus != "Other"
    ],
    "Other"
  )

  plot_data$Genus <- factor(
    plot_data$Genus,
    levels = genus_order
  )

  palette <- c(
    stats::setNames(
      genus_rank_colors,
      genus_order[
        genus_order != "Other"
      ]
    ),
    "Other" = "#A6CEE3"
  )

  ggplot(
    plot_data,
    aes(
      x = Season,
      y = RelativeAbundance,
      fill = Genus
    )
  ) +
    geom_col(
      width = 0.7,
      color = NA
    ) +
    scale_fill_manual(
      values = palette,
      breaks = rev(genus_order)
    ) +
    scale_y_continuous(
      labels = scales::percent_format(
        accuracy = 1
      ),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    labs(
      x = NULL,
      y = "Relative abundance",
      fill = "Genus"
    ) +
    theme_composition() +
    theme(
      legend.text = element_text(
        face = "italic",
        size = 11
      ),
      legend.position = "right",
      axis.text.x = element_text(
        face = "bold",
        color = "black"
      )
    )
}

# ------------------------------------------------------------
# 14. Figure 5A — Aoudad
# ------------------------------------------------------------

aoudad_spring <- make_figure5_panel(
  "Aoudad",
  "Spring"
)

aoudad_winter <- make_figure5_panel(
  "Aoudad",
  "Winter"
)

figure5_aoudad <- (
  aoudad_spring +
    aoudad_winter
) +
  patchwork::plot_layout(
    ncol = 2
  ) +
  patchwork::plot_annotation(
    title = "Aoudad",
    tag_levels = "A"
  )

# ------------------------------------------------------------
# 15. Figure 5B — Desert bighorn sheep
# ------------------------------------------------------------

bighorn_spring <- make_figure5_panel(
  "Desert bighorn sheep",
  "Spring"
)

bighorn_winter <- make_figure5_panel(
  "Desert bighorn sheep",
  "Winter"
)

figure5_bighorn <- (
  bighorn_spring +
    bighorn_winter
) +
  patchwork::plot_layout(
    ncol = 2
  ) +
  patchwork::plot_annotation(
    title = "Desert bighorn sheep"
  )

# ------------------------------------------------------------
# 16. Figure 5C — Mule deer
# ------------------------------------------------------------

mule_spring <- make_figure5_panel(
  "Mule deer",
  "Spring"
)

mule_winter <- make_figure5_panel(
  "Mule deer",
  "Winter"
)

figure5_mule <- (
  mule_spring +
    mule_winter
) +
  patchwork::plot_layout(
    ncol = 2
  ) +
  patchwork::plot_annotation(
    title = "Mule deer"
  )

# ------------------------------------------------------------
# 17. Combine Figure 5
# ------------------------------------------------------------

figure5 <- (
  (aoudad_spring | aoudad_winter) /
    (bighorn_spring | bighorn_winter) /
    (mule_spring | mule_winter)
) +
  patchwork::plot_layout(
    heights = c(1, 1, 1)
  )

save_figure(
  figure5,
  "Figure_5_genus_composition.png",
  width = 18,
  height = 20
)

# ------------------------------------------------------------
# 18. Export individual Figure 5 panels
# ------------------------------------------------------------

save_figure(
  aoudad_spring,
  "Figure_5A_aoudad_spring.png",
  width = 9,
  height = 6
)

save_figure(
  aoudad_winter,
  "Figure_5A_aoudad_winter.png",
  width = 9,
  height = 6
)

save_figure(
  bighorn_spring,
  "Figure_5B_bighorn_spring.png",
  width = 9,
  height = 6
)

save_figure(
  bighorn_winter,
  "Figure_5B_bighorn_winter.png",
  width = 9,
  height = 6
)

save_figure(
  mule_spring,
  "Figure_5C_mule_deer_spring.png",
  width = 9,
  height = 6
)

save_figure(
  mule_winter,
  "Figure_5C_mule_deer_winter.png",
  width = 9,
  height = 6
)

# ------------------------------------------------------------
# 19. Save objects
# ------------------------------------------------------------

saveRDS(
  list(
    functional = functional_summary,
    family = family_summary,
    genus = figure5_data,
    figure4 = figure4,
    figure5 = figure5
  ),
  "results/objects/composition_objects.rds"
)

# ------------------------------------------------------------
# 20. Completion
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "Composition analysis completed\n",
  "Figure 4: Functional groups and families\n",
  "Figure 5: Seasonal genus composition\n",
  "Resolution: 600 DPI\n",
  "Background: White\n",
  "========================================\n",
  sep = ""
)
