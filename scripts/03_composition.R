
source("scripts/00_setup.R")

# ------------------------------------------------------------
# 1. Load curated data
# ------------------------------------------------------------

obj <- readRDS("results/phyloseq_objects.rds")

ps_genus  <- obj$genus
ps_family <- obj$family

host_levels <- c(
  "Aoudad",
  "Desert bighorn sheep",
  "Mule deer"
)

season_levels <- c("Spring", "Winter")

# ------------------------------------------------------------
# 2. Convert phyloseq objects to relative abundance tables
# ------------------------------------------------------------

composition_table <- function(ps, rank) {

  mat <- get_mat(ps)
  meta <- get_meta(ps)

  if (any(rowSums(mat) <= 0)) {
    stop("Samples with zero total counts detected")
  }

  prop <- sweep(mat, 1, rowSums(mat), "/") * 100

  tax <- as.data.frame(
    phyloseq::tax_table(ps),
    stringsAsFactors = FALSE
  )

  if (!rank %in% names(tax)) {
    stop("Taxonomic rank not found: ", rank)
  }

  tax$TaxonID <- rownames(tax)

  long <- as.data.frame(prop) %>%
    tibble::rownames_to_column("SampleID") %>%
    tidyr::pivot_longer(
      cols = -SampleID,
      names_to = "TaxonID",
      values_to = "Percent"
    ) %>%
    dplyr::left_join(
      tax[, c("TaxonID", rank), drop = FALSE],
      by = "TaxonID"
    ) %>%
    dplyr::rename(Taxon = dplyr::all_of(rank)) %>%
    dplyr::left_join(
      meta[, c("SampleID", "Host", "Season")],
      by = "SampleID"
    )

  if (anyNA(long$Host) || anyNA(long$Season)) {
    stop("Sample metadata mismatch")
  }

  long %>%
    dplyr::mutate(
      Host = factor(
        as.character(Host),
        levels = host_levels
      ),
      Season = factor(
        as.character(Season),
        levels = season_levels
      )
    )
}

genus_long <- composition_table(
  ps_genus,
  "Genus"
)

family_long <- composition_table(
  ps_family,
  "Family"
)

write_table(
  genus_long,
  "genus_relative_abundance_per_sample.csv"
)

write_table(
  family_long,
  "family_relative_abundance_per_sample.csv"
)

# ------------------------------------------------------------
# 3. Mean abundance by host and season
# ------------------------------------------------------------

# First sum within samples, then average across samples.
# This avoids treating individual taxon records as replicates.

mean_composition <- function(df, grouping) {

  df %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(
        c("SampleID", grouping, "Taxon")
      ))
    ) %>%
    dplyr::summarise(
      Percent = sum(Percent),
      .groups = "drop"
    ) %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(
        c(grouping, "Taxon")
      ))
    ) %>%
    dplyr::summarise(
      Percent = mean(Percent),
      .groups = "drop"
    )
}

genus_season <- mean_composition(
  genus_long,
  c("Host", "Season")
)

write_table(
  genus_season,
  "genus_mean_abundance_host_season.csv"
)

# ------------------------------------------------------------
# 4. Overall family composition
# ------------------------------------------------------------

# Calculate mean family abundance within each season,
# then give Spring and Winter equal weight.
# This follows the season-averaging approach in the
# original plotting code, rather than pooling reads.

family_season <- mean_composition(
  family_long,
  c("Host", "Season")
)

family_overall <- family_season %>%
  dplyr::group_by(Host, Taxon) %>%
  dplyr::summarise(
    Percent = mean(Percent),
    .groups = "drop"
  )

write_table(
  family_overall,
  "family_overall_abundance.csv"
)

# Apply the original >=2% threshold separately by host.

family_plot <- family_overall %>%
  dplyr::mutate(
    Family = ifelse(
      Percent >= 2,
      Taxon,
      "Other"
    )
  ) %>%
  dplyr::group_by(Host, Family) %>%
  dplyr::summarise(
    Percent = sum(Percent),
    .groups = "drop"
  ) %>%
  dplyr::group_by(Host) %>%
  dplyr::mutate(
    Percent = 100 * Percent / sum(Percent)
  ) %>%
  dplyr::ungroup()

family_levels <- family_plot %>%
  dplyr::group_by(Family) %>%
  dplyr::summarise(
    Mean = mean(Percent),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    Family != "Other",
    Mean
  ) %>%
  dplyr::pull(Family)

family_plot$Family <- factor(
  family_plot$Family,
  levels = family_levels
)

p_family <- ggplot(
  family_plot,
  aes(
    x = Host,
    y = Percent,
    fill = Family
  )
) +
  geom_col(width = 0.8) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_brewer(
    palette = "Paired",
    drop = TRUE
  ) +
  labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Plant family"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(
      face = "bold",
      color = "black"
    ),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.background = element_rect(
      fill = "white",
      color = NA
    )
  )

ggsave(
  "results/figures/family_composition.png",
  plot = p_family,
  width = 10,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# 5. Functional-group composition
# ------------------------------------------------------------

# Functional groups must come from a verified mapping.
# Do not automatically classify genera by their names.

functional_file <- "config/functional_groups.csv"

if (!file.exists(functional_file)) {
  stop(
    "Missing verified functional-group mapping: ",
    functional_file,
    "\nExpected columns: Genus, FunctionalGroup"
  )
}

functional_map <- readr::read_csv(
  functional_file,
  show_col_types = FALSE
)

stopifnot(all(
  c("Genus", "FunctionalGroup") %in%
    names(functional_map)
))

functional_map <- functional_map %>%
  dplyr::distinct(Genus, FunctionalGroup)

if (anyDuplicated(functional_map$Genus)) {
  stop("Each genus must have one verified functional group")
}

missing_genera <- setdiff(
  unique(genus_long$Taxon),
  functional_map$Genus
)

if (length(missing_genera) > 0) {
  stop(
    "Functional groups missing for: ",
    paste(missing_genera, collapse = ", ")
  )
}

functional_long <- genus_long %>%
  dplyr::left_join(
    functional_map,
    by = c("Taxon" = "Genus")
  ) %>%
  dplyr::group_by(
    SampleID, Host, Season, FunctionalGroup
  ) %>%
  dplyr::summarise(
    Percent = sum(Percent),
    .groups = "drop"
  )

functional_season <- functional_long %>%
  dplyr::group_by(
    Host, Season, FunctionalGroup
  ) %>%
  dplyr::summarise(
    Percent = mean(Percent),
    .groups = "drop"
  )

functional_overall <- functional_season %>%
  dplyr::group_by(
    Host, FunctionalGroup
  ) %>%
  dplyr::summarise(
    Percent = mean(Percent),
    .groups = "drop"
  )

write_table(
  functional_overall,
  "functional_group_overall_abundance.csv"
)

p_functional <- ggplot(
  functional_overall,
  aes(
    x = Host,
    y = Percent,
    fill = FunctionalGroup
  )
) +
  geom_col(width = 0.8) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Functional group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(
      face = "bold",
      color = "black"
    ),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.background = element_rect(
      fill = "white",
      color = NA
    )
  )

ggsave(
  "results/figures/functional_composition.png",
  plot = p_functional,
  width = 10,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# 6. Genus-level seasonal composition
# ------------------------------------------------------------

# Select top 10 genera per host across both seasons.
# Remaining genera are combined as "Other".

top_genera <- genus_season %>%
  dplyr::group_by(Host, Taxon) %>%
  dplyr::summarise(
    Mean = mean(Percent),
    .groups = "drop"
  ) %>%
  dplyr::group_by(Host) %>%
  dplyr::slice_max(
    order_by = Mean,
    n = 10,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(Host, Taxon)

genus_plot <- genus_season %>%
  dplyr::left_join(
    top_genera %>%
      dplyr::mutate(Selected = TRUE),
    by = c("Host", "Taxon")
  ) %>%
  dplyr::mutate(
    Genus = ifelse(
      !is.na(Selected),
      Taxon,
      "Other"
    )
  ) %>%
  dplyr::group_by(
    Host, Season, Genus
  ) %>%
  dplyr::summarise(
    Percent = sum(Percent),
    .groups = "drop"
  )

write_table(
  genus_plot,
  "genus_top10_host_season.csv"
)

# This is a host-by-season summary plot.
# The original per-sample Figure 4 layout will be
# reconstructed separately using its verified palette.

p_genus <- ggplot(
  genus_plot,
  aes(
    x = Season,
    y = Percent,
    fill = Genus
  )
) +
  geom_col(width = 0.8) +
  facet_wrap(
    ~ Host,
    ncol = 1
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Plant genus"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(
      face = "bold",
      color = "black"
    ),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(
      face = "bold",
      size = 14
    ),
    legend.title = element_text(face = "bold"),
    plot.background = element_rect(
      fill = "white",
      color = NA
    )
  )

ggsave(
  "results/figures/genus_composition_host_season.png",
  plot = p_genus,
  width = 12,
  height = 12,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# 7. Save plot objects for combining figures
# ------------------------------------------------------------

saveRDS(
  list(
    functional = p_functional,
    family = p_family,
    genus = p_genus
  ),
  "results/composition_plots.rds"
)

cat(
  "\nComposition analysis completed.\n",
  "Figures saved as PNG at 600 DPI.\n"
)
