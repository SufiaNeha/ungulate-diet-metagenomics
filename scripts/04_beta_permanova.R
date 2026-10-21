source("scripts/00_setup.R")

# ------------------------------------------------------------
# 1. Load curated phyloseq objects
# ------------------------------------------------------------

objects <- readRDS(
  "results/phyloseq_objects.rds"
)

ps_genus <- objects$genus
ps_family <- objects$family

genus_counts <- get_mat(ps_genus)
family_counts <- get_mat(ps_family)

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

required_variables <- c(
  "SampleID",
  "Host",
  "Season",
  "Age",
  "Sex",
  "Body_condition"
)

missing_variables <- setdiff(
  required_variables,
  colnames(metadata)
)

if (length(missing_variables) > 0) {
  stop(
    "Missing metadata variables: ",
    paste(missing_variables, collapse = ", ")
  )
}

if (anyDuplicated(metadata$SampleID)) {
  stop("Duplicate sample IDs")
}

metadata$Host <- factor(
  metadata$Host,
  levels = host_order
)

metadata$Season <- factor(
  metadata$Season,
  levels = season_order
)

metadata$Sex <- factor(metadata$Sex)

if (
  anyNA(metadata$Host) ||
  anyNA(metadata$Season)
) {
  stop("Missing or invalid Host/Season values")
}

metadata$Group <- interaction(
  metadata$Host,
  metadata$Season,
  sep = " - ",
  drop = TRUE
)

# ------------------------------------------------------------
# 2. Match sample order
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
  any(rowSums(genus_counts) == 0) ||
  any(rowSums(family_counts) == 0)
) {
  stop("Zero-count samples detected")
}

# ------------------------------------------------------------
# 3. Calculate Bray-Curtis dissimilarity
# ------------------------------------------------------------

# Bray-Curtis is calculated from genus-level
# abundance data without rarefaction.

dist_genus <- vegan::vegdist(
  genus_counts,
  method = "bray"
)

dist_family <- vegan::vegdist(
  family_counts,
  method = "bray"
)

saveRDS(
  dist_genus,
  "results/objects/bray_genus.rds"
)

saveRDS(
  dist_family,
  "results/objects/bray_family.rds"
)

# ------------------------------------------------------------
# 4. Prepare complete-case model data
# ------------------------------------------------------------

model_variables <- c(
  "Host",
  "Season",
  "Age",
  "Sex",
  "Body_condition"
)

complete_ids <- metadata$SampleID[
  stats::complete.cases(
    metadata[, model_variables, drop = FALSE]
  )
]

model_metadata <- metadata[
  metadata$SampleID %in% complete_ids,
  ,
  drop = FALSE
]

model_metadata <- droplevels(
  model_metadata
)

model_genus <- genus_counts[
  model_metadata$SampleID,
  ,
  drop = FALSE
]

model_family <- family_counts[
  model_metadata$SampleID,
  ,
  drop = FALSE
]

if (nrow(model_metadata) < 8) {
  stop("Insufficient complete observations")
}

for (variable in model_variables) {

  values <- model_metadata[[variable]]

  if (
    is.factor(values) &&
    nlevels(values) < 2
  ) {
    stop(
      "Insufficient levels for predictor: ",
      variable
    )
  }

  if (
    is.numeric(values) &&
    length(unique(values)) < 2
  ) {
    stop(
      "No variation in predictor: ",
      variable
    )
  }
}

# Check the design matrix before fitting.

interaction_formula <- stats::as.formula(
  "~ Host * Season + Age + Sex + Body_condition"
)

design_matrix <- stats::model.matrix(
  interaction_formula,
  data = model_metadata
)

if (
  qr(design_matrix)$rank <
  ncol(design_matrix)
) {
  stop(
    "The interaction model is rank deficient. ",
    "Check demographic coding and group structure."
  )
}

write_table(
  data.frame(
    SampleID = model_metadata$SampleID,
    Host = model_metadata$Host,
    Season = model_metadata$Season,
    Age = model_metadata$Age,
    Sex = model_metadata$Sex,
    Body_condition = model_metadata$Body_condition
  ),
  "permanova_model_metadata.csv"
)

# ------------------------------------------------------------
# 5. Genus-level interaction PERMANOVA
# ------------------------------------------------------------

# Full model:
# Host * Season + Age + Sex + Body_condition
#
# Marginal tests evaluate terms conditional on
# the other terms retained in the model.

set.seed(345)

permanova_interaction <- vegan::adonis2(
  model_genus ~
    Host * Season +
    Age +
    Sex +
    Body_condition,
  data = model_metadata,
  method = "bray",
  permutations = 999,
  by = "margin"
)

interaction_table <- as.data.frame(
  permanova_interaction
)

interaction_table$Term <- rownames(
  interaction_table
)

write_table(
  interaction_table,
  "PERMANOVA_genus_interaction.csv"
)

# Overall interaction model significance.

set.seed(345)

permanova_interaction_overall <- vegan::adonis2(
  model_genus ~
    Host * Season +
    Age +
    Sex +
    Body_condition,
  data = model_metadata,
  method = "bray",
  permutations = 999
)

overall_table <- as.data.frame(
  permanova_interaction_overall
)

overall_table$Term <- rownames(
  overall_table
)

write_table(
  overall_table,
  "PERMANOVA_genus_interaction_overall.csv"
)

# ------------------------------------------------------------
# 6. Genus-level additive PERMANOVA
# ------------------------------------------------------------

set.seed(345)

permanova_additive <- vegan::adonis2(
  model_genus ~
    Host +
    Season +
    Age +
    Sex +
    Body_condition,
  data = model_metadata,
  method = "bray",
  permutations = 999,
  by = "margin"
)

additive_table <- as.data.frame(
  permanova_additive
)

additive_table$Term <- rownames(
  additive_table
)

write_table(
  additive_table,
  "PERMANOVA_genus_additive.csv"
)

# ------------------------------------------------------------
# 7. Family-level sensitivity analysis
# ------------------------------------------------------------

set.seed(345)

permanova_family <- vegan::adonis2(
  model_family ~
    Host * Season +
    Age +
    Sex +
    Body_condition,
  data = model_metadata,
  method = "bray",
  permutations = 999,
  by = "margin"
)

family_table <- as.data.frame(
  permanova_family
)

family_table$Term <- rownames(
  family_table
)

write_table(
  family_table,
  "PERMANOVA_family_interaction.csv"
)

set.seed(345)

permanova_family_overall <- vegan::adonis2(
  model_family ~
    Host * Season +
    Age +
    Sex +
    Body_condition,
  data = model_metadata,
  method = "bray",
  permutations = 999
)

family_overall_table <- as.data.frame(
  permanova_family_overall
)

family_overall_table$Term <- rownames(
  family_overall_table
)

write_table(
  family_overall_table,
  "PERMANOVA_family_interaction_overall.csv"
)

# ------------------------------------------------------------
# 8. Pairwise PERMANOVA
# ------------------------------------------------------------

# Pairwise tests compare the six host-season groups.
# Benjamini-Hochberg adjustment is applied across
# all pairwise comparisons.

group_levels <- levels(
  droplevels(metadata$Group)
)

group_pairs <- utils::combn(
  group_levels,
  2,
  simplify = FALSE
)

pairwise_results <- lapply(
  seq_along(group_pairs),
  function(i) {

    pair <- group_pairs[[i]]

    keep <- metadata$Group %in% pair

    pair_metadata <- metadata[
      keep,
      ,
      drop = FALSE
    ]

    pair_metadata$Group <- droplevels(
      pair_metadata$Group
    )

    pair_counts <- genus_counts[
      pair_metadata$SampleID,
      ,
      drop = FALSE
    ]

    if (
      nlevels(pair_metadata$Group) != 2 ||
      any(table(pair_metadata$Group) < 2)
    ) {
      return(NULL)
    }

    set.seed(345 + i)

    result <- vegan::adonis2(
      pair_counts ~ Group,
      data = pair_metadata,
      method = "bray",
      permutations = 999
    )

    result_df <- as.data.frame(
      result
    )

    data.frame(
      Group1 = pair[1],
      Group2 = pair[2],
      N1 = sum(
        pair_metadata$Group == pair[1]
      ),
      N2 = sum(
        pair_metadata$Group == pair[2]
      ),
      Df = result_df$Df[1],
      F = result_df$F[1],
      R2 = result_df$R2[1],
      P = result_df$`Pr(>F)`[1],
      stringsAsFactors = FALSE
    )
  }
)

pairwise_table <- dplyr::bind_rows(
  pairwise_results
)

if (nrow(pairwise_table) == 0) {
  stop("No valid pairwise comparisons")
}

pairwise_table$P_adjusted <- stats::p.adjust(
  pairwise_table$P,
  method = "BH"
)

pairwise_table$Significant <- (
  !is.na(pairwise_table$P_adjusted) &
    pairwise_table$P_adjusted < 0.05
)

# ------------------------------------------------------------
# 9. Classify pairwise comparisons
# ------------------------------------------------------------

group_info <- metadata %>%
  dplyr::mutate(
    Group = as.character(Group)
  ) %>%
  dplyr::distinct(
    Group,
    Host,
    Season
  )

pairwise_table <- pairwise_table %>%
  dplyr::left_join(
    group_info %>%
      dplyr::rename(
        Group1 = Group,
        Host1 = Host,
        Season1 = Season
      ),
    by = "Group1"
  ) %>%
  dplyr::left_join(
    group_info %>%
      dplyr::rename(
        Group2 = Group,
        Host2 = Host,
        Season2 = Season
      ),
    by = "Group2"
  ) %>%
  dplyr::mutate(
    ComparisonType = dplyr::case_when(
      Host1 == Host2 ~
        "Intraspecific seasonal",
      Season1 == Season2 ~
        "Interspecific within season",
      TRUE ~
        "Interspecific across seasons"
    )
  )

write_table(
  pairwise_table,
  "PERMANOVA_pairwise_all.csv"
)

significant_pairs <- pairwise_table %>%
  dplyr::filter(Significant)

write_table(
  significant_pairs,
  "PERMANOVA_pairwise_significant.csv"
)

# ------------------------------------------------------------
# 10. Homogeneity of multivariate dispersions
# ------------------------------------------------------------

# Test dispersion among the six host-season groups.

dispersion <- vegan::betadisper(
  dist_genus,
  group = metadata$Group,
  type = "median",
  bias.adjust = FALSE
)

dispersion_anova <- stats::anova(
  dispersion
)

set.seed(345)

dispersion_permutation <- vegan::permutest(
  dispersion,
  permutations = 999
)

write_table(
  data.frame(
    Term = rownames(dispersion_anova),
    as.data.frame(dispersion_anova),
    row.names = NULL,
    check.names = FALSE
  ),
  "PERMDISP_anova.csv"
)

capture.output(
  dispersion_permutation,
  file = file.path(
    "results/tables",
    "PERMDISP_permutation.txt"
  )
)

# ------------------------------------------------------------
# 11. Mean pairwise Bray-Curtis dissimilarity
# ------------------------------------------------------------

bray_matrix <- as.matrix(
  dist_genus
)

mean_bray_results <- lapply(
  group_pairs,
  function(pair) {

    ids1 <- metadata$SampleID[
      metadata$Group == pair[1]
    ]

    ids2 <- metadata$SampleID[
      metadata$Group == pair[2]
    ]

    pair_values <- as.vector(
      bray_matrix[
        ids1,
        ids2,
        drop = FALSE
      ]
    )

    data.frame(
      Group1 = pair[1],
      Group2 = pair[2],
      N = length(pair_values),
      MeanBray = mean(pair_values),
      SDBray = stats::sd(pair_values),
      stringsAsFactors = FALSE
    )
  }
)

mean_bray_table <- dplyr::bind_rows(
  mean_bray_results
)

mean_bray_table <- mean_bray_table %>%
  dplyr::left_join(
    group_info %>%
      dplyr::rename(
        Group1 = Group,
        Host1 = Host,
        Season1 = Season
      ),
    by = "Group1"
  ) %>%
  dplyr::left_join(
    group_info %>%
      dplyr::rename(
        Group2 = Group,
        Host2 = Host,
        Season2 = Season
      ),
    by = "Group2"
  ) %>%
  dplyr::mutate(
    ComparisonType = dplyr::case_when(
      Host1 == Host2 ~
        "Intraspecific seasonal",
      Season1 == Season2 ~
        "Interspecific within season",
      TRUE ~
        "Interspecific across seasons"
    )
  )

write_table(
  mean_bray_table,
  "mean_pairwise_bray_curtis.csv"
)

# ------------------------------------------------------------
# 12. NMDS ordination
# ------------------------------------------------------------

set.seed(345)

nmds <- vegan::metaMDS(
  genus_counts,
  distance = "bray",
  k = 2,
  trymax = 100,
  maxit = 200,
  autotransform = FALSE,
  trace = FALSE
)

saveRDS(
  nmds,
  "results/objects/NMDS_genus.rds"
)

# ------------------------------------------------------------
# 13. Extract NMDS sample coordinates
# ------------------------------------------------------------

site_scores <- as.data.frame(
  vegan::scores(
    nmds,
    display = "sites",
    choices = 1:2
  )
)

site_scores$SampleID <- rownames(
  site_scores
)

site_scores <- dplyr::left_join(
  site_scores,
  metadata,
  by = "SampleID"
)

write_table(
  site_scores,
  "NMDS_sample_coordinates.csv"
)

# ------------------------------------------------------------
# 14. Fit genus vectors to NMDS
# ------------------------------------------------------------

set.seed(345)

vector_fit <- vegan::envfit(
  nmds,
  genus_counts,
  permutations = 999,
  na.rm = FALSE
)

vector_scores <- as.data.frame(
  vegan::scores(
    vector_fit,
    display = "vectors",
    choices = 1:2
  )
)

vector_scores$GenusID <- rownames(
  vector_scores
)

vector_scores$R2 <- vector_fit$vectors$r

vector_scores$P <- vector_fit$vectors$pvals

vector_scores$P_adjusted <- stats::p.adjust(
  vector_scores$P,
  method = "BH"
)

# ------------------------------------------------------------
# 15. Add genus names to vector results
# ------------------------------------------------------------

taxonomy <- as(
  phyloseq::tax_table(ps_genus),
  "matrix"
)

genus_lookup <- data.frame(
  GenusID = rownames(taxonomy),
  Genus = as.character(
    taxonomy[, "Genus"]
  ),
  stringsAsFactors = FALSE
)

vector_scores <- dplyr::left_join(
  vector_scores,
  genus_lookup,
  by = "GenusID"
)

write_table(
  vector_scores,
  "NMDS_genus_vectors_all.csv"
)

# ------------------------------------------------------------
# 16. Select fitted vectors for visualization
# ------------------------------------------------------------

selected_vectors <- vector_scores %>%
  dplyr::filter(
    !is.na(P_adjusted),
    P_adjusted < 0.05
  ) %>%
  dplyr::arrange(
    P_adjusted,
    dplyr::desc(R2)
  ) %>%
  dplyr::slice_head(
    n = 15
  )

write_table(
  selected_vectors,
  "NMDS_genus_vectors_selected.csv"
)

# ------------------------------------------------------------
# 17. Scale NMDS vectors for plotting
# ------------------------------------------------------------

if (nrow(selected_vectors) > 0) {

  selected_vectors$VectorLength <- sqrt(
    selected_vectors$NMDS1^2 +
      selected_vectors$NMDS2^2
  )

  selected_vectors$NMDS1 <- (
    selected_vectors$NMDS1 *
      sqrt(selected_vectors$R2)
  )

  selected_vectors$NMDS2 <- (
    selected_vectors$NMDS2 *
      sqrt(selected_vectors$R2)
  )

  x_range <- diff(
    range(site_scores$NMDS1)
  )

  y_range <- diff(
    range(site_scores$NMDS2)
  )

  arrow_range <- max(
    abs(c(
      selected_vectors$NMDS1,
      selected_vectors$NMDS2
    ))
  )

  if (
    is.finite(arrow_range) &&
    arrow_range > 0
  ) {

    arrow_multiplier <- (
      0.35 * min(x_range, y_range)
    ) / arrow_range

    selected_vectors$NMDS1 <- (
      selected_vectors$NMDS1 *
        arrow_multiplier
    )

    selected_vectors$NMDS2 <- (
      selected_vectors$NMDS2 *
        arrow_multiplier
    )
  }
}

# ------------------------------------------------------------
# 18. NMDS colors and symbols
# ------------------------------------------------------------

host_colors <- c(
  "Aoudad" = "#E69F00",
  "Desert bighorn sheep" = "#0072B2",
  "Mule deer" = "#009E73"
)

season_shapes <- c(
  "Spring" = 16,
  "Winter" = 17
)

# ------------------------------------------------------------
# 19. Figure 7: NMDS
# ------------------------------------------------------------

p_nmds <- ggplot2::ggplot(
  site_scores,
  ggplot2::aes(
    x = NMDS1,
    y = NMDS2
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      color = Host,
      shape = Season
    ),
    size = 3.2,
    alpha = 0.9
  ) +
  ggplot2::scale_color_manual(
    values = host_colors
  ) +
  ggplot2::scale_shape_manual(
    values = season_shapes
  ) +
  ggplot2::labs(
    x = "NMDS1",
    y = "NMDS2",
    color = "Host species",
    shape = "Season",
    title = paste0(
      "Stress = ",
      round(nmds$stress, 3)
    )
  ) +
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
      color = "black"
    ),
    legend.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5
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

# ------------------------------------------------------------
# 20. Add 95% confidence ellipses
# ------------------------------------------------------------

ellipse_groups <- site_scores %>%
  dplyr::group_by(Group) %>%
  dplyr::filter(
    dplyr::n() >= 3
  ) %>%
  dplyr::ungroup()

if (nrow(ellipse_groups) > 0) {

  p_nmds <- p_nmds +
    ggplot2::stat_ellipse(
      data = ellipse_groups,
      ggplot2::aes(
        group = Group,
        color = Host,
        linetype = Season
      ),
      type = "t",
      level = 0.95,
      linewidth = 0.8,
      show.legend = FALSE
    )
}

# ------------------------------------------------------------
# 21. Add statistically selected genus vectors
# ------------------------------------------------------------

if (nrow(selected_vectors) > 0) {

  p_nmds <- p_nmds +
    ggplot2::geom_segment(
      data = selected_vectors,
      ggplot2::aes(
        x = 0,
        y = 0,
        xend = NMDS1,
        yend = NMDS2
      ),
      inherit.aes = FALSE,
      arrow = grid::arrow(
        length = grid::unit(
          0.18,
          "cm"
        )
      ),
      linewidth = 0.65,
      color = "black"
    ) +
    ggrepel::geom_text_repel(
      data = selected_vectors,
      ggplot2::aes(
        x = NMDS1,
        y = NMDS2,
        label = Genus
      ),
      inherit.aes = FALSE,
      size = 3.6,
      fontface = "italic",
      color = "black",
      max.overlaps = Inf,
      seed = 345
    )
}

# ------------------------------------------------------------
# 22. Export Figure 7
# ------------------------------------------------------------

save_figure(
  p_nmds,
  "Figure_7_NMDS.png",
  width = 12,
  height = 9
)

# ------------------------------------------------------------
# 23. Save analysis objects
# ------------------------------------------------------------

saveRDS(
  list(
    dist_genus = dist_genus,
    dist_family = dist_family,
    interaction = permanova_interaction,
    interaction_overall = permanova_interaction_overall,
    additive = permanova_additive,
    family = permanova_family,
    family_overall = permanova_family_overall,
    pairwise = pairwise_table,
    significant_pairs = significant_pairs,
    dispersion = dispersion,
    dispersion_anova = dispersion_anova,
    dispersion_permutation = dispersion_permutation,
    mean_bray = mean_bray_table,
    nmds = nmds,
    vector_fit = vector_fit,
    selected_vectors = selected_vectors,
    figure7 = p_nmds
  ),
  "results/objects/beta_diversity_objects.rds"
)

# ------------------------------------------------------------
# 24. Completion
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "Beta diversity analysis completed\n",
  "Genus-level PERMANOVA completed\n",
  "Family-level sensitivity analysis completed\n",
  "Pairwise PERMANOVA completed\n",
  "PERMDISP completed\n",
  "Mean Bray-Curtis dissimilarities completed\n",
  "NMDS stress: ",
  round(nmds$stress, 3), "\n",
  "Figure 7: NMDS\n",
  "Figures: PNG, 600 DPI, white background\n",
  "========================================\n",
  sep = ""
)
