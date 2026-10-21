source("scripts/00_setup.R")

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

dir.create(
  "results/objects",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Load genus-level data
# ------------------------------------------------------------

ps_objects <- readRDS(
  "results/phyloseq_objects.rds"
)

ps_genus <- ps_objects$genus

counts <- get_mat(ps_genus)
metadata <- get_meta(ps_genus)

if (!"SampleID" %in% names(metadata)) {
  stop("Metadata must contain SampleID")
}

if (
  anyNA(metadata$SampleID) ||
  anyDuplicated(metadata$SampleID)
) {
  stop("Missing or duplicated sample IDs")
}

if (!setequal(
  rownames(counts),
  metadata$SampleID
)) {
  stop("Count matrix and metadata do not match")
}

metadata <- metadata[
  match(rownames(counts), metadata$SampleID),
  ,
  drop = FALSE
]

rownames(metadata) <- metadata$SampleID

stopifnot(
  identical(
    rownames(counts),
    rownames(metadata)
  )
)

required_columns <- c(
  "Host",
  "Season",
  "Age",
  "Sex",
  "Body_condition"
)

missing_columns <- setdiff(
  required_columns,
  names(metadata)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_columns, collapse = ", ")
  )
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

metadata$Sex <- factor(metadata$Sex)

if (
  anyNA(metadata$Host) ||
  anyNA(metadata$Season)
) {
  stop("Host and Season must be complete")
}

if (
  any(!is.finite(counts)) ||
  any(counts < 0) ||
  any(abs(counts - round(counts)) > 1e-8)
) {
  stop("Invalid genus-level count values")
}

if (any(rowSums(counts) == 0)) {
  stop("Zero-count samples detected")
}

# Six host-season groups, ordered by host.

group_levels <- c(
  "Aoudad - Spring",
  "Aoudad - Winter",
  "Desert bighorn sheep - Spring",
  "Desert bighorn sheep - Winter",
  "Mule deer - Spring",
  "Mule deer - Winter"
)

metadata$Group <- factor(
  paste(
    metadata$Host,
    metadata$Season,
    sep = " - "
  ),
  levels = group_levels
)

if (anyNA(metadata$Group)) {
  stop("Invalid host-season group")
}

# Colors shared across Figures 2 and 6.

group_colors <- c(
  "Aoudad - Spring" = "#F8766D",
  "Aoudad - Winter" = "#C49A00",
  "Desert bighorn sheep - Spring" = "#00BA38",
  "Desert bighorn sheep - Winter" = "#00BFC4",
  "Mule deer - Spring" = "#619CFF",
  "Mule deer - Winter" = "#F564E3"
)

# ------------------------------------------------------------
# 2. Rarefy genus counts
# ------------------------------------------------------------

rarefaction_depth <- 20545L

sample_depths <- rowSums(counts)

if (any(sample_depths < rarefaction_depth)) {
  stop(
    "One or more samples contain fewer than ",
    rarefaction_depth,
    " genus-assigned reads"
  )
}

set.seed(345)

counts_rarefied <- vegan::rrarefy(
  counts,
  sample = rarefaction_depth
)

stopifnot(
  all(
    rowSums(counts_rarefied) ==
      rarefaction_depth
  )
)

ps_rarefied <- phyloseq::phyloseq(
  phyloseq::otu_table(
    counts_rarefied,
    taxa_are_rows = FALSE
  ),
  phyloseq::tax_table(
    phyloseq::tax_table(ps_genus)
  ),
  phyloseq::sample_data(metadata)
)

ps_rarefied <- phyloseq::prune_taxa(
  phyloseq::taxa_sums(ps_rarefied) > 0,
  ps_rarefied
)

saveRDS(
  ps_rarefied,
  "results/objects/phyloseq_genus_rarefied.rds"
)

write_table(
  data.frame(
    SampleID = rownames(counts_rarefied),
    ReadsBefore = sample_depths,
    ReadsAfter = rowSums(counts_rarefied)
  ),
  "rarefaction_depths.csv"
)

# ------------------------------------------------------------
# 3. Calculate alpha diversity
# ------------------------------------------------------------

alpha <- data.frame(
  SampleID = rownames(counts_rarefied),
  Observed = vegan::specnumber(
    counts_rarefied
  ),
  Shannon = vegan::diversity(
    counts_rarefied,
    index = "shannon"
  )
)

alpha <- dplyr::left_join(
  alpha,
  metadata,
  by = "SampleID"
)

alpha$Group <- factor(
  alpha$Group,
  levels = group_levels
)

write_table(
  alpha,
  "alpha_diversity_by_sample.csv"
)

alpha_summary <- alpha %>%
  dplyr::group_by(
    Host,
    Season,
    Group
  ) %>%
  dplyr::summarise(
    N = dplyr::n(),
    MeanObserved = mean(Observed),
    SDObserved = stats::sd(Observed),
    MeanShannon = mean(Shannon),
    SDShannon = stats::sd(Shannon),
    .groups = "drop"
  )

write_table(
  alpha_summary,
  "alpha_diversity_summary.csv"
)

# ------------------------------------------------------------
# 4. Alpha-diversity ANOVA
# ------------------------------------------------------------

# Retain the original measurement scale of Age
# and Body_condition. Categorical predictors
# should be stored as factors.

run_alpha_model <- function(response, data) {

  predictors <- c(
    "Host",
    "Season",
    "Age",
    "Sex",
    "Body_condition"
  )

  required <- c(
    response,
    predictors
  )

  dat <- data[
    stats::complete.cases(
      data[, required, drop = FALSE]
    ),
    ,
    drop = FALSE
  ]

  dat <- droplevels(dat)

  if (nrow(dat) < 8) {
    stop(
      "Insufficient complete observations for ",
      response
    )
  }

  for (variable in predictors) {

    values <- dat[[variable]]

    if (
      is.factor(values) &&
      nlevels(values) < 2
    ) {
      stop(
        "Predictor has fewer than two levels: ",
        variable
      )
    }

    if (
      is.numeric(values) &&
      length(unique(values)) < 2
    ) {
      stop(
        "Predictor has no variation: ",
        variable
      )
    }
  }

  model_formula <- stats::reformulate(
    predictors,
    response = response
  )

  design <- stats::model.matrix(
    model_formula,
    data = dat
  )

  if (
    qr(design)$rank < ncol(design)
  ) {
    stop(
      "Rank-deficient alpha-diversity model for ",
      response
    )
  }

  fit <- stats::aov(
    model_formula,
    data = dat
  )

  anova_table <- as.data.frame(
    stats::anova(fit)
  )

  anova_table$Term <- rownames(
    anova_table
  )

  write_table(
    anova_table,
    paste0(
      "ANOVA_",
      response,
      ".csv"
    )
  )

  # Shapiro-Wilk residual normality.

  residual_values <- stats::residuals(
    fit
  )

  if (
    length(residual_values) >= 3 &&
    length(residual_values) <= 5000
  ) {

    shapiro_result <- stats::shapiro.test(
      residual_values
    )

    capture.output(
      shapiro_result,
      file = file.path(
        "results/tables",
        paste0(
          "Shapiro_",
          response,
          ".txt"
        )
      )
    )
  }

  # Levene's test across host-season groups.

  levene_result <- car::leveneTest(
    dat[[response]],
    group = interaction(
      dat$Host,
      dat$Season,
      drop = TRUE
    ),
    center = median
  )

  capture.output(
    levene_result,
    file = file.path(
      "results/tables",
      paste0(
        "Levene_",
        response,
        ".txt"
      )
    )
  )

  # Tukey-adjusted host comparisons from
  # the fitted model.

  if (!requireNamespace(
    "emmeans",
    quietly = TRUE
  )) {
    stop(
      "Install emmeans before running ",
      "alpha-diversity models"
    )
  }

  host_emmeans <- emmeans::emmeans(
    fit,
    ~ Host
  )

  host_pairs <- as.data.frame(
    emmeans::pairs(
      host_emmeans,
      adjust = "tukey"
    )
  )

  write_table(
    host_pairs,
    paste0(
      "Tukey_Host_",
      response,
      ".csv"
    )
  )

  capture.output(
    summary(fit),
    file = file.path(
      "results/tables",
      paste0(
        "model_summary_",
        response,
        ".txt"
      )
    )
  )

  fit
}

model_observed <- run_alpha_model(
  "Observed",
  alpha
)

model_shannon <- run_alpha_model(
  "Shannon",
  alpha
)

saveRDS(
  list(
    Observed = model_observed,
    Shannon = model_shannon
  ),
  "results/objects/alpha_models.rds"
)

# ------------------------------------------------------------
# 5. Figure 6A: Observed genus richness
# ------------------------------------------------------------

alpha_theme <- ggplot2::theme_classic(
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
    axis.text.x = ggplot2::element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    legend.position = "none",
    panel.background = ggplot2::element_rect(
      fill = "white",
      color = NA
    ),
    plot.background = ggplot2::element_rect(
      fill = "white",
      color = NA
    )
  )

p_richness <- ggplot2::ggplot(
  alpha,
  ggplot2::aes(
    x = Group,
    y = Observed,
    fill = Group
  )
) +
  ggplot2::geom_violin(
    trim = FALSE,
    scale = "width",
    linewidth = 0.5,
    color = "black"
  ) +
  ggplot2::geom_jitter(
    width = 0.08,
    height = 0,
    size = 1.8,
    color = "black",
    alpha = 0.9
  ) +
  ggplot2::scale_fill_manual(
    values = group_colors,
    drop = FALSE
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Observed genus richness"
  ) +
  alpha_theme

# ------------------------------------------------------------
# 6. Figure 6B: Shannon diversity
# ------------------------------------------------------------

p_shannon <- ggplot2::ggplot(
  alpha,
  ggplot2::aes(
    x = Group,
    y = Shannon,
    fill = Group
  )
) +
  ggplot2::geom_violin(
    trim = FALSE,
    scale = "width",
    linewidth = 0.5,
    color = "black"
  ) +
  ggplot2::geom_jitter(
    width = 0.08,
    height = 0,
    size = 1.8,
    color = "black",
    alpha = 0.9
  ) +
  ggplot2::scale_fill_manual(
    values = group_colors,
    drop = FALSE
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Shannon diversity"
  ) +
  alpha_theme

# ------------------------------------------------------------
# 7. Combine Figure 6A and 6B
# ------------------------------------------------------------

figure6 <- (
  p_richness +
    p_shannon
) +
  patchwork::plot_layout(
    ncol = 2
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 16
    )
  )

save_figure(
  p_richness,
  "Figure_6A_observed_genus_richness.png",
  width = 9,
  height = 7
)

save_figure(
  p_shannon,
  "Figure_6B_shannon_diversity.png",
  width = 9,
  height = 7
)

save_figure(
  figure6,
  "Figure_6_alpha_diversity.png",
  width = 17,
  height = 8
)

# ------------------------------------------------------------
# 8. Prepare host-season sample groups
# ------------------------------------------------------------

group_ids <- split(
  rownames(counts),
  metadata$Group,
  drop = TRUE
)

# ------------------------------------------------------------
# 9. Genus accumulation curves
# ------------------------------------------------------------

accumulation_results <- lapply(
  names(group_ids),
  function(group_name) {

    ids <- group_ids[[group_name]]

    group_mat <- counts_rarefied[
      ids,
      ,
      drop = FALSE
    ]

    if (nrow(group_mat) < 2) {
      return(NULL)
    }

    set.seed(345)

    accumulation <- vegan::specaccum(
      group_mat,
      method = "random",
      permutations = 999
    )

    data.frame(
      Group = group_name,
      Samples = accumulation$sites,
      Richness = accumulation$richness,
      SD = accumulation$sd
    )
  }
)

accumulation_df <- dplyr::bind_rows(
  accumulation_results
)

accumulation_df$Group <- factor(
  accumulation_df$Group,
  levels = group_levels
)

write_table(
  accumulation_df,
  "genus_accumulation_curves.csv"
)

# ------------------------------------------------------------
# 10. Figure 2A: Genus accumulation
# ------------------------------------------------------------

p_accumulation <- ggplot2::ggplot(
  accumulation_df,
  ggplot2::aes(
    x = Samples,
    y = Richness,
    color = Group,
    fill = Group
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = pmax(0, Richness - SD),
      ymax = Richness + SD
    ),
    alpha = 0.12,
    color = NA
  ) +
  ggplot2::geom_line(
    linewidth = 1
  ) +
  ggplot2::scale_color_manual(
    values = group_colors
  ) +
  ggplot2::scale_fill_manual(
    values = group_colors
  ) +
  ggplot2::labs(
    x = "Number of samples",
    y = "Cumulative plant genera",
    color = "Host - Season",
    fill = "Host - Season"
  ) +
  ggplot2::theme_classic(
    base_size = 13
  ) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(
      face = "bold"
    ),
    axis.text = ggplot2::element_text(
      face = "bold",
      color = "black"
    ),
    legend.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.background = ggplot2::element_rect(
      fill = "white",
      color = NA
    )
  )

save_figure(
  p_accumulation,
  "Figure_2A_genus_accumulation.png",
  width = 10,
  height = 7
)

# ------------------------------------------------------------
# 11. Individual-sample rarefaction curves
# ------------------------------------------------------------

# Calculate expected genus richness at a series
# of sequencing depths from unrarefied counts.
#
# Depths are determined separately for each sample.
# No taxa or sample-specific results are predefined.

calculate_rarefaction <- function(
    sample_counts,
    sample_id,
    n_points = 100
) {

  total_reads <- sum(sample_counts)

  depths <- unique(
    as.integer(
      round(
        seq(
          1,
          total_reads,
          length.out = n_points
        )
      )
    )
  )

  richness <- vapply(
    depths,
    function(depth) {

      if (depth == total_reads) {
        return(
          sum(sample_counts > 0)
        )
      }

      as.numeric(
        vegan::rarefy(
          sample_counts,
          sample = depth
        )
      )
    },
    numeric(1)
  )

  data.frame(
    SampleID = sample_id,
    Reads = depths,
    Richness = richness
  )
}

rarefaction_results <- lapply(
  rownames(counts),
  function(sample_id) {

    calculate_rarefaction(
      sample_counts = counts[
        sample_id,
        ,
        drop = TRUE
      ],
      sample_id = sample_id
    )
  }
)

rare_df <- dplyr::bind_rows(
  rarefaction_results
)

rare_df <- dplyr::left_join(
  rare_df,
  metadata %>%
    dplyr::select(
      SampleID,
      Group
    ),
  by = "SampleID"
)

rare_df$Group <- factor(
  rare_df$Group,
  levels = group_levels
)

write_table(
  rare_df,
  "genus_rarefaction_curves.csv"
)

# ------------------------------------------------------------
# 12. Figure 2B: Rarefaction curves
# ------------------------------------------------------------

p_rarefaction <- ggplot2::ggplot(
  rare_df,
  ggplot2::aes(
    x = Reads,
    y = Richness,
    group = SampleID,
    color = Group
  )
) +
  ggplot2::geom_line(
    linewidth = 0.65,
    alpha = 0.8
  ) +
  ggplot2::geom_vline(
    xintercept = rarefaction_depth,
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.6
  ) +
  ggplot2::scale_color_manual(
    values = group_colors
  ) +
  ggplot2::labs(
    x = "Sequencing depth (reads)",
    y = "Expected genus richness",
    color = "Host - Season"
  ) +
  ggplot2::theme_classic(
    base_size = 13
  ) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(
      face = "bold"
    ),
    axis.text = ggplot2::element_text(
      face = "bold",
      color = "black"
    ),
    legend.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.background = ggplot2::element_rect(
      fill = "white",
      color = NA
    )
  )

save_figure(
  p_rarefaction,
  "Figure_2B_genus_rarefaction.png",
  width = 10,
  height = 7
)

# ------------------------------------------------------------
# 13. Combine Figure 2A and 2B
# ------------------------------------------------------------

figure2 <- (
  p_accumulation +
    p_rarefaction
) +
  patchwork::plot_layout(
    ncol = 2,
    guides = "collect"
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 16
    )
  )

save_figure(
  figure2,
  "Figure_2_accumulation_rarefaction.png",
  width = 17,
  height = 8
)

# ------------------------------------------------------------
# 14. iNEXT sample coverage
# ------------------------------------------------------------

# Pool original genus counts within each
# host-season group.

coverage_input <- lapply(
  group_ids,
  function(ids) {

    pooled <- colSums(
      counts[
        ids,
        ,
        drop = FALSE
      ]
    )

    as.numeric(
      pooled[pooled > 0]
    )
  }
)

coverage_results <- iNEXT::DataInfo(
  coverage_input,
  datatype = "abundance"
)

write_table(
  coverage_results,
  "iNEXT_sample_coverage.csv"
)

# ------------------------------------------------------------
# 15. Pianka's dietary niche overlap
# ------------------------------------------------------------

# Calculate genus proportions for each sample,
# then average within each host-season group.

proportions <- make_prop(
  counts_rarefied
)

group_profiles <- lapply(
  group_ids,
  function(ids) {

    colMeans(
      proportions[
        ids,
        ,
        drop = FALSE
      ]
    )
  }
)

profile_matrix <- do.call(
  rbind,
  group_profiles
)

# Pianka overlap:
# sum(p_i * q_i) /
# sqrt(sum(p_i^2) * sum(q_i^2))

pianka_index <- function(p, q) {

  denominator <- sqrt(
    sum(p^2) * sum(q^2)
  )

  if (denominator == 0) {
    return(NA_real_)
  }

  sum(p * q) / denominator
}

group_names <- rownames(
  profile_matrix
)

pianka_matrix <- matrix(
  NA_real_,
  nrow = length(group_names),
  ncol = length(group_names),
  dimnames = list(
    group_names,
    group_names
  )
)

for (i in seq_along(group_names)) {

  for (j in seq_along(group_names)) {

    pianka_matrix[i, j] <- pianka_index(
      profile_matrix[i, ],
      profile_matrix[j, ]
    )
  }
}

write_table(
  tibble::rownames_to_column(
    as.data.frame(
      pianka_matrix,
      check.names = FALSE
    ),
    var = "Group"
  ),
  "pianka_overlap_matrix.csv"
)

# ------------------------------------------------------------
# 16. Pairwise Pianka overlap summary
# ------------------------------------------------------------

pianka_pairs <- as.data.frame(
  as.table(pianka_matrix),
  stringsAsFactors = FALSE
)

names(pianka_pairs) <- c(
  "Group1",
  "Group2",
  "Pianka"
)

pianka_pairs <- pianka_pairs %>%
  dplyr::mutate(
    Group1 = as.character(Group1),
    Group2 = as.character(Group2)
  ) %>%
  dplyr::filter(
    Group1 < Group2
  )

group_info <- metadata %>%
  dplyr::distinct(
    Group,
    Host,
    Season
  ) %>%
  dplyr::mutate(
    Group = as.character(Group)
  )

pianka_pairs <- pianka_pairs %>%
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
    Comparison = dplyr::case_when(
      Host1 == Host2 &
        Season1 != Season2 ~
        "Within-host seasonal",
      Host1 != Host2 &
        Season1 == Season2 ~
        "Between-host same season",
      TRUE ~
        "Between-host different seasons"
    )
  )

write_table(
  pianka_pairs,
  "pianka_pairwise_comparisons.csv"
)

# ------------------------------------------------------------
# 17. Save analysis objects
# ------------------------------------------------------------

saveRDS(
  list(
    counts_rarefied = counts_rarefied,
    phyloseq_rarefied = ps_rarefied,
    alpha = alpha,
    alpha_summary = alpha_summary,
    model_observed = model_observed,
    model_shannon = model_shannon,
    accumulation = accumulation_df,
    rarefaction = rare_df,
    coverage = coverage_results,
    pianka_matrix = pianka_matrix,
    pianka_pairs = pianka_pairs,
    figure2 = figure2,
    figure6 = figure6
  ),
  "results/objects/diversity_overlap_objects.rds"
)

# ------------------------------------------------------------
# 18. Completion
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "Diversity and niche overlap completed\n",
  "Samples: ", nrow(counts_rarefied), "\n",
  "Genera after rarefaction: ",
  sum(colSums(counts_rarefied) > 0), "\n",
  "Rarefaction depth: ",
  rarefaction_depth, "\n",
  "Figure 2A: Genus accumulation\n",
  "Figure 2B: Rarefaction\n",
  "Figure 6A: Observed genus richness\n",
  "Figure 6B: Shannon diversity\n",
  "PNG: 600 DPI, white background\n",
  "========================================\n",
  sep = ""
)

