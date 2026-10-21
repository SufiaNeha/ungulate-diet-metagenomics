source("scripts/00_setup.R")

# ------------------------------------------------------------
# 1. Load genus-level data
# ------------------------------------------------------------

ps_objects <- readRDS("results/phyloseq_objects.rds")

ps_genus <- ps_objects$genus

counts <- get_mat(ps_genus)
metadata <- get_meta(ps_genus)

metadata <- metadata[
  rownames(counts),
  ,
  drop = FALSE
]

stopifnot(
  identical(rownames(counts), rownames(metadata))
)

if (anyNA(metadata$Host) || anyNA(metadata$Season)) {
  stop("Host and Season must be complete")
}

if (any(counts < 0) ||
    any(!is.finite(counts))) {
  stop("Invalid genus abundance values")
}

if (any(abs(counts - round(counts)) > 1e-8)) {
  stop("Genus abundances must be integer counts")
}

if (any(rowSums(counts) == 0)) {
  stop("Samples with zero genus-level reads detected")
}

metadata$Host <- factor(
  metadata$Host,
  levels = host_levels
)

metadata$Season <- factor(
  metadata$Season,
  levels = season_levels
)

metadata$Group <- interaction(
  metadata$Host,
  metadata$Season,
  sep = " - ",
  drop = TRUE
)

# ------------------------------------------------------------
# 2. Rarefy to standardized sequencing depth
# ------------------------------------------------------------

rarefaction_depth <- 20545L

sample_depths <- rowSums(counts)

if (any(sample_depths < rarefaction_depth)) {
  stop(
    "Some samples have fewer than ",
    rarefaction_depth,
    " genus-assigned reads. Check input curation."
  )
}

set.seed(345)

counts_rarefied <- vegan::rrarefy(
  counts,
  sample = rarefaction_depth
)

stopifnot(
  all(rowSums(counts_rarefied) == rarefaction_depth)
)

ps_rarefied <- phyloseq::phyloseq(
  phyloseq::otu_table(
    counts_rarefied,
    taxa_are_rows = FALSE
  ),
  phyloseq::tax_table(
    phyloseq::tax_table(ps_genus)
  ),
  phyloseq::sample_data(
    metadata
  )
)

# Remove taxa absent after rarefaction

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
# 3. Alpha diversity
# ------------------------------------------------------------

alpha <- data.frame(
  SampleID = rownames(counts_rarefied),
  Observed = vegan::specnumber(counts_rarefied),
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

write_table(
  alpha,
  "alpha_diversity_by_sample.csv"
)

alpha_summary <- alpha %>%
  dplyr::group_by(Host, Season) %>%
  dplyr::summarise(
    N = dplyr::n(),
    MeanObserved = mean(Observed),
    SDObserved = sd(Observed),
    MeanShannon = mean(Shannon),
    SDShannon = sd(Shannon),
    .groups = "drop"
  )

write_table(
  alpha_summary,
  "alpha_diversity_summary.csv"
)

# ------------------------------------------------------------
# 4. ANOVA and model diagnostics
# ------------------------------------------------------------

# Age and Body_condition retain their imported
# measurement scales. Missing predictor values are
# excluded only from the corresponding fitted model.

run_alpha_model <- function(response, data) {

  predictors <- c(
    "Host",
    "Season",
    "Age",
    "Sex",
    "Body_condition"
  )

  required <- c(response, predictors)

  dat <- data[
    stats::complete.cases(data[, required]),
    required,
    drop = FALSE
  ]

  dat <- droplevels(dat)

  if (nrow(dat) < 8) {
    warning(
      "Insufficient complete observations for ",
      response,
      " ANOVA"
    )
    return(NULL)
  }

  for (variable in predictors) {
    if (is.factor(dat[[variable]]) &&
        nlevels(dat[[variable]]) < 2) {
      warning(
        "Predictor ", variable,
        " has fewer than two levels"
      )
      return(NULL)
    }
  }

  formula <- stats::reformulate(
    predictors,
    response = response
  )

  fit <- stats::lm(
    formula,
    data = dat
  )

  if (fit$rank < length(stats::coef(fit))) {
    warning(
      "Rank-deficient model for ",
      response,
      "; review predictor coding"
    )
  }

  anova_table <- as.data.frame(
    stats::anova(fit)
  )

  anova_table$Term <- rownames(anova_table)

  write_table(
    anova_table,
    paste0("ANOVA_", response, ".csv")
  )

  # Shapiro-Wilk residual normality test

  residuals_model <- stats::residuals(fit)

  if (length(residuals_model) >= 3 &&
      length(residuals_model) <= 5000) {

    shapiro_result <- stats::shapiro.test(
      residuals_model
    )

    capture.output(
      shapiro_result,
      file = file.path(
        "results/tables",
        paste0("Shapiro_", response, ".txt")
      )
    )
  }

  # Levene's test for Host x Season groups

  levene_result <- car::leveneTest(
    dat[[response]],
    interaction(
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
      paste0("Levene_", response, ".txt")
    )
  )

  # Tukey-adjusted host comparisons
  # from the fitted additive model

  if (requireNamespace("emmeans", quietly = TRUE)) {

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
      paste0("Tukey_Host_", response, ".csv")
    )

  } else {

    warning(
      "Install emmeans for adjusted ",
      "host comparisons"
    )
  }

  capture.output(
    summary(fit),
    file = file.path(
      "results/tables",
      paste0("model_summary_", response, ".txt")
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
# 5. Alpha diversity figures
# ------------------------------------------------------------

p_richness <- ggplot(
  alpha,
  aes(
    x = Host,
    y = Observed,
    fill = Season
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    aes(color = Season),
    width = 0.12,
    height = 0,
    size = 2,
    alpha = 0.85
  ) +
  scale_fill_manual(
    values = c(
      Spring = "#A7D3F2",
      Winter = "#F9D291"
    )
  ) +
  scale_color_manual(
    values = c(
      Spring = "#4A90C2",
      Winter = "#D89A35"
    )
  ) +
  labs(
    x = "Host species",
    y = "Observed genus richness",
    fill = "Season",
    color = "Season"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(
      face = "bold",
      color = "black"
    ),
    legend.title = element_text(face = "bold")
  )

p_shannon <- ggplot(
  alpha,
  aes(
    x = Host,
    y = Shannon,
    fill = Season
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    aes(color = Season),
    width = 0.12,
    height = 0,
    size = 2,
    alpha = 0.85
  ) +
  scale_fill_manual(
    values = c(
      Spring = "#A7D3F2",
      Winter = "#F9D291"
    )
  ) +
  scale_color_manual(
    values = c(
      Spring = "#4A90C2",
      Winter = "#D89A35"
    )
  ) +
  labs(
    x = "Host species",
    y = "Shannon diversity",
    fill = "Season",
    color = "Season"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(
      face = "bold",
      color = "black"
    ),
    legend.title = element_text(face = "bold")
  )

save_figure(
  p_richness,
  "alpha_observed_richness.png",
  width = 10,
  height = 6
)

save_figure(
  p_shannon,
  "alpha_shannon_diversity.png",
  width = 10,
  height = 6
)

# ------------------------------------------------------------
# 6. Species accumulation curves
# ------------------------------------------------------------

# Accumulation is calculated separately for each
# host-season group using randomized sample order.

group_ids <- split(
  rownames(counts_rarefied),
  metadata$Group,
  drop = TRUE
)

accumulation_results <- lapply(
  names(group_ids),
  function(group_name) {

    sample_ids <- group_ids[[group_name]]

    group_mat <- counts_rarefied[
      sample_ids,
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

write_table(
  accumulation_df,
  "genus_accumulation_curves.csv"
)

p_accumulation <- ggplot(
  accumulation_df,
  aes(
    x = Samples,
    y = Richness,
    color = Group,
    fill = Group
  )
) +
  geom_ribbon(
    aes(
      ymin = pmax(0, Richness - SD),
      ymax = Richness + SD
    ),
    alpha = 0.12,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  labs(
    x = "Sample size",
    y = "Number of plant genera detected",
    color = "Host - Season",
    fill = "Host - Season"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(
      face = "bold",
      color = "black"
    ),
    legend.title = element_text(face = "bold")
  )

save_figure(
  p_accumulation,
  "genus_accumulation.png",
  width = 10,
  height = 7
)

# ------------------------------------------------------------
# 7. Individual-sample rarefaction curves
# ------------------------------------------------------------

# Curves are calculated from pre-rarefaction counts.

rare_list <- vegan::rarecurve(
  counts,
  step = 500,
  label = FALSE,
  tidy = FALSE
)

rare_df <- dplyr::bind_rows(
  lapply(
    seq_along(rare_list),
    function(i) {

      curve <- rare_list[[i]]

      data.frame(
        SampleID = rownames(counts)[i],
        Reads = attr(curve, "Subsample"),
        Richness = as.numeric(curve)
      )
    }
  )
)

rare_df <- dplyr::left_join(
  rare_df,
  metadata %>%
    dplyr::select(SampleID, Group),
  by = "SampleID"
)

write_table(
  rare_df,
  "genus_rarefaction_curves.csv"
)

p_rarefaction <- ggplot(
  rare_df,
  aes(
    x = Reads,
    y = Richness,
    group = SampleID,
    color = Group
  )
) +
  geom_line(
    linewidth = 0.7,
    alpha = 0.75
  ) +
  labs(
    x = "Sequencing depth (reads)",
    y = "Observed genus richness",
    color = "Host - Season"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(
      face = "bold",
      color = "black"
    ),
    legend.title = element_text(face = "bold")
  )

save_figure(
  p_rarefaction,
  "genus_rarefaction.png",
  width = 10,
  height = 7
)

# ------------------------------------------------------------
# 8. Combined accumulation and rarefaction figure
# ------------------------------------------------------------

p_sampling <- (
  p_accumulation +
    ggplot2::ggtitle("A")
) +
  (
    p_rarefaction +
      ggplot2::ggtitle("B")
  ) +
  patchwork::plot_layout(ncol = 2)

save_figure(
  p_sampling,
  "Figure_2_accumulation_rarefaction.png",
  width = 16,
  height = 7
)

# ------------------------------------------------------------
# 9. iNEXT sample coverage
# ------------------------------------------------------------

# Pool genus counts by host-season group.
# Incidence and abundance estimators are not mixed.

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

    as.numeric(pooled[pooled > 0])
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
# 10. Pianka's niche overlap
# ------------------------------------------------------------

# Mean sample-wise genus proportions for each
# host-season group.

prop <- make_prop(counts_rarefied)

group_profiles <- lapply(
  group_ids,
  function(ids) {

    colMeans(
      prop[
        ids,
        ,
        drop = FALSE
      ]
    )
  }
)

profile_mat <- do.call(
  rbind,
  group_profiles
)

# Explicit Pianka index:
# sum(p_i * p_j) /
# sqrt(sum(p_i^2) * sum(p_j^2))

pianka_index <- function(p, q) {

  denominator <- sqrt(
    sum(p^2) * sum(q^2)
  )

  if (denominator == 0) {
    return(NA_real_)
  }

  sum(p * q) / denominator
}

group_names <- rownames(profile_mat)

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
      profile_mat[i, ],
      profile_mat[j, ]
    )
  }
}

write_table(
  tibble::rownames_to_column(
    as.data.frame(pianka_matrix),
    "Group"
  ),
  "pianka_overlap_matrix.csv"
)

# ------------------------------------------------------------
# 11. Pairwise niche overlap summary
# ------------------------------------------------------------

pianka_pairs <- as.data.frame(
  as.table(pianka_matrix)
)

names(pianka_pairs) <- c(
  "Group1",
  "Group2",
  "Pianka"
)

pianka_pairs <- pianka_pairs %>%
  dplyr::filter(
    as.character(Group1) <
      as.character(Group2)
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
      Host1 == Host2 & Season1 != Season2 ~
        "Within-host seasonal",
      Host1 != Host2 & Season1 == Season2 ~
        "Between-host same season",
      TRUE ~ "Between-host different seasons"
    )
  )

write_table(
  pianka_pairs,
  "pianka_pairwise_comparisons.csv"
)

# ------------------------------------------------------------
# 12. Save diversity objects
# ------------------------------------------------------------

saveRDS(
  list(
    counts_rarefied = counts_rarefied,
    phyloseq_rarefied = ps_rarefied,
    alpha = alpha,
    accumulation = accumulation_df,
    rarefaction = rare_df,
    coverage = coverage_results,
    pianka_matrix = pianka_matrix,
    pianka_pairs = pianka_pairs
  ),
  "results/objects/diversity_overlap_objects.rds"
)

# ------------------------------------------------------------
# 13. Final summary
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "Diversity and niche overlap completed\n",
  "Samples: ", nrow(counts_rarefied), "\n",
  "Genera after rarefaction: ",
  sum(colSums(counts_rarefied) > 0), "\n",
  "Rarefaction depth: ", rarefaction_depth, "\n",
  "Figures: PNG, 600 DPI\n",
  "========================================\n",
  sep = ""
)
