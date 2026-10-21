# Dietary diversity and overlap among sympatric ungulates
source("scripts/00_setup.R")

obj <- readRDS("results/phyloseq_objects.rds")
gn <- obj$genus

# ------------------------------------------------------------
# 1. Rarefaction
# ------------------------------------------------------------

mat <- get_mat(gn)

if (any(abs(mat - round(mat)) > 1e-8)) {
  stop("Rarefaction requires integer counts")
}

mat <- round(mat)

depth <- min(rowSums(mat))

cat("Rarefaction depth:", depth, "\n")

# Check against the reported manuscript depth
if (depth != 20545) {
  warning(
    "Rarefaction depth differs from manuscript value of 20,545"
  )
}

set.seed(345)

rare <- vegan::rrarefy(
  mat,
  sample = depth
)

meta <- get_meta(gn)

meta <- meta[
  rownames(rare),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 2. Alpha diversity
# ------------------------------------------------------------

alpha <- meta %>%
  mutate(
    Observed = vegan::specnumber(rare),

    Shannon = vegan::diversity(
      rare,
      index = "shannon"
    )
  )

write_table(
  alpha,
  "alpha_per_sample.csv"
)

# Summary by host species and season

alpha_summary <- alpha %>%
  group_by(Host, Season) %>%
  summarise(
    N = n(),

    Mean_Richness = mean(Observed),
    SD_Richness = sd(Observed),

    Mean_Shannon = mean(Shannon),
    SD_Shannon = sd(Shannon),

    .groups = "drop"
  )

write_table(
  alpha_summary,
  "alpha_host_season_summary.csv"
)

# ------------------------------------------------------------
# 3. Observed genera by host and season
# ------------------------------------------------------------

groups <- split(
  seq_len(nrow(rare)),
  meta$Group
)

genus_counts <- lapply(
  names(groups),

  function(g) {

    ii <- groups[[g]]

    data.frame(
      Group = g,
      N = length(ii),

      DetectedGenera = sum(
        colSums(
          rare[ii, , drop = FALSE] > 0
        ) > 0
      )
    )
  }
)

write_table(
  dplyr::bind_rows(genus_counts),
  "genera_by_host_season.csv"
)

# ------------------------------------------------------------
# 4. ANOVA and assumption testing
# ------------------------------------------------------------

if ("Body condition" %in% names(alpha)) {

  names(alpha)[
    names(alpha) == "Body condition"
  ] <- "Body_condition"

}

needed <- c(
  "Host",
  "Season",
  "Age",
  "Sex",
  "Body_condition"
)

if (!all(needed %in% names(alpha))) {
  stop(
    "Missing demographic covariates required for ANOVA"
  )
}

alpha$Age <- as.numeric(alpha$Age)

alpha$Sex <- factor(alpha$Sex)

alpha$Body_condition <- factor(
  alpha$Body_condition
)

cc <- alpha[
  complete.cases(alpha[, needed]),
  ,
  drop = FALSE
]

for (metric in c("Observed", "Shannon")) {

  # Additive model

  formula_main <- as.formula(
    paste(
      metric,
      "~ Host + Season + Age + Sex + Body_condition"
    )
  )

  # Host × Season interaction model

  formula_inter <- as.formula(
    paste(
      metric,
      "~ Host * Season + Age + Sex + Body_condition"
    )
  )

  model_main <- aov(
    formula_main,
    data = cc
  )

  model_inter <- aov(
    formula_inter,
    data = cc
  )

  # Save additive ANOVA results

  capture.output(
    summary(model_main),

    file = paste0(
      "results/tables/",
      metric,
      "_ANOVA.txt"
    )
  )

  # Save interaction ANOVA results

  capture.output(
    summary(model_inter),

    file = paste0(
      "results/tables/",
      metric,
      "_interaction_ANOVA.txt"
    )
  )

  # Tukey post hoc comparison

  capture.output(
    TukeyHSD(
      model_main,
      which = "Host"
    ),

    file = paste0(
      "results/tables/",
      metric,
      "_Tukey_host.txt"
    )
  )

  # Shapiro-Wilk normality test

  capture.output(
    shapiro.test(
      residuals(model_main)
    ),

    file = paste0(
      "results/tables/",
      metric,
      "_Shapiro.txt"
    )
  )

  # Levene's test for homogeneity of variance

  capture.output(
    car::leveneTest(
      as.formula(
        paste(metric, "~ Host")
      ),
      data = cc
    ),

    file = paste0(
      "results/tables/",
      metric,
      "_Levene.txt"
    )
  )
}

# ------------------------------------------------------------
# 5. Genus accumulation curves
# ------------------------------------------------------------

acc <- lapply(
  groups,

  function(ii) {

    if (length(ii) < 2) {
      return(NULL)
    }

    a <- vegan::specaccum(
      rare[ii, , drop = FALSE],
      method = "random",
      permutations = 999
    )

    data.frame(
      N = a$sites,
      Richness = a$richness,
      SD = a$sd,

      Group = as.character(
        meta$Group[ii[1]]
      )
    )
  }
)

acc <- dplyr::bind_rows(acc)

write_table(
  acc,
  "genus_accumulation.csv"
)

p_acc <- ggplot(
  acc,
  aes(
    x = N,
    y = Richness,
    color = Group
  )
) +

  geom_line(
    linewidth = 0.9
  ) +

  geom_ribbon(
    aes(
      ymin = pmax(
        0,
        Richness - SD
      ),
      ymax = Richness + SD,
      fill = Group
    ),
    alpha = 0.12,
    color = NA
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(
    axis.title = element_text(
      face = "bold"
    ),

    axis.text = element_text(
      face = "bold"
    ),

    legend.title = element_text(
      face = "bold"
    )
  ) +

  labs(
    x = "Number of samples",
    y = "Accumulated plant genera",
    color = "Host and season",
    fill = "Host and season"
  )

ggsave(
  filename = "results/figures/accumulation.png",
  plot = p_acc,
  width = 8,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# 6. Individual-sample rarefaction curves
# ------------------------------------------------------------

rc <- vegan::rarecurve(
  mat,
  step = max(
    1,
    floor(depth / 30)
  ),
  sample = depth,
  label = FALSE,
  return = TRUE
)

rc_df <- dplyr::bind_rows(
  lapply(
    names(rc),

    function(s) {

      data.frame(
        SampleID = s,
        Reads = attr(
          rc[[s]],
          "Subsample"
        ),
        Richness = as.numeric(
          rc[[s]]
        )
      )
    }
  )
)

rc_df <- rc_df %>%
  left_join(
    meta[, c("SampleID", "Group")],
    by = "SampleID"
  )

p_rare <- ggplot(
  rc_df,
  aes(
    x = Reads,
    y = Richness,
    group = SampleID,
    color = Group
  )
) +

  geom_line(
    linewidth = 0.6,
    alpha = 0.7
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(
    axis.title = element_text(
      face = "bold"
    ),

    axis.text = element_text(
      face = "bold"
    ),

    legend.title = element_text(
      face = "bold"
    )
  ) +

  labs(
    x = "Plant sequence depth",
    y = "Observed genus richness",
    color = "Host and season"
  )

ggsave(
  filename = "results/figures/rarefaction.png",
  plot = p_rare,
  width = 8,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# 7. Sample coverage using iNEXT
# ------------------------------------------------------------

incidence <- lapply(
  groups,

  function(ii) {

    c(
      length(ii),

      colSums(
        rare[ii, , drop = FALSE] > 0
      )
    )
  }
)

coverage <- iNEXT::DataInfo(
  incidence,
  datatype = "incidence_freq"
)

write_table(
  coverage,
  "iNEXT_incidence_coverage.csv"
)

# ------------------------------------------------------------
# 8. Pianka's dietary overlap index
# ------------------------------------------------------------

prop <- make_prop(rare)

group_names <- names(groups)

profiles <- do.call(
  rbind,

  lapply(
    group_names,

    function(g) {

      colMeans(
        prop[
          as.character(meta$Group) == g,
          ,
          drop = FALSE
        ]
      )
    }
  )
)

rownames(profiles) <- group_names

# Pianka overlap function

pianka <- function(x, y) {

  sum(x * y) /
    sqrt(
      sum(x^2) *
        sum(y^2)
    )

}

P <- outer(
  seq_len(nrow(profiles)),
  seq_len(nrow(profiles)),

  Vectorize(
    function(i, j) {

      pianka(
        profiles[i, ],
        profiles[j, ]
      )

    }
  )
)

dimnames(P) <- list(
  group_names,
  group_names
)

write_table(
  data.frame(
    Group = rownames(P),
    P,
    check.names = FALSE
  ),
  "pianka_overlap_matrix.csv"
)

# ------------------------------------------------------------
# 9. Alpha diversity figures
# ------------------------------------------------------------

# Observed genus richness

p1 <- ggplot(
  alpha,
  aes(
    x = Group,
    y = Observed,
    fill = Host
  )
) +

  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.8
  ) +

  geom_jitter(
    width = 0.1,
    size = 1.5
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(
    axis.title = element_text(
      face = "bold"
    ),

    axis.text = element_text(
      face = "bold"
    ),

    axis.text.x = element_text(
      angle = 55,
      hjust = 1
    ),

    legend.title = element_text(
      face = "bold"
    )
  ) +

  labs(
    x = NULL,
    y = "Observed genus richness",
    fill = "Host species"
  )

# Shannon diversity

p2 <- ggplot(
  alpha,
  aes(
    x = Group,
    y = Shannon,
    fill = Host
  )
) +

  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.8
  ) +

  geom_jitter(
    width = 0.1,
    size = 1.5
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(
    axis.title = element_text(
      face = "bold"
    ),

    axis.text = element_text(
      face = "bold"
    ),

    axis.text.x = element_text(
      angle = 55,
      hjust = 1
    ),

    legend.title = element_text(
      face = "bold"
    )
  ) +

  labs(
    x = NULL,
    y = "Shannon diversity",
    fill = "Host species"
  )

# Combine alpha diversity panels

combined_alpha <- p1 + p2 +
  patchwork::plot_layout(
    ncol = 2
  )

ggsave(
  filename = "results/figures/alpha_diversity.png",
  plot = combined_alpha,
  width = 13,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# 10. Save analysis objects for downstream scripts
# ------------------------------------------------------------

saveRDS(
  list(
    rare = rare,
    raw = mat,
    metadata = meta,
    alpha = alpha,
    profiles = profiles
  ),
  "results/analysis_objects.rds"
)

cat(
  "\nDiversity and overlap analysis completed.\n"
)

cat(
  "Figures saved as PNG at 600 DPI.\n"
)

cat(
  "Results directory: results/\n"
)
