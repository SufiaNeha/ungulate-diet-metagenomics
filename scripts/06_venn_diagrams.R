
source("scripts/00_setup.R")

library(VennDiagram)
library(grid)
library(gridExtra)
library(dplyr)

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/objects",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Load rarefied genus-level data
# ------------------------------------------------------------

rarefied_file <- paste0(
  "results/objects/",
  "phyloseq_genus_rarefied.rds"
)

if (!file.exists(rarefied_file)) {
  stop(
    "Rarefied genus-level data not found. ",
    "Run scripts/02_diversity_overlap.R first."
  )
}

ps_rarefied <- readRDS(rarefied_file)

if (!inherits(ps_rarefied, "phyloseq")) {
  stop("Rarefied object must be a phyloseq object")
}

counts <- get_mat(ps_rarefied)
metadata <- get_meta(ps_rarefied)

required_columns <- c(
  "SampleID",
  "Host",
  "Season"
)

if (!all(required_columns %in% names(metadata))) {
  stop("Required metadata columns are missing")
}

if (
  anyNA(metadata$SampleID) ||
  anyDuplicated(metadata$SampleID)
) {
  stop("Missing or duplicate SampleID values")
}

if (!setequal(
  rownames(counts),
  metadata$SampleID
)) {
  stop("Counts and metadata do not match")
}

metadata <- metadata[
  match(rownames(counts), metadata$SampleID),
  ,
  drop = FALSE
]

rownames(metadata) <- metadata$SampleID

host_levels <- c(
  "Aoudad",
  "Desert bighorn sheep",
  "Mule deer"
)

season_levels <- c(
  "Spring",
  "Winter"
)

metadata$Host <- factor(
  metadata$Host,
  levels = host_levels
)

metadata$Season <- factor(
  metadata$Season,
  levels = season_levels
)

if (
  anyNA(metadata$Host) ||
  anyNA(metadata$Season)
) {
  stop("Invalid Host or Season values")
}

if (
  anyNA(counts) ||
  any(counts < 0)
) {
  stop("Invalid genus count values")
}

# ------------------------------------------------------------
# 2. Aggregate counts by genus name
# ------------------------------------------------------------

taxonomy <- as(
  phyloseq::tax_table(ps_rarefied),
  "matrix"
)

if (!"Genus" %in% colnames(taxonomy)) {
  stop("Genus rank missing from taxonomy")
}

genus_labels <- trimws(
  as.character(
    taxonomy[
      colnames(counts),
      "Genus"
    ]
  )
)

if (
  anyNA(genus_labels) ||
  any(genus_labels == "")
) {
  stop("Missing genus labels")
}

counts_by_genus <- t(
  rowsum(
    t(counts),
    group = genus_labels,
    reorder = TRUE
  )
)

# ------------------------------------------------------------
# 3. Get genera present in a sample group
# ------------------------------------------------------------

get_genera <- function(sample_ids) {

  if (length(sample_ids) == 0) {
    stop("No samples found for selected group")
  }

  selected <- counts_by_genus[
    sample_ids,
    ,
    drop = FALSE
  ]

  sort(
    colnames(selected)[
      colSums(selected) > 0
    ]
  )
}

get_host_season <- function(host, season) {

  ids <- metadata$SampleID[
    metadata$Host == host &
      metadata$Season == season
  ]

  get_genera(ids)
}

# ------------------------------------------------------------
# 4. Generate genus sets by host and season
# ------------------------------------------------------------

aoudad_spring <- get_host_season(
  "Aoudad",
  "Spring"
)

aoudad_winter <- get_host_season(
  "Aoudad",
  "Winter"
)

bighorn_spring <- get_host_season(
  "Desert bighorn sheep",
  "Spring"
)

bighorn_winter <- get_host_season(
  "Desert bighorn sheep",
  "Winter"
)

mule_spring <- get_host_season(
  "Mule deer",
  "Spring"
)

mule_winter <- get_host_season(
  "Mule deer",
  "Winter"
)

# ------------------------------------------------------------
# 5. Overall spring and winter genus sets
# ------------------------------------------------------------

overall_spring <- sort(
  unique(c(
    aoudad_spring,
    bighorn_spring,
    mule_spring
  ))
)

overall_winter <- sort(
  unique(c(
    aoudad_winter,
    bighorn_winter,
    mule_winter
  ))
)

# ------------------------------------------------------------
# 6. Define Figure 3 comparisons
# ------------------------------------------------------------

pairwise_sets <- list(

  A = list(
    title = "Overall",
    sets = list(
      Spring = overall_spring,
      Winter = overall_winter
    )
  ),

  B = list(
    title = "Aoudad",
    sets = list(
      Spring = aoudad_spring,
      Winter = aoudad_winter
    )
  ),

  C = list(
    title = "Desert bighorn sheep",
    sets = list(
      Spring = bighorn_spring,
      Winter = bighorn_winter
    )
  ),

  D = list(
    title = "Mule deer",
    sets = list(
      Spring = mule_spring,
      Winter = mule_winter
    )
  )
)

triple_sets <- list(

  E = list(
    title = "Spring",
    sets = list(
      Aoudad = aoudad_spring,
      "Mule deer" = mule_spring,
      "Desert bighorn sheep" = bighorn_spring
    )
  ),

  F = list(
    title = "Winter",
    sets = list(
      Aoudad = aoudad_winter,
      "Mule deer" = mule_winter,
      "Desert bighorn sheep" = bighorn_winter
    )
  )
)

# ------------------------------------------------------------
# 7. Pairwise overlap summaries
# ------------------------------------------------------------

pairwise_summary <- lapply(
  names(pairwise_sets),
  function(panel) {

    comparison <- pairwise_sets[[panel]]

    a <- comparison$sets[[1]]
    b <- comparison$sets[[2]]

    shared <- intersect(a, b)

    data.frame(
      Panel = panel,
      Comparison = comparison$title,
      Set1 = names(comparison$sets)[1],
      Set2 = names(comparison$sets)[2],
      Set1_Total = length(a),
      Set2_Total = length(b),
      Set1_Only = length(setdiff(a, b)),
      Shared = length(shared),
      Set2_Only = length(setdiff(b, a)),
      Union = length(union(a, b)),
      stringsAsFactors = FALSE
    )
  }
)

pairwise_summary <- dplyr::bind_rows(
  pairwise_summary
)

write_table(
  pairwise_summary,
  "Figure_3_pairwise_overlap_summary.csv"
)

# ------------------------------------------------------------
# 8. Three-host overlap summaries
# ------------------------------------------------------------

triple_overlap_stats <- function(
    sets,
    panel,
    season
) {

  A <- sets[[1]]
  B <- sets[[2]]
  C <- sets[[3]]

  ABC <- Reduce(
    intersect,
    list(A, B, C)
  )

  AB <- intersect(A, B)
  AC <- intersect(A, C)
  BC <- intersect(B, C)

  data.frame(
    Panel = panel,
    Season = season,
    Host1 = names(sets)[1],
    Host2 = names(sets)[2],
    Host3 = names(sets)[3],
    Host1_Total = length(A),
    Host2_Total = length(B),
    Host3_Total = length(C),
    Host1_Only = length(
      setdiff(A, union(B, C))
    ),
    Host2_Only = length(
      setdiff(B, union(A, C))
    ),
    Host3_Only = length(
      setdiff(C, union(A, B))
    ),
    Host1_Host2_Only = length(
      setdiff(AB, C)
    ),
    Host1_Host3_Only = length(
      setdiff(AC, B)
    ),
    Host2_Host3_Only = length(
      setdiff(BC, A)
    ),
    All_Three = length(ABC),
    Union = length(
      Reduce(union, list(A, B, C))
    ),
    stringsAsFactors = FALSE
  )
}

triple_summary <- dplyr::bind_rows(
  lapply(
    names(triple_sets),
    function(panel) {

      triple_overlap_stats(
        sets = triple_sets[[panel]]$sets,
        panel = panel,
        season = triple_sets[[panel]]$title
      )
    }
  )
)

write_table(
  triple_summary,
  "Figure_3_three_host_overlap_summary.csv"
)

# ------------------------------------------------------------
# 9. Export genus membership for all six panels
# ------------------------------------------------------------

all_comparisons <- c(
  pairwise_sets,
  triple_sets
)

membership_results <- lapply(
  names(all_comparisons),
  function(panel) {

    sets <- all_comparisons[[panel]]$sets

    all_genera <- sort(
      Reduce(union, sets)
    )

    membership <- data.frame(
      Panel = panel,
      Comparison = all_comparisons[[panel]]$title,
      Genus = all_genera,
      stringsAsFactors = FALSE
    )

    for (i in seq_along(sets)) {

      column_name <- paste0(
        "Set",
        i,
        "_Present"
      )

      membership[[column_name]] <- (
        all_genera %in% sets[[i]]
      )
    }

    membership
  }
)

membership_table <- dplyr::bind_rows(
  membership_results
)

write_table(
  membership_table,
  "Figure_3_genus_membership.csv"
)

# ------------------------------------------------------------
# 10. Pairwise Venn plotting function
# ------------------------------------------------------------

pair_venn_equal <- function(
    a,
    b,
    panel,
    title,
    rotation.degree = 0
) {

  venn_grob <- VennDiagram::draw.pairwise.venn(
    area1 = length(a),
    area2 = length(b),
    cross.area = length(intersect(a, b)),
    category = c("Spring", "Winter"),
    fill = c(
      "#A7D3F2",
      "#BFE8BF"
    ),
    alpha = c(0.7, 0.7),
    lwd = 0,
    cex = 1.2,
    cat.cex = 1.2,
    cat.col = c(
      "#0B5FA5",
      "#1B7F3A"
    ),
    cat.fontface = "bold",
    cat.pos = c(0, 0),
    cat.dist = c(-0.05, -0.05),
    scaled = FALSE,
    rotation.degree = rotation.degree,
    ind = FALSE
  )

  grid::grobTree(
    grid::textGrob(
      paste0(
        "(",
        panel,
        ") ",
        title
      ),
      x = 0.02,
      y = 0.98,
      just = c("left", "top"),
      gp = grid::gpar(
        fontsize = 16,
        fontface = "bold"
      )
    ),
    grid::gTree(
      children = venn_grob
    )
  )
}

# ------------------------------------------------------------
# 11. Three-host Venn plotting function
# ------------------------------------------------------------

triple_venn <- function(
    sets,
    panel,
    title
) {

  A <- sets[[1]]
  B <- sets[[2]]
  C <- sets[[3]]

  venn_grob <- VennDiagram::draw.triple.venn(
    area1 = length(A),
    area2 = length(B),
    area3 = length(C),
    n12 = length(intersect(A, B)),
    n13 = length(intersect(A, C)),
    n23 = length(intersect(B, C)),
    n123 = length(
      Reduce(
        intersect,
        list(A, B, C)
      )
    ),
    category = names(sets),
    fill = c(
      "#A7D3F2",
      "#BFE8BF",
      "#F9D291"
    ),
    alpha = c(0.7, 0.7, 0.7),
    lwd = 0,
    cex = 1.8,
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.col = c(
      "#0B5FA5",
      "#1B7F3A",
      "#9B6A1D"
    ),
    cat.pos = c(10, 30, 220),
    cat.dist = c(
      -0.08,
      -0.08,
      -0.10
    ),
    cat.just = list(
      c(1, 1),
      c(0, 1),
      c(0, 0.2)
    ),
    ind = FALSE
  )

  grid::grobTree(
    grid::textGrob(
      paste0(
        "(",
        panel,
        ") ",
        title
      ),
      x = 0.02,
      y = 0.98,
      just = c("left", "top"),
      gp = grid::gpar(
        fontsize = 16,
        fontface = "bold"
      )
    ),
    grid::gTree(
      children = venn_grob
    )
  )
}

# ------------------------------------------------------------
# 12. Create Figure 3A-D
# ------------------------------------------------------------

g_overall <- pair_venn_equal(
  overall_spring,
  overall_winter,
  panel = "A",
  title = "Overall"
)

g_aoudad <- pair_venn_equal(
  aoudad_spring,
  aoudad_winter,
  panel = "B",
  title = "Aoudad"
)

g_bighorn <- pair_venn_equal(
  bighorn_spring,
  bighorn_winter,
  panel = "C",
  title = "Desert bighorn sheep"
)

g_mule <- pair_venn_equal(
  mule_spring,
  mule_winter,
  panel = "D",
  title = "Mule deer",
  rotation.degree = 180
)

# ------------------------------------------------------------
# 13. Create Figure 3E-F
# ------------------------------------------------------------

g_spring <- triple_venn(
  sets = triple_sets$E$sets,
  panel = "E",
  title = "Spring"
)

g_winter <- triple_venn(
  sets = triple_sets$F$sets,
  panel = "F",
  title = "Winter"
)

# ------------------------------------------------------------
# 14. Combine Figure 3A-F
# ------------------------------------------------------------

figure3 <- gridExtra::arrangeGrob(
  g_overall,
  g_aoudad,
  g_bighorn,
  g_mule,
  g_spring,
  g_winter,
  ncol = 2,
  nrow = 3,
  widths = c(1, 1),
  heights = c(1, 1, 1)
)

# ------------------------------------------------------------
# 15. Save PNG files
# ------------------------------------------------------------

# Save grid graphics directly using a PNG device.

save_grid_png <- function(
    grob,
    filename,
    width,
    height
) {

  grDevices::png(
    filename = file.path(
      "results/figures",
      filename
    ),
    width = width,
    height = height,
    units = "in",
    res = 600,
    bg = "white",
    type = "cairo"
  )

  tryCatch(
    {
      grid::grid.newpage()
      grid::grid.draw(grob)
    },
    finally = {
      grDevices::dev.off()
    }
  )
}

individual_panels <- list(
  A = g_overall,
  B = g_aoudad,
  C = g_bighorn,
  D = g_mule,
  E = g_spring,
  F = g_winter
)

for (panel in names(individual_panels)) {

  save_grid_png(
    grob = individual_panels[[panel]],
    filename = paste0(
      "Figure_3",
      panel,
      "_venn.png"
    ),
    width = 7,
    height = 5
  )
}

save_grid_png(
  grob = figure3,
  filename = "Figure_3_venn_diagrams.png",
  width = 10,
  height = 14
)

# ------------------------------------------------------------
# 16. Save analysis objects
# ------------------------------------------------------------

saveRDS(
  list(
    pairwise_sets = pairwise_sets,
    triple_sets = triple_sets,
    pairwise_summary = pairwise_summary,
    triple_summary = triple_summary,
    membership = membership_table,
    figure3 = figure3
  ),
  "results/objects/venn_objects.rds"
)

# ------------------------------------------------------------
# 17. Completion
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "Figure 3 Venn diagrams completed\n",
  "A: Overall Spring vs Winter\n",
  "B: Aoudad Spring vs Winter\n",
  "C: Bighorn Spring vs Winter\n",
  "D: Mule deer Spring vs Winter\n",
  "E: Spring three-host overlap\n",
  "F: Winter three-host overlap\n",
  "PNG: 600 DPI, white background\n",
  "========================================\n",
  sep = ""
)
