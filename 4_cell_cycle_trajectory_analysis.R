# 1. Settings

H5AD_PATH <- file.path(
  "./1_scRNAseq_preprocessing_and_RGC_extraction_results",
  "subset_cluster3_for_tricycle.h5ad"
)

# The gene set used for the manuscript analysis was generated on April 24, 2026
# and is provided in this repository as
# 'GO_union_genes_expressed_in_cluster3_RGC_pct05_gene_symbols_only_2026-04-24.csv'.
GENESET_CSV <- "/path/to/GO_union_genes_expressed_in_cluster3_RGC_pct05_gene_symbols_only_2026-04-24.csv"

DATASET_LABEL <- "RGC_GO_pct05"
SPECIES <- "mouse"

LOESS_SPAN <- 0.25
N_GRID <- 300
N_PERM <- 1000
FDR_CUTOFF <- 0.05
RANDOM_SEED <- 123

OUTDIR <- "./4_cell_cycle_trajectory_analysis_results"

PLOT_DIR <- file.path(OUTDIR, "gene_theta_plots_all")
HIST_DIR <- file.path(OUTDIR, "permutation_histograms_all")
SIG_PLOT_DIR <- file.path(OUTDIR, "gene_theta_plots_significant")
SIG_HIST_DIR <- file.path(OUTDIR, "permutation_histograms_significant")
QC_DIR <- file.path(OUTDIR, "RGC_cell_cycle_marker_group_plots")
TWO_PANEL_DIR <- file.path(OUTDIR, "gene_theta_two_panel_plots_all")
SIG_TWO_PANEL_DIR <- file.path(OUTDIR, "gene_theta_two_panel_plots_significant")

for (d in c(
  OUTDIR,
  PLOT_DIR,
  HIST_DIR,
  SIG_PLOT_DIR,
  SIG_HIST_DIR,
  QC_DIR,
  TWO_PANEL_DIR,
  SIG_TWO_PANEL_DIR
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

set.seed(RANDOM_SEED)



# 2. Packages

required_pkgs <- c(
  "ggplot2",
  "scales",
  "egg",
  "gridExtra",
  "tricycle",
  "zellkonverter",
  "SingleCellExperiment",
  "SummarizedExperiment"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    "Required packages are not installed: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(egg)
  library(grid)
  library(gridExtra)
  library(tricycle)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})



# 3. Helper functions

safe_filename <- function(x) {
  x <- gsub("/", "_", x)
  x <- gsub(" ", "_", x)
  gsub("[^A-Za-z0-9_\\-\\.]", "_", x)
}


get_gene_actual <- function(gene_to_find, sce) {
  gene_lookup <- setNames(rownames(sce), tolower(rownames(sce)))
  gene_actual <- unname(gene_lookup[tolower(gene_to_find)])

  if (length(gene_actual) == 0 || is.na(gene_actual)) {
    return(NA_character_)
  }

  gene_actual
}


fit_periodic_loess <- function(theta, y, theta_grid, span = 0.25) {
  keep <- is.finite(theta) & is.finite(y)
  theta <- theta[keep]
  y <- y[keep]

  if (length(theta) < 10) {
    return(rep(NA_real_, length(theta_grid)))
  }

  if (stats::sd(y, na.rm = TRUE) == 0) {
    return(rep(mean(y, na.rm = TRUE), length(theta_grid)))
  }

  loess_df <- rbind(
    data.frame(theta = theta - 2*pi, y = y),
    data.frame(theta = theta,        y = y),
    data.frame(theta = theta + 2*pi, y = y)
  )

  fit_obj <- tryCatch(
    loess(
      y ~ theta,
      data = loess_df,
      span = span,
      degree = 2,
      surface = "direct",
      control = loess.control(surface = "direct")
    ),
    error = function(e) NULL
  )

  if (is.null(fit_obj)) {
    return(rep(NA_real_, length(theta_grid)))
  }

  pred <- tryCatch(
    predict(fit_obj, newdata = data.frame(theta = theta_grid)),
    error = function(e) rep(NA_real_, length(theta_grid))
  )

  as.numeric(pred)
}


calc_curve_stats <- function(theta, y, theta_grid, f_obs) {
  keep <- is.finite(theta) & is.finite(y)
  theta <- theta[keep]
  y <- y[keep]

  if (all(!is.finite(f_obs))) {
    return(list(
      T = NA_real_,
      peak_angle = NA_real_,
      trough_angle = NA_real_,
      amplitude = NA_real_,
      standardized_amplitude = NA_real_,
      pseudo_R2 = NA_real_,
      SSE_loess = NA_real_,
      SSE_flat = NA_real_
    ))
  }

  f_mean <- mean(f_obs, na.rm = TRUE)
  T_val <- mean((f_obs - f_mean)^2, na.rm = TRUE)

  peak_idx <- which.max(f_obs)
  trough_idx <- which.min(f_obs)

  peak_angle <- theta_grid[peak_idx]
  trough_angle <- theta_grid[trough_idx]
  amplitude <- max(f_obs, na.rm = TRUE) - min(f_obs, na.rm = TRUE)

  y_sd <- stats::sd(y, na.rm = TRUE)

  if (is.finite(y_sd) && y_sd > 0) {
    standardized_amplitude <- amplitude / y_sd
  } else {
    standardized_amplitude <- NA_real_
  }

  fitted_at_cells <- fit_periodic_loess(
    theta,
    y,
    theta,
    span = LOESS_SPAN
  )

  SSE_loess <- sum((y - fitted_at_cells)^2, na.rm = TRUE)
  SSE_flat <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)

  if (is.finite(SSE_flat) && SSE_flat > 0) {
    pseudo_R2 <- 1 - SSE_loess / SSE_flat
  } else {
    pseudo_R2 <- NA_real_
  }

  list(
    T = T_val,
    peak_angle = peak_angle,
    trough_angle = trough_angle,
    amplitude = amplitude,
    standardized_amplitude = standardized_amplitude,
    pseudo_R2 = pseudo_R2,
    SSE_loess = SSE_loess,
    SSE_flat = SSE_flat
  )
}


make_gene_plot <- function(plot_df, pred_df, gene_name, result_row, outfile) {
  theta_color_scale <- scale_color_gradientn(
    colours = grDevices::rainbow(100),
    limits = c(0, 2*pi),
    breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
    labels = c("0", "0.5pi", "pi", "1.5pi", "2pi"),
    oob = scales::squish
  )

  p <- ggplot(plot_df, aes(x = theta, y = expr, color = theta)) +
    geom_point(size = 0.85, alpha = 0.55) +
    geom_line(
      data = pred_df,
      aes(x = theta, y = fitted),
      inherit.aes = FALSE,
      linewidth = 1.1,
      color = "black"
    ) +
    scale_x_continuous(
      limits = c(0, 2*pi),
      breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
      labels = c("0", "0.5pi", "pi", "1.5pi", "2pi")
    ) +
    theta_color_scale +
    labs(
      title = paste0(DATASET_LABEL, ": ", gene_name, " vs tricycle theta"),
      subtitle = paste0(
        "T=", sprintf("%.5f", result_row$T),
        "   p=", sprintf("%.4g", result_row$p_value),
        "   FDR=", sprintf("%.4g", result_row$FDR),
        "   amp=", sprintf("%.4f", result_row$amplitude),
        "   peak=", sprintf("%.3f", result_row$peak_angle)
      ),
      x = "tricycle theta",
      y = paste0(gene_name, " expression"),
      color = "theta"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 8.7),
      plot.subtitle = element_text(size = 6.8)
    )

  ggsave(outfile, p, width = 6.0, height = 4.5, dpi = 300)
  p
}


make_hist_plot <- function(T_perm, T_obs, gene_name, result_row, outfile) {
  T_all <- c(T_obs, T_perm)
  T_all <- T_all[is.finite(T_all)]

  if (length(T_all) < 2) {
    return(NULL)
  }

  h <- hist(T_all, breaks = 40, plot = FALSE)

  hist_df <- data.frame(
    xmin = h$breaks[-length(h$breaks)],
    xmax = h$breaks[-1],
    count = h$counts
  )

  hist_df$contains_observed <- (
    hist_df$xmin <= T_obs &
      T_obs <= hist_df$xmax
  )

  if (!any(hist_df$contains_observed)) {
    closest_idx <- which.min(
      abs((hist_df$xmin + hist_df$xmax) / 2 - T_obs)
    )
    hist_df$contains_observed[closest_idx] <- TRUE
  }

  p <- ggplot(
    hist_df,
    aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = count)
  ) +
    geom_rect(aes(fill = contains_observed), color = "white") +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
    geom_vline(xintercept = T_obs, color = "red", linewidth = 0.9) +
    labs(
      title = paste0(gene_name, ": permutation distribution of T"),
      subtitle = paste0(
        "Observed T=", sprintf("%.5f", T_obs),
        "   p=", sprintf("%.4g", result_row$p_value),
        "   FDR=", sprintf("%.4g", result_row$FDR),
        "   B=", N_PERM
      ),
      x = "T = mean((f(theta_grid) - mean(f))^2)",
      y = "Count",
      fill = "Observed bin"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "right")

  ggsave(outfile, p, width = 6.0, height = 4.5, dpi = 300)
  p
}


format_axis_label_05 <- function(x) {
  ifelse(
    abs(x - round(x)) < 1e-9,
    sprintf("%d", as.integer(round(x))),
    sprintf("%.1f", x)
  )
}


make_axis_05_equal_max4 <- function(max_value) {
  if (!is.finite(max_value) || max_value <= 0) {
    max_value <- 0.5
  }

  step_candidates <- seq(0.5, 100, by = 0.5)

  for (step in step_candidates) {
    axis_max <- step * 3

    if (axis_max >= max_value) {
      breaks <- seq(0, axis_max, by = step)

      return(list(
        axis_min = 0,
        axis_max = axis_max,
        breaks = breaks,
        labels = format_axis_label_05(breaks),
        step = step
      ))
    }
  }

  stop("Failed to define axis breaks.")
}



# 4. Read RGC H5AD

if (!file.exists(H5AD_PATH)) {
  stop("H5AD file was not found: ", H5AD_PATH)
}

message("Reading H5AD: ", H5AD_PATH)
sce <- zellkonverter::readH5AD(H5AD_PATH)

message("Loaded RGC object:")
print(sce)



# 5. Detect gene ID type

gene_names <- rownames(sce)

if (is.null(gene_names)) {
  stop("rownames(sce) is NULL. Gene names are required.")
}

ensembl_fraction <- mean(grepl("^ENSMUSG|^ENSG", gene_names))

if (is.na(ensembl_fraction)) {
  ensembl_fraction <- 0
}

if (ensembl_fraction > 0.5) {
  GNAME_TYPE <- "ENSEMBL"
} else {
  GNAME_TYPE <- "SYMBOL"
}

message("Detected gene name type: ", GNAME_TYPE)



# 6. Prepare logcounts

assay_names_now <- SummarizedExperiment::assayNames(sce)

message(
  "Assays before processing: ",
  paste(assay_names_now, collapse = ", ")
)

if (!"logcounts" %in% assay_names_now) {
  if ("X" %in% assay_names_now) {
    SummarizedExperiment::assay(
      sce,
      "logcounts"
    ) <- SummarizedExperiment::assay(
      sce,
      "X"
    )
  } else {
    stop(
      "Expected 'X' assay containing the log-normalized expression matrix."
    )
  }
}

message(
  "Assays after processing: ",
  paste(SummarizedExperiment::assayNames(sce), collapse = ", ")
)



# 7. Run tricycle

message("Running project_cycle_space() ...")

sce <- tricycle::project_cycle_space(
  sce,
  species = SPECIES,
  gname.type = GNAME_TYPE,
  exprs_values = "logcounts"
)

message("Running estimate_cycle_position() ...")

sce <- tricycle::estimate_cycle_position(
  sce,
  species = SPECIES,
  gname.type = GNAME_TYPE,
  exprs_values = "logcounts",
  dimred = "tricycleEmbedding"
)



# 8. Collect metadata and define theta

meta_df <- as.data.frame(
  SummarizedExperiment::colData(sce)
)

meta_df$cell_id <- colnames(sce)

emb <- SingleCellExperiment::reducedDim(
  sce,
  "tricycleEmbedding"
)

meta_df$tricycle_PC1 <- emb[, 1]
meta_df$tricycle_PC2 <- emb[, 2]

meta_df$theta <- as.numeric(
  meta_df$tricyclePosition
)

meta_df$theta <- meta_df$theta %% (2*pi)

if ("total_counts" %in% colnames(meta_df)) {
  meta_df$total_counts <- as.numeric(meta_df$total_counts)
}

write.csv(
  meta_df,
  file = file.path(
    OUTDIR,
    paste0(DATASET_LABEL, "_tricycle_metadata.csv")
  ),
  row.names = FALSE
)

saveRDS(
  sce,
  file = file.path(
    OUTDIR,
    paste0(DATASET_LABEL, "_sce_with_tricycle.rds")
  )
)

message("RGC cells: ", ncol(sce))
message("Genes: ", nrow(sce))



# 9. Tricycle embedding plot

p_embed <- ggplot(
  meta_df,
  aes(x = tricycle_PC1, y = tricycle_PC2, color = theta)
) +
  geom_point(size = 1.2, alpha = 0.85) +
  coord_equal() +
  scale_color_gradientn(
    colours = grDevices::rainbow(100),
    limits = c(0, 2*pi),
    breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
    labels = c("0", "0.5pi", "pi", "1.5pi", "2pi"),
    oob = scales::squish
  ) +
  labs(
    title = paste0(DATASET_LABEL, ": tricycle embedding"),
    x = "Projected PC1",
    y = "Projected PC2",
    color = "theta"
  ) +
  theme_bw(base_size = 13)

ggsave(
  filename = file.path(
    OUTDIR,
    paste0(DATASET_LABEL, "_tricycle_embedding.png")
  ),
  plot = p_embed,
  width = 5.5,
  height = 4.8,
  dpi = 300
)



# 10. Cell-cycle marker QC

marker_groups <- list(
  "S marker genes" = c(
    "Pcna",
    "Mcm5",
    "Mcm6",
    "Mcm7",
    "Rrm2",
    "Tyms"
  ),
  "G2/M marker genes" = c(
    "Top2a",
    "Smc4",
    "Cdk1",
    "Ccnb1",
    "Ube2c",
    "Aurka",
    "Birc5"
  )
)

QC_PDF_WIDTH <- 6.2
QC_PDF_HEIGHT <- 4.8
QC_POINT_SIZE <- 0.75
QC_POINT_ALPHA <- 0.30
QC_LINE_WIDTH <- 0.70

build_marker_group_data <- function(
  gene_vec,
  sce,
  meta_df,
  loess_span = 0.25,
  n_grid = 300
) {
  point_df_list <- list()
  curve_df_list <- list()
  summary_list <- list()
  missing_genes <- character(0)

  for (g in gene_vec) {
    gene_actual <- get_gene_actual(g, sce)

    if (is.na(gene_actual)) {
      missing_genes <- c(missing_genes, g)
      next
    }

    expr_vec <- as.numeric(
      SummarizedExperiment::assay(
        sce,
        "logcounts"
      )[gene_actual, ]
    )

    plot_df <- data.frame(
      theta = meta_df$theta,
      expr = expr_vec,
      gene = gene_actual,
      requested_gene = g,
      stringsAsFactors = FALSE
    )

    plot_df <- plot_df[
      is.finite(plot_df$theta) &
        is.finite(plot_df$expr),
      ,
      drop = FALSE
    ]

    if (nrow(plot_df) < 10) {
      missing_genes <- c(missing_genes, g)
      next
    }

    loess_df <- rbind(
      data.frame(theta = plot_df$theta - 2*pi, y = plot_df$expr),
      data.frame(theta = plot_df$theta,        y = plot_df$expr),
      data.frame(theta = plot_df$theta + 2*pi, y = plot_df$expr)
    )

    fit_obj <- tryCatch(
      loess(
        y ~ theta,
        data = loess_df,
        span = loess_span,
        degree = 2,
        surface = "direct"
      ),
      error = function(e) NULL
    )

    if (is.null(fit_obj)) {
      missing_genes <- c(missing_genes, g)
      next
    }

    pred_df <- data.frame(
      theta = seq(0, 2*pi, length.out = n_grid)
    )

    pred_df$fitted <- as.numeric(
      predict(fit_obj, newdata = pred_df)
    )

    pred_df$gene <- gene_actual
    pred_df$requested_gene <- g

    point_df_list[[gene_actual]] <- plot_df
    curve_df_list[[gene_actual]] <- pred_df

    summary_list[[gene_actual]] <- data.frame(
      requested_gene = g,
      matched_gene = gene_actual,
      n_cells = length(expr_vec),
      n_finite_cells = nrow(plot_df),
      n_detected_cells = sum(expr_vec > 0, na.rm = TRUE),
      detected_fraction = mean(expr_vec > 0, na.rm = TRUE),
      mean_expr = mean(expr_vec, na.rm = TRUE),
      max_expr = max(expr_vec, na.rm = TRUE),
      peak_angle = pred_df$theta[which.max(pred_df$fitted)],
      trough_angle = pred_df$theta[which.min(pred_df$fitted)],
      amplitude = (
        max(pred_df$fitted, na.rm = TRUE) -
          min(pred_df$fitted, na.rm = TRUE)
      ),
      stringsAsFactors = FALSE
    )
  }

  list(
    point_df = if (length(point_df_list) > 0) {
      do.call(rbind, point_df_list)
    } else {
      NULL
    },
    curve_df = if (length(curve_df_list) > 0) {
      do.call(rbind, curve_df_list)
    } else {
      NULL
    },
    summary_df = if (length(summary_list) > 0) {
      do.call(rbind, summary_list)
    } else {
      NULL
    },
    missing_genes = unique(missing_genes)
  )
}


make_marker_group_plot <- function(
  group_title,
  point_df,
  curve_df
) {
  ggplot() +
    geom_point(
      data = point_df,
      aes(x = theta, y = expr, color = gene),
      size = QC_POINT_SIZE,
      alpha = QC_POINT_ALPHA,
      stroke = 0
    ) +
    geom_line(
      data = curve_df,
      aes(x = theta, y = fitted, color = gene),
      linewidth = QC_LINE_WIDTH,
      lineend = "round"
    ) +
    scale_x_continuous(
      breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
      labels = c("0", "0.5pi", "pi", "1.5pi", "2pi"),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    coord_cartesian(xlim = c(0, 2*pi)) +
    labs(
      title = group_title,
      x = "tricycle theta",
      y = "expression",
      color = "Gene"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "plain"),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      panel.grid.major = element_line(linewidth = 0.35, color = "grey88"),
      panel.grid.minor = element_line(linewidth = 0.25, color = "grey93"),
      panel.border = element_rect(linewidth = 0.7, color = "black")
    )
}


qc_plot_list <- list()
qc_summary_list <- list()
qc_missing_list <- list()

for (group_name in names(marker_groups)) {
  message("QC group: ", group_name)

  res <- build_marker_group_data(
    gene_vec = marker_groups[[group_name]],
    sce = sce,
    meta_df = meta_df,
    loess_span = LOESS_SPAN,
    n_grid = N_GRID
  )

  if (is.null(res$point_df) || is.null(res$curve_df)) {
    next
  }

  p <- make_marker_group_plot(
    group_title = group_name,
    point_df = res$point_df,
    curve_df = res$curve_df
  )

  qc_plot_list[[group_name]] <- p

  ggsave(
    filename = file.path(
      QC_DIR,
      paste0("RGC_", safe_filename(group_name), ".pdf")
    ),
    plot = p,
    device = "pdf",
    width = QC_PDF_WIDTH,
    height = QC_PDF_HEIGHT,
    units = "in"
  )

  if (!is.null(res$summary_df)) {
    tmp_summary <- res$summary_df
    tmp_summary$group <- group_name
    qc_summary_list[[group_name]] <- tmp_summary
  }

  if (length(res$missing_genes) > 0) {
    qc_missing_list[[group_name]] <- data.frame(
      group = group_name,
      missing_gene = res$missing_genes,
      stringsAsFactors = FALSE
    )
  }
}

if (length(qc_plot_list) > 0) {
  grDevices::pdf(
    file = file.path(
      QC_DIR,
      "RGC_cell_cycle_marker_groups_all.pdf"
    ),
    width = QC_PDF_WIDTH,
    height = QC_PDF_HEIGHT,
    onefile = TRUE
  )

  for (nm in names(qc_plot_list)) {
    print(qc_plot_list[[nm]])
  }

  dev.off()
}

if (length(qc_summary_list) > 0) {
  write.csv(
    do.call(rbind, qc_summary_list),
    file = file.path(
      QC_DIR,
      "RGC_cell_cycle_marker_groups_summary.csv"
    ),
    row.names = FALSE
  )
}

if (length(qc_missing_list) > 0) {
  write.csv(
    do.call(rbind, qc_missing_list),
    file = file.path(
      QC_DIR,
      "RGC_cell_cycle_marker_groups_missing_genes.csv"
    ),
    row.names = FALSE
  )
}



# 11. Read input gene set

if (!file.exists(GENESET_CSV)) {
  stop("Gene-set file was not found: ", GENESET_CSV)
}

message("Reading input gene set: ", GENESET_CSV)

gene_set_df <- read.csv(
  GENESET_CSV,
  stringsAsFactors = FALSE
)

if (!"gene_symbol" %in% colnames(gene_set_df)) {
  stop("gene_symbol column was not found in GENESET_CSV.")
}

GENES_OF_INTEREST <- unique(
  as.character(gene_set_df$gene_symbol)
)

GENES_OF_INTEREST <- GENES_OF_INTEREST[
  !is.na(GENES_OF_INTEREST) &
    GENES_OF_INTEREST != ""
]

if (length(GENES_OF_INTEREST) == 0) {
  stop("No valid gene symbols were found in GENESET_CSV.")
}

message(
  "Input gene set: ",
  length(GENES_OF_INTEREST),
  " genes"
)



# 12. Gene lookup

gene_lookup <- setNames(
  rownames(sce),
  tolower(rownames(sce))
)

requested_genes_unique <- unique(
  GENES_OF_INTEREST
)

requested_lower <- tolower(
  requested_genes_unique
)

present_idx <- requested_lower %in% names(gene_lookup)

present_requested <- requested_genes_unique[present_idx]
missing_requested <- requested_genes_unique[!present_idx]

present_actual <- unname(
  gene_lookup[tolower(present_requested)]
)

write.csv(
  data.frame(
    requested_gene = present_requested,
    matched_gene = present_actual
  ),
  file = file.path(
    OUTDIR,
    paste0(DATASET_LABEL, "_present_genes.csv")
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    missing_gene = missing_requested
  ),
  file = file.path(
    OUTDIR,
    paste0(DATASET_LABEL, "_missing_genes.csv")
  ),
  row.names = FALSE
)

message("Present genes: ", length(present_actual))
message("Missing genes: ", length(missing_requested))



# 13. Periodic LOESS permutation test

theta_vec <- meta_df$theta
theta_grid <- seq(0, 2*pi, length.out = N_GRID)
expr_mat <- SummarizedExperiment::assay(sce, "logcounts")

results_list <- vector("list", length(present_actual))
names(results_list) <- present_actual

curve_list <- vector("list", length(present_actual))
names(curve_list) <- present_actual

perm_list <- vector("list", length(present_actual))
names(perm_list) <- present_actual

message("Starting permutation test...")
message("Number of genes: ", length(present_actual))
message("Permutation per gene: ", N_PERM)
message("Loess span: ", LOESS_SPAN)
message("Theta grid: ", N_GRID)

for (i in seq_along(present_actual)) {
  gene_req <- present_requested[i]
  gene_act <- present_actual[i]

  message(
    "[",
    i,
    "/",
    length(present_actual),
    "] ",
    gene_act
  )

  y <- as.numeric(expr_mat[gene_act, ])

  f_obs <- fit_periodic_loess(
    theta_vec,
    y,
    theta_grid,
    span = LOESS_SPAN
  )

  obs_stats <- calc_curve_stats(
    theta_vec,
    y,
    theta_grid,
    f_obs
  )

  T_obs <- obs_stats$T
  T_perm <- rep(NA_real_, N_PERM)

  for (b in seq_len(N_PERM)) {
    theta_perm <- sample(
      theta_vec,
      length(theta_vec),
      replace = FALSE
    )

    f_perm <- fit_periodic_loess(
      theta_perm,
      y,
      theta_grid,
      span = LOESS_SPAN
    )

    T_perm[b] <- mean(
      (
        f_perm -
          mean(f_perm, na.rm = TRUE)
      )^2,
      na.rm = TRUE
    )
  }

  p_val <- (
    1 +
      sum(T_perm >= T_obs, na.rm = TRUE)
  ) / (
    1 +
      sum(is.finite(T_perm))
  )

  det_rate <- mean(y > 0, na.rm = TRUE)
  mean_expr <- mean(y, na.rm = TRUE)
  max_expr <- max(y, na.rm = TRUE)

  results_list[[i]] <- data.frame(
    requested_gene = gene_req,
    matched_gene = gene_act,
    detection_rate = det_rate,
    mean_expr = mean_expr,
    max_expr = max_expr,
    T = T_obs,
    p_value = p_val,
    peak_angle = obs_stats$peak_angle,
    trough_angle = obs_stats$trough_angle,
    amplitude = obs_stats$amplitude,
    maximum_amplitude = obs_stats$amplitude,
    standardized_amplitude = obs_stats$standardized_amplitude,
    pseudo_R2 = obs_stats$pseudo_R2,
    SSE_loess = obs_stats$SSE_loess,
    SSE_flat = obs_stats$SSE_flat,
    stringsAsFactors = FALSE
  )

  curve_list[[i]] <- data.frame(
    gene = gene_act,
    theta = theta_grid,
    fitted = f_obs
  )

  perm_list[[i]] <- data.frame(
    gene = gene_act,
    permutation = seq_len(N_PERM),
    T_perm = T_perm
  )
}



# 14. BH-FDR and result tables

results_df <- do.call(rbind, results_list)

results_df$FDR <- p.adjust(
  results_df$p_value,
  method = "BH"
)

results_df$significant_FDR_0.05 <- (
  results_df$FDR < FDR_CUTOFF
)

results_df <- results_df[
  order(
    results_df$FDR,
    results_df$p_value,
    -results_df$T
  ),
]

expected_result_columns <- c(
  "requested_gene",
  "matched_gene",
  "detection_rate",
  "mean_expr",
  "max_expr",
  "T",
  "p_value",
  "peak_angle",
  "trough_angle",
  "amplitude",
  "maximum_amplitude",
  "standardized_amplitude",
  "pseudo_R2",
  "SSE_loess",
  "SSE_flat",
  "FDR",
  "significant_FDR_0.05"
)

if (!identical(colnames(results_df), expected_result_columns)) {
  stop(
    "Result CSV schema differs from the original analysis."
  )
}

write.csv(
  results_df,
  file = file.path(
    OUTDIR,
    paste0(
      DATASET_LABEL,
      "_periodic_loess_permutation_results_all.csv"
    )
  ),
  row.names = FALSE
)

sig_results_df <- results_df[
  results_df$significant_FDR_0.05,
]

write.csv(
  sig_results_df,
  file = file.path(
    OUTDIR,
    paste0(
      DATASET_LABEL,
      "_periodic_loess_permutation_results_significant_FDR005.csv"
    )
  ),
  row.names = FALSE
)

curve_df <- do.call(rbind, curve_list)

if (!identical(
  colnames(curve_df),
  c("gene", "theta", "fitted")
)) {
  stop(
    "Fitted-curve CSV schema differs from the original analysis."
  )
}

write.csv(
  curve_df,
  file = file.path(
    OUTDIR,
    paste0(
      DATASET_LABEL,
      "_periodic_loess_fitted_curves_all.csv"
    )
  ),
  row.names = FALSE
)

perm_df <- do.call(rbind, perm_list)

if (!identical(
  colnames(perm_df),
  c("gene", "permutation", "T_perm")
)) {
  stop(
    "Permutation CSV schema differs from the original analysis."
  )
}

write.csv(
  perm_df,
  file = file.path(
    OUTDIR,
    paste0(
      DATASET_LABEL,
      "_permutation_T_values_all.csv"
    )
  ),
  row.names = FALSE
)

saveRDS(
  list(
    results_df = results_df,
    curve_list = curve_list,
    perm_list = perm_list,
    theta_grid = theta_grid,
    settings = list(
      LOESS_SPAN = LOESS_SPAN,
      N_GRID = N_GRID,
      N_PERM = N_PERM,
      FDR_CUTOFF = FDR_CUTOFF,
      RANDOM_SEED = RANDOM_SEED
    )
  ),
  file = file.path(
    OUTDIR,
    paste0(
      DATASET_LABEL,
      "_periodic_test_full_results.rds"
    )
  )
)

message(
  "Significant genes at FDR < ",
  FDR_CUTOFF,
  ": ",
  nrow(sig_results_df)
)



# 15. Original single-panel plots and permutation histograms

pdf(
  file = file.path(
    OUTDIR,
    paste0(DATASET_LABEL, "_all_genes_vs_theta.pdf")
  ),
  width = 6.0,
  height = 4.5
)

for (i in seq_along(present_actual)) {
  gene_act <- present_actual[i]

  res_row <- results_df[
    results_df$matched_gene == gene_act,
  ][1, ]

  y <- as.numeric(expr_mat[gene_act, ])

  plot_df <- data.frame(
    theta = theta_vec,
    expr = y
  )

  plot_df <- plot_df[
    is.finite(plot_df$theta) &
      is.finite(plot_df$expr),
    ,
    drop = FALSE
  ]

  pred_df <- curve_list[[gene_act]]
  gene_file <- safe_filename(gene_act)

  p_gene <- make_gene_plot(
    plot_df = plot_df,
    pred_df = pred_df,
    gene_name = gene_act,
    result_row = res_row,
    outfile = file.path(
      PLOT_DIR,
      paste0(
        DATASET_LABEL,
        "_",
        gene_file,
        "_vs_theta.png"
      )
    )
  )

  print(p_gene)

  make_hist_plot(
    T_perm = perm_list[[gene_act]]$T_perm,
    T_obs = res_row$T,
    gene_name = gene_act,
    result_row = res_row,
    outfile = file.path(
      HIST_DIR,
      paste0(
        DATASET_LABEL,
        "_",
        gene_file,
        "_T_histogram.png"
      )
    )
  )
}

dev.off()


if (nrow(sig_results_df) > 0) {
  pdf(
    file = file.path(
      OUTDIR,
      paste0(
        DATASET_LABEL,
        "_significant_genes_vs_theta_FDR005.pdf"
      )
    ),
    width = 6.0,
    height = 4.5
  )

  for (i in seq_len(nrow(sig_results_df))) {
    gene_act <- sig_results_df$matched_gene[i]
    res_row <- sig_results_df[i, ]

    y <- as.numeric(expr_mat[gene_act, ])

    plot_df <- data.frame(
      theta = theta_vec,
      expr = y
    )

    plot_df <- plot_df[
      is.finite(plot_df$theta) &
        is.finite(plot_df$expr),
      ,
      drop = FALSE
    ]

    pred_df <- curve_list[[gene_act]]
    gene_file <- safe_filename(gene_act)

    p_gene <- make_gene_plot(
      plot_df = plot_df,
      pred_df = pred_df,
      gene_name = gene_act,
      result_row = res_row,
      outfile = file.path(
        SIG_PLOT_DIR,
        paste0(
          DATASET_LABEL,
          "_",
          gene_file,
          "_vs_theta.png"
        )
      )
    )

    print(p_gene)

    make_hist_plot(
      T_perm = perm_list[[gene_act]]$T_perm,
      T_obs = res_row$T,
      gene_name = gene_act,
      result_row = res_row,
      outfile = file.path(
        SIG_HIST_DIR,
        paste0(
          DATASET_LABEL,
          "_",
          gene_file,
          "_T_histogram.png"
        )
      )
    )
  }

  dev.off()
}



# 16. Two-panel publication figures

TWO_PANEL_FIG_WIDTH <- 7
TWO_PANEL_FIG_HEIGHT <- 7.2

FIXED_PANEL_WIDTH <- 4.80
FIXED_PANEL_HEIGHT <- 2.05
PANEL_PADDING <- 0.20

TOP_POINT_SIZE <- 3.00
TOP_POINT_ALPHA <- 1.00
TOP_LINE_WIDTH <- 1.80

BOTTOM_LINE_WIDTH <- 1.80

HIST_FILL <- "#727171"
HIST_ALPHA <- 0.70

BASE_AXIS_TEXT_SIZE <- 12
AXIS_TEXT_SIZE <- BASE_AXIS_TEXT_SIZE * 1.9

BASE_TITLE_SIZE <- 20
TITLE_TEXT_SIZE <- BASE_TITLE_SIZE * 2.3

TWO_PANEL_SUFFIX <- paste0(
  "_two_panel_tricycle_theta_",
  "axis_05_equal_max4_grayhist_",
  "fixed_panel_v4_dummy_right_axis"
)

theta_color_scale <- scale_color_gradientn(
  colours = grDevices::rainbow(100),
  limits = c(0, 2*pi),
  breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
  labels = c("0", "0.5pi", "pi", "1.5pi", "2pi"),
  oob = scales::squish,
  guide = "none"
)

two_panel_summary_list <- list()
two_panel_missing_genes <- character(0)

for (i in seq_along(present_actual)) {
  gene_req <- present_requested[i]
  gene_act <- present_actual[i]

  message("Two-panel plot: ", gene_act)

  expr_vec <- as.numeric(expr_mat[gene_act, ])

  plot_df <- data.frame(
    theta = theta_vec,
    expr = expr_vec,
    stringsAsFactors = FALSE
  )

  plot_df <- plot_df[
    is.finite(plot_df$theta) &
      is.finite(plot_df$expr),
    ,
    drop = FALSE
  ]

  if (nrow(plot_df) < 10) {
    two_panel_missing_genes <- c(
      two_panel_missing_genes,
      gene_req
    )
    next
  }

  pred_df <- curve_list[[gene_act]]

  if (
    is.null(pred_df) ||
      all(!is.finite(pred_df$fitted))
  ) {
    two_panel_missing_genes <- c(
      two_panel_missing_genes,
      gene_req
    )
    next
  }

  A <- mean(pred_df$fitted, na.rm = TRUE)

  if (!is.finite(A) || A == 0) {
    two_panel_missing_genes <- c(
      two_panel_missing_genes,
      gene_req
    )
    next
  }

  pred_df$fitted_div_A <- pred_df$fitted / A

  hist_theta <- plot_df$theta[
    is.finite(plot_df$theta) &
      is.finite(plot_df$expr) &
      plot_df$expr > 0
  ]

  hist_obj <- hist(
    hist_theta,
    breaks = theta_grid,
    plot = FALSE,
    include.lowest = TRUE,
    right = FALSE
  )

  hist_df <- data.frame(
    xmin = hist_obj$breaks[-length(hist_obj$breaks)],
    xmax = hist_obj$breaks[-1],
    xmid = (
      hist_obj$breaks[-length(hist_obj$breaks)] +
        hist_obj$breaks[-1]
    ) / 2,
    count = hist_obj$counts,
    stringsAsFactors = FALSE
  )

  top_axis <- make_axis_05_equal_max4(
    max(
      c(plot_df$expr, pred_df$fitted),
      na.rm = TRUE
    )
  )

  bottom_curve_axis <- make_axis_05_equal_max4(
    max(pred_df$fitted_div_A, na.rm = TRUE)
  )

  hist_axis <- make_axis_05_equal_max4(
    max(hist_df$count, na.rm = TRUE)
  )

  hist_scale_factor <- (
    bottom_curve_axis$axis_max /
      hist_axis$axis_max
  )

  hist_df$count_scaled <- (
    hist_df$count *
      hist_scale_factor
  )

  gene_title_expr <- bquote(
    italic(.(gene_act))
  )

  p_top <- ggplot() +
    geom_line(
      data = pred_df,
      aes(x = theta, y = fitted),
      inherit.aes = FALSE,
      color = "black",
      linewidth = TOP_LINE_WIDTH,
      lineend = "round"
    ) +
    geom_point(
      data = plot_df,
      aes(x = theta, y = expr, color = theta),
      size = TOP_POINT_SIZE,
      alpha = TOP_POINT_ALPHA,
      stroke = 0
    ) +
    scale_x_continuous(
      breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
      labels = NULL,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      breaks = top_axis$breaks,
      labels = top_axis$labels,
      sec.axis = sec_axis(
        trans = ~ .,
        name = NULL,
        breaks = top_axis$breaks,
        labels = rep("", length(top_axis$breaks))
      )
    ) +
    coord_cartesian(
      xlim = c(0, 2*pi),
      ylim = c(top_axis$axis_min, top_axis$axis_max)
    ) +
    theta_color_scale +
    labs(
      title = gene_title_expr,
      x = NULL,
      y = "expression"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = TITLE_TEXT_SIZE,
        face = "plain",
        color = "black"
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_text(
        size = AXIS_TEXT_SIZE,
        color = "black"
      ),
      axis.title.y.right = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(
        size = AXIS_TEXT_SIZE,
        color = "black"
      ),
      axis.text.y.right = element_text(
        size = AXIS_TEXT_SIZE,
        color = "white"
      ),
      axis.ticks.y = element_line(
        color = "black",
        linewidth = 0.5
      ),
      axis.ticks.y.right = element_line(
        color = "white",
        linewidth = 0.5
      ),
      panel.grid.major = element_line(
        linewidth = 0.35,
        color = "grey88"
      ),
      panel.grid.minor = element_line(
        linewidth = 0.25,
        color = "grey93"
      ),
      panel.border = element_rect(
        linewidth = 0.7,
        color = "black"
      ),
      plot.margin = margin(8, 8, 4, 8)
    )

  p_bottom <- ggplot() +
    geom_rect(
      data = hist_df,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = 0,
        ymax = count_scaled
      ),
      fill = HIST_FILL,
      alpha = HIST_ALPHA,
      color = NA
    ) +
    geom_line(
      data = pred_df,
      aes(x = theta, y = fitted_div_A),
      color = "black",
      linewidth = BOTTOM_LINE_WIDTH,
      lineend = "round"
    ) +
    scale_x_continuous(
      breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
      labels = c("0", "0.5pi", "pi", "1.5pi", "2pi"),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      name = "relative exp.",
      breaks = bottom_curve_axis$breaks,
      labels = bottom_curve_axis$labels,
      sec.axis = sec_axis(
        trans = ~ . / hist_scale_factor,
        name = "cell count",
        breaks = hist_axis$breaks,
        labels = hist_axis$labels
      )
    ) +
    coord_cartesian(
      xlim = c(0, 2*pi),
      ylim = c(
        bottom_curve_axis$axis_min,
        bottom_curve_axis$axis_max
      )
    ) +
    labs(
      x = "tricycle theta"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.title.x = element_text(
        size = AXIS_TEXT_SIZE,
        color = "black"
      ),
      axis.title.y = element_text(
        size = AXIS_TEXT_SIZE,
        color = "black"
      ),
      axis.title.y.right = element_text(
        size = AXIS_TEXT_SIZE,
        color = HIST_FILL
      ),
      axis.text.x = element_text(
        size = AXIS_TEXT_SIZE,
        color = "black"
      ),
      axis.text.y = element_text(
        size = AXIS_TEXT_SIZE,
        color = "black"
      ),
      axis.text.y.right = element_text(
        size = AXIS_TEXT_SIZE,
        color = HIST_FILL
      ),
      axis.ticks.x = element_line(
        color = "black",
        linewidth = 0.5
      ),
      axis.ticks.y = element_line(
        color = "black",
        linewidth = 0.5
      ),
      axis.ticks.y.right = element_line(
        color = "black",
        linewidth = 0.5
      ),
      panel.grid.major = element_line(
        linewidth = 0.35,
        color = "grey88"
      ),
      panel.grid.minor = element_line(
        linewidth = 0.25,
        color = "grey93"
      ),
      panel.border = element_rect(
        linewidth = 0.7,
        color = "black"
      ),
      plot.margin = margin(4, 8, 8, 8)
    )

  p_top_fixed <- egg::set_panel_size(
    p_top,
    width = grid::unit(FIXED_PANEL_WIDTH, "in"),
    height = grid::unit(FIXED_PANEL_HEIGHT, "in")
  )

  p_bottom_fixed <- egg::set_panel_size(
    p_bottom,
    width = grid::unit(FIXED_PANEL_WIDTH, "in"),
    height = grid::unit(FIXED_PANEL_HEIGHT, "in")
  )

  p_combined_fixed <- gridExtra::arrangeGrob(
    p_top_fixed,
    p_bottom_fixed,
    ncol = 1,
    padding = grid::unit(PANEL_PADDING, "in")
  )

  gene_file <- safe_filename(gene_act)

  pdf_file <- file.path(
    TWO_PANEL_DIR,
    paste0(gene_file, TWO_PANEL_SUFFIX, ".pdf")
  )

  png_file <- file.path(
    TWO_PANEL_DIR,
    paste0(gene_file, TWO_PANEL_SUFFIX, ".png")
  )

  summary_file <- file.path(
    TWO_PANEL_DIR,
    paste0(gene_file, TWO_PANEL_SUFFIX, "_summary.csv")
  )

  grDevices::pdf(
    file = pdf_file,
    width = TWO_PANEL_FIG_WIDTH,
    height = TWO_PANEL_FIG_HEIGHT,
    onefile = FALSE,
    useDingbats = FALSE
  )

  grid::grid.newpage()
  grid::grid.draw(p_combined_fixed)
  grDevices::dev.off()

  png(
    filename = png_file,
    width = TWO_PANEL_FIG_WIDTH,
    height = TWO_PANEL_FIG_HEIGHT,
    units = "in",
    res = 300
  )

  grid::grid.newpage()
  grid::grid.draw(p_combined_fixed)
  dev.off()

  summary_df <- data.frame(
    requested_gene = gene_req,
    matched_gene = gene_act,
    n_cells_total = length(expr_vec),
    n_cells_finite_for_plot = nrow(plot_df),
    n_cells_expr_positive_for_hist = sum(
      plot_df$expr > 0,
      na.rm = TRUE
    ),
    loess_span = LOESS_SPAN,
    n_theta_grid_points = length(theta_grid),
    n_hist_intervals = length(theta_grid) - 1,
    A_mean_of_300_fitted_values = A,
    mean_expr = mean(plot_df$expr, na.rm = TRUE),
    max_expr = max(plot_df$expr, na.rm = TRUE),
    top_y_axis_max = top_axis$axis_max,
    top_y_breaks = paste(
      top_axis$labels,
      collapse = ";"
    ),
    bottom_curve_y_axis_max = (
      bottom_curve_axis$axis_max
    ),
    bottom_curve_y_breaks = paste(
      bottom_curve_axis$labels,
      collapse = ";"
    ),
    hist_y_axis_max = hist_axis$axis_max,
    hist_y_breaks = paste(
      hist_axis$labels,
      collapse = ";"
    ),
    histogram_fill = HIST_FILL,
    histogram_alpha = HIST_ALPHA,
    fixed_panel_width_in = FIXED_PANEL_WIDTH,
    fixed_panel_height_in = FIXED_PANEL_HEIGHT,
    panel_padding_in = PANEL_PADDING,
    title_text_size = TITLE_TEXT_SIZE,
    dummy_right_axis_in_top_panel = TRUE,
    peak_angle = pred_df$theta[
      which.max(pred_df$fitted)
    ],
    trough_angle = pred_df$theta[
      which.min(pred_df$fitted)
    ],
    amplitude = (
      max(pred_df$fitted, na.rm = TRUE) -
        min(pred_df$fitted, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    summary_df,
    file = summary_file,
    row.names = FALSE
  )

  two_panel_summary_list[[gene_act]] <- summary_df

  res_row <- results_df[
    results_df$matched_gene == gene_act,
  ][1, ]

  if (isTRUE(res_row$significant_FDR_0.05)) {
    file.copy(
      pdf_file,
      file.path(
        SIG_TWO_PANEL_DIR,
        basename(pdf_file)
      ),
      overwrite = TRUE
    )

    file.copy(
      png_file,
      file.path(
        SIG_TWO_PANEL_DIR,
        basename(png_file)
      ),
      overwrite = TRUE
    )

    file.copy(
      summary_file,
      file.path(
        SIG_TWO_PANEL_DIR,
        basename(summary_file)
      ),
      overwrite = TRUE
    )
  }
}

if (length(two_panel_summary_list) > 0) {
  two_panel_summary_df <- do.call(
    rbind,
    two_panel_summary_list
  )

  write.csv(
    two_panel_summary_df,
    file = file.path(
      TWO_PANEL_DIR,
      paste0(
        DATASET_LABEL,
        "_two_panel_tricycle_theta_summary_all.csv"
      )
    ),
    row.names = FALSE
  )

  if (nrow(sig_results_df) > 0) {
    sig_two_panel_summary_df <- (
      two_panel_summary_df[
        two_panel_summary_df$matched_gene %in%
          sig_results_df$matched_gene,
        ,
        drop = FALSE
      ]
    )

    write.csv(
      sig_two_panel_summary_df,
      file = file.path(
        SIG_TWO_PANEL_DIR,
        paste0(
          DATASET_LABEL,
          "_two_panel_tricycle_theta_summary_significant_FDR005.csv"
        )
      ),
      row.names = FALSE
    )
  }
}

if (length(two_panel_missing_genes) > 0) {
  write.csv(
    data.frame(
      missing_gene = unique(two_panel_missing_genes)
    ),
    file = file.path(
      TWO_PANEL_DIR,
      paste0(
        DATASET_LABEL,
        "_two_panel_tricycle_theta_missing_genes.csv"
      )
    ),
    row.names = FALSE
  )
}



# 17. Summary plots

p_volcano <- ggplot(
  results_df,
  aes(
    x = T,
    y = -log10(FDR),
    color = significant_FDR_0.05
  )
) +
  geom_point(alpha = 0.75, size = 1.8) +
  scale_color_manual(
    values = c(
      "FALSE" = "grey60",
      "TRUE" = "red"
    )
  ) +
  labs(
    title = paste0(
      DATASET_LABEL,
      ": periodic loess permutation test"
    ),
    x = "T = mean((f(theta_grid) - mean(f))^2)",
    y = "-log10(FDR)",
    color = "FDR < 0.05"
  ) +
  theme_bw(base_size = 13)

ggsave(
  filename = file.path(
    OUTDIR,
    paste0(
      DATASET_LABEL,
      "_T_vs_minuslog10FDR.png"
    )
  ),
  plot = p_volcano,
  width = 6.0,
  height = 4.8,
  dpi = 300
)


p_amp <- ggplot(
  results_df,
  aes(
    x = peak_angle,
    y = amplitude,
    color = significant_FDR_0.05
  )
) +
  geom_point(alpha = 0.75, size = 1.8) +
  scale_color_manual(
    values = c(
      "FALSE" = "grey60",
      "TRUE" = "red"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 2*pi),
    breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
    labels = c("0", "0.5pi", "pi", "1.5pi", "2pi")
  ) +
  labs(
    title = paste0(
      DATASET_LABEL,
      ": peak angle and amplitude"
    ),
    x = "peak angle",
    y = "maximum amplitude = max(f_obs) - min(f_obs)",
    color = "FDR < 0.05"
  ) +
  theme_bw(base_size = 13)

ggsave(
  filename = file.path(
    OUTDIR,
    paste0(
      DATASET_LABEL,
      "_peak_angle_vs_amplitude.png"
    )
  ),
  plot = p_amp,
  width = 6.0,
  height = 4.8,
  dpi = 300
)



# 18. Session information and summary

capture.output(
  sessionInfo(),
  file = file.path(
    OUTDIR,
    paste0(DATASET_LABEL, "_sessionInfo.txt")
  )
)

summary_df <- data.frame(
  n_cells = ncol(sce),
  n_genes_total = nrow(sce),
  n_genes_in_csv = length(GENES_OF_INTEREST),
  n_genes_present_in_sce = length(present_actual),
  n_genes_missing_in_sce = length(missing_requested),
  n_grid = N_GRID,
  n_permutation = N_PERM,
  loess_span = LOESS_SPAN,
  fdr_cutoff = FDR_CUTOFF,
  n_significant_FDR005 = nrow(sig_results_df),
  theta_min = min(meta_df$theta, na.rm = TRUE),
  theta_max = max(meta_df$theta, na.rm = TRUE)
)

write.csv(
  summary_df,
  file = file.path(
    OUTDIR,
    paste0(DATASET_LABEL, "_summary.csv")
  ),
  row.names = FALSE
)



# 19. Finish

message("")
message("Done.")
message("Output directory: ", OUTDIR)
message("Input gene set: ", length(GENES_OF_INTEREST))
message("RGC cells: ", ncol(sce))
message("LOESS span: ", LOESS_SPAN)
message("Theta grid points: ", N_GRID)
message("Permutations per gene: ", N_PERM)
message("FDR cutoff: ", FDR_CUTOFF)

message("")
message("Main outputs:")
message(
  "  ",
  file.path(
    OUTDIR,
    paste0(
      DATASET_LABEL,
      "_periodic_loess_permutation_results_all.csv"
    )
  )
)
message(
  "  ",
  file.path(
    OUTDIR,
    paste0(
      DATASET_LABEL,
      "_periodic_loess_permutation_results_significant_FDR005.csv"
    )
  )
)
message(
  "  ",
  file.path(
    OUTDIR,
    paste0(
      DATASET_LABEL,
      "_periodic_loess_fitted_curves_all.csv"
    )
  )
)
message(
  "  ",
  file.path(
    OUTDIR,
    paste0(
      DATASET_LABEL,
      "_permutation_T_values_all.csv"
    )
  )
)
message(
  "  ",
  file.path(
    OUTDIR,
    paste0(DATASET_LABEL, "_tricycle_metadata.csv")
  )
)
