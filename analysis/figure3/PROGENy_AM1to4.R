# PROGENy AM1to4
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(Matrix); library(pheatmap)
})

output_dir <- am_data("HTAN/am/PROGENy")
input_rds <- am_data("HTAN/am/object/Rename_am.rds")
model_file <- file.path(output_dir, "model_human_full.rda")
top_genes <- 500
force_recompute <- FALSE
am_order <- c("AM-1", "AM-2", "AM-3", "AM-4")
pathway_order <- c(
  "NFkB", "JAK-STAT", "TNFa", "EGFR", "MAPK", "PI3K", "p53",
  "Hypoxia", "VEGF", "TGFb", "WNT", "Trail", "Androgen", "Estrogen"
)
font_size <- 12
font_family <- "sans"
figure_width <- 8
figure_height <- 3.4
column_angle <- "90"
heat_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
heat_limits <- c(-2, 2)

row_standardize <- TRUE
show_numbers <- FALSE
column_gaps <- NULL
export_cell_csv <- TRUE
model_commit <- "cad6be0514c3248b9465e48f1cfd2f6a4c3dfb6f"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
f <- function(n) file.path(output_dir, n)
started <- Sys.time()
log_phase <- function(s) {
  line <- paste(format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), s,
    round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2),
    sep = "\t"
  )
  message(line); cat(line, "\n", file = f("PROGENy_run.log"), append = TRUE)
}
stopifnot(file.exists(input_rds), file.exists(model_file))
signature <- list(
  input = normalizePath(input_rds), input_size = file.info(input_rds)$size,
  input_mtime = as.numeric(file.info(input_rds)$mtime),
  model_md5 = unname(tools::md5sum(model_file)), top = top_genes,
  algorithm = "official_weighted_sum_all_original_cells_column_scale_v1"
)
cache_file <- f("PROGENy_scores_cache.rds")
cached <- if (!force_recompute && file.exists(cache_file)) readRDS(cache_file) else NULL
if (is.null(cached) || !identical(cached$signature, signature)) {
  log_phase("reading_original_Seurat_object")
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("首次计算需要Seurat包")
  obj <- readRDS(input_rds)
  expr <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")
  metadata <- obj@meta.data
  stopifnot(identical(rownames(metadata), colnames(expr)), "cell_type" %in% names(metadata))
  values <- if (inherits(expr, "sparseMatrix")) expr@x else as.vector(expr)
  stopifnot(all(is.finite(values)), all(values >= 0))
  audit <- data.frame(
    input = input_rds, assay = "RNA", layer = "data", n_genes = nrow(expr),
    n_cells = ncol(expr), minimum = min(0, values), maximum = max(values),
    duplicated_genes_removed = sum(duplicated(rownames(expr)))
  )
  write.csv(audit, f("PROGENy_input_audit.csv"), row.names = FALSE)
  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  rm(values)
  log_phase("input_checked_nonnegative_RNA_data")

  env <- new.env(parent = emptyenv()); load(model_file, envir = env)
  full <- env$model_human_full
  stopifnot(all(c("gene", "pathway", "weight", "p.value") %in% names(full)))
  model_long <- do.call(rbind, lapply(split(full, full$pathway), function(d) {
    cutoff <- sort(d$p.value)[min(top_genes, sum(is.finite(d$p.value)))]
    d[is.finite(d$p.value) & d$p.value <= cutoff, ]
  }))
  rownames(model_long) <- NULL
  stopifnot(!anyDuplicated(paste(model_long$gene, model_long$pathway, sep = "|")))
  model <- matrix(0,
    nrow = length(unique(model_long$gene)), ncol = length(unique(model_long$pathway)),
    dimnames = list(sort(unique(model_long$gene)), sort(unique(model_long$pathway)))
  )
  model[cbind(match(model_long$gene, rownames(model)), match(model_long$pathway, colnames(model)))] <- model_long$weight
  common <- intersect(rownames(expr), rownames(model))
  stopifnot(length(common) > 100)
  model_long$present_in_expression <- model_long$gene %in% common
  write.csv(model_long, f("PROGENy_model_top500_genes.csv"), row.names = FALSE)
  write.csv(model, f("PROGENy_model_top500_weights.csv"))
  coverage <- data.frame(
    pathway = colnames(model), model_nonzero_genes = colSums(model != 0),
    matched_nonzero_genes = colSums(model[common, , drop = FALSE] != 0)
  )
  write.csv(coverage, f("PROGENy_model_gene_coverage.csv"), row.names = FALSE)

  ref <- dplyr::group_by(full, pathway)
  ref <- dplyr::top_n(ref, top_genes, wt = -p.value)
  ref <- dplyr::ungroup(ref)
  ref <- dplyr::select(ref, -p.value)
  ref <- tidyr::spread(ref, pathway, weight, fill = 0)
  ref <- data.frame(ref, check.names = FALSE)
  rownames(ref) <- ref$gene; ref$gene <- NULL
  stopifnot(isTRUE(all.equal(model, as.matrix(ref)[rownames(model), colnames(model)], check.attributes = TRUE)))
  log_phase(paste0("model_checked_", length(common), "_matched_genes"))

  input_model <- expr[common, , drop = FALSE]
  raw_scores <- as.matrix(Matrix::t(input_model) %*% model[common, , drop = FALSE])

  idx <- seq_len(min(100, ncol(expr)))
  reference <- t(as.matrix(input_model[, idx, drop = FALSE])) %*% model[common, , drop = FALSE]
  maximum_error <- max(abs(reference - raw_scores[idx, , drop = FALSE]))
  stopifnot(maximum_error < 1e-8)
  scores <- scale(raw_scores, center = TRUE, scale = TRUE)
  stopifnot(all(is.finite(scores)), identical(rownames(scores), rownames(metadata)))

  metadata$AM_group <- sub("^(AM-[0-9]+).*", "\\1", metadata$cell_type)
  cached <- list(
    signature = signature, metadata = metadata, raw_scores = raw_scores,
    scaled_scores = scores, model_commit = model_commit,
    matched_genes = common, validation_maximum_error = maximum_error
  )
  saveRDS(cached, cache_file, compress = "gzip")
  write.csv(
    data.frame(
      check = c("model_matches_official_getModel", "sparse_vs_dense_100_cells", "all_scores_finite"),
      status = "PASS", maximum_error = c(0, maximum_error, 0)
    ),
    f("PROGENy_validation.csv"),
    row.names = FALSE
  )
  rm(obj, expr, input_model); gc()
  log_phase("scores_computed_and_cached")
} else log_phase("loaded_score_cache")

metadata <- cached$metadata
scores <- cached$scaled_scores
keep <- metadata$AM_group %in% am_order
stopifnot(all(am_order %in% metadata$AM_group), all(pathway_order %in% colnames(scores)))
meta_cols <- intersect(c("sample", "cell_type", "AM_group", "Dataset2", "Group"), names(metadata))
cell_data <- cbind(
  data.frame(cell_id = rownames(metadata)[keep]), metadata[keep, meta_cols, drop = FALSE],
  as.data.frame(scores[keep, , drop = FALSE], check.names = FALSE)
)
if (export_cell_csv && (!file.exists(f("PROGENy_AM1to4_cell_scores.csv.gz")) || force_recompute || is.null(cached$export_done))) {
  con <- gzfile(f("PROGENy_AM1to4_cell_scores.csv.gz"), "wt")
  write.csv(cell_data, con, row.names = FALSE); close(con)
}
means <- do.call(rbind, lapply(am_order, function(am) colMeans(scores[metadata$AM_group == am, , drop = FALSE])))
rownames(means) <- am_order
means <- means[, pathway_order, drop = FALSE]
write.csv(means, f("PROGENy_AM1to4_mean_scores.csv"))
write.csv(data.frame(AM_group = am_order, n_cells = vapply(am_order, function(a) sum(metadata$AM_group == a), integer(1))),
  f("PROGENy_AM1to4_cell_counts.csv"),
  row.names = FALSE
)
sample_means <- aggregate(as.data.frame(scores[keep, , drop = FALSE], check.names = FALSE),
  by = list(sample = metadata$sample[keep], AM_group = metadata$AM_group[keep]), FUN = mean
)
write.csv(sample_means, f("PROGENy_AM1to4_sample_mean_scores.csv"), row.names = FALSE)

plot_matrix <- if (row_standardize) t(scale(t(means))) else means
stopifnot(identical(rownames(plot_matrix), am_order), !"AM-5" %in% rownames(plot_matrix))
write.csv(plot_matrix, f("PROGENy_AM1to4_heatmap_values.csv"))
breaks <- seq(heat_limits[1], heat_limits[2], length.out = length(heat_colors) + 1)
p <- pheatmap(plot_matrix,
  cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
  color = heat_colors, breaks = breaks, border_color = NA,
  fontsize = font_size, fontfamily = font_family,
  angle_col = column_angle, gaps_col = column_gaps,
  display_numbers = show_numbers, main = "PROGENy", silent = TRUE
)

if (interactive()) {
  grid::grid.newpage(); grid::grid.draw(p$gtable)
}
grDevices::cairo_pdf(f("PROGENy_AM1to4.pdf"), width = figure_width, height = figure_height)
grid::grid.newpage(); grid::grid.draw(p$gtable); dev.off()
png(f("PROGENy_AM1to4.png"), width = figure_width, height = figure_height, units = "in", res = 300, type = "cairo")
grid::grid.newpage(); grid::grid.draw(p$gtable); dev.off()
capture.output(sessionInfo(), file = f("PROGENy_sessionInfo.txt"))
saveRDS(
  list(
    signature = signature, displayed_AM = am_order, pathway_order = pathway_order,
    row_standardize = row_standardize, standardization_reference = "all original AM1-5 cells",
    aggregation = "mean of cell-level scaled scores, original cell-weighted convention",
    input_normalization = "existing nonnegative RNA data, unchanged",
    official_source = paste0("https://github.com/saezlab/progeny/tree/", model_commit)
  ),
  f("PROGENy_parameters.rds")
)
log_phase("DONE_AM1to4_heatmap_and_intermediate_data")
