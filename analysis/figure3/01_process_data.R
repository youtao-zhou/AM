# 01 process data
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(Seurat)
library(Matrix)
library(data.table)
library(uwot)

setwd(am_data("HTAN/am/tcell"))
tcell <- readRDS("data/Tcell_full100k_input.rds")
anno <- fread("data/Tcell_full100k_annotation.csv.gz")
stopifnot(
  ncol(tcell) == 100000L, nrow(anno) == 100000L,
  !anyDuplicated(anno$cell), setequal(colnames(tcell), anno$cell)
)
anno <- anno[match(colnames(tcell), cell)]

groups <- c(
  "CD4 Tn/Tcm", "CD4 Tem", "CD4 Tfh", "CD8 Tm", "CD8 CTL",
  "CD8 Tex", "Treg", "T-IFN", "Cycling T"
)
cols <- setNames(c(
  "#409B96", "#83B76F", "#CE9850", "#5286BF", "#6B60A4",
  "#BF6371", "#A86196", "#A2AA41", "#D87946"
), groups)
stopifnot(!anyNA(anno$T_short), setequal(anno$T_short, groups))

for (k in c("Group", "analysis_eligible", "annotation_status")) {
  tcell[[paste0(k, "_previous54")]] <- tcell[[k]][, 1]
}
tcell$Group <- anno$original_Group
tcell$analysis_eligible <- TRUE
tcell$T_primary_eligible <- TRUE
tcell$T_input_scope <- "Original 100000-cell atlas; user-approved full scope"
annotation_cols <- c(grep("^T_", names(anno), value = TRUE), "previous_analysis_eligible")
for (k in annotation_cols) tcell[[k]] <- anno[[k]]
tcell$T_short <- factor(anno$T_short, levels = groups)
tcell$sub_type <- tcell$T_short
tcell$T_major <- factor(ifelse(anno$T_short %in% groups[1:3], "CD4 T",
  ifelse(anno$T_short %in% groups[4:6], "CD8 T",
    ifelse(anno$T_short == "Treg", "Treg", "Cross-lineage state")
  )
))
Idents(tcell) <- "T_short"

old <- Embeddings(tcell, "umap")
colnames(old) <- c("UMAPOLD_1", "UMAPOLD_2")
tcell[["umap_previous54"]] <- CreateDimReducObject(old, key = "UMAPOLD_", assay = "RNA")
h <- Embeddings(tcell, "harmony")[, 1:30, drop = FALSE]
set.seed(58)
u <- uwot::umap(h,
  n_neighbors = 30, min_dist = 0.04, spread = 1,
  n_epochs = 400, metric = "euclidean", init = "spectral",
  n_threads = 16, n_sgd_threads = 1, repulsion_strength = 1,
  negative_sample_rate = 5, ret_nn = TRUE, ret_model = FALSE,
  verbose = TRUE, seed = 58
)
emb <- u$embedding
rownames(emb) <- colnames(tcell)
colnames(emb) <- c("UMAP_1", "UMAP_2")
tcell[["umap"]] <- CreateDimReducObject(emb, key = "UMAP_", assay = "RNA")
saveRDS(u$nn, "data/Tcell_full100k_Harmony_neighbors30.rds")
coords <- data.table(
  cell = colnames(tcell), UMAP1 = emb[, 1], UMAP2 = emb[, 2],
  T_short = as.character(tcell$T_short)
)
fwrite(coords, "data/Tcell_full100k_plot_coordinates.csv.gz")

genes <- c(
  "CD4", "IL7R", "CCR7", "TCF7", "CCL5", "GZMA", "CXCL13", "ICOS",
  "IL21", "CD8B", "GZMK", "FGFBP2", "PRF1", "HAVCR2", "ENTPD1",
  "LAG3", "FOXP3", "IL2RA", "ISG15", "MX1", "MKI67", "TOP2A"
)
expr <- GetAssayData(tcell, assay = "RNA", layer = "data")[genes, ]
stopifnot(all(is.finite(expr@x)), all(expr@x >= 0))
markers <- rbindlist(lapply(groups, function(g) {
  x <- expr[, tcell$T_short == g, drop = FALSE]
  data.table(
    group = g, gene = genes, n_cells = ncol(x),
    mean = Matrix::rowMeans(x), pct = Matrix::rowMeans(x > 0) * 100
  )
}))
markers[, z := if (sd(mean) > 0) as.numeric(scale(mean)) else 0, by = gene]
fwrite(markers, "data/Tcell_full100k_marker_plot_data.csv")
counts <- anno[, .(
  n_cells = .N, n_retained_reference = sum(T_in_previous_balanced),
  n_restored = sum(!T_in_previous_balanced),
  n_identity_review = sum(T_identity_review)
), by = T_short]
fwrite(counts, "data/Tcell_full100k_group_counts.csv")

parameters <- data.table(
  n_cells = ncol(tcell), dimensions = 30L, n_neighbors = 30L,
  min_dist = 0.04, spread = 1, n_epochs = 400L, seed = 58L,
  supervised = FALSE, manual_coordinate_offsets = FALSE
)
fwrite(parameters, "data/Tcell_full100k_UMAP_parameters.csv")
atlas_info <- list(
  name = "Full original 100k balanced T atlas", n_cells = ncol(tcell),
  groups = groups, colors = cols, group_column = "T_short", umap = parameters,
  input = "data/Tcell_full100k_input.rds",
  annotation = "Frozen balanced labels; provisional labels retained with review flags",
  expression = "Nonnegative log-normalized RNA/data; not raw counts",
  scope = "All original 100000 cells; no tissue-origin or AM-overlap restriction"
)
tcell@misc$full100k_balanced <- atlas_info
stopifnot(validObject(tcell), ncol(tcell) == 100000L, all(is.finite(emb)))
save(tcell, atlas_info,
  file = "data/Tcell_full100k_balanced.RData",
  compress = "gzip", compression_level = 1
)
print(counts)
message("已保存完整 100,000 细胞对象；02 和 03 可分别直接读入出图。")
