# 01 GRSA pseudobulk tumour vs normal
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ReporterScore)
})
base <- am_data("HTAN/am/GRSA"); dd <- file.path(base, "data")
t0 <- Sys.time()
log_phase <- function(p) {
  z <- sprintf("%s\t184_%s\t%.1f", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), p, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  cat(z, "\n"); cat(z, "\n", file = am_data("vibe/184_fig1e_grsa_redraw/run.log"), append = TRUE)
}

lu <- readRDS(am_data("HTAN/am/object/Rename_am.rds"))
no <- readRDS(am_data("HTAN/am_normal/Rename_am_normal.rds"))
log_phase("loaded")
expr <- function(o) LayerData(o, assay = "RNA", layer = "data")
smp <- function(o) {
  k <- intersect(c("sample", "Sample", "orig.ident"), colnames(o@meta.data))[1]; as.character(o@meta.data[[k]])
}
Xl <- expr(lu); Xn <- expr(no); sl <- smp(lu); sn <- smp(no)
genes <- intersect(rownames(Xl), rownames(Xn))
pseudobulk <- function(X, s) {
  n <- table(s); keep <- names(n)[n >= 20]
  list(mat = sapply(keep, function(k) Matrix::rowMeans(X[genes, s == k, drop = FALSE])), n = n[keep])
}
pl <- pseudobulk(Xl, sl); pn <- pseudobulk(Xn, sn)
colnames(pl$mat) <- paste0("LUAD|", colnames(pl$mat)); colnames(pn$mat) <- paste0("Normal|", colnames(pn$mat))
M <- cbind(pl$mat, pn$mat)
M <- M[rowSums(M > 0) >= 0.5 * ncol(M), ]
samples <- data.table(
  pseudobulk_id = colnames(M), sample = sub("^[^|]+\\|", "", colnames(M)),
  group = c(rep("LUAD", ncol(pl$mat)), rep("Normal", ncol(pn$mat))),
  n_AM = as.integer(c(pl$n, pn$n))
)
fwrite(samples, file.path(dd, "Fig1e_pseudobulk_samples.csv"))
cat(sprintf("pseudobulk samples: LUAD %d, normal %d; genes %d\n", ncol(pl$mat), ncol(pn$mat), nrow(M)))
log_phase("pseudobulk")

meta <- data.frame(Group = ifelse(samples$group == "LUAD", "am_tumor", "am_normal"), row.names = samples$pseudobulk_id)
set.seed(1)
res <- GRSA(as.data.frame(M), "Group", meta,
  mode = "directed", method = "wilcox.test", perm = 999,
  type = "hsa", feature = "gene", threads = 1, verbose = FALSE
)
rs <- as.data.table(res$reporter_s)
rs[, Description := sub(" - Homo sapiens \\(human\\)", "", Description)]
setorder(rs, -ReporterScore)
fwrite(rs, file.path(dd, "Fig1e_GRSA_pseudobulk_all_pathways.csv"))
writeLines(capture.output(sessionInfo()), file.path(dd, "Fig1e_sessionInfo.txt"))
cat(sprintf("pathways tested: %d; vs_group attribute: %s\n", nrow(rs), paste(attr(res$reporter_s, "vs_group"), collapse = " vs ")))
log_phase("grsa_done")
