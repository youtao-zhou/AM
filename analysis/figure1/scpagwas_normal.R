# scpagwas normal
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(Seurat); library(scPagwas); library(org.Hs.eg.db); library(Matrix)
})
out <- am_data("vibe/191_scpagwas_hg19/normal_allcelltype"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time(); log_phase <- function(p) {
  z <- sprintf("%s\t%s\t%.1f", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), p, as.numeric(difftime(Sys.time(), t0, units = "secs"))); cat(z, "\n"); cat(z, "\n", file = file.path(out, "run.log"), append = TRUE)
}
src <- readLines(am_code("analysis/figure1/scPagwas_cycle.R"))
ex <- parse(text = src); for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% c("RenameGenesSeurat", "RenameGenesInSeuratList")) eval(e)
stopifnot(exists("RenameGenesSeurat"), exists("RenameGenesInSeuratList"))
E <- am_data("vibe/191_scpagwas_hg19/normal_allcelltype/export")
X <- as(Matrix::readMM(file.path(E, "matrix.mtx")), "CsparseMatrix"); genes <- readLines(file.path(E, "genes.tsv")); cells <- data.table::fread(file.path(E, "cells.csv"))
dimnames(X) <- list(genes, cells[[1]])
options(Seurat.object.assay.version = "v3")
so <- CreateSeuratObject(counts = X, meta.data = data.frame(cell_type = cells$cell_type, row.names = cells[[1]]), min.cells = 0, min.features = 0)
so[["RNA"]]@data <- X
stopifnot(inherits(so[["RNA"]], "Assay"))
log_phase(sprintf("normal_subsample_loaded cells=%d", ncol(so)))
data.table::fwrite(data.table::data.table(cell = colnames(so)), file.path(out, "subsample_cells.csv"))
lst <- SplitObject(so, split.by = "cell_type"); rm(so); gc()
lst <- RenameGenesInSeuratList(seurat_list = lst)
m <- lst[[1]]; if (length(lst) > 1) for (j in 2:length(lst)) m <- merge(m, y = lst[[j]]); rm(lst); gc()
Idents(m) <- m@meta.data$cell_type
log_phase(sprintf("subsample_ready cells=%d types=%d", ncol(m), length(unique(Idents(m)))))
options(
  scPagwas.bootstrap.cores = 32L, scPagwas.pathway.cores = 16L, scPagwas.pathway.inner.cores = 1L,
  scPagwas.pathway.checkpoint.dir = file.path(out, "checkpoints"), scPagwas.pathway.backingpath = file.path(tempdir(), "scp191n")
)
dir.create(file.path(out, "checkpoints"), showWarnings = FALSE); dir.create(file.path(tempdir(), "scp191n"), showWarnings = FALSE)
set.seed(20260729)
Pagwas <- scPagwas_main(
  Pagwas = NULL, gwas_data = am_data("vibe/191_scpagwas_hg19/data/lcgwas_Pagwas_hg19.txt"), Single_data = m,
  output.prefix = "subset1", output.dirs = "vibe/191_scpagwas_hg19/normal_allcelltype/subset1", Pathway_list = Genes_by_pathway_kegg,
  n.cores = 32, assay = "RNA", singlecell = TRUE, iters_singlecell = 100, celltype = TRUE, iters_celltype = 200,
  block_annotation = block_annotation_hg37, chrom_ld = chrom_ld, seurat_return = FALSE
)
saveRDS(Pagwas, file.path(out, "Normal_scPagwas_hg19.rds"), compress = FALSE)
write.csv(as.data.frame(Pagwas$bootstrap_results), file.path(out, "subset1_celltypes_bootstrap_results_hg19.csv"))
log_phase("normal_allcelltype_complete")
