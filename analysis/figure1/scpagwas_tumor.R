# scpagwas tumor
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(Seurat); library(scPagwas)
})
out <- am_data("vibe/191_scpagwas_hg19/tumor_allcelltype"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time(); log_phase <- function(p) {
  z <- sprintf("%s\t%s\t%.1f", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), p, as.numeric(difftime(Sys.time(), t0, units = "secs"))); cat(z, "\n"); cat(z, "\n", file = file.path(out, "run.log"), append = TRUE)
}
obj <- readRDS(am_data("HTAN/scPagwas/scPagwas_lc/Tumor_allcelltype/subset1_ENSEMBL.rds"))
log_phase(sprintf("object_loaded cells=%d types=%d", ncol(obj), length(unique(Idents(obj)))))
options(
  scPagwas.bootstrap.cores = 16L, scPagwas.pathway.cores = 6L, scPagwas.pathway.inner.cores = 1L,
  scPagwas.pathway.checkpoint.dir = file.path(out, "checkpoints"), scPagwas.pathway.backingpath = file.path(tempdir(), "scp191")
)
dir.create(file.path(out, "checkpoints"), showWarnings = FALSE); dir.create(file.path(tempdir(), "scp191"), showWarnings = FALSE)
set.seed(20260729)
Pagwas <- scPagwas_main(
  Pagwas = NULL, gwas_data = am_data("vibe/191_scpagwas_hg19/data/lcgwas_Pagwas_hg19.txt"), Single_data = obj,
  output.prefix = "subset1", output.dirs = "vibe/191_scpagwas_hg19/tumor_allcelltype/subset1", Pathway_list = Genes_by_pathway_kegg,
  n.cores = 32, assay = "RNA", singlecell = TRUE, iters_singlecell = 100, celltype = TRUE, iters_celltype = 200,
  block_annotation = block_annotation_hg37, chrom_ld = chrom_ld, seurat_return = FALSE
)
saveRDS(Pagwas, file.path(out, "Tumor_scPagwas_hg19.rds"), compress = FALSE)
write.csv(as.data.frame(Pagwas$bootstrap_results), file.path(out, "subset1_celltypes_bootstrap_results_hg19.csv"))
log_phase("tumor_allcelltype_complete")
