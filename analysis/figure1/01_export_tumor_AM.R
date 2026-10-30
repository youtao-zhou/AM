# 01 export tumor AM
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table)
})

root <- am_data("HTAN")
out_dir <- file.path(root, "am/cNMF/data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

t0 <- Sys.time()
log_phase <- function(phase, msg) {
  el <- sprintf("%.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs")))
  cat(sprintf(
    "%s\t%s\t%s\t%s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    phase, el, msg
  ), file = file.path(out_dir, "run.log"), append = TRUE)
  message(sprintf("[%s] %s :: %s", el, phase, msg))
}

log_phase("run_started", "export full tumor AM object (user grouping)")

obj <- readRDS(file.path(root, "am/object/Rename_am.rds"))
log_phase("phase1_object_loaded", sprintf("%d cells x %d genes", ncol(obj), nrow(obj)))

meta <- obj@meta.data
meta$barcode <- rownames(meta)

ph <- unique(fread(file.path(root, "cellhint/TRN_phenotype.csv"),
  select = c("sample", "origin")
), by = "sample")
meta$salcher_origin <- ph$origin[match(as.character(meta$sample), ph$sample)]
meta$salcher_origin[is.na(meta$salcher_origin)] <- "HTAN_not_in_Salcher"

audit <- as.data.frame(table(
  cell_type = meta$cell_type,
  cohort = meta$Dataset2,
  origin = meta$salcher_origin
))
audit <- audit[audit$Freq > 0, ]
fwrite(audit, file.path(out_dir, "tumorobject_input_audit.csv"))
log_phase("phase2_audit_written", sprintf(
  "%d non-empty state x cohort x origin cells",
  nrow(audit)
))

cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
gs <- Matrix::rowSums(cnt)
cnt <- cnt[gs > 0, , drop = FALSE]
log_phase("phase3_genes_filtered", sprintf("%d genes retained", nrow(cnt)))

obs <- meta[, c(
  "barcode", "cell_type", "sample", "Dataset2", "salcher_origin",
  "n_genes_by_counts", "total_counts", "pct_counts_mt"
)]
obs$macro_AM <- sub("\\..*$", "", obs$cell_type)
fwrite(obs, file.path(out_dir, "tumorobject_cell_metadata.csv"))

writeMM(cnt, file.path(out_dir, "tumorobject_counts.mtx"))
writeLines(rownames(cnt), file.path(out_dir, "tumorobject_genes.txt"))
writeLines(colnames(cnt), file.path(out_dir, "tumorobject_barcodes.txt"))

log_phase("phase4_export_done", sprintf("wrote %d x %d", nrow(cnt), ncol(cnt)))
cat("R_ALL_DONE\n")
