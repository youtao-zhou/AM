# scPagwas cycle
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(Seurat)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(scCustomize)
library(scPagwas)
library(dior)

RenameGenesSeurat <- function(obj, newnames) {
  print("在整合之前运行此函数。它仅更改obj@assays$RNA@counts、@data 和 @scale.data。")
  if ("RNA" %in% names(obj@assays)) {
    RNA <- obj@assays$RNA
    if (nrow(RNA@counts) == length(newnames)) {
      dimnames(RNA@counts)[[1]] <- newnames
      if (!is.null(RNA@data) && nrow(RNA@data) == length(newnames)) {
        dimnames(RNA@data)[[1]] <- newnames
      }
      if (!is.null(RNA@scale.data) && nrow(RNA@scale.data) == length(newnames)) {
        dimnames(RNA@scale.data)[[1]] <- newnames
      }
      obj@assays$RNA <- RNA
    } else {
      stop("基因集不匹配：nrow(RNA) != length(newnames)")
    }
  } else {
    stop("提供的 Seurat 对象没有 RNA 测序。")
  }
  return(obj)
}

RenameGenesInSeuratList <- function(seurat_list) {
  for (i in seq_along(seurat_list)) {
    sce <- seurat_list[[i]]

    ids = AnnotationDbi::select(org.Hs.eg.db,
      keys = rownames(sce),
      columns = c('ENSEMBL', 'SYMBOL'),
      keytype = 'ENSEMBL'
    )
    ids = na.omit(ids)

    ids = ids[!duplicated(ids$SYMBOL), ]
    ids = ids[!duplicated(ids$ENSEMBL), ]

    ids <- ids[!duplicated(ids$ENSEMBL) & !duplicated(ids$SYMBOL), ]

    pos = match(ids$ENSEMBL, rownames(sce))
    sce = sce[pos, ]

    seurat_list[[i]] = RenameGenesSeurat(
      obj = sce,
      newnames = ids$SYMBOL
    )
  }
  return(seurat_list)
}

for (i in 1:1) {
  subset_name <- sprintf("subset_%d.h5", i)
  print(sprintf("正在处理 %s...", subset_name))

  tryCatch(
    {
      seurat_object <- dior::read_h5(subset_name)
      ifnb.list <- SplitObject(seurat_object, split.by = "cell_type")

      ifnb.list <- RenameGenesInSeuratList(seurat_list = ifnb.list)
      rm(seurat_object)
      gc()

      merged.seurat <- ifnb.list[[1]]
      if (length(ifnb.list) > 1) {
        for (j in 2:length(ifnb.list)) {
          merged.seurat <- merge(merged.seurat, y = ifnb.list[[j]])
        }
      }
      Idents(merged.seurat) <- merged.seurat@meta.data$cell_type
      saveRDS(merged.seurat, sprintf("subset%d_ENSEMBL.rds", i))

      Pagwas <- scPagwas_main(
        Pagwas = NULL,
        gwas_data = "lcgwas_Pagwas.txt",
        Single_data = sprintf("subset%d_ENSEMBL.rds", i),
        output.prefix = sprintf("subset%d", i),
        output.dirs = sprintf("subset%d", i),
        Pathway_list = Genes_by_pathway_kegg,
        n.cores = 1,
        assay = "RNA",
        singlecell = TRUE,
        iters_singlecell = 100,
        celltype = TRUE,
        block_annotation = block_annotation_hg37,
        chrom_ld = chrom_ld
      )
      save(Pagwas, file = sprintf("./subset%d_scPagwas.RData", i))
    },
    error = function(e) {
      print(sprintf("错误处理 %s: %s", subset_name, e$message))
    }
  )
}
