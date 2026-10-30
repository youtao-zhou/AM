# am
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
setwd(am_data("HTAN/am"))

pbmc <- readRDS("am_cluster12.rds")
markers <- readRDS("marker12.rds")

library(Seurat)
library(dplyr)
library(tidyr)
library(data.table)
library(readxl)
library(org.Hs.eg.db)
library(clusterProfiler)

markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

genes <- top10_markers$gene

expression_matrix <- FetchData(seurat_object, vars = genes)

cluster_info <- Idents(seurat_object)

average_expression <- sapply(unique(cluster_info), function(cluster) {
  colMeans(expression_matrix[cluster_info == cluster, ])
})

colnames(average_expression) <- paste0("Cluster_", unique(cluster_info))

ordered_clusters <- paste0("Cluster_", 0:15)
average_expression_df <- average_expression_df[, c(ordered_clusters)]

average_expression_df <- as.data.frame(average_expression)
average_expression_df$Gene <- genes
average_expression_df <- average_expression_df[, c("Gene", paste0("Cluster_", unique(cluster_info)))]

top10_markers_list <- top10_markers %>%
  group_by(cluster) %>%
  summarise(genes = list(gene))
marker_list <- split(top10_markers_list$genes, top10_markers_list$cluster)

markers_res_0_5 <- list(
  `0` = c("DLX5", "ALDOA", "FKBP2", "CD9", "FABP4", "CES1", "CXCL8", "HLA-DQB2", "FOLR3", "DEFB1"),
  `1` = c("ALDOA", "DLX5", "FKBP2", "IGFBP2", "AKR1C3", "RND3", "SFTPC", "RPS4Y1", "MT1E", "PTP4A1"),
  `2` = c("RNASEK", "GABARAP", "EEF1G", "NME2", "ATP6V0C", "CRIP1", "IFI30", "RPL17", "GNG10", "POLD4"),
  `3` = c("C5AR1", "ADAMTSL4", "GNB1", "TMBIM1", "SEPTIN6", "SIRPB1", "CS", "CD4", "MME", "HDGF"),
  `4` = c("RNASEK", "GABARAP", "NME2", "EEF1G", "ATP6V0C", "IFI30", "APOC2", "RPL17", "GNG10", "CRIP1"),
  `5` = c("ALDOA", "DLX5", "CCL20", "CCL4", "MT2A", "FKBP2", "CCL3", "CCL18", "IL1B", "CXCL2"),
  `6` = c("SPP1", "RNASE1", "MARCKS", "RGS1", "CCL3", "IER3", "ZFP36L1", "IGKC", "LGMN", "CCL4L2"),
  `7` = c("ARHGAP26", "LRMDA", "FCHSD2", "DLEU2", "SPIDR", "MAML2", "SLC16A10", "PRKCE", "TPRG1", "DPYD"),
  `8` = c("IGHA1", "IGKC", "COL1A1", "COL1A2", "MGP", "COL3A1", "IGHGP", "IGFBP7", "SPARCL1", "RGS1"),
  `9` = c("MCF2L2", "LSAMP", "SLC27A4", "MAML2", "RABGEF1", "ATG7", "SSH2", "SNTB1", "SPIDR", "RAD51B"),
  `10` = c("LSAMP", "TTC39B", "SIK3", "ATP10A", "ZSWIM6", "KLHL29", "TANC2", "ZFAND3", "RAD51B", "NFAT5"),
  `11` = c("CCL21", "IGHA1", "BEST1", "HLA-DRB5", "SCGB3A1", "SFTPC", "HLA-DQA2", "NBEAL1", "HLA-DQB2", "ALDOA"),
  `12` = c("VCAN", "BASP1", "CD84", "PMP22", "PEA15", "A2M", "ZFP36L1", "PLTP", "RNASE1", "C5AR1")
)

long_df <- do.call(rbind, lapply(names(markers_res_0_5), function(cluster) {
  data.frame(Gene = markers_res_0_5[[cluster]], Group = cluster, stringsAsFactors = FALSE)
}))

long_df$Expression <- NA

for (cluster in unique(long_df$Group)) {
  cluster_cells <- WhichCells(pbmc, idents = cluster)
  cluster_data <- GetAssayData(pbmc, slot = "data")[, cluster_cells]
  for (gene in markers_res_0_5[[cluster]]) {
    if (gene %in% rownames(cluster_data)) {
      long_df$Expression[long_df$Gene == gene & long_df$Group == cluster] <- sum(cluster_data[gene, ])
    }
  }
}

wide_df <- long_df %>%
  pivot_wider(names_from = Group, values_from = Expression, names_prefix = "Cluster_", values_fill = list(Expression = 0))

DimPlot(pbmc, reduction = "umap")

unique_clusters <- unique(markers$cluster)

kegg_results_list <- list()
combined_kegg_results <- data.frame()

for (cluster_id in unique_clusters) {
  cluster_genes <- markers[markers$cluster == cluster_id, "gene"]

  entrez_ids <- bitr(cluster_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  if (nrow(entrez_ids) > 0) {
    kegg_result <- enrichKEGG(
      gene = entrez_ids$ENTREZID,
      organism = 'hsa',
      keyType = 'kegg',
      pvalueCutoff = 1
    )

    kegg_results_list[[as.character(cluster_id)]] <- kegg_result

    kegg_df <- as.data.frame(kegg_result@result)
    kegg_df$cluster <- cluster_id

    combined_kegg_results <- rbind(combined_kegg_results, kegg_df)
  }
}

top10_kegg_results <- combined_kegg_results %>%
  group_by(cluster) %>%
  arrange(pvalue) %>%
  slice_head(n = 10) %>%
  ungroup()

head(top10_kegg_results)

write.csv(top10_kegg_results, "kegg_top10.csv", row.names = FALSE)

filtered_kegg_results <- combined_kegg_results %>%
  filter(pvalue < 0.05)

head(filtered_kegg_results)

write.csv(filtered_kegg_results, "kegg_p005.csv", row.names = FALSE)

description_counts <- filtered_kegg_results %>%
  group_by(cluster) %>%
  summarise(Description_Count = n_distinct(Description))

unique_clusters <- unique(markers$cluster)

go_results_list <- list()
combined_go_results <- data.frame()

for (cluster_id in unique_clusters) {
  cluster_genes <- markers[markers$cluster == cluster_id, "gene"]

  entrez_ids <- bitr(cluster_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  if (nrow(entrez_ids) > 0) {
    go_result <- enrichGO(
      gene = entrez_ids$ENTREZID,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pvalueCutoff = 1,
      keyType = 'ENTREZID'
    )

    go_results_list[[as.character(cluster_id)]] <- go_result

    go_df <- as.data.frame(go_result@result)
    go_df$cluster <- cluster_id

    combined_go_results <- rbind(combined_go_results, go_df)
  }
}

top10_go_results <- combined_go_results %>%
  group_by(cluster) %>%
  arrange(pvalue) %>%
  slice_head(n = 10) %>%
  ungroup()

filtered_go_results <- combined_go_results %>%
  filter(p.adjust < 0.05)

head(filtered_go_results)

write.csv(filtered_go_results, "go_p005.csv", row.names = FALSE)

description_counts <- filtered_go_results %>%
  group_by(cluster) %>%
  summarise(Description_Count = n_distinct(Description))

description_markers <- markers %>%
  group_by(cluster) %>%
  summarise(gene = n_distinct(gene))

library(Seurat)

new_cluster_names <- c(
  "AM-1.FABP4",
  "AM-1.MT1E",
  "AM-1.APOC2- RNASEK+",
  "AM-2.CD4",
  "AM-1.APOC2+ RNASEK+",
  "AM-2.CK",
  "AM-3.SPP1",
  "AM-4.SPIDR+ RAD51B-",
  "AM-5.IGHA1",
  "AM-4.SPIDR+ RAD51B+",
  "AM-4.SPIDR- RAD51B+",
  "AM-1.HLA",
  "AM-5.VCAN"
)

old_cluster_ids <- 0:12
names(new_cluster_names) <- old_cluster_ids
pbmc <- RenameIdents(pbmc, new_cluster_names)
pbmc@meta.data$cell_type <- Idents(pbmc)
umap_data <- Embeddings(pbmc, "umap")
x_range <- range(umap_data[, 1])
y_range <- range(umap_data[, 2])

x_range_extended <- c(x_range[1] - 3, x_range[2] + 1)
y_range_extended <- c(y_range[1] - 1, y_range[2] + 1)

umap_plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) +
  NoLegend() +
  coord_cartesian(xlim = x_range_extended, ylim = y_range_extended)

print(umap_plot)
saveRDS(pbmc, "Rename_am.rds")
library(dior)
dior::write_h5(pbmc, "Rename_am.h5")
