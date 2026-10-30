# STARTRAC am
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)

am <- readRDS(am_data("HTAN/am/object/Rename_am.rds"))
am_normal <- readRDS(am_data("HTAN/am_normal/Rename_am_normal.rds"))
am_metadata <- am@meta.data
am_normal_metadata <- am_normal@meta.data
colnames(am_normal_metadata)[11] <- "RNA_snn_res.0.5"
am_metadata$Group <- "Tumor"
am_normal_metadata$Group <- "Normal"
metadata <- rbind(am_metadata, am_normal_metadata)
metadata$Cell_Type <- substring(metadata$cell_type, 1, 4)

in.dat <- metadata
R_oe <- calTissueDist(in.dat,
  byPatient = F,
  colname.cluster = "Cell_Type",
  colname.tissue = "Group",
  method = "chisq",
  min.rowSum = 0
)
R_oe
col_fun <- colorRamp2(c(0, 0.75, 1.5), c("#FFFEE7", "#FFF0CF", "#E65100"))
row_order <- c("AM-1", "AM-2", "AM-3", "AM-4", "AM-5")
pdf("STARTRAC_heatmap_Group2.pdf")
Heatmap(as.matrix(R_oe),
  show_heatmap_legend = TRUE,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  row_names_side = 'right',
  show_column_names = TRUE,
  show_row_names = TRUE,
  col = col_fun,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    title = "R_o/e",
    at = seq(0, 1.5, by = 0.5),
    labels = seq(0, 1.5, by = 0.5),
    legend_gp = gpar(fill = col_fun(seq(0, 1.5, by = 0.5)))
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", R_oe[i, j]), x, y, gp = gpar(fontsize = 8))
  },
  row_order = row_order
)
dev.off()
