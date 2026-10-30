# tumor umap
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(dplyr)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(SeuratObject)
library(scales)
library(grid)
library(ComplexHeatmap)
library(grDevices)
library(colorRamp2)
library(MetBrewer)
library(RColorBrewer)
source(am_code("analysis/figure1/BetterVisDotPlot_Heatmap_Seurat.R"))

sc_dataset <- readRDS(am_data("HTAN/am/object/Rename_am.rds"))
sc_dataset$cell_type2 <- substring(sc_dataset$cell_type, 1, 4)

tumor_colors <- readRDS(am_data("HTAN/am/am_tumor_palette.rds"))

emb <- as.data.frame(Embeddings(sc_dataset, "umap"))
colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
df_plot <- emb %>%
  mutate(
    ident      = as.character(Idents(sc_dataset)),
    cell_type2 = as.character(sc_dataset$cell_type2)
  )

am_levels <- c("AM-1", "AM-2", "AM-3", "AM-4", "AM-5")
df_am <- df_plot %>% filter(cell_type2 %in% am_levels)

am_labels <- df_am %>%
  group_by(cell_type2) %>%
  summarise(cx = mean(UMAP_1), cy = mean(UMAP_2), .groups = "drop")

xr <- range(df_plot$UMAP_1, na.rm = TRUE)
yr <- range(df_plot$UMAP_2, na.rm = TRUE)
dx <- diff(xr); dy <- diff(yr)
x0 <- xr[1]; y0 <- yr[1]
ax_len <- 0.12 * min(dx, dy)
hx <- ax_len; hy <- ax_len

p <- ggplot(df_plot, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = ident), size = 0.6, alpha = 0.9, show.legend = FALSE) +
  scale_color_manual(values = tumor_colors, guide = "none") +
  stat_ellipse(
    data = df_am,
    aes(group = cell_type2),
    type = "norm", level = 0.9,
    linetype = "dashed", linewidth = 0.35, color = "black",
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = am_labels,
    aes(x = cx, y = cy, label = cell_type2),
    inherit.aes = FALSE,
    size = 6, fontface = "bold",
    box.padding = 0.3, point.padding = 0.2,
    min.segment.length = 0, segment.size = 0.3
  ) +
  geom_segment(
    inherit.aes = FALSE,
    aes(x = x0, y = y0, xend = x0 + hx, yend = y0),
    linewidth = 0.8,
    arrow = arrow(length = unit(0.2, "cm"), type = "closed")
  ) +
  geom_segment(
    inherit.aes = FALSE,
    aes(x = x0, y = y0, xend = x0, yend = y0 + hy),
    linewidth = 0.8,
    arrow = arrow(length = unit(0.2, "cm"), type = "closed")
  ) +
  annotate("text", x = x0 + hx * 0.55, y = y0 - dy * 0.03, label = "UMAP1", size = 5) +
  annotate("text", x = x0 - dx * 0.03, y = y0 + hy * 0.55, label = "UMAP2", size = 5, angle = 90, vjust = 0.5) +
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(legend.position = "none")

pdf("am_tumor_UMAP.pdf")
p
dev.off()

df_legend <- data.frame(ident = unique(df_plot$ident))
p_legend <- ggplot(df_legend, aes(x = 1, color = ident)) +
  geom_point(aes(y = 1), size = 6) +
  scale_color_manual(values = tumor_colors, name = "") +
  theme_void() +
  theme(legend.position = "right", legend.text = element_text(size = 13))
legend_only <- cowplot::get_legend(p_legend)
pdf("am_tumor_umap_legend.pdf")
grid::grid.newpage(); grid::grid.draw(legend_only)
dev.off()

library(Seurat)
library(ggplot2)

genes_of_interest <- c(
  "FABP4", "MT1E", "APOC2", "RNASEK", "CD4", "CCL20", "CXCL2", "CCL4", "CCL3", "CCL18",
  "SPP1", "CD80", "IGHA1", "HLA-DRB5", "HLA-DQA2", "VCAN"
)
sc_dataset$cell_type <- factor(
  sub("^(AM-\\d+).*", "\\1", sc_dataset$cell_type),
  levels = paste0("AM-", 1:5)
)

custom_markers <- c(
  "AM-1.FABP4" = "AM-1",
  "AM-1.MT1E" = "AM-1",
  "AM-1.APOC2- RNASEK+" = "AM-1",
  "AM-1.APOC2+ RNASEK+" = "AM-1",
  "AM-1.HLA" = "AM-1",
  "AM-2.CD4" = "AM-2",
  "AM-2.CK" = "AM-2",
  "AM-3.SPP1" = "AM-3",
  "AM-4.CD80" = "AM-4",
  "AM-5.IGHA1" = "AM-5",
  "AM-5.VCAN" = "AM-5"
)

custom_colors <- c(
  "AM-1" = "#7AC7C9",
  "AM-2" = "#FBBE09",
  "AM-3" = "#E3A7C5",
  "AM-4" = "#91C79E",
  "AM-5" = "#DD7572"
)

BetterVis_DotPlot_Heatmap_Seurat(sc_dataset,
  features = genes_of_interest, cluster.features = F, cluster.idents = F, rotate.axis = F, idents.kmeans = 1, features.kmeans = 1,
  show.data.legend = F, show.idents.legend = F, show.split.legend = F, row.names.side = "right", column.names.side = "top",
  gene_markers = custom_markers, marker_colors = custom_colors, annotation_marker_legend = F,
  row.split = c(1, 1, 1, 1, 1, 2, 2, 3, 4, 5, 5),
  col.split = c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 4, 5, 1, 1, 5),
  data.colors = c("#1289B2", "#FCF6AF", "#E53849"), dots.size = 4, dots.type = "square root", inner.border = F, border_gp_control = F,
  features.names.size = 12, idents.names.size = 12, keep_marker_order = F
)

genes_of_interest <- c("PSORS1C1", "NOTCH4", "FTCD", "CH25H", "IGFBP1")
pdf(am_data("HTAN/am/umap_proportion_dotplot/DotPlot_am_IGHA1_target_Dotplot.pdf"), width = 12, height = 8)
dot_plot <- DotPlot(sc_dataset, features = genes_of_interest, group.by = "cell_type") +
  scale_color_gradient(low = "lightgrey", high = "red") +
  scale_size(range = c(4, 10)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 18),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

print(dot_plot)
dev.off()
