## Figrue 2
library(Seurat)
library(scRNAseq)
library(dplyr)
library(tidyr)
library(tidyverse)
library(data.table)
library(stringr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(scales)
library(ggnewscale)
library(RColorBrewer)
library(MetBrewer)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(slingshot)
library(tradeSeq)
library(KernSmooth)
library(SCENIC)
library(AUCell)
library(SCopeLoomR)
library(CellChat)
library(Startrac)
library(future)
library(BiocParallel)
library(tictoc)
library(sscVis)



data1<-read_excel("am_enrichment.xlsx",sheet = 1)
data2<-read_excel("am_enrichment.xlsx",sheet = 2)

common_descriptions <- intersect(data1$Description, data2$Description)

data1_common <- data1 %>% filter(Description %in% common_descriptions)
data2_common <- data2 %>% filter(Description %in% common_descriptions)

data1_common$Group <- "Tumor"
data2_common$Group <- "Normal"

combined_data <- bind_rows(data1_common, data2_common)



padj_cut    <- 0.05
top_per_AM  <- 20          

df_sig <- combined_data %>%
  filter(!is.na(p.adjust), p.adjust < padj_cut) %>%
  mutate(
    AM = sub(".*(AM-\\d+).*", "\\1", cell_type),
    score = -log10(p.adjust),
    term_id = if ("ID" %in% names(.)) ID else str_extract(Description, "GO:\\d+")
  ) %>%
  group_by(Group, AM) %>%
  slice_max(order_by = score, n = top_per_AM, with_ties = FALSE) %>%
  ungroup()


shared_ids <- df_sig %>%
  distinct(term_id, Group) %>%
  count(term_id) %>%
  filter(n >= 2) %>%
  pull(term_id)

df_shared <- df_sig %>%
  filter(!is.na(term_id)) %>%                 
  filter(term_id %in% shared_ids) %>%        
  mutate(
    Tissue = dplyr::recode(Group, "Normal" = "Normal AM", "Tumor" = "Tumor AM"),
    score  = -log10(p.adjust),             
    AM     = factor(AM, levels = c("AM-1","AM-2","AM-3","AM-4","AM-5"))
  ) %>%
  select(AM, term_id, Tissue, score, Count, p.adjust, Group, Description) %>%
  arrange(AM, desc(score))

df_shared2 <- df_shared %>%
  dplyr::mutate(
    Type = dplyr::recode(Group, "Normal" = "Normal AM", "Tumor" = "Tumor AM"),
    score  = -log10(p.adjust)
  )

BetterVis_GO_Dotplot_Facet <- function(
    df,                                
    classified_var = "Type",       
    facet_var      = "AM",            
    y_var          = "term_id",       
    x_var          = "score",         
    size_var       = "Count",         
    colors         = c("Normal AM" = "#4C97F3", "Tumor AM" = "#E15759"),
    ncol           = 3,
    strokeSize     = 0.25,
    point_size_min = 1.5,           
    point_size_max = 7,                
    strokeColor    = "grey30",
    alpha          = 0.85,
    legend_title   = "Cell Type",
    legend_position = "right",
    legend_title_size = 14,
    legend_text_size  = 12,
    legend_key_size   = 10,
    legend_point_size = 6,             
    base_size      = 13,
    x_title        = "-log10(FDR)",
    y_title        = "GO ID",
    x_limits       = NULL,           
    theme_custom   = NULL,            
    guides_custom  = NULL             
){
  stopifnot(all(c(classified_var, facet_var, y_var, x_var) %in% names(df)))
  suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(forcats); library(scales) })
  
  df_plot <- df %>%
    mutate(
      !!classified_var := .data[[classified_var]],
      !!facet_var      := .data[[facet_var]],
      !!y_var          := .data[[y_var]],
      !!x_var          := as.numeric(.data[[x_var]])
    ) %>%
    group_by(.data[[facet_var]]) %>%
    mutate(y_fac = forcats::fct_reorder(.data[[y_var]], .data[[x_var]], .fun = max)) %>%
    ungroup()
  if (!size_var %in% names(df_plot)) df_plot[[size_var]] <- 1
  rng <- range(df_plot[[size_var]], na.rm = TRUE)
  if (diff(rng) == 0) {
    df_plot$.size_vis <- (point_size_min + point_size_max) / 2
  } else {
    df_plot$.size_vis <- scales::rescale(df_plot[[size_var]],
                                         to = c(point_size_min, point_size_max),
                                         from = rng)
  }
  
  p <- ggplot(df_plot, aes(x = .data[[x_var]], y = y_fac)) +
    geom_point(
      aes(fill = .data[[classified_var]]),
      shape = 21, color = strokeColor, stroke = strokeSize,
      alpha = alpha, size = df_plot$.size_vis, show.legend = FALSE
    ) +
    geom_point(
      aes(fill = .data[[classified_var]], size = .data[[size_var]]),
      shape = 21, color = strokeColor, stroke = strokeSize,
      alpha = 0, show.legend = TRUE
    ) +
    facet_wrap(reformulate(facet_var), scales = "free", ncol = ncol) +
    scale_fill_manual(values = colors, name = legend_title,
                      guide = guide_legend(override.aes = list(size = legend_point_size, alpha = 1))) +
    scale_size_area(max_size = 7, limits = c(0, NA),
                    guide = guide_legend(title = "Gene Count",
                                         override.aes = list(alpha = 1))) +
    scale_x_continuous(limits = c(0,NA),breaks = scales::pretty_breaks(n = 4),
                       expand = expansion(mult = c(0.02, 0.08))) +
    scale_y_discrete(drop = TRUE) +
    labs(x = x_title, y = y_title) +
    theme_bw(base_size = base_size) +
    theme(
      panel.grid.major.x = element_line(color = "grey90"),
      panel.grid.minor   = element_blank(),
      axis.title.x       = element_text(face = "bold"),
      axis.title.y       = element_text(face = "bold"),
      strip.text         = element_text(face = "bold"),
      legend.position    = legend_position,
      legend.title       = element_text(size = legend_title_size, face = "bold"),
      legend.text        = element_text(size = legend_text_size),
      legend.key.size    = unit(legend_key_size, "pt")
    )
  
  if (!is.null(x_limits))     p <- p + coord_cartesian(xlim = x_limits)
  if (!is.null(theme_custom)) p <- p + theme_custom
  if (!is.null(guides_custom))p <- p + guides_custom
  return(p)
}



my_colors <- c("Normal AM" = "#80B1D3", "Tumor AM" = "#FB8072")

p <- BetterVis_GO_Dotplot_Facet(
  df = df_shared2,
  classified_var = "Type",
  facet_var      = "AM",
  y_var          = "term_id",
  x_var          = "score",
  size_var       = "Count",
  point_size_min = 2.5,
  point_size_max = 6,
  colors         = my_colors,
  ncol           = 3,
  legend_title   = "Type",
  legend_position = c(0.85,0.25),       
  legend_title_size = 12,
  legend_text_size  = 12,
  legend_key_size   = 15,
  legend_point_size = 4,
  x_title        = "-log10(FDR)",  
  y_title        = "",
  theme_custom   = theme(axis.text.y = element_text(size = 10))
)
print(p)

###################

am<-readRDS("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am/Rename_am.rds")
am_normal<-readRDS("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am_normal/Rename_am_normal.rds")
am_metadata<-am@meta.data
am_normal_metadata<-am_normal@meta.data
colnames(am_normal_metadata)[11]<-"RNA_snn_res.0.5"
am_metadata$Group<-"Tumor"
am_normal_metadata$Group<-"Normal"
metadata<-rbind(am_metadata,am_normal_metadata)
metadata$Cell_Type<-substring(metadata$cell_type,1,4)

in.dat<-metadata
R_oe <- calTissueDist(in.dat,
                      byPatient = F,
                      colname.cluster = "Cell_Type",
                      colname.tissue = "Group",
                      method = "chisq",
                      min.rowSum = 0)
R_oe


R_oe<-read.csv("R_oe.csv",row.names = 1)
col_fun <- colorRamp2(c(0, 0.75, 1.5), c("#04a3ff", "#ffffff","#ff349c"))
row_order <- c("AM-1", "AM-2", "AM-3", "AM-4", "AM-5")
Heatmap(as.matrix(R_oe),
        show_heatmap_legend = TRUE,
        cluster_rows = F,  
        cluster_columns = FALSE, 
        row_names_side = 'left',
        column_names_side = "top",
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        column_names_rot = 0,
        row_names_gp = gpar(fontsize = 18),
        column_names_gp = gpar(fontsize = 18),
        heatmap_legend_param = list(
          title = "Ro/e",
          title_gp = gpar(fontsize = 16, fontface = "bold"),
          at = seq(0, 1.5, by = 0.5),
          labels = seq(0, 1.5, by = 0.5),
          labels_gp = gpar(fontsize = 15),
          legend_gp = gpar(fill = col_fun(seq(0, 1.5, by = 0.5)))
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", R_oe[i, j]), x, y, gp = gpar(fontsize = 16))
        },
        row_order = row_order 
)

#######################

pbmc<-readRDS("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am/Rename_am.rds")
sce <- as.SingleCellExperiment(pbmc, assay = "RNA")

sce_slingshot2 <- slingshot(sce,
                            reducedDim = 'UMAP', 
                            clusterLabels = sce$cell_type,start.clus =c("AM-4.SPIDR+ RAD51B-","AM-2.CD4"),omega=TRUE) 

sim<-sce_slingshot2


color <- c(brewer.pal(12, 'Set3'), brewer.pal(8, 'Dark2'))[1:13]
color[color == "#FFED6F"] <- "#E41A1C" 
umap_data <- reducedDims(sim)$UMAP
cell_types <- sim$cell_type

lin1 <- getLineages(sim, clusterLabels = "cell_type", reducedDim = "UMAP") 
pdf("slingshot_point.pdf")
plot(umap_data, col= color[cell_types], pch=16, asp = 1, xlab = "UMAP1", ylab = "UMAP2")
lines(SlingshotDataSet(lin1), lwd=2, col = 'black', type = 'lineages')
unique_cell_types <- unique(cell_types)
for (cell_type in unique_cell_types) {
  cells_of_type <- which(cell_types == cell_type)
  x <- mean(umap_data[cells_of_type, 1])
  y <- mean(umap_data[cells_of_type, 2])
  text(x, y, labels = cell_type, col = color[cell_type], pos = 3, cex = 0.8, font = 2)
}
dev.off()

umap_data <- reducedDims(sim)$UMAP                
cell_types <- sim$cell_type
barcodes <- colnames(sim)


cyto_tumor<-readRDS("result.rds")
phenot<-readRDS("phenot_tumor.rds")
ct_all <- cyto_tumor$CytoTRACE                   
ct <- ct_all[barcodes]                            
rng <- range(ct, na.rm = TRUE)

pal <- colorRampPalette(c("#1289B2","#FCF6AF","#E53849"))(256)
ct01 <- (ct - rng[1]) / diff(rng)
cols <- pal[pmax(1, pmin(256, 1 + floor(ct01 * 255)))]

ord <- order(ct, na.last = TRUE)

lin1 <- getLineages(sim, clusterLabels = "cell_type", reducedDim = "UMAP")

pdf("slingshot_point_cytotrace_TumorAM.pdf", width = 6, height = 6)
op <- par(mar = c(4,4,3,6)) 

plot(umap_data[ord, 1], umap_data[ord, 2],
     col = ifelse(is.na(cols[ord]), "#DDDDDD", cols[ord]),
     pch = 16, cex = 0.6, asp = 1,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n",
     main = "")

lines(SlingshotDataSet(lin1), lwd = 2, col = "black", type = "lineages")

unique_cell_types <- unique(cell_types)
for (cell_type in unique_cell_types) {
  idx <- which(cell_types == cell_type)
  x <- mean(umap_data[idx, 1])
  y <- mean(umap_data[idx, 2])
  points(x, y, pch = 16, cex = 1.3, col = "black")
  text(x, y, labels = cell_type, pos = 3, cex = 0.8, font = 2)
}

usr <- par("usr")
dx <- diff(usr[1:2]); dy <- diff(usr[3:4])
x0 <- usr[1] + 0.06 * dx; y0 <- usr[3] + 0.06 * dy
ax_len <- 0.12 * min(dx, dy)
arrows(x0, y0, x0 + ax_len, y0, angle = 25, length = 0.08, lwd = 2, xpd = NA)
arrows(x0, y0, x0, y0 + ax_len, angle = 25, length = 0.08, lwd = 2, xpd = NA)
text(x0 + ax_len*0.55, y0 - 0.03*dy, "UMAP1", cex = 1.2)
text(x0 - 0.03*dx,     y0 + ax_len*0.55, "UMAP2", cex = 1.2, srt = 90)

leg_xmin <- usr[2] + 0.02 * dx
leg_xmax <- leg_xmin + 0.035 * dx
leg_ymin <- usr[3] + 0.10 * dy
leg_ymax <- usr[3] + 0.90 * dy

rimg <- as.raster(matrix(pal, nrow = 256, ncol = 1))  
rasterImage(rimg, leg_xmin, leg_ymin, leg_xmax, leg_ymax, xpd = NA)
rect(leg_xmin, leg_ymin, leg_xmax, leg_ymax, border = "grey30", xpd = NA)

ticks <- pretty(rng, n = 5)
ty <- leg_ymin + (ticks - rng[1]) / diff(rng) * (leg_ymax - leg_ymin)
segments(leg_xmax, ty, leg_xmax + 0.01*dx, ty, xpd = NA)
text(leg_xmax + 0.015*dx, ty, labels = format(ticks, digits = 3),
     adj = c(0, 0.5), cex = 0.9, xpd = NA)
text((leg_xmin + leg_xmax)/2, leg_ymax + 0.035*dy, "CytoTRACE",
     font = 2, cex = 1.0, xpd = NA)

par(op)
dev.off()


pdf("slingshot_curve_group.pdf")
plot(umap_data, col= color[cell_groups], pch=16, asp = 1, xlab = "UMAP1", ylab = "UMAP2")
lines(SlingshotDataSet(sim), lwd=2, col='black')
unique_cell_groups <- unique(cell_groups)
for (cell_group in unique_cell_groups) {
  cells_of_group <- which(cell_groups == cell_group)
  x <- mean(umap_data[cells_of_group, 1])
  y <- mean(umap_data[cells_of_group, 2])
  text(x, y, labels = cell_group, col = color[cell_group], pos = 3, cex = 0.8, font = 2)
}
dev.off()

##############
safelog <- function(x) ifelse(is.finite(x), x, NA_real_)

.get_mat <- function(x){
  if (inherits(x, "SummarizedExperiment")) {
    tryCatch(SummarizedExperiment::assay(x), error = function(e) as.matrix(x))
  } else {
    as.matrix(x)
  }
}

am_prefix <- function(v) sub("^(AM-[0-9]+).*", "\\1", v)

diff_regulon_AM <- function(x_normal, anno_normal, x_tumor, anno_tumor,
                            am_levels = NULL, min_cells = 15,
                            exact_wilcox = FALSE, fdr_method = "fdr") {
  
  M_n <- .get_mat(x_normal)
  M_t <- .get_mat(x_tumor)

  common_reg <- intersect(rownames(M_n), rownames(M_t))
  M_n <- M_n[common_reg, , drop = FALSE]
  M_t <- M_t[common_reg, , drop = FALSE]
  anno_normal$AM_group <- am_prefix(anno_normal$celltype)
  anno_tumor$AM_group  <- am_prefix(anno_tumor$celltype)
  
  if (is.null(am_levels)) {
    am_levels <- intersect(unique(anno_normal$AM_group), unique(anno_tumor$AM_group))
    am_levels <- sort(am_levels)
  }

  out_list <- vector("list", length(am_levels))
  names(out_list) <- am_levels
  
  for (am in am_levels) {
    cells_n <- rownames(anno_normal)[anno_normal$AM_group == am]
    cells_t <- rownames(anno_tumor )[anno_tumor$AM_group  == am]
    if (length(cells_n) < min_cells || length(cells_t) < min_cells) {
      warning(sprintf("Skip %s: too few cells (Normal=%d, Tumor=%d)", am, length(cells_n), length(cells_t)))
      next
    }
    Xn <- M_n[, intersect(cells_n, colnames(M_n)), drop = FALSE]
    Xt <- M_t[, intersect(cells_t, colnames(M_t)), drop = FALSE]
    pvec   <- numeric(nrow(M_n))
    mN     <- numeric(nrow(M_n))
    mT     <- numeric(nrow(M_n))
    
    for (i in seq_len(nrow(M_n))) {
      vn <- Xn[i, ]; vt <- Xt[i, ]
      mN[i] <- mean(vn, na.rm = TRUE)
      mT[i] <- mean(vt, na.rm = TRUE)
      pvec[i] <- tryCatch(
        wilcox.test(vt, vn, exact = exact_wilcox)$p.value,
        error = function(e) NA_real_
      )
    }
    
    df <- data.frame(
      AM        = am,
      regulon   = rownames(M_n),
      mean_tumor  = mT,
      mean_normal = mN,
      delta       = mT - mN,               
      p_value     = pvec,
      n_tumor     = ncol(Xt),
      n_normal    = ncol(Xn),
      stringsAsFactors = FALSE
    )
    df$FDR <- p.adjust(df$p_value, method = fdr_method)
    df <- df[order(df$FDR, -abs(df$delta)), ]
    out_list[[am]] <- df
  }

  res_tbl <- dplyr::bind_rows(out_list)
  ams_present <- unique(res_tbl$AM)
  regs_all    <- unique(res_tbl$regulon)
  
  to_mat <- function(value_col){
    M <- matrix(NA_real_, nrow = length(regs_all), ncol = length(ams_present),
                dimnames = list(regs_all, ams_present))
    for (am in ams_present) {
      dd <- out_list[[am]]
      M[dd$regulon, am] <- dd[[value_col]]
    }
    M
  }
  delta_mat   <- to_mat("delta")
  normal_mat  <- to_mat("mean_normal")
  tumor_mat   <- to_mat("mean_tumor")
  fdr_mat     <- to_mat("FDR")
  
  list(
    table     = res_tbl,  
    delta_mat = delta_mat,  
    normal_mat = normal_mat, 
    tumor_mat  = tumor_mat, 
    fdr_mat    = fdr_mat    
  )
}


load("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am_normal/SCENIC/for_rss_and_visual.Rdata")
auc_normal  <- sub_regulonAUC
anno_normal <- cellTypes

load("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am/SCENIC/for_rss_and_visual.Rdata")
auc_tumor   <- sub_regulonAUC
anno_tumor  <- cellTypes

res <- diff_regulon_AM(
  x_normal = auc_normal,  anno_normal = anno_normal,
  x_tumor  = auc_tumor,   anno_tumor  = anno_tumor,
  min_cells = 20,        
  exact_wilcox = FALSE,   
  fdr_method = "fdr"
)

head(res$table)

library(pheatmap)

alpha <- 0.05   
top_k <- 8      

pick <- res$table %>%
  dplyr::filter(FDR < alpha) %>%
  dplyr::group_by(AM) %>%
  dplyr::slice_max(order_by = abs(delta), n = top_k, with_ties = FALSE) %>%
  dplyr::ungroup()

regs_use <- unique(pick$regulon)
ams_use  <- unique(pick$AM)

pair_mat <- cbind(res$normal_mat[regs_use, ams_use, drop=FALSE],
                  res$tumor_mat [regs_use, ams_use, drop=FALSE])
colnames(pair_mat) <- as.vector(rbind(paste0(ams_use, "_Normal"),
                                      paste0(ams_use, "_Tumor")))

pair_scaled <- t(scale(t(pair_mat)))
pair_scaled[!is.finite(pair_scaled)] <- 0

ann_col <- data.frame(
  AM   = rep(ams_use, each = 2),
  Type = rep(c("Normal","Tumor"), times = length(ams_use)),
  row.names = colnames(pair_scaled)
)

hp <- pheatmap(pair_scaled, cluster_rows = TRUE, cluster_cols = FALSE)
tf_order <- rownames(pair_scaled)[hp$tree_row$order]  


pheatmap(pair_scaled,
         cluster_rows = TRUE, cluster_cols = FALSE,
         annotation_col = ann_col,
         color = colorRampPalette(c("#3B76AF","white","#D84A3A"))(300),
         breaks = seq(-1, 2, length.out = 300),  
         main = "Top regulon differences per AM (ΔAUC with FDR<0.05)",
         fontsize_row = 6, fontsize_col = 9, border_color = NA)

######## Visualization
BetterVis_Heatmap_Multi <- function(df, x_var, number_vars = NULL, annotation_var = NULL,
                                    color_heatmap = met.brewer("Cassatt1"), 
                                    color_annotation = NULL,
                                    legend_heatmap_position = "right",
                                    legend_annotation_position = c(1.3, -4),
                                    value_limits = c(-1, 2),          
                                    heatmap_theme = theme(),  
                                    annotation_theme = theme(),  
                                    heatmap_plot_margin = unit(c(0, 2, -20, 0.5), "cm"),
                                    row_order = NULL) {                 # ★ 新增
  x_var <- as.character(x_var)
  
  if (is.null(number_vars) || length(number_vars) == 0) {
    number_vars <- df %>% select(where(is.numeric)) %>% colnames()
  }
  
  df_char <- df %>%
    select(all_of(x_var), where(~ is.character(.) | is.factor(.))) %>%
    pivot_longer(-all_of(x_var), names_to = "name", values_to = "value")
  
  df_char$sample <- factor(df_char[[x_var]], levels = unique(df_char[[x_var]]))
  df_char$name   <- factor(df_char$name, levels = unique(df_char$name) %>% rev())
  
  if (!is.null(annotation_var)) {
    missing_cols <- setdiff(annotation_var, colnames(df))
    if (length(missing_cols) > 0) {
      stop(paste("Error: The following columns in annotation_var do not exist in the data frame:", 
                 paste(missing_cols, collapse = ", ")))
    }
    legend_data <- lapply(annotation_var, function(var) {
      subset_df <- df_char %>%
        filter(name == var) %>%
        dplyr::rename(!!var := value) %>%
        select(all_of(x_var), name, !!var)
      return(subset_df)
    })
  }
  
  df_long <- df %>%
    pivot_longer(cols = all_of(number_vars), names_to = "variable", values_to = "value")
  
  if (!is.null(row_order)) {
    keep <- intersect(row_order, unique(df_long$variable))
    df_long$variable <- factor(df_long$variable, levels = rev(keep)) 
  } else {
    df_long$variable <- factor(df_long$variable, levels = unique(df_long$variable) %>% rev())
  }
  
  p2 <- ggplot(df_long, aes(x = !!sym(x_var), y = variable)) +
    geom_tile(aes(fill = value), color = "black") +
    scale_fill_gradientn(
      colors = color_heatmap,
      limits = value_limits,
      oob = scales::squish,
      na.value = NA
    ) +
    labs(x = NULL, y = NULL) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = "black"),
      legend.position = legend_heatmap_position,
      legend.direction = "vertical",
      legend.spacing.x = unit(0.1, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.5, "cm"),
      plot.margin = heatmap_plot_margin,
      plot.background = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.text = element_text(color = "black"),
      legend.box.spacing = unit(0.5, "cm")
    ) +
    heatmap_theme
  
  if (!is.null(annotation_var)) {
    p1 <- ggplot(df_char, aes(x = !!sym(x_var), y = name))
    for (i in seq_along(legend_data)) {
      unique_values <- unique(legend_data[[i]][[annotation_var[i]]])
      num_colors <- length(unique_values)
      color_values <- if (!is.null(color_annotation) && length(color_annotation) >= i) {
        color_annotation[[i]]
      } else {
        scales::brewer_pal(palette = "Set3")(num_colors)
      }
      p1 <- p1 + 
        geom_tile(data = legend_data[[i]], color = "black", aes(fill = !!sym(annotation_var[i]))) +
        scale_fill_manual(values = color_values) +
        ggnewscale::new_scale_fill()
    }
    
    p1 <- p1 + 
      scale_y_discrete(position = "right") +
      labs(x = NULL, y = NULL) +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        legend.position = legend_annotation_position,
        legend.direction = "vertical",
        legend.spacing.x = unit(0.1, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.margin = unit(c(0, 1, -20, 0.5), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color = "black")
      ) + annotation_theme
    
    return(p1 / p2 + patchwork::plot_layout(heights = c(0.3, 2)))
  } else {
    return(p2)
  }
}
mat <- pair_scaled
mat[mat < -1] <- -1; mat[mat > 2] <- 2  

df_heat <- cbind(
  case.ID = colnames(mat),
  as.data.frame(t(mat), check.names = FALSE)
)

df_heat$Type <- ifelse(grepl("_Tumor$", df_heat$case.ID), "Tumor", "Normal")
df_heat$AM   <- sub("_(Normal|Tumor)$", "", df_heat$case.ID)

num_cols <- tf_order[tf_order %in% colnames(df_heat)]
x_levels <- as.vector(rbind(paste0("AM-1","_Normal"), paste0("AM-1","_Tumor"),
                            paste0("AM-2","_Normal"), paste0("AM-2","_Tumor"),
                            paste0("AM-3","_Normal"), paste0("AM-3","_Tumor"),
                            paste0("AM-4","_Normal"), paste0("AM-4","_Tumor"),
                            paste0("AM-5","_Normal"), paste0("AM-5","_Tumor")))
df_heat$case.ID <- factor(df_heat$case.ID, levels = intersect(x_levels, df_heat$case.ID))
p<-BetterVis_Heatmap_Multi(
  df = df_heat,
  x_var = "case.ID",
  number_vars = num_cols,                       
  annotation_var = c("AM","Type"),
  color_heatmap = rev(met.brewer("Cassatt1")),
  color_annotation = list(
    c("#adc9ca","#a7a5d7","#ffac5f","#7aa5d2","#e49898"), 
    c("#80B1D3","#FB8072")                                 
  ),
  value_limits = c(-1, 2),                                 
  legend_heatmap_position = c(1.08, 0.7),
  legend_annotation_position = c(1.1, -5),
  heatmap_theme = theme(axis.text.y = element_text(face="bold",size = 12),
                        legend.title = element_text(size=14),
                        legend.text = element_text(size=12)),
  annotation_theme = theme(axis.text.y = element_text(face="bold",size=12),
                           legend.title = element_text(size=12),
                           legend.text = element_text(size=12)),
  row_order = tf_order                            
)
p
ggsave("SCENIC_vs_Normal.pdf",p,width=6.5,height=4.72)

############################
dir_tumor  <- "/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am/cytosig"
dir_normal <- "/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am_normal/cytosig"
rds_tumor  <- "/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am/Rename_am.rds"
rds_normal <- "/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am_normal/Rename_am_normal.rds"
dir.create("results", showWarnings = FALSE); dir.create("plots", showWarnings = FALSE)

FDR_THR <- 0.05
D_THR   <- 0.60      

read_coef <- function(dir_path, nThread = max(1, parallel::detectCores()-1)) {
  df <- fread(file.path(dir_path, "cytosig_scores.Coef"), header = FALSE, sep = "\t",
              data.table = FALSE, nThread = nThread, showProgress = TRUE)
  rn <- df[[1]]; df <- df[, -1, drop = FALSE]; rownames(df) <- rn
  as.matrix(df)
}

coef_tum  <- read_coef(dir_tumor)
coef_norm <- read_coef(dir_normal)

seu_tum  <- readRDS(rds_tumor)
seu_norm <- readRDS(rds_normal)
stopifnot(ncol(coef_tum)  == ncol(seu_tum),
          ncol(coef_norm) == ncol(seu_norm))
colnames(coef_tum)  <- colnames(seu_tum)
colnames(coef_norm) <- colnames(seu_norm)


pick_field <- function(meta, candidates){
  hit <- intersect(candidates, colnames(meta))
  if (!length(hit)) stop("not found", paste(candidates, collapse=", "))
  hit[1]
}

pseudo_by_donor_AM <- function(mat_coef, seurat_obj,
                               subtype_field = "cell_type", donor_field = NULL) {
  meta <- seurat_obj@meta.data
  if (is.null(donor_field)) donor_field <- pick_field(meta, c("donor","patient","case.ID","sample","orig.ident"))
  AM  <- sub("^(AM-[1-5]).*", "\\1", as.character(meta[[subtype_field]]))
  dnr <- as.character(meta[[donor_field]])
  keep <- !is.na(AM) & !is.na(dnr)
  mat_coef <- mat_coef[, keep, drop=FALSE]; AM <- AM[keep]; dnr <- dnr[keep]
  
  grp <- factor(paste(dnr, AM, sep="|"))
  G <- sparseMatrix(i = seq_along(grp), j = as.integer(grp), x = 1,
                    dims = c(length(grp), nlevels(grp)),
                    dimnames = list(NULL, levels(grp)))
  sum_mat <- mat_coef %*% G
  n_cells <- colSums(G)
  mean_mat <- sweep(sum_mat, 2, n_cells, "/")
  mean_mat <- as.matrix(mean_mat); storage.mode(mean_mat) <- "double"
  
  info <- data.frame(do.call(rbind, strsplit(colnames(mean_mat), "\\|")),
                     stringsAsFactors = FALSE)
  colnames(info) <- c("donor","AM"); rownames(info) <- colnames(mean_mat)
  list(mat = mean_mat, info = info)
}

pb_tum  <- pseudo_by_donor_AM(coef_tum,  seu_tum)
pb_norm <- pseudo_by_donor_AM(coef_norm, seu_norm)

compute_effect_size <- function(expr, samp, paired = TRUE){
  idxN <- which(samp$Tissue == "Normal")
  idxT <- which(samp$Tissue == "Tumor")
  donorsN <- samp$donor[idxN]; donorsT <- samp$donor[idxT]
  common  <- intersect(donorsN, donorsT)
  
  if (paired && length(common) >= 3){
    n_idx <- idxN[match(common, donorsN)]
    t_idx <- idxT[match(common, donorsT)]
    D <- expr[, t_idx, drop=FALSE] - expr[, n_idx, drop=FALSE]  
    dz <- rowMeans(D, na.rm=TRUE) / apply(D, 1, sd, na.rm=TRUE)
    return(list(x = dz, label = "dz", n_pair = length(common)))
  } else {
    xN <- expr[, idxN, drop=FALSE]; xT <- expr[, idxT, drop=FALSE]
    mN <- rowMeans(xN, na.rm=TRUE);  mT <- rowMeans(xT, na.rm=TRUE)
    sN <- apply(xN, 1, sd, na.rm=TRUE); sT <- apply(xT, 1, sd, na.rm=TRUE)
    nN <- ncol(xN); nT <- ncol(xT)
    sp <- sqrt(((nN-1)*sN^2 + (nT-1)*sT^2) / pmax(nN+nT-2, 1))
    d  <- (mT - mN) / sp
    return(list(x = d, label = "d", n_pair = 0))
  }
}

run_one_AM <- function(AMx, fdr_thr = FDR_THR, d_thr = D_THR,
                       out_res = "results", out_fig = "plots") {
  n_col <- which(pb_norm$info$AM == AMx)
  t_col <- which(pb_tum$info$AM  == AMx)
  if (length(n_col) + length(t_col) < 6) {
    message(AMx, "< 6"); return(NULL)
  }
  expr <- cbind(pb_norm$mat[, n_col, drop=FALSE],
                pb_tum$mat [, t_col, drop=FALSE])
  samp <- data.frame(
    sample = colnames(expr),
    donor  = c(pb_norm$info$donor[n_col], pb_tum$info$donor[t_col]),
    Tissue = factor(c(rep("Normal", length(n_col)), rep("Tumor", length(t_col))),
                    levels = c("Normal", "Tumor"))
  )
  rownames(samp) <- samp$sample
  
  design <- model.matrix(~ Tissue, data = samp)
  corfit <- duplicateCorrelation(expr, design, block = samp$donor)
  fit <- lmFit(expr, design, block = samp$donor, correlation = corfit$consensus)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "TissueTumor", number = Inf, sort.by = "P")
  tt$Effect   <- tt$logFC
  tt$FDR      <- tt$adj.P.Val
  tt$Cytokine <- rownames(tt)
  tt$AM       <- AMx
  
  eff  <- compute_effect_size(expr, samp, paired = TRUE)
  xvec <- eff$x[match(tt$Cytokine, rownames(expr))]
  xvar <- eff$label
  tt[[xvar]] <- xvec
  
  tt$Sig <- (tt$FDR < fdr_thr) & (abs(tt[[xvar]]) >= d_thr)
  fwrite(tt[, c("Cytokine","AM","Effect", xvar, "P.Value","FDR","Sig")],
         file = file.path(out_res, paste0(AMx, "_limma_", xvar, ".csv")))

  xlim_val <- 1.1 * quantile(abs(tt[[xvar]]), 0.99, na.rm = TRUE)
  ylim_val <- 1.1 * quantile(-log10(tt$FDR), 0.99, na.rm = TRUE)
  
  p <- ggplot(tt, aes_string(x = xvar, y = "-log10(FDR)")) +
    geom_point(aes(alpha = Sig), size = 2) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
    geom_vline(xintercept = c(-d_thr, d_thr), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed", linewidth = 0.4) +
    coord_cartesian(xlim = c(-xlim_val, xlim_val), ylim = c(0, ylim_val)) +
    labs(title = paste0(AMx, " Volcano (", xvar, ", paired n=", eff$n_pair, ")"),
         x = paste0("Cohen’s ", xvar, " (Tumor − Normal)"),
         y = "-log10(FDR)") +
    theme_classic(base_size = 12)
  
  lab_df <- tt %>% filter(Sig) %>% arrange(FDR) %>% slice_head(n = 10)
  if (nrow(lab_df) > 0) {
    p <- p + ggrepel::geom_text_repel(data = lab_df, aes(label = Cytokine),
                                      size = 3, max.overlaps = 20, box.padding = 0.3)
  }
  ggsave(file.path(out_fig, paste0(AMx, "_volcano_", xvar, ".pdf")), p, width = 5.3, height = 4.3)
  ggsave(file.path(out_fig, paste0(AMx, "_volcano_", xvar, ".png")), p, width = 5.3, height = 4.3, dpi = 300)
  invisible(tt)
}

AMs <- intersect(paste0("AM-", 1:5),
                 intersect(unique(pb_norm$info$AM), unique(pb_tum$info$AM)))
res_list <- lapply(AMs, run_one_AM)  
names(res_list) <- AMs

df <- bind_rows(res_list) %>% as.data.frame()

xvar <- if ("dz" %in% names(df)) "dz" else if ("d" %in% names(df)) "d" else "Effect"
fdr_thr <- FDR_THR
d_thr   <- if (xvar %in% c("d","dz")) D_THR else max(0.01, 0.5*quantile(abs(df$Effect), 0.75, na.rm=TRUE))

df <- df %>% mutate(SigStat = FDR < fdr_thr & abs(.data[[xvar]]) >= d_thr)

sig_only <- df %>% filter(SigStat) %>%
  mutate(Direction = ifelse(.data[[xvar]] > 0, "Up", "Down")) %>%
  group_by(Cytokine) %>%
  summarise(n_AM = n_distinct(AM),
            dir_same = n_distinct(Direction) == 1,
            .groups="drop")

df <- df %>% left_join(sig_only, by = "Cytokine") %>%
  mutate(n_AM = ifelse(is.na(n_AM), 0L, n_AM),
         OverlapCat = case_when(
           !SigStat ~ "NS",
           n_AM >= 3 & dir_same ~ "Same Direction≥3",
           n_AM >= 2 & dir_same ~ "Same Direction≥2",
           n_AM >= 2 & !dir_same ~ "Diff Direction",
           TRUE ~ "Unique"
         ))

xlim_val <- 1.1 * quantile(abs(df[[xvar]]), 0.99, na.rm = TRUE)

cols <- c("NS"="grey85","Unique"="#6ea6ff",
          "Same Direction≥2"="#ff9f43","Same Direction≥3"="#ff5d5d",
          "Diff Direction"="#9b59b6")

p_volcano <- ggplot(df, aes(x = .data[[xvar]], y = -log10(FDR))) +
  facet_wrap(~ AM, ncol = 3, scales = "free_y") +
  geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed", linewidth = 0.35) +
  geom_vline(xintercept = c(-d_thr, d_thr), linetype = "dashed", linewidth = 0.35) +
  geom_point(aes(fill = OverlapCat), shape = 21, size = 2.2, alpha = 0.95,
             colour = "black", stroke = 0.15) +
  scale_fill_manual(values = cols, drop = FALSE) +
  coord_cartesian(xlim = c(-xlim_val, xlim_val)) +
  labs(title = "",
       x = if (xvar %in% c("d","dz")) paste0("Cohen's ", xvar, " (Tumor − Normal)") else "Effect (Tumor − Normal)",
       y = "-log10(FDR)", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.8,0.28), 
        legend.text = element_text(size=12),
        legend.title =  element_text(size=12),
        axis.title = element_text(size=12),
        panel.grid.minor = element_blank())
p_volcano

lab_uni <- df %>% filter(SigStat, OverlapCat == "Unique") %>%
  group_by(AM) %>% slice_max(order_by = abs(.data[[xvar]]), n = 8, with_ties = FALSE) %>% ungroup()
lab_mix <- df %>% filter(SigStat, OverlapCat == "Shared (mixed dir)") %>%
  group_by(AM) %>% slice_max(order_by = abs(.data[[xvar]]), n = 8, with_ties = FALSE) %>% ungroup()

if (nrow(lab_uni) > 0)
  p_volcano <- p_volcano + ggrepel::geom_text_repel(
    data = lab_uni, aes(label = Cytokine),
    size = 3, box.padding = 0.3, max.overlaps = Inf, min.segment.length = 0.03)


if (nrow(lab_mix) > 0)
  p_volcano <- p_volcano + ggrepel::geom_text_repel(
    data = lab_mix, aes(label = Cytokine),
    size = 3, box.padding = 0.25, max.overlaps = 80, min.segment.length = 0.03, alpha = 0.9)
p_volcano

text_adjust    <- 1.25    
y_jitter       <- 0.2   
label_box_size <- 0.22     
label_topN_per_am <- Inf  

cols <- c(
  "NS"                 = "grey85",
  "Unique"             = "#6ea6ff",
  "Same Direction≥2"   = "#ff9f43",
  "Same Direction≥3"   = "#ff5d5d",
  "Diff Direction"     = "#9b59b6"
)

xvar      <- if ("dz" %in% names(df)) "dz" else if ("d" %in% names(df)) "d" else "Effect"
am_levels <- levels(factor(df$AM))
am_to_y   <- setNames(seq_along(am_levels), am_levels)

if (!"sizeFDR" %in% names(df)) df <- df %>% mutate(sizeFDR = pmin(-log10(FDR), 10))

label_cats <- c("Unique", "Diff Direction")
lab_df <- df %>%
  filter(SigStat, OverlapCat %in% label_cats) %>%
  group_by(AM) %>%
  { if (is.finite(label_topN_per_am)) slice_max(., order_by = abs(.data[[xvar]]),
                                                n = label_topN_per_am, with_ties = FALSE) else . } %>%
  ungroup()

df$y_num <- am_to_y[df$AM]
lab_df$y_num <- am_to_y[lab_df$AM]

set.seed(7)

p_dotplot<- ggplot() +
  geom_vline(xintercept = c(-d_thr, d_thr), linetype = "dashed", linewidth = 0.35) +
  geom_point(
    data = df %>% filter(!SigStat),
    aes(x = .data[[xvar]], y = y_num, size = sizeFDR),
    shape = 21, stroke = 0.15, fill = "grey85", colour = "black", alpha = 0.45,
    position = position_jitter(height = y_jitter, width = 0)
  ) +

  geom_point(
    data = df %>% filter(SigStat),
    aes(x = .data[[xvar]], y = y_num, fill = OverlapCat, size = sizeFDR),
    shape = 21, stroke = 0.2, colour = "black",
    position = position_jitter(height = y_jitter, width = 0)
  ) +

  ggrepel::geom_label_repel(
    data = lab_df,
    aes(x = .data[[xvar]], y = y_num, label = Cytokine, colour = OverlapCat),
    fill = "white",
    label.size = label_box_size,
    segment.color = NA,
    size = 2.6 * text_adjust,
    box.padding = 0.4, point.padding = 0.25,
    show.legend = FALSE
  ) +

  scale_fill_manual(
    values = cols, drop = FALSE, name = NULL,
    guide = guide_legend(
      override.aes = list(size = 4, shape = 21, colour = "black", stroke = 0.2)
    )
  ) +

  scale_size(
    range = c(1.2, 4.6),
    breaks = c(2.5, 5, 7.5, 10),
    labels = c("2.5", "5.0", "7.5", "10.0"),
    guide = guide_legend(
      title = "-log10(FDR)",
      override.aes = list(shape = 21, fill = "grey70", colour = "black", stroke = 0.2)
    )
  ) +
  
  scale_colour_manual(values = cols, guide = "none") +
  scale_y_continuous(breaks = seq_along(am_levels), labels = am_levels,
                     expand = expansion(mult = c(0.03, 0.06))) +
  labs(
    x = if (xvar %in% c("d","dz")) paste0("Cohen's ", xvar, " (Tumor vs Normal)") else "Effect (Tumor − Normal)",
    y = NULL
  ) +
  theme_bw(base_size = 11 * text_adjust) +
  theme(
    legend.position = "right",
    legend.text  = element_text(size = 10.5 * text_adjust, colour = "black"),
    legend.title = element_text(size = 10.5 * text_adjust, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y  = element_text(size = 10.5 * text_adjust, colour = "black"),
    axis.text.x  = element_text(size = 10.5 * text_adjust, colour = "black")
  )
p_dotplot
ggsave("Cohensd_dotplot_compare.pdf", p_dotplot, width = 8, height = 4.5, dpi = 300)

###################

cellchat_tumor <- readRDS("cellchat_am_tumor_computeCommunProb.rds")
cellchat_norm  <- readRDS("cellchat_am_normal_computeCommunProb.rds")

for (objname in c("cellchat_tumor","cellchat_norm")){
  obj <- get(objname)
  obj <- filterCommunication(obj, min.cells = 10)
  obj <- computeCommunProbPathway(obj)
  obj <- aggregateNet(obj)
  assign(objname, obj)
}

ensure_patterns <- function(obj, mode = c("outgoing","incoming"), k = 2){
  mode <- match.arg(mode)
  if (is.null(obj@netP$pattern[[mode]])) obj <- identifyCommunicationPatterns(obj, pattern = mode, k = k)
  obj
}

left_sankey_df <- function(obj, pattern = c("outgoing","incoming"), cutoff = 0.55){
  pattern <- match.arg(pattern)
  df <- obj@netP$pattern[[pattern]]$pattern$cell
  if (is.null(df)) stop("No pattern results for ", pattern, ". Run identifyCommunicationPatterns() first.")
  df$Contribution[df$Contribution < cutoff] <- 0
  df %>%
    mutate(AM = stringr::str_extract(as.character(CellGroup), "AM-\\d+")) %>%
    group_by(AM, Pattern) %>%
    summarise(Contribution = sum(Contribution), .groups = "drop") %>%
    filter(Contribution > 0)
}

plot_left_sankey <- function(df, am_cols, pat_palette, title = "", font.size = 4.8) {
  pats <- sort(unique(df$Pattern), key = function(x) as.integer(sub("Pattern ", "", x)))
  pat_cols <- setNames(pat_palette[seq_along(pats)], pats)
  plot.data <- df %>% transmute(CellGroup = AM, Pattern, Contribution)
  plot.data.long <- ggalluvial::to_lodes_form(plot.data, axes = 1:2, id = "connection")
  am_used  <- as.character(unique(plot.data.long$stratum[plot.data.long$x == "CellGroup"]))
  pat_used <- as.character(unique(plot.data.long$stratum[plot.data.long$x == "Pattern"]))
  col_map  <- c(setNames(am_cols[am_used], am_used),
                setNames(pat_cols[pat_used], pat_used))
  ggplot(plot.data.long,
         aes(x = factor(x, levels = c("CellGroup","Pattern")),
             y = Contribution, stratum = stratum, alluvium = connection,
             fill = stratum, label = stratum)) +
    geom_flow(width = 1/4, aes.flow = "backward") +
    geom_stratum(width = 1/4, size = 0.1, color = NA, alpha = 0.95) +
    geom_text(stat = "stratum", size = font.size) +
    scale_x_discrete(limits = c(), labels = c("", "")) +
    scale_fill_manual(values = scales::alpha(col_map, 0.90), drop = FALSE) +  
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title  = element_text(face = "bold",hjust = 0.5,vjust = -1,size = 16),
      axis.title = element_blank(), axis.text.y = element_blank(),
      panel.grid = element_blank(), panel.border = element_blank(),
      axis.ticks = element_blank(), axis.text.x = element_text(size = 10)
    ) +
    ggtitle(title)
}

cellchat_tumor <- ensure_patterns(cellchat_tumor, "incoming", k = 5)
cellchat_tumor <- ensure_patterns(cellchat_tumor, "outgoing", k = 2)
cellchat_norm  <- ensure_patterns(cellchat_norm,  "incoming", k = 2)
cellchat_norm  <- ensure_patterns(cellchat_norm,  "outgoing", k = 3)

cutoff <- 0.55
df_in_T  <- left_sankey_df(cellchat_tumor, "incoming", cutoff)
df_out_T <- left_sankey_df(cellchat_tumor, "outgoing", cutoff)
df_in_N  <- left_sankey_df(cellchat_norm,  "incoming", cutoff)
df_out_N <- left_sankey_df(cellchat_norm,  "outgoing", cutoff)

am_cols <- c("AM-1"="#ADC9CA","AM-2"="#A7A5D7","AM-3"="#FFAC5F","AM-4"="#7AA5D2","AM-5"="#E49898")
pat_palette <- c("#B8E3FC","#F7E6C3","#FBC1C1","#B2FDB5","#D7AFFA")

pdf("Tumor Outgoing.pdf",p1,w=6,height = 4.5)
p1 <- plot_left_sankey(df_out_T, am_cols, pat_palette,
                       sprintf("Tumor Outgoing",  length(unique(df_out_T$Pattern))))
p1
dev.off()

pdf("Normal Outgoing.pdf",p1,w=6,height = 4.5)
p2 <- plot_left_sankey(df_out_N, am_cols, pat_palette,
                       sprintf("Normal Outgoing", length(unique(df_out_N$Pattern))))
p2
dev.off()

pdf("Tumor Incoming.pdf",p1,w=6,height = 4.5)
p3 <- plot_left_sankey(df_in_T,  am_cols, pat_palette,
                       sprintf("Tumor Incoming",  length(unique(df_in_T$Pattern))))
p3
dev.off()

pdf("Normal Incoming.pdf",p1,w=6,height = 4.5)
p4 <- plot_left_sankey(df_in_N,  am_cols, pat_palette,
                       sprintf("Normal Incoming", length(unique(df_in_N$Pattern))))
p4
dev.off()


