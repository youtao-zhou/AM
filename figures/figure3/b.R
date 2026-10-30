# PROGENy activity standardized across 14 pathways within each AM state.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({library(Matrix); library(pheatmap)})

output_dir <- am_output("Figure3", "b")
am_order <- c("AM-1", "AM-2", "AM-3", "AM-4")
pathway_order <- c("NFkB", "JAK-STAT", "TNFa", "EGFR", "MAPK", "PI3K", "p53",
                   "Hypoxia", "VEGF", "TGFb", "WNT", "Trail", "Androgen", "Estrogen")
font_size <- 12
font_family <- "sans"
figure_width <- 8
figure_height <- 3.4
column_angle <- "90"
heat_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
heat_limits <- c(-2, 2)

row_standardize <- TRUE
show_numbers <- FALSE
column_gaps <- NULL

dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
f <- function(n) file.path(output_dir,n)
means <- as.matrix(read.csv(am_data("HTAN/am/PROGENy/PROGENy_AM1to4_mean_scores.csv"),row.names=1,check.names=FALSE))
means <- means[am_order,pathway_order,drop=FALSE]

plot_matrix <- if (row_standardize) t(scale(t(means))) else means
stopifnot(identical(rownames(plot_matrix), am_order), !"AM-5" %in% rownames(plot_matrix))
write.csv(plot_matrix, f("PROGENy_AM1to4_heatmap_values.csv"))
breaks <- seq(heat_limits[1], heat_limits[2], length.out = length(heat_colors) + 1)
p <- pheatmap(plot_matrix, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
               color = heat_colors, breaks = breaks, border_color = NA,
               fontsize = font_size, fontfamily = font_family,
               angle_col = column_angle, gaps_col = column_gaps,
               display_numbers = show_numbers, main = "PROGENy", silent = TRUE)

if (interactive()) {grid::grid.newpage(); grid::grid.draw(p$gtable)}
grDevices::cairo_pdf(f("PROGENy_AM1to4.pdf"), width = figure_width, height = figure_height)
grid::grid.newpage(); grid::grid.draw(p$gtable); dev.off()
png(f("PROGENy_AM1to4.png"), width = figure_width, height = figure_height, units = "in", res = 300, type = "cairo")
grid::grid.newpage(); grid::grid.draw(p$gtable); dev.off()
