# Marker expression across AM states, using the saved dot-plot matrix.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
source(am_code("R/figure1_data.R"))
library(ggplot2)
library(patchwork)
groups <- paste0("AM-", 1:5)

marker_def <- data.frame(
  block = c(
    rep("Identity", 4), rep("AM-1", 2), rep("AM-2", 2),
    "AM-3", rep("AM-4", 3), "AM-5"
  ),
  gene = c(
    "CD68", "PPARG", "MARCO", "C1QA",
    "FABP4", "LDHB", "CD4", "C5AR1", "SPP1",
    "CD80", "SPIDR", "RAD51B", "IGHA1"
  ),
  stringsAsFactors = FALSE
)
marker_def$target_group <- ifelse(marker_def$block == "Identity", "All AM", marker_def$block)

gene_order <- marker_def$gene
n_genes <- length(gene_order)
gene_y <- setNames(n_genes - seq_along(gene_order) + 1, gene_order)

d <- expand.grid(gene = gene_order, group = groups, stringsAsFactors = FALSE)
src <- read_source("F1d_markers")
ix <- match(paste(d$gene, d$group), paste(src$gene, src$group))
stopifnot(!anyNA(ix))
d$pct_expressing <- src$pct_expressing[ix]
d$relative_expression <- src$relative_expression[ix]
d$target_group <- marker_def$target_group[match(d$gene, marker_def$gene)]
d$is_target <- d$target_group == "All AM" | d$target_group == d$group
d$y <- unname(gene_y[d$gene])
d$group <- factor(d$group, levels = groups)

annotation_colors <- c(
  "Identity" = "#4C78A8", "AM-1" = "#00A6A6", "AM-2" = "#F2B701",
  "AM-3" = "#B279A2", "AM-4" = "#59A14F", "AM-5" = "#E45756"
)
bettervis_colors <- c("#1289B2", "#FCF6AF", "#E53849")

p_dot <- ggplot(d, aes(group, y)) +
  geom_tile(width = 0.94, height = 0.94, fill = "#F5F6F7", color = "white", linewidth = 0.55) +
  geom_tile(
    data = subset(d, is_target & target_group != "All AM"),
    width = 0.94, height = 0.94, fill = "#FFF1CF", color = "white", linewidth = 0.55
  ) +
  geom_point(
    data = subset(d, !is_target),
    aes(size = pct_expressing, fill = relative_expression),
    shape = 21, color = "#777777", stroke = 0.24
  ) +
  geom_point(
    data = subset(d, is_target),
    aes(size = pct_expressing, fill = relative_expression),
    shape = 21, color = "#202020", stroke = 0.2
  ) +
  scale_x_discrete(position = "top", drop = FALSE) +
  scale_y_continuous(
    breaks = unname(gene_y[gene_order]), labels = gene_order,
    limits = c(0.5, n_genes + 0.5), expand = c(0, 0)
  ) +
  scale_size_area(
    name = "Percent", max_size = 8.5, limits = c(0, 100), breaks = c(10, 25, 50, 75, 100),
    guide = guide_legend(
      order = 1, title.position = "top", nrow = 1,
      override.aes = list(fill = "#BDBDBD", color = "#404040", stroke = 0.3)
    )
  ) +
  scale_fill_gradientn(
    name = "Expression", colours = bettervis_colors,
    values = scales::rescale(c(0, 0.5, 1)), limits = c(0, 1), breaks = c(0, 0.5, 1),
    oob = scales::squish,
    guide = guide_colorbar(
      order = 2, title.position = "top", direction = "horizontal",
      barwidth = grid::unit(30, "mm"), barheight = grid::unit(3.8, "mm")
    )
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 10.5) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(face = "bold", size = 11.2, color = "black", margin = margin(b = 5)),
    axis.text.y = element_text(face = "italic", size = 10.3, color = "black"),
    axis.ticks = element_blank(),
    legend.position = "bottom", legend.box = "horizontal", legend.box.just = "center",
    legend.title = element_text(face = "bold", size = 8.8),
    legend.text = element_text(size = 8.2),
    legend.spacing.x = grid::unit(2.5, "mm"),
    legend.box.spacing = grid::unit(1, "mm"),
    plot.margin = margin(5, 7, 4, 1)
  )

anno_df <- marker_def
anno_df$y <- unname(gene_y[anno_df$gene])
p_anno <- ggplot(anno_df, aes(1, y, fill = block)) +
  geom_tile(width = 0.92, height = 0.94, color = "white", linewidth = 0.55) +
  scale_fill_manual(
    name = "AM", values = annotation_colors, breaks = names(annotation_colors),
    guide = guide_legend(
      order = 3, title.position = "top", nrow = 2, byrow = TRUE,
      override.aes = list(color = NA)
    )
  ) +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, n_genes + 0.5), expand = c(0, 0)) +
  theme_void(base_size = 10.5) +
  theme(
    legend.position = "bottom", legend.box = "horizontal", legend.box.just = "center",
    legend.title = element_text(face = "bold", size = 8.8),
    legend.text = element_text(size = 8.2),
    legend.spacing.x = grid::unit(2.5, "mm"),
    legend.box.spacing = grid::unit(1, "mm"),
    plot.margin = margin(5, 0, 4, 4)
  )

final_plot <- (p_anno + p_dot +
  patchwork::plot_layout(widths = c(0.18, 6.3), guides = "collect")) &
  theme(
    legend.position = "bottom", legend.box = "horizontal", legend.box.just = "center",
    legend.spacing.x = grid::unit(2.5, "mm"), legend.box.spacing = grid::unit(1, "mm")
  )
if (interactive()) print(final_plot)
ggplot2::ggsave(file.path(OUT_ROOT, "Figure1d.pdf"), final_plot, width = 6.69444444, height = 5.75, device = grDevices::cairo_pdf, bg = "white")
ggplot2::ggsave(file.path(OUT_ROOT, "Figure1d.png"), final_plot, width = 6.69444444, height = 5.75, dpi = 180, bg = "white")
write.csv(d, file.path(DERIVED_ROOT, "d_plotted_values.csv"), row.names = FALSE)
