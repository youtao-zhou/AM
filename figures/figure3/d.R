# T-cell UMAP using the saved coordinates and labels for 100,000 cells.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(ggplot2)
library(ggrepel)

setwd(am_data("HTAN/am/tcell"))

groups <- c(
  "CD4 Tn/Tcm", "CD4 Tem", "CD4 Tfh", "CD8 Tm", "CD8 CTL",
  "CD8 Tex", "Treg", "T-IFN", "Cycling T"
)
cols <- setNames(c(
  "#409B96", "#83B76F", "#CE9850", "#5286BF", "#6B60A4",
  "#BF6371", "#A86196", "#A2AA41", "#D87946"
), groups)

coords <- read.csv(gzfile("data/Tcell_full100k_plot_coordinates.csv.gz"))
stopifnot(nrow(coords) == 100000L, !anyNA(coords$T_short))
d <- data.frame(UMAP1 = coords$UMAP1, UMAP2 = coords$UMAP2, group = factor(coords$T_short, groups))

set.seed(58)
d <- d[sample(nrow(d)), ]

lab <- do.call(rbind, lapply(split(d, d$group), function(x) {
  x[which.min((x$UMAP1 - median(x$UMAP1))^2 + (x$UMAP2 - median(x$UMAP2))^2), ]
}))

r <- max(diff(range(d$UMAP1)), diff(range(d$UMAP2)))
x0 <- -12
y0 <- -7.5
len <- 0.15 * r
arw <- arrow(length = unit(5, "pt"), type = "closed")

p <- ggplot(d, aes(UMAP1, UMAP2, color = group)) +
  geom_point(size = 0.28, alpha = 0.82, stroke = 0) +
  geom_label_repel(
    data = lab, aes(label = group), seed = 58,
    size = 3.3, family = "Noto Sans", fontface = "bold",
    fill = alpha("white", 0.9), label.size = 0,
    label.padding = unit(0.13, "lines"), box.padding = 0.35,
    min.segment.length = 0.1, segment.color = "#888888",
    segment.size = 0.25, max.overlaps = Inf, show.legend = FALSE
  ) +
  annotate("segment", x = x0, xend = x0 + len, y = y0, yend = y0, linewidth = 0.4, arrow = arw) +
  annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + len, linewidth = 0.4, arrow = arw) +
  annotate("text",
    x = x0 + len / 2, y = y0 - 0.025 * r, label = "UMAP1",
    vjust = 1, size = 3.5, family = "Noto Sans"
  ) +
  annotate("text",
    x = x0 - 0.025 * r, y = y0 + len / 2, label = "UMAP2",
    angle = 90, vjust = 0, size = 3.5, family = "Noto Sans"
  ) +
  scale_color_manual(values = cols, name = NULL) +
  coord_equal(
    xlim = c(x0 - 0.06 * r, max(d$UMAP1)),
    ylim = c(y0 - 0.06 * r, max(d$UMAP2)), clip = "off"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_void(base_size = 11, base_family = "Noto Sans") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.key.height = unit(16, "pt"),
    plot.margin = margin(10, 10, 16, 16)
  )

if (interactive()) print(p)
ggsave(am_output("Figure3", "Tcell_umap.pdf"), p, width = 6.33, height = 4.51, device = cairo_pdf, bg = "white")
