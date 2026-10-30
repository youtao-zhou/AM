# PCoA centroids from the historical 13-subcluster analysis.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
source(am_code("R/figure1_data.R"))
library(ggplot2)
library(ggrepel)
AM_ORDER <- paste0("AM-", 1:5)
AM2_COLOR <- "#326FAC"
AM4_COLOR <- "#C44D3F"
GREY_COLOR <- "#9CA3AA"

centroids <- read_source("F1h_MDS_centroids")
rownames(centroids) <- centroids$AM
plot_data <- centroids
plot_data$AM <- factor(plot_data$AM, levels = AM_ORDER)
plot_data$focus <- ifelse(plot_data$AM == "AM-2", "AM-2",
  ifelse(plot_data$AM == "AM-4", "AM-4", "Other")
)

p <- ggplot(plot_data, aes(Coordinate1, Coordinate2)) +
  geom_hline(yintercept = 0, color = "#D8DCE0", linewidth = 0.42) +
  geom_vline(xintercept = 0, color = "#D8DCE0", linewidth = 0.42) +
  geom_point(
    data = subset(plot_data, AM %in% c("AM-2", "AM-4")),
    aes(color = focus), size = 9, alpha = 0.12, stroke = 0, show.legend = FALSE
  ) +
  geom_point(aes(fill = focus, size = focus),
    shape = 21, color = "white",
    stroke = 0.8, show.legend = FALSE
  ) +
  geom_text_repel(
    aes(label = AM),
    size = 4,
    fontface = ifelse(as.character(plot_data$AM) %in% c("AM-2", "AM-4"), "bold", "plain"),
    color = "#252A30",
    box.padding = 0.45, point.padding = 0.35,
    min.segment.length = Inf, seed = 1, max.overlaps = Inf
  ) +
  scale_color_manual(values = c("AM-2" = AM2_COLOR, "AM-4" = AM4_COLOR, "Other" = GREY_COLOR), guide = "none") +
  scale_fill_manual(values = c("AM-2" = AM2_COLOR, "AM-4" = AM4_COLOR, "Other" = GREY_COLOR), guide = "none") +
  scale_size_manual(values = c("AM-2" = 4.7, "AM-4" = 4.7, "Other" = 3.6), guide = "none") +
  labs(x = "Coordinate 1", y = "Coordinate 2") +
  theme_classic(base_size = 11) +
  theme(
    axis.title  = element_text(size = 12, color = "black"),
    axis.text   = element_text(size = 12, color = "black"),
    axis.line   = element_line(color = "black", linewidth = 0.2),
    plot.margin = margin(9, 14, 8, 9)
  )
if (interactive()) print(p)
ggplot2::ggsave(file.path(OUT_ROOT, "Figure1h.pdf"), p, width = 5.41666667, height = 3.77777778, device = grDevices::cairo_pdf, bg = "white")
ggplot2::ggsave(file.path(OUT_ROOT, "Figure1h.png"), p, width = 5.41666667, height = 3.77777778, dpi = 180, bg = "white")
write.csv(plot_data, file.path(DERIVED_ROOT, "h_plotted_values.csv"), row.names = FALSE)
