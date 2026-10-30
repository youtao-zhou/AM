# Representative AM-4 spatial windows; scale bar: 200 pixels.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(ggplot2)
library(patchwork)

setwd(am_data("tsa"))

font <- "Noto Sans"
ink <- "#20292F"
cols <- c("AM-4" = "#2E44B8", "Other macrophages" = "#5DBFA0", "CD8" = "#B23A8A", "Tumor" = "#D9A43A")
labels <- c("AM-4" = "AM-4", "Other macrophages" = "Other macrophages", "CD8" = "CD8 T cells", "Tumor" = "Tumor cells (HALO label)")
sizes <- c("AM-4" = 2, "Other macrophages" = 1.9, "CD8" = 1.9, "Tumor" = 1.05)
win <- 1000
scale_bar_px <- 200
fig_width <- 4.97
fig_height <- 3.50

cells <- fread("data/AM4_pCR_NR_zoom_cells.csv", colClasses = c(patient_id = "character"))
wins <- fread("data/AM4_pCR_NR_zoom_windows.csv", colClasses = c(patient_id = "character"))

panel <- function(pid, title, show_legend) {
  z <- cells[patient_id == pid]
  dummy <- data.table(xp = -1e4, yp = -1e4, cell_type = names(cols))
  p <- ggplot() +
    geom_point(data = z[cell_type == "Other cells"], aes(xp, yp), colour = "#D3D8DC", size = 0.6, shape = 16) +
    geom_point(
      data = z[cell_type != "Other cells"][order(match(cell_type, c("Tumor", "CD8", "Other macrophages", "AM-4")))],
      aes(xp, yp, colour = cell_type, size = cell_type), shape = 16
    ) +
    geom_point(data = dummy, aes(xp, yp, colour = cell_type, size = cell_type), shape = 16) +
    scale_colour_manual(values = cols, breaks = names(cols), labels = labels[names(cols)], drop = FALSE) +
    scale_size_manual(values = sizes, guide = "none") +
    annotate("segment", x = 40, xend = 40 + scale_bar_px, y = -35, yend = -35, linewidth = 0.7, colour = ink) +
    annotate("text", x = 40 + scale_bar_px / 2, y = -85, label = paste(scale_bar_px, "px"), size = 2.9, family = font, colour = ink) +
    annotate("rect", xmin = 0, xmax = win, ymin = 0, ymax = win, fill = NA, colour = "grey70", linewidth = 0.3) +
    coord_equal(xlim = c(0, win), ylim = c(-110, win), expand = FALSE, clip = "off") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2.8))) +
    labs(title = title, colour = NULL) +
    theme_void(base_family = font, base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5), legend.position = "bottom", legend.text = element_text(size = 8.5))
  if (!show_legend) p <- p + guides(colour = "none")
  p
}

p <- (panel("2007957", "pCR", FALSE) + theme(plot.margin = margin(2, 8, 2, 2)) |
  panel("2009455", "NR", TRUE) + theme(plot.margin = margin(2, 2, 2, 8))) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

if (interactive()) print(p)
ggsave(am_output("Figure4", "AM4_pCR_NR_zoom.pdf"), p, width = 3.65, height = 2.91, device = cairo_pdf, bg = "white")
