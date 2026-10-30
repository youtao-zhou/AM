# Representative AM-2/CD8 spatial windows; scale bar: 200 pixels.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)

library(ggplot2)
library(patchwork)

setwd(am_data("tsa"))

font <- "Noto Sans"
ink <- "#20292F"

cols <- c(
  "AM-2" = "#E25563",
  "Other macrophages" = "#3E97D1",
  "CD8" = "#6E4FB0",
  "Tumor" = "#D9A43A"
)

labels <- c(
  "AM-2" = "AM-2",
  "Other macrophages" = "Other macrophages",
  "CD8" = "CD8 T cells",
  "Tumor" = "Tumor cells"
)

sizes <- c(
  "AM-2" = 2,
  "Other macrophages" = 1.9,
  "CD8" = 1.9,
  "Tumor" = 1.05
)

scale_bar_px <- 200
fig_width <- 4.97
fig_height <- 3.50

cells <- fread(
  "data/AM2_CD8_pCR_NR_zoom_cells.csv",
  colClasses = list(character = "patient_id")
)

wins <- fread(
  "data/AM2_CD8_pCR_NR_zoom_windows.csv",
  colClasses = list(character = "patient_id")
)

stopifnot(
  nrow(wins) == 2,
  setequal(wins$response, c("pCR", "NR"))
)

expected <- melt(
  wins[, .(
    patient_id,
    `AM-2` = n_AM2,
    `Other macrophages` = n_other_macrophages,
    CD8 = n_CD8,
    Tumor = n_tumor_label
  )],
  id.vars = "patient_id",
  variable.name = "cell_type",
  value.name = "expected",
  variable.factor = FALSE
)

observed <- cells[
  cell_type != "Other cells",
  .N,
  by = .(patient_id, cell_type)
]

check <- merge(
  expected,
  observed,
  by = c("patient_id", "cell_type"),
  all.x = TRUE
)

check[is.na(N), N := 0L]
stopifnot(all(check$N == check$expected))

cells[, cell_type := factor(
  cell_type,
  levels = c("Other cells", "Tumor", "Other macrophages", "CD8", "AM-2")
)]

panel <- function(response_name) {
  w <- wins[response == response_name]
  stopifnot(nrow(w) == 1)

  pid <- w$patient_id
  win <- w$width_px
  z <- copy(cells[patient_id == pid])

  z[, `:=`(
    xp = x_in_window,
    yp = win - y_in_window
  )]

  ggplot() +
    geom_point(
      data = z[cell_type == "Other cells"],
      aes(xp, yp),
      colour = "#C3CAD0",
      size = 0.75,
      shape = 16,
      show.legend = FALSE
    ) +
    geom_point(
      data = z[cell_type != "Other cells"][order(cell_type)],
      aes(xp, yp, colour = cell_type, size = cell_type),
      shape = 16,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = cols,
      guide = "none"
    ) +
    scale_size_manual(
      values = sizes,
      guide = "none"
    ) +
    annotate(
      "segment",
      x = win * .04,
      xend = win * .04 + scale_bar_px,
      y = -win * .035,
      yend = -win * .035,
      linewidth = .35,
      colour = ink
    ) +
    annotate(
      "text",
      x = win * .04 + scale_bar_px / 2,
      y = -win * .085,
      label = paste(scale_bar_px, "px"),
      size = 2.9,
      family = font,
      colour = ink
    ) +
    annotate(
      "rect",
      xmin = 0,
      xmax = win,
      ymin = 0,
      ymax = win,
      fill = NA,
      colour = "grey70",
      linewidth = .3
    ) +
    coord_equal(
      xlim = c(0, win),
      ylim = c(-win * .11, win),
      expand = FALSE,
      clip = "off"
    ) +
    labs(title = response_name) +
    theme_void(base_family = font, base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 11.5,
        colour = ink,
        hjust = .5,
        margin = margin(b = 3)
      ),
      legend.position = "none",
      plot.margin = margin(4, 6, 2, 6)
    )
}

legend_types <- c("AM-2", "Other macrophages", "CD8", "Tumor")

stopifnot(
  all(legend_types %in% names(cols)),
  all(legend_types %in% names(labels))
)

legend_x <- c(0.06, 0.53, 0.06, 0.53)
legend_y <- c(0.73, 0.73, 0.27, 0.27)

legend_items <- lapply(seq_along(legend_types), function(i) {
  type <- legend_types[i]

  grid::grobTree(
    grid::circleGrob(
      x = grid::unit(legend_x[i], "npc"),
      y = grid::unit(legend_y[i], "npc"),
      r = grid::unit(1.8, "mm"),
      gp = grid::gpar(
        fill = unname(cols[type]),
        col = NA
      )
    ),
    grid::textGrob(
      label = unname(labels[type]),
      x = grid::unit(legend_x[i], "npc") + grid::unit(3, "mm"),
      y = grid::unit(legend_y[i], "npc"),
      just = c("left", "centre"),
      gp = grid::gpar(
        fontfamily = font,
        fontsize = 11.5,
        col = ink
      )
    )
  )
})

legend_grob <- do.call(grid::grobTree, legend_items)

legend_panel <- patchwork::wrap_elements(
  full = legend_grob,
  clip = FALSE
)

panels <- panel("pCR") | panel("NR")

p <- (panels / legend_panel) +
  plot_layout(heights = c(1, 0.20))

if (interactive()) print(p)

ggsave(
  am_output("Figure3", "AM2_CD8_pCR_NR_zoom.pdf"),
  plot = p,
  width = fig_width,
  height = fig_height,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)
