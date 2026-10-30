# AM-2/T-cell associations using within-study ranks and saved permutation statistics.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(ggplot2)
library(patchwork)

base_dir <- am_data("HTAN/am/cytosig")
fig_dir <- file.path(base_dir, "figure", "AM2CD4_killing_T")

font_family <- "Noto Sans"
cols <- c(Low = "#70CDBE", High = "#AC99D2")

n_col <- 2
fig_width <- 6.18
panel_height <- 4.81

stats <- fread(file.path(base_dir, "AM2CD4_killing_T_eight_panel_statistics.csv"))
points <- fread(file.path(base_dir, "AM2CD4_killing_T_eight_panel_points.csv"))
points[, group := factor(group, levels = c("Low", "High"))]

plot_panel <- function(k) {
  s <- stats[panel == k]
  d <- points[panel == k]

  stopifnot(nrow(s) == 1, nrow(d) == s$n_within)

  p_label <- if (s$p_within < 0.001) {
    "p < 0.001"
  } else {
    sprintf("p = %.3f", s$p_within)
  }

  label <- sprintf(
    "r = %.3f\n%s",
    s$rho_within, p_label
  )

  scatter <- ggplot(d, aes(x_rank_in_study, y_rank_in_study)) +
    geom_point(aes(color = group), size = 2) +
    geom_vline(xintercept = 0.5, linetype = "dashed") +
    geom_smooth(
      method = "lm", formula = y ~ x,
      color = "grey40", se = FALSE
    ) +
    annotate(
      "text",
      x = 0.03, y = 0.93, label = label,
      hjust = 0, vjust = 1, size = 3.1, family = font_family
    ) +
    annotate(
      "text",
      x = c(0.05, 0.95), y = 1.04,
      label = names(cols), color = unname(cols),
      hjust = c(0, 1), size = 4, family = font_family
    ) +
    scale_color_manual(values = cols) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_y_continuous(limits = c(0, 1.07), breaks = c(0, 0.5, 1)) +
    labs(
      x = paste0("AM-2.CD4 ", s$response, " score"),
      y = s$ylab
    ) +
    theme_test(base_family = font_family) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10)
    )

  density <- ggplot(d, aes(y_rank_in_study, fill = group, color = group)) +
    geom_density(alpha = 0.35, linewidth = 0.4) +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    scale_x_continuous(limits = c(0, 1.07)) +
    coord_flip() +
    theme_void() +
    theme(legend.position = "none")

  scatter + density + plot_layout(widths = c(5, 1))
}

panels_ctl <- stats[grepl("^CD8 CTL", axis)][order(panel), panel]
panels_other <- stats[!grepl("^CD8 CTL", axis)][order(panel), panel]

save_figure <- function(panels, name) {
  figure <- wrap_plots(lapply(panels, plot_panel), ncol = n_col)
  ggsave(
    file.path(am_output("Figure3", ""), paste0(name, ".pdf")), figure,
    width = fig_width, height = panel_height,
    device = cairo_pdf, bg = "white"
  )
  figure
}

figure_ctl <- save_figure(panels_ctl, "AM2CD4_CytoSig_T_correlation_CD8_CTL")
figure_other <- save_figure(panels_other, "AM2CD4_CytoSig_T_correlation_other_T")
if (interactive()) print(figure_ctl)
if (interactive()) print(figure_other)
