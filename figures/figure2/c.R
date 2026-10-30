# AM distribution across spatial domains in 26 sections; saved RCTD estimates.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork); library(scales); library(grid)
})
setwd(am_data("sp/BANKSY"))
base_family <- "sans"
if (requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE)) {
  tryCatch(
    {
      ff <- sysfonts::font_files()
      sysfonts::font_add("Helvetica", ff$path[grepl("Helvetica|Arial", ff$family)][1])
      showtext::showtext_auto(); base_family <- "Helvetica"
    },
    error = function(e) NULL
  )
}

D <- read.csv(gzfile("data_domain_am_cd8_26.csv.gz"), stringsAsFactors = FALSE)
SAMPLES <- c("P7_LUAD2", "P10_LUAD", "P16_LUAD", "P21_LUAD", "P24_LUAD", "P25_LUAD")
D <- D[D$sample %in% SAMPLES, ]
D$sample <- factor(D$sample, levels = SAMPLES)
D$Y <- -D$array_row; D$X <- D$array_col

COMP <- c("Tumor", "AlvEpi", "Airway", "Stroma", "Vascular", "Lymphoid", "Myeloid")
CCOL <- c(
  Tumor = "#d65a5a", AlvEpi = "#8fbf7f", Airway = "#6cc2bb", Stroma = "#d9c27e",
  Vascular = "#7d9fc9", Lymphoid = "#a87cc0", Myeloid = "#e0a266"
)
D$compartment <- factor(D$compartment, levels = COMP)

norm01 <- function(v) {
  hi <- quantile(v, .99, na.rm = TRUE); pmin(v, hi) / (hi + 1e-9)
}
D <- D %>%
  group_by(sample) %>%
  mutate(AM4n = norm01(AM4), AM2n = norm01(AM2), CD8n = norm01(CD8)) %>%
  ungroup()

PT <- 1.15
GRADV <- c(0, 0.30, 0.62, 1)
RAMP <- list(
  AM2 = c("#ffffff", "#c6dbef", "#4292c6", "#08306b"),
  AM4 = c("#ffffff", "#fcbba1", "#ef3b2c", "#67000d"),
  CD8 = c("#ffffff", "#c7e9c0", "#41ab5d", "#00441b")
)

sq_lim <- function(d) {
  cx <- mean(range(d$X)); cy <- mean(range(d$Y))
  sp <- max(diff(range(d$X)), diff(range(d$Y))) * 1.04
  list(x = c(cx - sp / 2, cx + sp / 2), y = c(cy - sp / 2, cy + sp / 2))
}
base_theme <- theme_void(base_family = base_family) +
  theme(
    plot.title = element_text(hjust = .5, size = 11, face = "bold", margin = margin(b = 3)),
    plot.margin = margin(2, 2, 2, 2), legend.position = "none",
    panel.border = element_rect(color = "grey78", fill = NA, linewidth = .5)
  )

mk <- function(s, ch) {
  d <- D[D$sample == s, ]; L <- sq_lim(d)
  ttl <- if (s == SAMPLES[1]) switch(ch,
    dom = "Spatial domains",
    AM2 = "AM-2",
    AM4 = "AM-4",
    CD8 = "CD8 T"
  ) else ""
  if (ch == "dom") {
    g <- ggplot(d, aes(X, Y)) +
      geom_point(aes(color = compartment), size = PT, shape = 16) +
      scale_color_manual(values = CCOL, drop = FALSE, name = NULL)
  } else {
    cfg <- list(AM2 = list("AM2n", RAMP$AM2), AM4 = list("AM4n", RAMP$AM4), CD8 = list("CD8n", RAMP$CD8))[[ch]]
    d <- d[order(d[[cfg[[1]]]]), ]
    g <- ggplot(d, aes(X, Y)) +
      geom_point(aes(color = .data[[cfg[[1]]]]), size = PT, shape = 16) +
      scale_color_gradientn(
        colours = cfg[[2]], values = GRADV, limits = c(0, 1),
        name = switch(ch,
          AM2 = "AM-2",
          AM4 = "AM-4",
          CD8 = "CD8 T"
        )
      )
  }
  g + coord_cartesian(expand = FALSE) + labs(title = ttl) + base_theme +
    {
      if (ch == "dom") theme(axis.title.y = element_text(
        angle = 90, size = 10, face = "bold",
        color = "grey20", margin = margin(r = 4)
      )) else NULL
    } + {
      if (ch == "dom") ylab(sub("_.*", "", s)) else NULL
    }
}

plots <- list()
for (s in SAMPLES) for (ch in c("dom", "AM2", "AM4", "CD8")) plots[[paste(s, ch)]] <- mk(s, ch)
fig <- wrap_plots(plots, ncol = 4, byrow = TRUE)

DOMLAB <- c(
  Tumor = "Tumour", AlvEpi = "Alveolar epithelium", Airway = "Airway", Stroma = "Stroma",
  Vascular = "Vascular", Lymphoid = "Lymphoid", Myeloid = "Myeloid"
)

leg_dom <- cowplot::get_legend(
  ggplot(D, aes(X, Y, fill = compartment)) +
    geom_tile() +
    scale_fill_manual(values = CCOL, labels = DOMLAB[COMP], name = NULL) +
    guides(fill = guide_legend(
      nrow = 1, keywidth = unit(12, "pt"), keyheight = unit(12, "pt"),
      override.aes = list(color = "white")
    )) +
    theme(
      legend.position = "bottom", legend.text = element_text(size = 9.5),
      legend.key = element_blank(), legend.spacing.x = unit(3, "pt")
    )
)

mk_cbar <- function(cols, nm) cowplot::get_legend(
  ggplot(data.frame(x = 1, v = c(0, 1)), aes(x, x, color = v)) +
    geom_point() +
    scale_color_gradientn(
      colours = cols, values = GRADV, limits = c(0, 1),
      breaks = c(0, 1), labels = c("low", "high"), name = nm
    ) +
    guides(color = guide_colorbar(
      barwidth = unit(42, "pt"), barheight = unit(7, "pt"),
      title.position = "top", title.hjust = .5, ticks = FALSE, frame.colour = "grey60",
      frame.linewidth = .3
    )) +
    theme(
      legend.position = "bottom", legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8.5), legend.margin = margin(0, 4, 0, 4)
    )
)
mk_dlegend <- function(cols, nm, nrow = 1) {
  sp <- if (utils::packageVersion("ggplot2") >= "3.5.0")
    theme(legend.key.spacing.x = unit(8, "pt"))
  else
    theme(legend.spacing.x = unit(8, "pt"))

  cowplot::get_legend(
    ggplot(
      data.frame(x = 1, v = factor(names(cols), levels = names(cols))),
      aes(x, x, fill = v)
    ) +
      geom_tile() +
      scale_fill_manual(values = cols, name = nm) +
      guides(fill = guide_legend(
        nrow = nrow, byrow = TRUE,
        title.position = "top", title.hjust = .5,
        keywidth = unit(9, "pt"), keyheight = unit(9, "pt")
      )) +
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10, margin = margin(l = 3, r = 12)),
        legend.margin = margin(0, 4, 0, 4)
      ) +
      sp
  )
}

leg_dom <- mk_dlegend(CCOL, "Compartment", nrow = 1)
leg_am2 <- mk_cbar(RAMP$AM2, "AM-2"); leg_am4 <- mk_cbar(RAMP$AM4, "AM-4"); leg_cd8 <- mk_cbar(RAMP$CD8, "CD8 T")

legend_row <- patchwork::wrap_elements(leg_dom) |
  patchwork::wrap_elements(leg_am2) |
  patchwork::wrap_elements(leg_am4) |
  patchwork::wrap_elements(leg_cd8)
legend_row <- legend_row + plot_layout(widths = c(7, 1, 1, 1))
final <- fig / legend_row + plot_layout(heights = c(1, .07))
if (interactive()) print(final)
CELL <- 2.5; W <- 4 * CELL; H <- 6 * CELL * 1.09 + 0.2
if (interactive()) print(final)
ggsave(am_output("Figure2", "Figure2c.pdf"), final, width = W, height = H, device = cairo_pdf, limitsize = FALSE)
cat("saved figs/banksy_domain_am_cd8.{png,pdf}\n")
