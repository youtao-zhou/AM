# CellChat probabilities normalized to the maximum within each of nine T-cell groups.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(ggplot2)
library(patchwork)
library(MetBrewer)

setwd(am_data("HTAN/am/cellchat"))

d <- fread("data/AM2_CD4_T_all_pairs.csv.gz")

main_T <- c("CD4 Tn/Tcm", "CD4 Tem", "CD8 Tm", "CD8 CTL", "Treg")
other_T <- c("CD4 Tfh", "CD8 Tex", "T-IFN", "Cycling T")

sel <- data.table(
  direction = c(rep("T to AM", 3), rep("AM to T", 4)),
  interaction_name = c(
    "CD40LG_CD40", "IFNG_IFNGR1_IFNGR2", "CCL5_CCR1",
    "HLA-B_CD8A", "HLA-DRA_CD4", "CD86_CD28", "CXCL16_CXCR6"
  )
)

d[, rel := {
  m <- max(prob); if (m > 0) prob / m else 0
}, by = .(direction, interaction_name)]
d <- merge(d, sel, by = c("direction", "interaction_name"))
rows <- lapply(
  split(sel$interaction_name, sel$direction),
  function(z) unique(d[interaction_name %in% z][match(z, interaction_name), interaction_name_2])
)
lim <- max(d$prob)

pal <- met.brewer("Hiroshige", 100, direction = -1)
p_lv <- c("P ≤ 0.01", "0.01 < P ≤ 0.05", "P > 0.05")
fnt <- "Noto Sans"
lw <- 0.3

dot_panel <- function(dir, states, relative = FALSE) {
  pairs <- rev(rows[[dir]])
  k <- length(states); n <- length(pairs)

  x <- d[direction == dir & T_state %in% states]
  x[, xi := match(T_state, states)]
  x[, yi := match(interaction_name_2, pairs)]
  x[, pgrp := factor(fcase(p_plus1 <= 0.01, p_lv[1], p_plus1 <= 0.05, p_lv[2], default = p_lv[3]), p_lv)]
  x[, a := ifelse(paired_support, 1, 0.3)]
  x[, val := if (relative) rel else prob]

  if (dir == "T to AM") {
    strip <- data.table(xmin = 1:k - 0.5, xmax = 1:k + 0.5, lab = states)
    bottom <- rep("AM-2.CD4", k)
  } else {
    strip <- data.table(xmin = 0.5, xmax = k + 0.5, lab = "AM-2.CD4")
    bottom <- states
  }

  top <- n + 0.5; sh <- 0.8
  yb <- 0.5 - 0.4
  cl <- if (relative) 1 else lim

  ggplot(x, aes(xi, yi)) +
    geom_vline(xintercept = 1:k, colour = "grey75", linewidth = 0.25, linetype = "dotted") +
    geom_point(data = x[prob > 0], aes(size = pgrp, colour = val, alpha = a)) +
    geom_point(data = x[prob == 0], shape = 4, size = 1.6, colour = "grey75") +
    geom_rect(
      data = strip, aes(xmin = xmin, xmax = xmax, ymin = top, ymax = top + sh),
      inherit.aes = FALSE, fill = "white", colour = "black", linewidth = lw
    ) +
    geom_text(
      data = strip, aes(x = (xmin + xmax) / 2, y = top + sh / 2, label = lab),
      inherit.aes = FALSE, size = 10 / .pt, family = fnt
    ) +
    annotate("text", x = 1:k, y = yb, label = bottom, size = 9 / .pt, family = fnt) +
    annotate("text",
      x = 0.4, y = c(top + sh / 2, yb), label = c("Sender", "Receiver"),
      hjust = 1, size = 10 / .pt, family = fnt
    ) +
    coord_cartesian(xlim = c(0.5, k + 0.5), ylim = c(0.5, n + 0.5), expand = FALSE, clip = "off") +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = 1:n, labels = pairs) +
    scale_size_manual(values = setNames(c(4.8, 3.5, 2), p_lv), drop = FALSE, name = "Permutation P") +
    scale_colour_gradientn(
      colours = pal, limits = c(0, cl), breaks = c(0, cl / 2, cl),
      labels = function(z) formatC(z, format = "f", digits = if (relative) 1 else 3),
      name = if (relative) "Relative Prob" else "Commun Prob"
    ) +
    scale_alpha_identity() +
    labs(x = NULL, y = NULL) +
    guides(
      size = guide_legend(
        order = 1, nrow = 1,
        override.aes = list(size = c(4.8, 3.5, 2), alpha = 1, colour = "black")
      ),
      colour = guide_colourbar(
        order = 2, barwidth = unit(3.5, "cm"), barheight = unit(0.3, "cm"),
        theme = theme(
          legend.title = element_text(size = 9, margin = margin(r = 10, b = 6)),
          legend.text = element_text(size = 8, margin = margin(t = 3))
        )
      )
    ) +
    theme_minimal(base_size = 10, base_family = fnt) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = lw),
      panel.grid.major.y = element_line(colour = "grey75", linewidth = 0.25, linetype = "dotted"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(colour = "black", size = 9),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      plot.margin = margin(30, 10, 22, 8)
    )
}

two_way <- function(states, relative = FALSE) {
  (dot_panel("T to AM", states, relative) / dot_panel("AM to T", states, relative)) +
    plot_layout(
      guides = "collect",
      heights = c(length(rows[["T to AM"]]), length(rows[["AM to T"]]))
    ) &
    theme(
      legend.position = "bottom", legend.box = "horizontal",
      legend.box.spacing = unit(24, "pt")
    )
}

p <- two_way(main_T, relative = TRUE)
p2 <- two_way(other_T, relative = TRUE)

if (interactive()) print(p)
if (interactive()) print(p2)
ggsave(am_output("Figure3", "AM2_CD4_T_cellchat_main.pdf"), p, width = 8.30, height = 4.71, device = cairo_pdf, bg = "white")
ggsave(am_output("Figure3", "AM2_CD4_T_cellchat_supp.pdf"), p2, width = 8.30, height = 4.71, device = cairo_pdf, bg = "white")
