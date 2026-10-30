# Functional modules from the historical 13-subcluster analysis.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
source(am_code("R/figure1_data.R"))
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
AM_ORDER <- paste0("AM-", 1:5)
src <- read_source("F1i_modules")
module_z <- t(as.matrix(src[, AM_ORDER]))
colnames(module_z) <- src[[1]]
MODULE_ORDER <- colnames(module_z)

MODULE_LABELS <- c(
  "MHC-II antigen presentation"       = "MHC-II antigen presentation",
  "T-cell regulation"                 = "T-cell regulation",
  "Immune-cell recruitment"           = "Immune-cell recruitment",
  "Cell-matrix adhesion"              = "Cell-matrix adhesion",
  "Tissue morphogenesis / remodeling" = "Tissue remodeling",
  "Vesicle trafficking / autophagy"   = "Intracellular trafficking",
  "Stress adaptation / wound repair"  = "Tissue repair"
)

heatmap_long <- as.data.frame(module_z, check.names = FALSE) %>%
  mutate(AM = rownames(.)) %>%
  pivot_longer(all_of(MODULE_ORDER), names_to = "module", values_to = "z") %>%
  mutate(
    AM = factor(AM, levels = AM_ORDER),
    module = factor(module,
      levels = rev(MODULE_ORDER),
      labels = unname(MODULE_LABELS[rev(MODULE_ORDER)])
    )
  )

HEAT_COLORS <- c(
  "#4575B4", "#74ADD1", "#ABD9E9", "#FFFFBF",
  "#FFFFBF", "#FEE090", "#F46D43", "#D73027"
)

HEAT_VALUES <- c(0.000, 0.300, 0.375, 0.450, 0.550, 0.625, 0.700, 1.000)

p_heat <- ggplot(heatmap_long, aes(AM, module, fill = z)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_gradientn(
    colors = HEAT_COLORS, values = HEAT_VALUES,
    limits = c(-2, 2), oob = squish, breaks = c(-2, 0, 2),
    name = "Zscore"
  ) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), labels = scales::label_wrap(16)) +
  guides(fill = guide_colorbar(
    title.position = "top", title.hjust = 0.5,
    barwidth = grid::unit(4, "mm"), barheight = grid::unit(32, "mm"),
    ticks = FALSE, frame.colour = NA
  )) +
  labs(x = NULL, y = NULL) +
  theme(
    text = element_text(color = "black"),
    axis.text.x = element_text(
      size = 12, color = "black",
      margin = margin(b = 14)
    ),
    axis.text.y = element_text(
      size = 12, color = "black",
      lineheight = 0.9, margin = margin(r = 6)
    ),
    axis.ticks = element_blank(), axis.line = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    plot.margin = margin(t = 10, r = 8, b = 0, l = 2)
  )
if (interactive()) print(p_heat)
ggplot2::ggsave(file.path(OUT_ROOT, "Figure1i.pdf"), p_heat, width = 5.41666667, height = 3.77777778, device = grDevices::cairo_pdf, bg = "white")
ggplot2::ggsave(file.path(OUT_ROOT, "Figure1i.png"), p_heat, width = 5.41666667, height = 3.77777778, dpi = 180, bg = "white")
write.csv(data.frame(module = colnames(module_z), t(module_z), check.names = FALSE), file.path(DERIVED_ROOT, "i_plotted_values.csv"), row.names = FALSE)
