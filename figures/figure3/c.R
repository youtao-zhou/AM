# CytoSig response profiles from the saved activity matrix.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
out <- am_output("Figure3", "")
library(ggplot2)

base_dir <- am_data("HTAN/am/cytosig")

fine <- as.matrix(read.csv(
  file.path(base_dir, "B_CD4_CK_heatmap_values_shared_scale.csv"),
  row.names = 1, check.names = FALSE
))

groups <- c("AM-1", "AM-2.CD4", "AM-2.CK", "AM-3", "AM-4")
responses <- c("IFN1", "IFNG", "IL12", "CD40L", "CXCL12", "GCSF")

fine <- fine[, groups, drop = FALSE]
main <- fine[responses, , drop = FALSE]

xlab <- c(
  "AM-1" = "AM-1",
  "AM-2.CD4" = "AM-2 CD4",
  "AM-2.CK" = "AM-2 CK",
  "AM-3" = "AM-3",
  "AM-4" = "AM-4"
)

pal <- c(
  "#006D77", "#4FA3A5", "#A8D5D1", "#F7F7F7",
  "#E9B8CF", "#C9579A", "#7A0F4F"
)

to_long <- function(m) {
  data.frame(
    response = rep(rownames(m), ncol(m)),
    group = rep(colnames(m), each = nrow(m)),
    value = as.vector(m)
  )
}

dm <- to_long(main)
dm$response <- factor(dm$response, levels = rev(responses))
dm$group <- factor(dm$group, levels = groups)

p_main <- ggplot(dm, aes(group, response)) +
  geom_point(
    aes(size = abs(value), fill = value),
    shape = 21, color = "grey30", stroke = 0.25
  ) +
  scale_fill_gradientn(
    colors = pal, limits = c(-2.5, 2.5), breaks = c(-2, 0, 2),
    oob = scales::squish, name = "Zscore"
  ) +
  scale_size_area(max_size = 9, name = NULL) +
  scale_x_discrete(labels = xlab) +
  labs(x = NULL, y = NULL) +
  guides(
    fill = guide_colorbar(order = 1),
    size = guide_legend(
      order = 2,
      override.aes = list(fill = "grey70", color = "grey30")
    )
  ) +
  theme_minimal(base_size = 11, base_family = "DejaVu Sans") +
  theme(
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 10),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.title = element_text(color = "black", size = 10),
    legend.text = element_text(color = "black", size = 10),
    plot.margin = margin(8, 10, 8, 8)
  )

ds <- to_long(fine)
ds$response <- factor(ds$response, levels = rev(rownames(fine)))
ds$group <- factor(ds$group, levels = groups)
k <- match("AM-2.CD4", groups)

p_supp <- ggplot(ds, aes(group, response, fill = value)) +
  geom_tile(color = "white", linewidth = 0.18) +
  annotate(
    "rect",
    xmin = k - 0.5, xmax = k + 0.5,
    ymin = 0.5, ymax = nrow(fine) + 0.5,
    fill = NA, color = "grey15", linewidth = 0.4
  ) +
  scale_fill_gradientn(
    colors = pal, limits = c(-2.5, 2.5), breaks = c(-2, 0, 2),
    oob = scales::squish, name = "Zscore"
  ) +
  scale_x_discrete(labels = xlab, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11, base_family = "DejaVu Sans") +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 8.5),
    legend.title = element_text(color = "black", size = 10),
    legend.text = element_text(color = "black", size = 10),
    plot.margin = margin(8, 10, 8, 8)
  )

if (interactive()) print(p_main)
if (interactive()) print(p_supp)
ggsave(
  file.path(am_output("Figure3", "Figure3c.pdf")), p_main,
  width = 4.97, height = 3.50,
  device = cairo_pdf, bg = "white"
)

ggsave(
  file.path(out, "CytoSig_AM2_heatmap.pdf"), p_supp,
  width = 5.31, height = 5.74,
  device = cairo_pdf, bg = "white"
)
