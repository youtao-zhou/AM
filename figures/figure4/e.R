# T-cell signatures across AM-4 score deciles in 319 OAK patients.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(ggplot2)

base_dir <- am_data("HTAN/EGA/OAK/AM4CD4")
fig_dir <- file.path(base_dir, "figure")

font_family <- "Noto Sans"
z_limit <- 0.6
axis_order <- c("All T cells", "CD8 effector", "Memory", "Helper", "Treg")

d <- fread(file.path(base_dir, "data", "AM4_heatmap_plot_data.csv"))
stopifnot(setequal(unique(d$axis), axis_order))

d[, axis := factor(axis, levels = axis_order)]
d[, gene := factor(gene, levels = rev(unique(d[order(gene_order), gene])))]

groups <- sort(unique(d$group))
group_labels <- rep("", length(groups))
group_labels[c(1, length(groups))] <- c("Low", "High")
d[, group := factor(group, levels = groups)]

p <- ggplot(d, aes(group, gene, fill = mean_z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  facet_grid(axis ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#C8423B",
    limits = c(-z_limit, z_limit),
    oob = scales::squish,
    name = "Zscore"
  ) +
  scale_x_discrete(labels = group_labels) +
  labs(x = "AM-4 Signature", y = NULL) +
  theme_test(base_size = 11, base_family = font_family) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black", size = 11),
    axis.text.y = element_text(face = "italic"),
    axis.title.x = element_text(size = 11),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(
      colour = "black", angle = 0, hjust = 1, size = 11
    ),
    panel.spacing = grid::unit(1.5, "pt"),
    legend.title = element_text(colour = "black", size = 11),
    legend.text = element_text(colour = "black", size = 11),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3)
  )

if (interactive()) print(p)

ggsave(
  file.path(am_output("Figure4", "Figure4e.pdf")), p,
  width = 4.67, height = 4.78, device = cairo_pdf, bg = "white"
)
