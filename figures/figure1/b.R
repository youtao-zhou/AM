# scPagwas cell-type effects and bootstrap 95% confidence intervals.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
source(am_code("R/figure1_data.R"))
library(dplyr)
library(ggplot2)

df <- read.csv(file.path(INPUT_ROOT, "b", "Fig1b_plot_data.csv"))

df$annotation <- recode(df$annotation, "CD1c-positive myeloid dendritic cell" = "CD1c+ myeloid DC")
df <- transmute(df, side, annotation, beta = bias_corrected_estimate, lo = CI_lo, hi = CI_hi, neglogp = neglog10_p)

df <- df %>%
  group_by(side) %>%
  arrange(desc(abs(beta)), .by_group = TRUE) %>%
  ungroup()
df$y <- factor(df$annotation, levels = rev(unique(df$annotation)))
df$side <- factor(df$side, levels = c("Tumor", "Normal"))

font <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "Liberation Sans"
fs_text <- 16
fs_title <- 14
fs_strip <- 15
fs_leg <- 13
box_lw <- 0.2

p1b <- ggplot(df, aes(x = beta, y = y, colour = neglogp)) +
  facet_grid(rows = vars(side), scales = "free_y", space = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_segment(aes(x = lo, xend = hi, yend = y), linewidth = 2.1, lineend = "round") +
  geom_point(shape = 23, size = 3, stroke = 0.3, colour = "black", fill = "#FFC75F") +
  scale_colour_gradientn(
    colours = c("#FF9671", "#FF6F91", "#845EC2"),
    limits = c(1.5, 5.5), breaks = c(1.5, 3.5, 5.5), oob = scales::squish,
    name = expression(-log[10](p))
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "Bias Corrected Estimate (95% CI)", y = NULL) +
  theme_bw(base_size = fs_text, base_family = font) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = box_lw),
    axis.ticks = element_line(colour = "black", linewidth = box_lw),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(6, "pt"),
    axis.text.y = element_text(hjust = 1, margin = margin(r = 5)),
    axis.title.x = element_text(size = fs_title),
    strip.background = element_blank(),
    strip.text.y = element_text(face = "bold", size = fs_strip),
    legend.position = "inside", legend.position.inside = c(0.80, 0.22),
    legend.background = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    legend.title = element_text(size = fs_title),
    legend.text = element_text(size = fs_leg),
    plot.margin = margin(5, 5, 5, 12)
  )

if (interactive()) print(p1b)
ggsave(file.path(OUT_ROOT, "Figure1b.pdf"), p1b, width = 5.92, height = 4.03, device = cairo_pdf, bg = "white")

write.csv(df, file.path(DERIVED_ROOT, "b_plotted_values.csv"), row.names = FALSE)
