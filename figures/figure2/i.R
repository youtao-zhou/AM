# AM proximity to tumor cells in 10 untreated patients; paired Wilcoxon test.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})
tsa <- am_data("tsa")
WIDTH <- 7.2; HEIGHT <- 3.4
FONT <- "Noto Sans"; BS <- 11
cols <- c("AM-2" = "#3A8FB7", "AM-4" = "#D6604D")

P <- fread(file.path(tsa, "data/AM2_AM4_tumor_proximity_untreated_patient_data.csv"))
stopifnot(uniqueN(P$patient_id) == 10)
th <- theme_bw(base_family = FONT, base_size = BS) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#E7E7E7", linewidth = 0.35),
    panel.border = element_rect(colour = "#282828", fill = NA, linewidth = 0.55), axis.text = element_text(colour = "#252525", size = BS - 1),
    plot.title = element_text(hjust = 0.5, size = BS), plot.margin = margin(7, 8, 6, 7), legend.position = "none"
  )
panel <- function(r) {
  w <- dcast(P[radius_px == r], patient_id ~ panel, value.var = "AM_minus_other_pp")
  pv <- wilcox.test(w$`AM-4`, w$`AM-2`, paired = TRUE, exact = TRUE)$p.value
  d <- melt(w, id.vars = "patient_id", variable.name = "am", value.name = "v")[, x := as.numeric(factor(am, c("AM-2", "AM-4")))]
  top <- max(d$v, 0); rng <- diff(range(c(d$v, 0))); yb <- top + 0.10 * rng
  ggplot(d, aes(x, v)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "#858585", linewidth = 0.55) +
    geom_line(aes(group = patient_id), colour = "grey70", linewidth = 0.45) +
    geom_point(aes(colour = am), shape = 21, fill = "white", size = 2.4, stroke = 1.0) +
    annotate("segment", x = 1, xend = 2, y = yb, yend = yb, linewidth = 0.35) +
    annotate("text", x = 1.5, y = yb + 0.04 * rng, label = sprintf('italic(P) == "%.3f"', pv), parse = TRUE, vjust = 0, size = 3.4, family = FONT) +
    scale_colour_manual(values = cols) +
    scale_x_continuous(breaks = 1:2, labels = c("AM-2", "AM-4"), limits = c(0.6, 2.4)) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.16))) +
    labs(x = NULL, y = "Near tumor cells vs other Mφ (%)", title = sprintf("Within %d px of tumor cells", r)) +
    th
}
fig <- wrap_plots(list(panel(25), panel(50)), nrow = 1)
if (interactive()) print(fig)

ggsave(am_output("Figure2", "Figure2i.pdf"), fig, width = 7.40, height = 3.54, device = cairo_pdf, bg = "white")
