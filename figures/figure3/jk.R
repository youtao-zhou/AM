# Patient-level spatial curves and paired exact Wilcoxon comparison (21 patients).
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(ggplot2)
library(patchwork)

setwd(am_data("tsa"))

curves <- fread("data/AM2_CD8_post_clustering_patient_curves.csv")[selected21 == TRUE]
roster <- fread("data/AM2_CD8_post_clustering_roster.csv")[selected21 == TRUE]
metrics <- fread("data/AM2_CD8_post_clustering_original_family_metrics.csv")
units <- fread("data/AM2_CD8_post_clustering_units.csv")
u <- units[field == "display_radius_unit", value]

cv <- curves[, .(
  median = median(clustering) / 1000,
  q25 = quantile(clustering, 0.25) / 1000,
  q75 = quantile(clustering, 0.75) / 1000
), by = radius]
setorder(cv, radius)

pt <- metrics[patient_id %in% roster$patient_id]
pv <- wilcox.test(pt$AM2_100px_percent, pt$other_100px_percent,
  paired = TRUE, exact = TRUE
)$p.value

bx <- melt(pt,
  id.vars = "patient_id",
  measure.vars = c("other_100px_percent", "AM2_100px_percent"),
  variable.name = "pop", value.name = "percent"
)
bx[, pop := factor(pop, c("other_100px_percent", "AM2_100px_percent"), c("Other macrophages", "AM-2"))]

cols <- c("Other macrophages" = "#91CAE8", "AM-2" = "#F48892")
ink <- "#20292F"
fnt <- "Noto Sans"
ylab_k <- bquote("Cross-K enrichment (" * 10^3 * " " * .(u)^2 * ")")

th <- theme_classic(base_size = 11, base_family = fnt) +
  theme(
    axis.text = element_text(colour = ink, size = 10),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.line = element_line(colour = ink, linewidth = 0.25),
    axis.ticks = element_line(colour = ink, linewidth = 0.4),
    plot.title = element_text(face = "bold", size = 12, margin = margin(b = 9)),
    plot.margin = margin(10, 14, 10, 10),
    legend.position = "none"
  )

A <- ggplot(cv, aes(radius, median)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.45, linetype = "22") +
  geom_line(colour = cols[["AM-2"]], linewidth = 1.2) +
  geom_point(shape = 21, fill = "white", colour = cols[["AM-2"]], size = 2.7, stroke = 1) +
  scale_x_continuous(breaks = seq(0, 200, 50), limits = c(0, 207), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 40, 10), expand = c(0, 0)) +
  coord_cartesian(ylim = c(-4, 40)) +
  labs(title = "AM-2 around CD8 T cells", x = paste0("r (", u, ")"), y = ylab_k) +
  th

B <- ggplot(bx, aes(pop, percent, fill = pop)) +
  stat_boxplot(geom = "errorbar", width = 0.17, linewidth = 0.5, colour = ink) +
  geom_boxplot(width = 0.48, outlier.shape = NA, colour = ink, alpha = 0.6, linewidth = 0.55) +
  geom_point(
    shape = 21, colour = ink, stroke = 0.3, size = 2,
    position = position_jitter(width = 0.1, height = 0, seed = 104)
  ) +
  annotate("segment", x = 1, xend = 2, y = 105, yend = 105, linewidth = 0.45, colour = ink) +
  annotate("segment", x = c(1, 2), xend = c(1, 2), y = 103, yend = 105, linewidth = 0.45, colour = ink) +
  annotate("text",
    x = 1.5, y = 109, label = sprintf("P = %.4f", pv),
    size = 3.6, family = fnt, colour = ink
  ) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(labels = c("Other\nmacrophages", "AM-2"), expand = expansion(add = 0.6)) +
  scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 114), clip = "off") +
  labs(title = "", x = NULL, y = "Near CD8 T cells (%)") +
  th

S <- ggplot(cv, aes(radius, median)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = cols[["AM-2"]], alpha = 0.3) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4, linetype = "22") +
  geom_line(colour = cols[["AM-2"]], linewidth = 0.9) +
  geom_point(shape = 21, fill = "white", colour = cols[["AM-2"]], size = 2, stroke = 0.7) +
  scale_x_continuous(breaks = seq(0, 200, 50), limits = c(0, 207), expand = c(0, 0)) +
  labs(title = "AM-2 around CD8 T cells", x = paste0("r (", u, ")"), y = ylab_k) +
  th

p <- A + B + plot_layout(widths = c(1.55, 1))

ggsave(am_output("Figure3", "AM2_CD8_CrossK.pdf"), A, width = 4.65, height = 3.38, device = cairo_pdf, bg = "white")
ggsave(am_output("Figure3", "AM2_otherAM_CD8_Crossk.pdf"), B, width = 4.65, height = 3.38, device = cairo_pdf, bg = "white")

if (interactive()) print(p)
fwrite(cv, am_output("Figure3", "Figure3j_values.csv"))
write.csv(data.frame(P = pv, n = nrow(pt)), am_output("Figure3", "Figure3k_test.csv"), row.names = FALSE)
