# AM-4 proximity to tumor and CD8 cells; paired comparisons in 21 patients.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(ggplot2)

setwd(am_data("tsa"))

font <- "Noto Sans"
ink <- "#20292F"
axis_lwd <- 0.25
cols <- c("Other macrophages" = "#5DBFA0", "AM-4" = "#2E44B8")
fig_width <- 3.6
fig_height <- 2.8

d <- fread("data/AM4_post21_tumor_CD8_proximity_patient_data.csv", colClasses = c(patient_id = "character"))
d[, population := factor(population, c("Other macrophages", "AM-4"))]
d[, anchor := factor(anchor, c("Tumor", "CD8"), c("Near tumor cells", "Near CD8 T cells"))]

w <- dcast(d, patient_id + anchor ~ population, value.var = "percent_within_100px")
st <- w[, .(
  n = .N,
  median_am4 = median(`AM-4`),
  median_other = median(`Other macrophages`),
  n_am4_higher = sum(`AM-4` > `Other macrophages`),
  p_wilcoxon_paired = wilcox.test(`AM-4`, `Other macrophages`, paired = TRUE, exact = .N <= 30)$p.value
), by = anchor]
st[, label := ifelse(p_wilcoxon_paired >= 0.05, "ns", sprintf("P = %.3f", p_wilcoxon_paired))]
fwrite(st, am_output("Figure4", "Figure4c_statistics.csv"))

ref <- fread(am_data("reference/figure4c_tests.csv"))
ref[, anchor := factor(anchor, c("Tumor", "CD8"), c("Near tumor cells", "Near CD8 T cells"))]
chk <- merge(st[, .(anchor, p_wilcoxon_paired)], ref[, .(anchor, p_ref = p_wilcoxon_paired)], by = "anchor")
stopifnot(all(abs(chk$p_wilcoxon_paired - chk$p_ref) < 1e-10))

p <- ggplot(d, aes(population, percent_within_100px)) +
  geom_line(aes(group = patient_id), colour = "grey82", linewidth = 0.25) +
  geom_boxplot(aes(fill = population), width = 0.55, outlier.shape = NA, alpha = 0.55, linewidth = 0.35) +
  geom_point(aes(fill = population), shape = 21, size = 1.3, stroke = 0.25, colour = ink) +
  facet_wrap(~anchor, scales = "free_y") +
  geom_text(data = st, aes(x = 1.5, y = Inf, label = label), vjust = 1.3, size = 2.6, family = font, inherit.aes = FALSE) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_discrete(labels = c("Other\nmacrophages", "AM-4")) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.15))) +
  labs(x = NULL, y = "Cells within 100 px (%)") +
  theme_classic(base_size = 9, base_family = font) +
  theme(
    axis.text = element_text(colour = "black", size = 8),
    axis.line = element_line(linewidth = axis_lwd),
    axis.ticks = element_line(linewidth = axis_lwd),
    strip.background = element_blank(),
    strip.text = element_text(size = 8.5, face = "bold")
  )
if (interactive()) print(p)
ggsave(am_output("Figure4", "AM4_tsa_tumor_CD8_distance.pdf"), p, width = 4.89, height = 3.28, device = cairo_pdf, bg = "white")
