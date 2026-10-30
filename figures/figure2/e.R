# Tumor enrichment in 20 paired patients, using the atlas reference.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(ggplot2)

root <- am_data("sp")
REFERENCE <- "atlas"
COLORS <- c('AM-2' = '#237F96', 'AM-4' = '#C74940')
POINT_SIZE <- 2.1
FONT <- 'Arial'
LINE <- 0.3
WIDTH <- 5.36
HEIGHT <- 3.47

sfx <- if (REFERENCE == "atlas") "_atlas" else ""

patient <- read.csv(file.path(root, sprintf('AM24_tumor_fraction_patient_data%s.csv', sfx)), check.names = FALSE)
long <- read.csv(file.path(root, sprintf('AM24_tumor_fraction_plot_data%s.csv', sfx)), check.names = FALSE)
long$AM <- factor(long$AM, levels = c('AM-2', 'AM-4'))

p <- wilcox.test(patient$AM4_tumor_pct, patient$AM2_tumor_pct,
  paired = TRUE, exact = TRUE, correct = FALSE
)$p.value
p_label <- if (p < 0.001) {
  e <- floor(log10(p)); sprintf('italic(P) == %.1f %%*%% 10^%d', p / 10^e, e)
} else sprintf('italic(P) == "%.3f"', p)

y_top <- max(long$tumor_fraction_percent, na.rm = TRUE)
Y_MAX <- ceiling((y_top + 6) / 10) * 10
y_brk <- y_top + 2.5

set.seed(2201)
long$x_jitter <- as.numeric(long$AM) + runif(nrow(long), -0.105, 0.105)

g <- ggplot(long, aes(AM, tumor_fraction_percent)) +
  geom_boxplot(aes(fill = AM),
    width = 0.45, outlier.shape = NA,
    alpha = 0.22, colour = '#444444', linewidth = LINE, show.legend = FALSE
  ) +
  geom_point(aes(x = x_jitter, fill = AM),
    shape = 21, size = POINT_SIZE,
    colour = 'white', stroke = 0.2, alpha = 0.8, show.legend = FALSE
  ) +
  scale_fill_manual(values = COLORS) +
  annotate('segment',
    x = c(1, 1, 2), xend = c(1, 2, 2),
    y = c(y_brk - 1, y_brk, y_brk), yend = c(y_brk, y_brk, y_brk - 1), linewidth = LINE
  ) +
  annotate('text', x = 1.5, y = y_brk + 2.2, label = p_label, parse = TRUE, size = 4, family = FONT) +
  scale_y_continuous(
    limits = c(0, Y_MAX), breaks = seq(0, Y_MAX, 10),
    expand = expansion(mult = 0)
  ) +
  scale_x_discrete(expand = expansion(add = 0.55)) +
  labs(x = NULL, y = 'AM in tumor region (%)') +
  theme_classic(base_size = 11, base_family = FONT) +
  theme(
    axis.text = element_text(colour = 'black'),
    axis.text.x = element_text(size = 12, margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 7)),
    axis.ticks = element_line(linewidth = LINE),
    axis.ticks.x = element_blank(),
    axis.line = element_line(linewidth = LINE),
    plot.margin = margin(8, 8, 5, 5)
  )

if (interactive()) print(g)

outdir <- am_output("Figure2", "e")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
stem <- file.path(outdir, sprintf('AM24_tumor_fraction_statistics_sp%s', sfx))
ggsave(paste0(stem, '.pdf'), g, width = WIDTH, height = HEIGHT, units = 'in', device = cairo_pdf)
ggsave(paste0(stem, '.png'), g, width = WIDTH, height = HEIGHT, units = 'in', dpi = 300, bg = 'white')
cat(sprintf(
  'Fig. 2e (%s reference): n = %d, median paired difference %.2f pp, %d/%d AM-4 higher, P = %.3g\n', REFERENCE,
  nrow(patient), median(patient$delta_pct), sum(patient$delta_pct > 0), nrow(patient), p
))
