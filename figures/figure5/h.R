# Tissue AM score = z(AM2 signature) - z(log2(CD80 O/E + 0.01)).
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(ggplot2)
library(patchwork)

setwd(am_data("HTAN/am/GSE243013"))

tis <- fread("data/Fig5g_patients.csv")
tis$response <- factor(tis$response, levels = c("non-MPR", "MPR", "pCR"))

zs <- function(x) (x - mean(x)) / sd(x)
tis$score <- zs(tis$AM2_signature) - zs(log2(tis$CD80_OE + 0.01))
stopifnot(max(abs(tis$score - tis$tissue_AM_score)) < 1e-10)

fmtP <- function(p) {
  if (p < 0.001) sprintf("P = %.1e", p)
  else if (p < 0.01) sprintf("P = %.3f", p)
  else sprintf("P = %.2f", p)
}

tis$ord <- as.integer(tis$response)
P_trend <- cor.test(tis$score, tis$ord, method = "spearman", exact = FALSE)$p.value

rho_reg <- cor(tis$score, tis$pathological_regression, method = "spearman", use = "complete.obs")
P_reg <- cor.test(tis$score, tis$pathological_regression, method = "spearman", exact = FALSE)$p.value
st <- data.table(
  n = nrow(tis), P_trend = P_trend, rho_reg = rho_reg, P_reg = P_reg,
  n_reg = sum(!is.na(tis$pathological_regression))
)

for (r in levels(tis$response)) st[[paste0("median_", r)]] <- median(tis[response == r]$score)
fwrite(st, am_output("Figure5", "Fig5g_stats.csv"))
print(st)
cat(sprintf("GSE243013 LUAD, post-treatment (n = %d)\n", st$n))

font <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "Liberation Sans"
fs <- 11
ts <- fs / .pt
axis_lw <- 0.2

cols <- c("non-MPR" = "#FBA72A", "MPR" = "#78B7C5", "pCR" = "#7294D4")
darker <- function(col, f = 0.72) rgb(t(col2rgb(col) * f), maxColorValue = 255)
cols_dk <- setNames(darker(cols), names(cols))
ylim_all <- c(-3.6, 3.6)
my_theme <- theme_classic(base_size = fs, base_family = font) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black", size = fs),
    axis.title = element_text(size = fs),
    axis.line = element_line(linewidth = axis_lw),
    axis.ticks = element_line(linewidth = axis_lw)
  )

set.seed(1)
p_left <- ggplot(tis, aes(response, score, fill = response, colour = response)) +
  geom_boxplot(outlier.shape = NA, width = 0.62, linewidth = 0.3, alpha = 0.4) +
  geom_jitter(shape = 21, width = 0.13, height = 0, size = 1.6, alpha = 0.55, stroke = 0.25) +
  annotate("text", x = 2, y = 3.35, label = fmtP(P_trend), size = ts, family = font) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_colour_manual(values = cols_dk, guide = "none") +
  scale_x_discrete(labels = c("Non-MPR", "MPR", "pCR")) +
  scale_y_continuous(breaks = c(-3, 0, 3), expand = c(0, 0)) +
  coord_cartesian(ylim = ylim_all) +
  labs(x = NULL, y = "Tissue AM score") +
  my_theme +
  theme(plot.margin = margin(2, 10, 2, 2))

tr <- tis[!is.na(pathological_regression)]
tr$reg_pct <- 100 * tr$pathological_regression
p_right <- ggplot(tr, aes(reg_pct, score)) +
  geom_smooth(
    method = "lm", formula = y ~ x, se = TRUE, colour = "#3A3F4B", fill = "#C9CED6",
    alpha = 0.45, linewidth = 0.5, lineend = "round"
  ) +
  geom_point(aes(fill = response, colour = response), shape = 21, size = 2, alpha = 0.7, stroke = 0.3) +
  annotate("text",
    x = -2, y = 3.35, label = sprintf("ρ = %.2f\n%s", rho_reg, fmtP(P_reg)),
    hjust = 0, vjust = 1, lineheight = 0.9, size = ts, family = font
  ) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_colour_manual(values = cols_dk, guide = "none") +
  scale_x_continuous(breaks = c(0, 50, 100)) +
  scale_y_continuous(breaks = c(-3, 0, 3)) +
  coord_cartesian(xlim = c(-6, 106), ylim = ylim_all, expand = FALSE) +
  labs(x = "Pathological regression (%)", y = NULL) +
  my_theme +
  theme(plot.margin = margin(2, 6, 2, 4))

p5g <- (p_left | p_right) + plot_layout(widths = c(1.25, 1.75))

if (interactive()) print(p5g)
ggsave(am_output("Figure5", "Tissue_AM_score_LUAD.pdf"), p5g, width = 5.13, height = 2.69, device = cairo_pdf, bg = "white")
