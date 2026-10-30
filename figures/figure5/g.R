# Treatment-by-score interactions, fitted jointly for AM and TWAS scores.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(survival)
library(ggplot2)
library(patchwork)

setwd(am_data("HTAN/am/TWAS_signature"))

d <- fread("data/Fig5f_patients.csv")
d[, TWz := TWAS_z]
d[, AMz := AM_z]
d$io <- ifelse(d$arm == "Atezolizumab", 1, 0)
d$pop <- ifelse(d$cohort == "POPLAR", 1, 0)
terms <- c("TWAS score" = "io:TWz", "AM score" = "io:AMz")

fmtP <- function(p) ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.2g", p))

get_int <- function(fit, trial) {
  s <- summary(fit)
  data.table(
    score = names(terms), trial = trial,
    HR = s$conf.int[terms, 1], lo = s$conf.int[terms, 3], hi = s$conf.int[terms, 4],
    P = s$coefficients[terms, 5], n = fit$n, events = fit$nevent
  )
}
res <- rbind(
  get_int(coxph(Surv(pfs_months, pfs_event) ~ io * TWz + io * AMz + strata(cohort), data = d), "OAK + POPLAR"),
  get_int(coxph(Surv(pfs_months, pfs_event) ~ io * TWz + io * AMz, data = d[cohort == "OAK"]), "OAK"),
  get_int(coxph(Surv(pfs_months, pfs_event) ~ io * TWz + io * AMz, data = d[cohort == "POPLAR"]), "POPLAR")
)

fh <- coxph(Surv(pfs_months, pfs_event) ~ io * TWz + io * AMz + io:pop + TWz:pop + AMz:pop + io:TWz:pop + io:AMz:pop + strata(cohort), data = d)
P_het <- summary(fh)$coefficients[c("io:TWz:pop", "io:AMz:pop"), 5]
res$P_het <- ifelse(res$score == "TWAS score", P_het[1], P_het[2])

reference <- fread(am_data("reference/figure5g_tests.csv"))
stopifnot(max(abs(res$HR - reference$HR)) < 1e-9, max(abs(res$P - reference$P)) < 1e-9)
fwrite(res, am_output("Figure5", "Fig5g_interaction.csv"))
print(res)

font <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "Liberation Sans"
fs <- 11
ts <- fs / .pt
axis_lw <- 0.3
col_term <- c("TWAS score" = "#C26A1B", "AM score" = "#E25563")

rows <- data.table()
for (s in names(terms)) {
  rows <- rbind(rows,
    data.table(score = s, kind = "head", label = s),
    res[score == s, .(score, kind = "row", label = trial, n, HR, lo, hi, P)],
    fill = TRUE
  )
}
rows$y <- rev(seq_len(nrow(rows)))
top <- max(rows$y) + 1
vars <- rows[kind == "row"]
vars$pooled <- vars$label == "OAK + POPLAR"

y_scale <- scale_y_continuous(limits = c(0.45, top + 0.55), expand = c(0, 0))

p_tab <- ggplot(rows, aes(y = y)) +
  geom_text(
    data = rows[kind == "head"], aes(x = 0, label = label, colour = score),
    hjust = 0, fontface = "bold", size = ts, family = font
  ) +
  geom_text(data = vars, aes(x = 0.1, label = label), hjust = 0, size = ts, family = font) +
  geom_text(data = vars, aes(x = 1.22, label = n), hjust = 1, size = ts, family = font) +
  geom_text(
    data = vars, aes(x = 1.32, label = sprintf("%.2f (%.2f–%.2f)", HR, lo, hi)),
    hjust = 0, size = ts, family = font
  ) +
  annotate("text",
    x = c(0, 1.22, 1.32), y = top, label = c("Score / trial", "n", "HR (95% CI)"),
    hjust = c(0, 1, 0), fontface = "bold", size = ts, family = font
  ) +
  scale_colour_manual(values = col_term, guide = "none") +
  scale_x_continuous(limits = c(0, 2.2), expand = c(0, 0)) +
  y_scale +
  theme_void(base_family = font)

ref <- data.table(y0 = min(vars$y) - 0.45, y1 = max(vars$y) + 0.45)
p_for <- ggplot(vars, aes(HR, y, colour = score)) +
  geom_segment(
    data = ref, aes(x = 1, xend = 1, y = y0, yend = y1), inherit.aes = FALSE,
    linetype = "22", colour = "grey45", linewidth = 0.3
  ) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.5) +
  geom_point(data = vars[pooled == FALSE], shape = 15, size = 2.6) +
  geom_point(data = vars[pooled == TRUE], shape = 18, size = 4) +
  scale_colour_manual(values = col_term, guide = "none") +
  scale_x_log10(
    breaks = c(0.5, 1, 1.5), labels = c("0.5", "1", "1.5"),
    limits = c(0.28, 1.85), expand = c(0, 0)
  ) +
  y_scale +
  coord_cartesian(clip = "off") +
  labs(x = "Interaction HR per SD\n(atezolizumab vs docetaxel)", y = NULL) +
  theme_classic(base_size = fs, base_family = font) +
  theme(
    axis.text.x = element_text(colour = "black", size = fs),
    axis.title.x = element_text(size = fs),
    axis.line.x = element_line(linewidth = axis_lw),
    axis.ticks.x = element_line(linewidth = axis_lw),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()
  )

p_p <- ggplot(vars, aes(y = y)) +
  geom_text(aes(x = 0.15, label = fmtP(P)), hjust = 0, size = ts, family = font) +
  annotate("text", x = 0.15, y = top, label = "P", hjust = 0, fontface = "bold", size = ts, family = font) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  y_scale +
  theme_void(base_family = font)

p5g <- (p_tab | p_for | p_p) +
  plot_layout(widths = c(2.2, 1.6, 0.45)) &
  theme(plot.margin = margin(2, 2, 2, 2))

if (interactive()) print(p5g)
ggsave(am_output("Figure5", "multivariable_forest.pdf"), p5g, width = 6.35, height = 2.99, device = cairo_pdf, bg = "white")
ggsave(am_output("Figure5", "multivariable_forest.png"), p5g, width = 6.35, height = 2.99, dpi = 300, bg = "white")
