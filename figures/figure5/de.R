# AM score survival and response: 421 randomized patients, 390 response-evaluable.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(survival)
library(ggplot2)
library(patchwork)

setwd(am_data("HTAN/am/TWAS_signature"))

d <- fread("data/Fig5cd_patients.csv")
d$AM_group <- factor(d$AM_group, levels = c("AM score-high", "AM score-low"))
d$arm <- factor(d$arm, levels = c("Atezolizumab", "Docetaxel"))
d$io <- ifelse(d$arm == "Atezolizumab", 1, 0)
d$high <- ifelse(d$AM_group == "AM score-high", 1, 0)

groups <- c("AM score-high", "AM score-low")

fmtP <- function(p) {
  if (p < 0.001) sprintf("P = %.1e", p)
  else if (p < 0.01) sprintf("P = %.3f", p)
  else sprintf("P = %.2f", p)
}

fit_int <- coxph(Surv(pfs_months, pfs_event) ~ io * high + strata(cohort), data = d)
P_int <- summary(fit_int)$coefficients["io:high", "Pr(>|z|)"]
cat("Interaction", fmtP(P_int), "\n")

hr <- data.table()
for (g in groups) {
  sub <- d[AM_group == g]
  s <- summary(coxph(Surv(pfs_months, pfs_event) ~ io + strata(cohort), data = sub))
  hr <- rbind(hr, data.table(
    group = g, n = nrow(sub),
    HR = s$conf.int["io", 1], lo = s$conf.int["io", 3], hi = s$conf.int["io", 4],
    P = s$coefficients["io", 5], P_int = P_int
  ))
}

e <- d[!is.na(objective_response)]
orr <- e[, .(responders = sum(objective_response), n = .N), by = .(group = AM_group, arm)]
orr$ORR <- 100 * orr$responders / orr$n
orr$P_fisher <- NA_real_
for (g in groups) {
  tab <- table(e[AM_group == g]$arm, e[AM_group == g]$objective_response)
  orr[group == g, P_fisher := fisher.test(tab)$p.value]
}
setorder(orr, group, arm)

fwrite(hr, am_output("Figure5", "Fig5c_pfs_hr.csv"))
fwrite(orr, am_output("Figure5", "Fig5d_orr.csv"))

stopifnot(
  abs(hr$HR - c(0.811169484404248, 1.59956629977625)) < 1e-6, abs(hr$P_int[1] - 0.000778302520132103) < 1e-9,
  abs(orr[arm == "Atezolizumab"]$P_fisher - c(0.859492681602191, 0.00731354485306428)) < 1e-9
)
print(hr)
print(orr)

font <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "Liberation Sans"
fs_c <- 11
fs_d <- 13
axis_lw <- 0.3

cols <- list(
  "AM score-high" = c(Atezolizumab = "#E25563", Docetaxel = "#F1AAB1"),
  "AM score-low"  = c(Atezolizumab = "#2E44B8", Docetaxel = "#97A2DC")
)

my_theme <- function(fs) {
  theme_classic(base_size = fs, base_family = font) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = fs),
      axis.title = element_text(size = fs),
      axis.line = element_line(linewidth = axis_lw),
      axis.ticks = element_line(linewidth = axis_lw),
      legend.text = element_text(size = fs),
      plot.title = element_text(face = "bold", hjust = 0.5, size = fs)
    )
}

ts_c <- fs_c / .pt

plot_km <- function(g, show_y) {
  sub <- d[AM_group == g]
  fit <- survfit(Surv(pfs_months, pfs_event) ~ arm, data = sub)

  km <- data.table(
    time = fit$time, surv = fit$surv, n.censor = fit$n.censor,
    arm = rep(sub("arm=", "", names(fit$strata)), fit$strata)
  )
  km <- rbind(data.table(time = 0, surv = 1, n.censor = 0, arm = c("Atezolizumab", "Docetaxel")), km)

  h <- hr[group == g]
  label <- sprintf("HR %.2f (%.2f–%.2f)\n%s", h$HR, h$lo, h$hi, fmtP(h$P))

  ggplot(km, aes(time, surv * 100, colour = arm)) +
    geom_step(linewidth = 0.45) +
    geom_point(data = km[n.censor > 0], shape = 124, size = 1.3) +
    annotate("text",
      x = 24, y = 100, label = label, hjust = 1, vjust = 1,
      size = ts_c, family = font, lineheight = 0.95
    ) +
    scale_colour_manual(values = cols[[g]], name = NULL) +
    scale_x_continuous(breaks = seq(0, 24, 6)) +
    scale_y_continuous(breaks = seq(0, 100, 25)) +
    coord_cartesian(xlim = c(0, 24), ylim = c(0, 100)) +
    labs(
      title = sprintf("%s (n = %d)", g, h$n), x = "Months",
      y = if (show_y) "Progression-free survival (%)" else NULL
    ) +
    my_theme(fs_c) +
    theme(
      plot.title = element_text(colour = cols[[g]][1]),
      legend.position = "inside", legend.position.inside = c(0.98, 0.62),
      legend.justification = c(1, 1), legend.background = element_blank(),
      legend.key.height = unit(12, "pt")
    )
}

plot_risk <- function(g, show_label) {
  fit <- survfit(Surv(pfs_months, pfs_event) ~ arm, data = d[AM_group == g])
  s <- summary(fit, times = seq(0, 24, 6), extend = TRUE)
  rk <- data.table(time = s$time, n = s$n.risk, arm = sub("arm=", "", s$strata))
  rk$arm <- factor(rk$arm, levels = c("Docetaxel", "Atezolizumab"))

  ggplot(rk, aes(time, arm, label = n, colour = arm)) +
    geom_text(size = ts_c, family = font) +
    scale_colour_manual(values = cols[[g]], guide = "none") +
    scale_y_discrete(labels = c(Atezolizumab = "Atezo", Docetaxel = "Doce")) +
    coord_cartesian(xlim = c(0, 24)) +
    labs(title = if (show_label) "No. at risk" else NULL, x = NULL, y = NULL) +
    theme_void(base_size = fs_c, base_family = font) +
    theme(
      axis.text.y = if (show_label) element_text(size = fs_c, hjust = 1) else element_blank(),
      plot.title = element_text(size = fs_c), plot.title.position = "plot"
    )
}

p5c <- (plot_km("AM score-high", TRUE) | plot_km("AM score-low", FALSE)) /
  (plot_risk("AM score-high", TRUE) | plot_risk("AM score-low", FALSE)) +
  plot_layout(heights = c(4, 1))

ts_d <- fs_d / .pt

orr$fill <- ifelse(orr$group == "AM score-high",
  ifelse(orr$arm == "Atezolizumab", "#E25563", "#F1AAB1"),
  ifelse(orr$arm == "Atezolizumab", "#2E44B8", "#97A2DC")
)

bar <- orr[, .(y = max(ORR) + 14, P = P_fisher[1]), by = group]
bar$label <- sapply(bar$P, fmtP)

p5d <- ggplot(orr, aes(arm, ORR)) +
  geom_col(aes(fill = I(fill)), width = 0.65) +
  geom_text(aes(y = ORR + 1, label = sprintf("%d/%d", responders, n)),
    vjust = 0,
    colour = "grey40", size = ts_d, family = font
  ) +
  geom_text(aes(y = ORR + 6.8, label = sprintf("%.0f%%", ORR)),
    vjust = 0,
    fontface = "bold", size = ts_d, family = font
  ) +
  geom_segment(data = bar, aes(x = 1, xend = 2, y = y, yend = y), linewidth = axis_lw) +
  geom_text(data = bar, aes(x = 1.5, y = y + 1.4, label = label), vjust = 0, size = ts_d, family = font) +
  facet_wrap(~group) +
  scale_x_discrete(labels = c(Atezolizumab = "Atezo", Docetaxel = "Doce")) +
  scale_y_continuous(breaks = seq(0, 100, 10), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, max(bar$y) + 9)) +
  labs(x = NULL, y = "Objective response rate (%)") +
  my_theme(fs_d) +
  theme(
    strip.background = element_blank(), strip.text = element_text(face = "bold", size = fs_d),
    axis.ticks.x = element_blank(), panel.spacing.x = unit(14, "pt")
  )

if (interactive()) print(p5c)
if (interactive()) print(p5d)
ggsave(am_output("Figure5", "KM_AMScore.pdf"), p5c, width = 6.32, height = 3.59, device = cairo_pdf, bg = "white")
ggsave(am_output("Figure5", "AMScore_ORR.pdf"), p5d, width = 6.01, height = 3.69, device = cairo_pdf, bg = "white")
ggsave(am_output("Figure5", "KM_AMScore.png"), p5c, width = 6.32, height = 3.59, dpi = 300, bg = "white")
ggsave(am_output("Figure5", "AMScore_ORR.png"), p5d, width = 6.01, height = 3.69, dpi = 300, bg = "white")
