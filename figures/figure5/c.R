# TWAS score and survival in 421 patients; OAK-reference z > 0 defines High.
# Pooled Cox models are stratified by trial.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(survival)
library(ggplot2)
library(patchwork)

setwd(am_data("HTAN/am/TWAS_signature"))

d <- fread("data/Fig5c_patients.csv")
d$TWAS_group <- factor(d$TWAS_group, levels = c("High", "Low"))
d$io <- ifelse(d$arm == "Atezolizumab", 1, 0)

eps <- list(PFS = c("pfs_months", "pfs_event"), OS = c("os_months", "os_event"))

fmt_p <- function(p) if (p < 1e-4) "p < 0.0001" else paste("p =", signif(p, 2))

stat <- data.table()
for (ep in names(eps)) {
  d$time <- d[[eps[[ep]][1]]]; d$event <- d[[eps[[ep]][2]]]

  fi <- summary(coxph(Surv(time, event) ~ io * TWAS_z + strata(cohort), data = d))
  P_int <- fi$coefficients["io:TWAS_z", 5]
  for (a in c("Atezolizumab", "Docetaxel")) {
    sub <- d[arm == a]
    s <- summary(coxph(Surv(time, event) ~ I(TWAS_group == "High") + strata(cohort), data = sub))
    lr <- survdiff(Surv(time, event) ~ TWAS_group + strata(cohort), data = sub)
    stat <- rbind(stat, data.table(
      endpoint = ep, arm = a, n = nrow(sub), events = sum(sub$event), n_high = sum(sub$TWAS_group == "High"),
      HR = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4],
      P_logrank = 1 - pchisq(lr$chisq, 1), P_int = P_int
    ))
  }
}

reference <- fread(am_data("reference/figure5c_tests.csv"))
pfs <- stat[endpoint == "PFS"]
stopifnot(max(abs(pfs$HR - reference$HR)) < 1e-9, max(abs(pfs$P_logrank - reference$P)) < 1e-9, abs(pfs$P_int[1] - 0.0127639842546729) < 1e-9)

os <- stat[endpoint == "OS"]
stopifnot(all(round(os$HR, 2) == c(1.38, 1.28)), all(signif(os$P_logrank, 2) == c(0.059, 0.14)), signif(os$P_int[1], 2) == 0.63)
fwrite(stat, am_output("Figure5", "Fig5c_twas_km_hr.csv"))
print(stat)

font <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "Liberation Sans"
pal <- c(High = "#D55E00", Low = "#0072B2")

stepxy <- function(t, y) {
  n <- length(t); data.table(time = c(t[1], rep(t[-1], each = 2)), y = c(rep(y[-n], each = 2), y[n]))
}

plot_arm <- function(a, ep) {
  sub <- d[arm == a]; st <- stat[arm == a & endpoint == ep]
  sub$time <- sub[[eps[[ep]][1]]]; sub$event <- sub[[eps[[ep]][2]]]
  km <- rbindlist(lapply(levels(d$TWAS_group), function(g) {
    f <- survfit(Surv(time, event) ~ 1, data = sub[TWAS_group == g])
    data.table(
      group = g, time = c(0, f$time), surv = c(1, f$surv), lower = nafill(c(1, f$lower), "locf"),
      upper = nafill(c(1, f$upper), "locf"), n.censor = c(0, f$n.censor)
    )
  }))
  band <- km[,
    {
      lo <- stepxy(time, lower); up <- stepxy(time, upper); .(time = lo$time, lower = lo$y, upper = up$y)
    },
    by = group
  ]
  xmax <- max(sub$time)

  lab <- sprintf(
    "HR = %.2f (%.2f–%.2f)%s", st$HR, st$lo, st$hi,
    if (a == "Atezolizumab") sprintf("\n%s", fmt_p(st$P_int)) else ""
  )
  ggplot(km, aes(time, surv, colour = group)) +
    geom_ribbon(data = band, aes(x = time, ymin = lower, ymax = upper, fill = group), inherit.aes = FALSE, alpha = 0.3) +
    geom_step(linewidth = 1) +
    geom_point(data = km[n.censor > 0], shape = 3, size = 2.5, show.legend = FALSE) +
    annotate("text", x = xmax / 50, y = 0.08, label = fmt_p(st$P_logrank), hjust = 0, size = 7, family = font) +
    annotate("text", x = xmax, y = 1, label = "OAK + POPLAR", hjust = 1, vjust = 1, size = 6, family = font, fontface = "bold") +
    annotate("text", x = xmax, y = 0.9, label = lab, hjust = 1, vjust = 1, size = 6, family = font, lineheight = 0.95) +
    scale_colour_manual(values = pal, name = "TWAS score") +
    scale_fill_manual(values = pal, name = "TWAS score") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = sprintf("%.2f", seq(0, 1, 0.25))) +
    labs(x = sprintf("%s - %s", a, ep), y = "Survival probability") +
    theme_bw(base_size = 25, base_family = font) +
    theme(
      legend.position = "top",
      axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
      axis.text = element_text(colour = "black"),
      panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.4),
      legend.title = element_text(size = 18), legend.text = element_text(size = 18)
    )
}

plots <- list()
for (ep in names(eps)) for (a in c("Atezolizumab", "Docetaxel")) {
  p <- plot_arm(a, ep)
  plots[[paste(ep, a)]] <- p
  if (interactive()) print(p)
  f <- sprintf(am_output("Figure5", "TWAS_signature_pooled_%s_%s_survivalcurve"), a, ep)
  ggsave(paste0(f, ".pdf"), p, width = 6, height = 4.5, device = cairo_pdf, bg = "white")
  ggsave(paste0(f, ".png"), p, width = 6, height = 4.5, dpi = 300, bg = "white")
}

p4 <- (plots[["PFS Docetaxel"]] | plots[["PFS Atezolizumab"]]) /
  (plots[["OS Docetaxel"]] | plots[["OS Atezolizumab"]]) +
  plot_layout(guides = "collect") & theme(legend.position = "top")
if (interactive()) print(p4)
ggsave(am_output("Figure5", "TWAS_signature_pooled_PFS_OS_2x2.pdf"), p4, width = 12, height = 8.6, device = cairo_pdf, bg = "white")
ggsave(am_output("Figure5", "TWAS_signature_pooled_PFS_OS_2x2.png"), p4, width = 12, height = 8.6, dpi = 300, bg = "white")
