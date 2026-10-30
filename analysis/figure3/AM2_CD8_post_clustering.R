# AM2 CD8 post clustering
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(ggplot2)
library(patchwork)

setwd(am_data("tsa"))
dir.create("figure", showWarnings = FALSE)

curves <- fread("data/AM2_CD8_post_clustering_patient_curves.csv")
metrics <- fread("data/AM2_CD8_post_clustering_original_family_metrics.csv")
units <- fread("data/AM2_CD8_post_clustering_units.csv")
radius_unit <- units[field == "display_radius_unit", value]
if (!radius_unit %in% c("px", "μm")) stop("Check the documented coordinate unit.")
stopifnot(
  uniqueN(curves$patient_id) == 30L, nrow(curves) == 150L,
  all(curves$counted == "am2_query"), all(curves$anchor == "clean_cd8"),
  !anyDuplicated(curves[, .(patient_id, radius)]),
  max(abs(curves$observed_K - curves$reference_K - curves$clustering)) < 1e-7,
  all(curves[, .(ok = identical(sort(radius), c(0L, 50L, 100L, 150L, 200L))), by = patient_id]$ok)
)
setorder(curves, patient_id, radius)

patient <- curves[radius >= 50, .(
  mean_clustering_50_200 = sum(diff(radius) * (head(clustering, -1) + tail(clustering, -1)) / 2) / 150,
  response = unique(response), selected21 = unique(selected21)
), by = patient_id]
patient[, clustering_thousands := mean_clustering_50_200 / 1000]

analyse_scope <- function(ids, scope_name) {
  cv <- curves[patient_id %in% ids]
  pt <- patient[patient_id %in% ids]
  x <- pt$mean_clustering_50_200
  fit <- wilcox.test(x, mu = 0, alternative = "two.sided", exact = TRUE, correct = FALSE)
  curve_summary <- cv[, .(
    n = .N, median = median(clustering) / 1000,
    q25 = quantile(clustering, .25) / 1000,
    q75 = quantile(clustering, .75) / 1000
  ), by = radius]
  setorder(curve_summary, radius)

  family <- merge(metrics[patient_id %in% ids], pt[, .(patient_id, CrossK_AM2_absolute = mean_clustering_50_200)],
    by = "patient_id"
  )
  variables <- c(
    "CrossK_AM2_absolute", "difference_100px_pp", "difference_50px_pp",
    "difference_200px_pp", "nn_meanlog_affinity"
  )
  family_tests <- rbindlist(lapply(variables, function(v) {
    xx <- family[[v]]
    f <- wilcox.test(xx, mu = 0, alternative = "two.sided", exact = TRUE, correct = FALSE)
    data.table(
      scope = scope_name, metric = v, n = length(xx), n_positive = sum(xx > 0),
      median = median(xx), p_two_sided = f$p.value
    )
  }))
  family_tests[, p_Holm_original_five := p.adjust(p_two_sided, "holm")]
  s <- data.table(
    scope = scope_name, n = length(x), n_positive = sum(x > 0),
    n_negative = sum(x < 0), median_mean_clustering = median(x),
    median_mean_clustering_thousands = median(x) / 1000,
    p_two_sided_signed_rank = fit$p.value,
    p_Holm_original_five = family_tests[metric == "CrossK_AM2_absolute", p_Holm_original_five],
    test = "Exact two-sided Wilcoxon signed-rank vs zero; patient-level 50-200 mean integral"
  )
  list(curve = curve_summary, patient = pt, stat = s, family = family_tests)
}
main <- analyse_scope(patient[selected21 == TRUE, patient_id], "Post21_excluding_PR")
full <- analyse_scope(patient$patient_id, "All30_reference")
stopifnot(
  main$stat$n == 21, full$stat$n == 30, full$stat$n_positive == 22,
  abs(full$stat$p_two_sided_signed_rank - 0.000418396666646004) < 1e-12
)

fwrite(patient, "data/AM2_CD8_post_clustering_patient_endpoints.csv")
fwrite(
  rbindlist(list(copy(main$curve)[, scope := main$stat$scope], copy(full$curve)[, scope := full$stat$scope])),
  "data/AM2_CD8_post_clustering_curve_summary.csv"
)
fwrite(rbindlist(list(main$stat, full$stat)), "data/AM2_CD8_post_clustering_statistics.csv")
fwrite(rbindlist(list(main$family, full$family)), "data/AM2_CD8_post_clustering_original_family_statistics.csv")

fill <- "#85C3DC"
line <- "#357F9F"
ink <- "#20292F"
reference <- "#899399"
font <- "Noto Sans"
fmt_p <- function(x) if (x < .001) formatC(x, format = "e", digits = 2) else formatC(x, format = "f", digits = 4)
base_theme <- theme_classic(base_size = 11, base_family = font) +
  theme(
    axis.text = element_text(colour = ink, size = 10),
    axis.title = element_text(colour = ink, size = 11),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.title.x = element_text(margin = margin(t = 7)),
    axis.line = element_line(colour = ink, linewidth = .45),
    axis.ticks = element_line(colour = ink, linewidth = .4),
    plot.title = element_text(face = "bold", size = 12, colour = ink, margin = margin(b = 9)),
    plot.margin = margin(10, 14, 10, 10), legend.position = "none",
    plot.tag = element_text(face = "bold", size = 14), plot.tag.position = c(0, 1)
  )
y_curve <- bquote("Cross-K enrichment (" * 10^3 * " " * .(radius_unit)^2 * ")")
y_patient <- bquote("Mean enrichment (" * 10^3 * " " * .(radius_unit)^2 * ")")

make_curve <- function(x, annotate_test = FALSE) {
  z <- x$curve
  yr <- range(c(z$q25, z$q75, 0))
  height <- diff(yr)
  ymin <- min(0, yr[1]) - .05 * height
  ymax <- yr[2] + .12 * height
  p <- ggplot(z, aes(radius, median)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = fill, alpha = .30, colour = NA) +
    geom_hline(yintercept = 0, linetype = "22", linewidth = .45, colour = reference) +
    geom_line(colour = line, linewidth = 1.05) +
    geom_point(colour = line, fill = "white", shape = 21, size = 2.4, stroke = .85) +
    scale_x_continuous(breaks = seq(0, 200, 50), limits = c(0, 207), expand = expansion(mult = 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    labs(title = "AM-2 around CD8 T cells", x = paste0("r (", radius_unit, ")"), y = y_curve) +
    base_theme +
    annotate("text",
      x = 5, y = ymax - .03 * height, hjust = 0, vjust = 1,
      label = "Median and interquartile range", colour = "#617078", family = font, size = 3
    ) +
    annotate("text",
      x = 205, y = -.025 * height, hjust = 1, vjust = 1,
      label = "Tissue reference", colour = reference, family = font, size = 2.8
    )
  if (annotate_test) {
    p <- p + annotate("text",
      x = 5, y = ymax - .14 * height, hjust = 0, vjust = 1,
      label = paste0(
        "P = ", fmt_p(x$stat$p_two_sided_signed_rank), "\n",
        x$stat$n_positive, "/", x$stat$n, " patients above reference"
      ),
      colour = ink, family = font, size = 3.25, lineheight = 1.15
    )
  }
  p
}

make_patient <- function(x) {
  z <- copy(x$patient)
  z[, category := "Post-neoadjuvant"]
  yr <- range(c(z$clustering_thousands, 0))
  height <- diff(yr)
  top <- max(yr) + .23 * height
  ggplot(z, aes(category, clustering_thousands)) +
    geom_hline(yintercept = 0, linetype = "22", linewidth = .45, colour = reference) +
    stat_boxplot(geom = "errorbar", width = .13, colour = line, linewidth = .5) +
    geom_boxplot(width = .36, outlier.shape = NA, fill = fill, alpha = .30, colour = line, linewidth = .55) +
    geom_point(
      position = position_jitter(width = .115, height = 0, seed = 103),
      size = 2.0, colour = line, alpha = .88
    ) +
    scale_x_discrete(labels = paste0("Post-neoadjuvant\n(n = ", x$stat$n, ")"), expand = expansion(add = .65)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    coord_cartesian(ylim = c(min(yr) - .08 * height, top)) +
    annotate("text",
      x = 1, y = top - .02 * height,
      label = paste0("P = ", fmt_p(x$stat$p_two_sided_signed_rank)),
      vjust = 1, family = font, size = 3.6, colour = ink
    ) +
    annotate("text",
      x = 1, y = top - .105 * height,
      label = paste0(x$stat$n_positive, "/", x$stat$n, " above reference"),
      vjust = 1, family = font, size = 3.0, colour = "#617078"
    ) +
    labs(title = "Patient-level summary", x = NULL, y = y_patient) +
    base_theme
}

save_plot <- function(p, name, width, height = 4.1) {
  ggsave(file.path("figure", paste0(name, ".pdf")), p,
    width = width, height = height,
    device = cairo_pdf, bg = "white"
  )
  ggsave(file.path("figure", paste0(name, ".png")), p,
    width = width, height = height,
    dpi = 320, bg = "white"
  )
  saveRDS(p, file.path("data", paste0(name, "_ggplot.rds")))
}

A <- make_curve(main)
B <- make_patient(main)
combined <- (A + labs(tag = "a")) + (B + labs(tag = "b")) + plot_layout(widths = c(1.55, 1))
save_plot(combined, "AM2_CD8_post21_clustering", 8.0)
save_plot(make_curve(main, TRUE), "AM2_CD8_post21_clustering_curve", 5.25)
save_plot(B, "AM2_CD8_post21_clustering_patient_summary", 3.5)
full_combined <- (make_curve(full) + labs(tag = "a")) + (make_patient(full) + labs(tag = "b")) + plot_layout(widths = c(1.55, 1))
save_plot(full_combined, "AM2_CD8_post30_clustering_reference", 8.0)
if (interactive()) print(combined)
print(rbindlist(list(main$stat, full$stat)))
