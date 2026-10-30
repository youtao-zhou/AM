# AM signatures by PD-L1 and TMB category in 319 OAK patients.
# Select the top 50 marker rows before removing duplicate genes.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggpubr)
})

arg <- function(key, default) {
  x <- grep(sprintf("^--%s=", key), commandArgs(TRUE), value = TRUE)
  if (length(x)) sub(sprintf("^--%s=", key), "", x[1]) else default
}

BASE <- normalizePath(arg("base", am_data("HTAN/EGA/OAK/TMB_PDL1")), mustWork = TRUE)
OUT <- arg("out", am_output("Figure2", "ab"))
MODE <- arg("mode", "plot"); stopifnot(MODE %in% c("full", "plot"))

DATA <- file.path(BASE, "data")
RESULTS <- file.path(OUT, "results")
FIGURES <- file.path(OUT, "figures")
invisible(lapply(c(RESULTS, FIGURES), dir.create, recursive = TRUE, showWarnings = FALSE))
stopifnot("data/ not found next to the script — pass --base=DIR" = dir.exists(DATA))

START <- Sys.time()
log_phase <- function(s) {
  z <- paste(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), MODE, s,
    round(as.numeric(difftime(Sys.time(), START, units = "secs")), 1),
    sep = "\t"
  )
  cat(z, "\n"); cat(z, "\n", file = file.path(RESULTS, "run.log"), append = TRUE)
}
save_table <- function(x, nm) fwrite(x, file.path(RESULTS, paste0(nm, ".csv")))

S <- c("AM2", "AM4", "AM2_CD4", "AM2_CK", "Teff")

checks <- data.table(check = character(), passed = logical(), max_difference = numeric())
check <- function(nm, ok, diff = NA_real_) {
  checks <<- rbind(checks, data.table(check = nm, passed = isTRUE(ok), max_difference = diff))
  if (!isTRUE(ok)) {
    save_table(checks, "validation_checks"); stop("Validation failed: ", nm)
  }
}

log_phase("start")

mk <- fread(file.path(DATA, "marker_candidates_AM24.csv"))

pick <- function(st) {
  z <- copy(mk[cluster == st & p_val_adj < .05 & avg_log2FC > .25])
  setorder(z, -avg_log2FC)
  unique(head(z$gene, 50))
}

sets <- list(
  AM2     = union(pick("AM-2.CD4"), pick("AM-2.CK")),
  AM4     = pick("AM-4.CD80"),
  AM2_CD4 = pick("AM-2.CD4"),
  AM2_CK  = pick("AM-2.CK"),
  Teff    = c("CD8A", "EOMES", "PRF1", "IFNG", "CD274")
)

frozen_sets <- fread(file.path(DATA, "signature_gene_sets_frozen.csv"))
for (s in S) check(paste0("frozen_gene_set_", s), setequal(sets[[s]], frozen_sets[signature == s, gene]))
check("AM2_100_AM4_50", length(sets$AM2) == 100L && length(sets$AM4) == 50L)
save_table(
  rbindlist(lapply(S, function(s) data.table(signature = s, gene = sets[[s]]))),
  "signature_gene_sets_rebuilt"
)

d <- fread(file.path(DATA, "OAK_patient_data_frozen.csv"))
check("319_unique_patients", nrow(d) == 319L && !anyDuplicated(d$sample) && !anyDuplicated(d$UNI_ID))
check("five_finite_frozen_scores", all(is.finite(as.matrix(d[, ..S]))))

rank_enrichment <- function(expr, sets) {
  N <- nrow(expr); r <- seq.int(N, 1); w <- r^.25
  masks <- lapply(sets, function(g) rownames(expr) %in% g)
  ans <- matrix(NA_real_, ncol(expr), length(sets),
    dimnames = list(colnames(expr), names(sets))
  )
  for (j in seq_len(ncol(expr))) {
    ord <- order(expr[, j], decreasing = TRUE)
    for (k in seq_along(masks)) {
      I <- masks[[k]][ord]; stopifnot(sum(I) > 0, sum(!I) > 0)
      ans[j, k] <- sum(w[I] * r[I]) / sum(w[I]) - mean(r[!I])
      if (j == 1L) {
        running <- sum(cumsum(I * w / sum(w[I])) - cumsum((!I) / sum(!I)))
        check(
          paste0("closed_form_matches_running_sum_", names(sets)[k]),
          abs(running - ans[j, k]) < 1e-7, abs(running - ans[j, k])
        )
      }
    }
    if (j %% 100 == 0) log_phase(sprintf("rank_scores_%d_of_%d", j, ncol(expr)))
  }
  ans
}

if (MODE == "full") {
  expr <- as.matrix(readRDS(file.path(DATA, "OAK_expression_319.rds")))
  check(
    "expression_integrity",
    all(is.finite(expr)) && !anyDuplicated(rownames(expr)) &&
      !anyDuplicated(colnames(expr)) && setequal(colnames(expr), d$sample)
  )
  expr <- expr[, d$sample, drop = FALSE]

  raw <- rank_enrichment(expr, sets)

  means <- colMeans(raw); sds <- apply(raw, 2, sd)
  check("nonzero_signature_sd", all(is.finite(sds) & sds > 0))
  z <- sweep(sweep(raw, 2, means, "-"), 2, sds, "/")

  for (s in S) {
    dif <- max(abs(z[, s] - d[[s]]))
    check(paste0("rebuild_matches_frozen_", s), dif < 1e-8, dif)
    d[, (s) := z[, s]]
  }

  save_table(as.data.table(raw, keep.rownames = "sample"), "rank_enrichment_raw_319")
  save_table(
    data.table(signature = S, raw_mean = means[S], raw_SD = sds[S], n = 319L),
    "cohort_z_scaling"
  )

  coverage <- rbindlist(lapply(S, function(s) data.table(
    signature = s,
    n_defined = length(sets[[s]]),
    n_present = sum(sets[[s]] %in% rownames(expr)),
    missing   = paste(setdiff(sets[[s]], rownames(expr)), collapse = ";")
  )))
  expected_cov <- fread(file.path(DATA, "signature_coverage_expected.csv"))
  for (s in S) check(
    paste0("gene_coverage_", s),
    coverage[signature == s, n_present] == expected_cov[signature == s, n_present]
  )
  save_table(coverage, "signature_coverage")

  rm(expr); log_phase("signatures_rebuilt_and_verified")
} else {
  log_phase("frozen_scores_loaded_no_expression_recomputation")
}

check(
  "known_PDL1_categories",
  all(na.omit(d$PDL1_TC_22C3) %in% c("[0,1)", "[1,50)", "[50,100]", "NOT EVALUATED"))
)
check(
  "known_tTMB_categories",
  all(na.omit(d$tTMB) %in% c("<16", ">=16", "NOT EVALUATED"))
)

d[, pdl1_original := factor(PDL1_TC_22C3, levels = c("[0,1)", "[1,50)", "[50,100]"))]
d[, pdl1 := factor(
  fcase(PDL1_TC_22C3 %in% c("[0,1)", "[1,50)"), "<50%",
    PDL1_TC_22C3 == "[50,100]", ">=50%",
    default = NA_character_
  ),
  levels = c("<50%", ">=50%")
)]
d[, tmb := factor(tTMB, levels = c("<16", ">=16"))]
check("153_PDL1_182_tTMB", sum(!is.na(d$pdl1)) == 153L && sum(!is.na(d$tmb)) == 182L)
save_table(d, "OAK_patient_scores_and_biomarkers")

BIOMARKERS <- c(pdl1 = "PD-L1 TC (22C3)", tmb = "Recorded tTMB category")

cnt <- rbindlist(lapply(names(BIOMARKERS), function(v) {
  z <- d[!is.na(get(v)), .(n = .N, atezolizumab = sum(IO == 1), docetaxel = sum(IO == 0)), by = c(v)]
  setnames(z, v, "group")
  z[, `:=`(group = as.character(group), biomarker = v)][]
}))
save_table(cnt, "group_counts")

TESTS <- list(); DESC <- list(); LONG <- list()

for (v in names(BIOMARKERS)) {
  z <- d[!is.na(get(v))]

  LONG[[v]] <- melt(z[, c("sample", "UNI_ID", "IO", v, "AM2", "AM4"), with = FALSE],
    id.vars = c("sample", "UNI_ID", "IO", v),
    measure.vars = c("AM2", "AM4"),
    variable.name = "signature", value.name = "score"
  )
  setnames(LONG[[v]], v, "group")
  LONG[[v]][, `:=`(group = as.character(group), biomarker = v)]

  lev <- levels(z[[v]])
  for (s in S) {
    f <- as.formula(paste(s, "~", v))
    kw <- kruskal.test(f, data = z)
    wx <- wilcox.test(f, data = z, exact = FALSE, correct = FALSE)
    check(paste0("KW_Wilcoxon_", v, "_", s), abs(kw$p.value - wx$p.value) < 1e-10)

    low <- z[[s]][z[[v]] == lev[1]]
    high <- z[[s]][z[[v]] == lev[2]]
    TESTS[[length(TESTS) + 1L]] <- data.table(
      biomarker = v, signature = s, n = nrow(z), n_low = length(low), n_high = length(high),
      statistic = unname(kw$statistic), df = unname(kw$parameter),
      P = kw$p.value, wilcoxon_P = wx$p.value,
      rank_biserial_high_vs_low = mean(outer(high, low, function(a, b) sign(a - b)))
    )

    DESC[[length(DESC) + 1L]] <- z[, .(
      n = .N,
      median = median(get(s)),
      Q1 = unname(quantile(get(s), .25)),
      Q3 = unname(quantile(get(s), .75))
    ),
    by = .(group = as.character(get(v)))
    ][, `:=`(biomarker = v, signature = s)]
  }
}

res <- rbindlist(TESTS)
res[, BH_q := p.adjust(P, "BH")]
res[, BH_family := "10 tests: 2 biomarkers x AM2, AM4, AM2_CD4, AM2_CK, Teff"]

expected <- fread(file.path(DATA, "expected_binary_biomarker_tests.csv"))
rr <- merge(res, expected, by = c("biomarker", "signature"), suffixes = c("_new", "_old"))
check(
  "all_10_P_values_reproduced", nrow(rr) == 10L && max(abs(rr$P_new - rr$P_old)) < 1e-12,
  max(abs(rr$P_new - rr$P_old))
)
check(
  "all_10_q_values_reproduced", max(abs(rr$BH_q_new - rr$BH_q_old)) < 1e-12,
  max(abs(rr$BH_q_new - rr$BH_q_old))
)
check("group_counts_119_34_141_41", identical(
  as.integer(c(
    cnt[biomarker == "pdl1" & group == "<50%", n], cnt[biomarker == "pdl1" & group == ">=50%", n],
    cnt[biomarker == "tmb" & group == "<16", n], cnt[biomarker == "tmb" & group == ">=16", n]
  )),
  c(119L, 34L, 141L, 41L)
))

save_table(res, "binary_biomarker_tests")
save_table(rbindlist(DESC), "group_summary")
save_table(rbindlist(LONG), "TMB_PDL1_AM24_plot_data_long")
log_phase("all_original_statistics_reproduced")

COL <- c(AM2 = "#4DBBD5", AM4 = "#E64B35")
LAB <- c(AM2 = "AM-2", AM4 = "AM-4")
XLAB <- c(pdl1 = "PD-L1", tmb = "TMB")

FONT <- local({
  cand <- c("Arial", "Helvetica", "Noto Sans", "DejaVu Sans")
  have <- tryCatch(systemfonts::system_fonts()$family, error = function(e) character())
  hit <- cand[cand %in% have]
  if (length(hit)) hit[1] else ""
})

theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size, base_family = FONT) +
    theme(
      axis.line = element_line(linewidth = .25, colour = "black"),
      axis.ticks = element_line(linewidth = .25, colour = "black"),
      axis.ticks.length = unit(3, "pt"),
      axis.text = element_text(colour = "black", size = rel(1)),
      axis.title = element_text(colour = "black", size = rel(1.1)),
      axis.title.x = element_text(margin = margin(t = 7)),
      axis.title.y = element_text(margin = margin(r = 7)),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold", size = rel(1.2),
        colour = "black", margin = margin(b = 6)
      ),
      panel.spacing = unit(14, "pt"),
      plot.margin = margin(8, 10, 6, 8),
      legend.position = "none"
    )
}

fmt_p <- function(p) {
  if (!is.finite(p)) return("P = NA")
  if (p >= 1e-3) return(sprintf("P = %.3f", p))
  e <- floor(log10(p))
  sprintf("P = %.1fe%d", signif(p / 10^e, 2), e)
}
nice <- function(g) gsub(">=", "\u2265", g, fixed = TRUE)

panel <- function(v) {
  z <- copy(LONG[[v]])
  lev <- if (v == "pdl1") c("<50%", ">=50%") else c("<16", ">=16")
  z[, `:=`(
    group = factor(group, levels = lev),
    signature = factor(signature, levels = c("AM2", "AM4"))
  )]

  ns <- cnt[biomarker == v]
  axis_lab <- setNames(
    vapply(lev, function(g) sprintf("%s\n(n = %d)", nice(g), ns[group == g, n]), character(1)), lev
  )

  ann <- merge(res[biomarker == v & signature %in% c("AM2", "AM4"), .(signature, P)],
    z[, .(hi = max(score), span = diff(range(score))), by = signature],
    by = "signature"
  )
  ann[, `:=`(
    group1 = lev[1], group2 = lev[2],
    p.label = vapply(P, fmt_p, character(1)),
    y.position = hi + .10 * span,
    y.top = hi + .30 * span,
    signature = factor(signature, levels = c("AM2", "AM4"))
  )]

  ggplot(z, aes(group, score)) +
    geom_hline(yintercept = 0, colour = "grey88", linewidth = .4) +
    geom_boxplot(aes(fill = signature),
      width = .55, colour = "grey20",
      linewidth = .5, fatten = 1.6,
      outlier.size = .9, outlier.colour = "grey35", outlier.stroke = 0
    ) +
    stat_pvalue_manual(ann,
      label = "p.label", tip.length = .015,
      bracket.size = .25, size = 4.2, colour = "black", family = FONT
    ) +
    geom_blank(aes(x = 1, y = y.top), data = ann, inherit.aes = FALSE) +
    facet_wrap(~signature, nrow = 1, labeller = as_labeller(LAB)) +
    scale_fill_manual(values = COL) +
    scale_x_discrete(labels = axis_lab) +
    scale_y_continuous(expand = expansion(mult = c(.05, .02))) +
    labs(x = unname(XLAB[v]), y = "AM Score") +
    theme_pub()
}

save_plot <- function(gg, nm, w, h) {
  ggsave(file.path(FIGURES, paste0(nm, ".pdf")), gg, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIGURES, paste0(nm, ".png")), gg, width = w, height = h, dpi = 600, bg = "white")
}

fig_pdl1 <- panel("pdl1")
fig_tmb <- panel("tmb")
save_plot(fig_pdl1, "PDL1_AM24_Signature", 6.42, 4.21)
save_plot(fig_tmb, "TMB_AM24_Signature", 6.42, 4.21)

combo <- (fig_pdl1 | fig_tmb) +
  plot_annotation(
    tag_levels = "a",
    caption = paste0(
      "Boxes, median and IQR; whiskers, 1.5x IQR; points, outliers. ",
      "Two-sided Wilcoxon rank-sum test within each signature.\n",
      "PD-L1 groups merged at the 50% boundary; TMB groups, patients and raw P values unchanged."
    ),
    theme = theme(
      plot.caption = element_text(
        size = 9, colour = "grey35", hjust = 0,
        family = FONT, margin = margin(t = 8)
      ),
      plot.margin = margin(8, 12, 8, 10)
    )
  ) &
  theme(plot.tag = element_text(face = "bold", size = 15, family = FONT))

save_plot(combo, "TMB_PDL1_AM24_Signature", 8.8, 3.9)

for (nm in c("validation_checks", paste0("validation_checks_", MODE))) save_table(checks, nm)
inputs <- list.files(DATA, full.names = TRUE)
save_table(
  data.table(file = basename(inputs), md5 = unname(tools::md5sum(inputs))),
  "input_MD5_manifest"
)
capture.output(sessionInfo(), file = file.path(RESULTS, paste0("session_info_", MODE, ".txt")))

print(res[signature %in% c("AM2", "AM4"), .(biomarker, signature, n_low, n_high, P, BH_q)])
log_phase("complete")
