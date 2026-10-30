# signature functions
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(data.table); library(survival)
})
OUT141 <- am_data("vibe/141_figure5_s6_clinical_prediction")
HZ141 <- file.path(OUT141, 'data/harmonized')
.t0_141 <- Sys.time()
log_phase <- function(s) {
  z <- paste(format(Sys.time(), '%Y-%m-%dT%H:%M:%S%z'), s, round(as.numeric(difftime(Sys.time(), .t0_141, units = 'secs')), 1), sep = '\t')
  cat(z, '\n'); cat(z, '\n', file = file.path(OUT141, 'run.log'), append = TRUE)
}

rank_enrichment <- function(expr, sets) {
  N <- nrow(expr); r <- seq.int(N, 1); w <- r^.25
  masks <- lapply(sets, function(g) rownames(expr) %in% g)
  ans <- matrix(NA_real_, ncol(expr), length(sets), dimnames = list(colnames(expr), names(sets)))
  for (j in seq_len(ncol(expr))) {
    ord <- order(expr[, j], decreasing = TRUE)
    for (k in seq_along(masks)) {
      I <- masks[[k]][ord]; ans[j, k] <- sum(w[I] * r[I]) / sum(w[I]) - mean(r[!I])
    }
  }
  ans
}
zs <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

score_cohort <- function(o, sets = readRDS(file.path(HZ141, 'gene_sets.rds'))) {
  E <- o$expr; cov <- vapply(sets, function(g) mean(g %in% rownames(E)), numeric(1))
  ok <- cov >= 0.6
  rs <- c('AM2', 'AM4', 'AM2_CD4')
  raw <- rank_enrichment(E, sets[rs])
  d <- data.table(sample = colnames(E))
  for (s in rs) d[[s]] <- if (ok[[s]]) zs(raw[, s]) else NA_real_
  gz <- function(g) {
    g <- intersect(g, rownames(E)); X <- E[g, , drop = FALSE]; X <- (X - rowMeans(X)) / apply(X, 1, sd); colMeans(X, na.rm = TRUE)
  }
  d$Teff8 <- if (ok[['Teff8']]) zs(gz(sets$Teff8)) else NA_real_
  d$TIS18 <- if (ok[['TIS18']]) zs(gz(sets$TIS18)) else NA_real_
  d$CXCL9_SPP1 <- if (ok[['CXCL9_SPP1']]) zs(E['CXCL9', ] - E['SPP1', ]) else NA_real_
  d$NI_raw <- d$AM2 - d$AM4
  d$NI <- zs(d$NI_raw)
  attr(d, 'coverage') <- cov
  d
}

cox_row <- function(fit, term, extra = list()) {
  s <- summary(fit)$coefficients; ci <- summary(fit)$conf.int
  if (!term %in% rownames(s)) return(NULL)
  data.table(
    term = term, HR = ci[term, 'exp(coef)'], lo = ci[term, 'lower .95'], hi = ci[term, 'upper .95'],
    logHR = s[term, "coef"], se = s[term, "se(coef)"], P = s[term, "Pr(>|z|)"], n = fit$n, events = fit$nevent
  )
}
cox_one <- function(d, tt, ee, x, extra = list()) {
  dd <- d[is.finite(get(tt)) & is.finite(get(ee)) & get(tt) > 0 & is.finite(get(x))]
  if (nrow(dd) < 10 || sum(dd[[ee]]) < 5) return(NULL)
  f <- coxph(as.formula(sprintf('Surv(%s, %s) ~ %s', tt, ee, x)), data = dd)
  r <- cox_row(f, x); for (k in names(extra)) r[[k]] <- extra[[k]]; r
}
