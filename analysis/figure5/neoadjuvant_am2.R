# neoadjuvant am2
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
source(am_code("analysis/figure5/signature_functions.R"))
D4 <- file.path(OUT141, 'data/m4'); dir.create(D4, recursive = TRUE, showWarnings = FALSE)
pr <- fread(am_data("vibe/140_figure4_v2_am4_niche/data/a/gse243013_macrophage_predictions.csv.gz"))
cl <- fread(am_data("vibe/140_figure4_v2_am4_niche/data/gse243013/sample_clinical.csv"))
pt <- pr[, .(
  n_mac = .N, AM2_like_pct = 100 * mean(pred_macro == 'AM-2'), AM4_like_pct = 100 * mean(pred_macro == 'AM-4'),
  sigAM2 = mean(sig_AM2), sigAM4 = mean(sig_AM4)
), by = sampleID]
pt <- merge(pt, cl[, .(sampleID, histology = cancer_type, response = pathological_response, regression = pathological_response_rate)], by = 'sampleID')
pt <- pt[n_mac >= 30 & response %in% c('pCR', 'MPR', 'non-MPR')]
pt[, NI_mac := zs(zs(sigAM2) - zs(sigAM4))][, `:=`(pCR = as.integer(response == 'pCR'), MPR = as.integer(response %in% c('pCR', 'MPR')))]
pt[, regression := suppressWarnings(as.numeric(regression))]
fwrite(pt, file.path(D4, 'gse243013_patient_readouts.csv'))
cat(sprintf('patients=%d  LUAD=%d  pCR=%d  MPR=%d\n', nrow(pt), sum(pt$histology == 'LUAD'), sum(pt$pCR), sum(pt$MPR)))
V <- c('AM2_like_pct', 'AM4_like_pct', 'NI_mac')
res <- list()
for (v in V) for (ep in c('pCR', 'MPR')) for (h in c('All (histology-adjusted)', 'LUAD', 'LUSC')) {
  d <- if (h == 'All (histology-adjusted)') pt else pt[histology == h]
  d <- copy(d)[, xs := zs(get(v))]
  f <- if (h == 'All (histology-adjusted)') glm(as.formula(sprintf('%s ~ xs + histology', ep)), binomial, d) else glm(as.formula(sprintf('%s ~ xs', ep)), binomial, d)
  b <- summary(f)$coefficients['xs', ]
  rs <- if (h == 'All (histology-adjusted)') NA_real_ else suppressWarnings(cor.test(d[[v]], d$regression, method = 'spearman', exact = FALSE))$estimate
  rp <- if (h == 'All (histology-adjusted)') NA_real_ else suppressWarnings(cor.test(d[[v]], d$regression, method = 'spearman', exact = FALSE))$p.value
  res[[length(res) + 1]] <- data.table(
    readout = v, endpoint = ep, subset = h, n = nrow(d), n_pos = sum(d[[ep]]), OR_perSD = exp(b[1]),
    lo = exp(b[1] - 1.96 * b[2]), hi = exp(b[1] + 1.96 * b[2]), P = b[4], rho_regression = rs, P_rho = rp
  )
}
res <- rbindlist(res)[, q_BH := p.adjust(P, 'BH'), by = subset]
fwrite(res, file.path(D4, 'gse243013_response_tests.csv'))
print(res[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])
log_phase('m4_neoadj_gse243013')
