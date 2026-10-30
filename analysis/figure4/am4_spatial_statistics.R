# am4 spatial statistics
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(data.table);library(RANN);library(parallel)
})
root <- am_data(""); tsa <- file.path(root, 'tsa')
out <- file.path(root, 'vibe/139_figure4_s5_am4_multimodal')
dir.create(file.path(out, 'data/m5_cache'), recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()
log_phase <- function(s) {
  z <- paste(format(Sys.time(), '%Y-%m-%dT%H:%M:%S%z'), s, round(as.numeric(difftime(Sys.time(), t0, units = 'secs')), 1), sep = '\t');cat(z, '\n');cat(z, '\n', file = file.path(out, 'run.log'), append = TRUE)
}

post <- fread(file.path(tsa, 'data/AM4_post_neoadjuvant_patient_data.csv'), colClasses = list(character = 'patient_id'))
ros_post <- post[, .(patient_id, cohort = 'Post', response = response_group, main = selected21, delivery_source, file = source_file, orange_rule = 'exported')]
old <- fread(file.path(tsa, 'data/AM2_AM4_CD8_pre_post_patient_data.csv'))[cohort == 'Untreated' & endpoint == 'AM-4', .(patient_id, file = file_path)]
new_files <- list.files(file.path(tsa, 'Halo_archive_20261009_new/ObjectData'), pattern = 'object_results\\.csv$', full.names = TRUE)
ros_unt <- rbind(
  old[, .(patient_id, file, orange_rule = 'exported', batch = 'historical_2025')],
  data.table(patient_id = sub('-AM4-.*', '', basename(new_files)), file = new_files, orange_rule = 'harmonized_gt20', batch = 'restained_2026')
)
ros_unt[, `:=`(cohort = 'Untreated', response = 'Untreated', main = patient_id != 'D2407937', delivery_source = batch)]
ros <- rbind(ros_post[, batch := 'post_2026'], ros_unt, use.names = TRUE, fill = TRUE)
stopifnot(nrow(ros_post) == 30, nrow(ros_unt) == 11, all(file.exists(ros$file)))
fwrite(ros, file.path(out, 'data/m5_roster.csv'))
log_phase('phase5_m5_roster')

radii <- c(25, 50, 75, 100, 150, 200)
one <- function(i) {
  r <- ros[i]
  hdr <- names(fread(r$file, nrows = 0))
  if (!'RED Positive Classification' %in% hdr) {
    need <- c('XMin', 'XMax', 'YMin', 'YMax', 'all macrophage', 'tumor cells', 'AM4 marophagee', 'CD8a positive cells')
    z <- fread(r$file, select = need, showProgress = FALSE)
    mac <- z[['all macrophage']] > 0; am4 <- z[['AM4 marophagee']] > 0; tum <- z[['tumor cells']] > 0 & !mac; green <- z[['CD8a positive cells']] > 0
    orange <- am4; cd8c <- green & !mac & !tum
    chk <- data.table(patient_id = r$patient_id, rule_mac = NA, rule_tumor = NA, rule_am4 = NA, rule_cd8 = NA, mode = 'HALO phenotype columns (2025 panel)')
  } else {
    need <- c(
      'XMin', 'XMax', 'YMin', 'YMax', 'DAPI Positive Classification', 'GREEN Positive Classification', 'ORANGE Positive Classification',
      'RED Positive Classification', 'CY-7 Positive Classification', 'all macrophage', 'tumor cells', 'AM4 marophagee', 'CD8A'
    )
    opt <- intersect('ORANGE Cytoplasm Intensity', hdr)
    miss <- setdiff(need, hdr); if (length(miss)) stop(r$patient_id, ': missing ', paste(miss, collapse = ';'))
    z <- fread(r$file, select = c(need, opt), showProgress = FALSE)
    dapi <- z[['DAPI Positive Classification']] > 0; green <- z[['GREEN Positive Classification']] > 0; red <- z[['RED Positive Classification']] > 0
    cy7 <- z[['CY-7 Positive Classification']] > 0; or_exp <- z[['ORANGE Positive Classification']] > 0
    chk <- data.table(
      patient_id = r$patient_id, rule_mac = identical(dapi & cy7, z[['all macrophage']] > 0), rule_tumor = identical(dapi & red & !cy7, z[['tumor cells']] > 0),
      rule_am4 = identical(dapi & cy7 & or_exp, z[['AM4 marophagee']] > 0), rule_cd8 = identical(dapi & green, z$CD8A > 0), mode = 'channel rules'
    )
    orange <- if (r$orange_rule == 'harmonized_gt20') {
      stopifnot(length(opt) == 1); z[['ORANGE Cytoplasm Intensity']] > 20
    } else or_exp
    mac <- dapi & cy7; am4 <- mac & orange; tum <- dapi & red & !cy7; cd8c <- dapi & green & !mac & !tum
  }
  if (!'AM4 marophagee' %in% names(z)) stop('no AM4 column')
  x <- (z$XMin + z$XMax) / 2; y <- (z$YMin + z$YMax) / 2
  q_am4 <- am4 & !green; q_oth <- mac & !orange & !green
  typ <- fifelse(cd8c, 'CD8', fifelse(q_am4, 'AM-4', fifelse(q_oth, 'Other macrophages', fifelse(tum, 'Tumor', 'Other cells'))))
  cache <- data.table(x, y, cell_type = factor(typ, c('Other cells', 'Tumor', 'Other macrophages', 'CD8', 'AM-4')))
  qm <- q_am4 | q_oth
  nnd <- function(a) if (sum(a) < 1) rep(NA_real_, sum(qm)) else RANN::nn2(cbind(x[a], y[a]), cbind(x[qm], y[qm]), k = 1)$nn.dists[, 1]
  dt <- nnd(tum); dc <- nnd(cd8c)
  ia <- q_am4[qm]
  cache[qm, `:=`(d_tumor = dt, d_cd8 = dc)]
  saveRDS(cache, file.path(out, 'data/m5_cache', paste0(r$patient_id, '_', r$cohort, '.rds')), compress = FALSE)
  gm <- function(d, a) 100 * exp(mean(log1p(d[a])) - mean(log1p(d)))
  prox <- rbindlist(lapply(radii, function(rr) data.table(
    radius = rr, anchor = c('Tumor', 'Tumor', 'CD8', 'CD8'), population = c('AM-4', 'Other macrophages', 'AM-4', 'Other macrophages'),
    percent = 100 * c(mean(dt[ia] <= rr), mean(dt[!ia] <= rr), mean(dc[ia] <= rr), mean(dc[!ia] <= rr))
  )))
  prox[, patient_id := r$patient_id]
  summ <- data.table(
    patient_id = r$patient_id, n_total = nrow(z), n_mac = sum(mac), n_am4 = sum(am4), n_am4_halo_exported = sum(z[['AM4 marophagee']] > 0),
    n_tumor = sum(tum), n_cd8 = sum(cd8c), n_q_am4 = sum(q_am4), n_q_oth = sum(q_oth),
    am4_pct_mac = 100 * sum(am4) / sum(mac), tumor_pct = 100 * sum(tum) / nrow(z), cd8_pct = 100 * sum(cd8c) / nrow(z),
    rel_dist_tumor = gm(dt, ia), rel_dist_cd8 = gm(dc, ia),
    med_d_tumor_am4 = median(dt[ia]), med_d_tumor_oth = median(dt[!ia]), med_d_cd8_am4 = median(dc[ia]), med_d_cd8_oth = median(dc[!ia])
  )
  list(summ = summ, prox = prox, chk = chk)
}
res <- mclapply(seq_len(nrow(ros)), function(i) tryCatch(one(i), error = function(e) e), mc.cores = 12, mc.preschedule = FALSE)
bad <- which(vapply(res, inherits, logical(1), 'error'))
if (length(bad)) {
  for (b in bad) cat('ERROR', ros$patient_id[b], conditionMessage(res[[b]]), '\n', file = file.path(out, 'audit/failures.md'), append = TRUE); res <- res[-bad]; ros_ok <- ros[-bad]
} else ros_ok <- ros
summ <- merge(rbindlist(lapply(res, `[[`, 'summ')), ros_ok[, .(patient_id, cohort, response, main, delivery_source, batch, orange_rule)], by = 'patient_id')
prox <- merge(rbindlist(lapply(res, `[[`, 'prox')), ros_ok[, .(patient_id, cohort, response, main, delivery_source)], by = 'patient_id')
chk <- rbindlist(lapply(res, `[[`, 'chk'))
fwrite(chk, file.path(out, 'audit/m5_halo_rule_check.csv'))

v <- merge(summ[cohort == 'Post', .(patient_id, n_am4, n_mac)], post[, .(patient_id, am4_count, allmac_count)], by = 'patient_id')
fwrite(v[, match := n_am4 == am4_count & n_mac == allmac_count], file.path(out, 'audit/m5_post_count_reproduction.csv'))
fwrite(summ, file.path(out, 'data/m5_patient_summary.csv')); fwrite(prox, file.path(out, 'data/m5_patient_proximity.csv'))
cat('rule checks all TRUE:', all(unlist(chk[mode == 'channel rules', .(rule_mac, rule_tumor, rule_am4, rule_cd8)])), ' post counts reproduce:', all(v$match), ' failures:', length(bad), '\n')
log_phase('phase5_m5_extracted')

strat_p <- function(d, y, g1, g2) {
  x <- d[response %in% c(g1, g2)]; yy <- x[[y]]; sc <- numeric(nrow(x)); al <- list()
  for (ids in split(seq_len(nrow(x)), x$delivery_source)) {
    n <- length(ids);k <- sum(x$response[ids] == g1);sc[ids] <- (rank(yy[ids]) - (n + 1) / 2) / (n + 1)
    al[[length(al) + 1]] <- if (k == 0) 0 else if (k == n) sum(sc[ids]) else combn(ids, k, function(j) sum(sc[j]))
  }
  nl <- Reduce(function(a, b) as.vector(outer(a, b, '+')), al); o <- sum(sc[x$response == g1]); mean(abs(nl) >= abs(o) - 1e-12)
}
w <- dcast(prox, patient_id + cohort + main + response + radius + anchor ~ population, value.var = 'percent')
setnames(w, c('AM-4', 'Other macrophages'), c('am4', 'oth'))
paired <- rbindlist(lapply(list(c('Post', TRUE, 'post21'), c('Post', NA, 'post30'), c('Untreated', TRUE, 'untreated10'), c('Untreated', NA, 'untreated11')), function(s) {
  z <- w[cohort == s[1] & !is.na(am4) & !is.na(oth)]; if (!is.na(s[2])) z <- z[main == TRUE]
  z[, .(
    scope = s[3], n = .N, median_am4 = median(am4), median_oth = median(oth), median_diff_pp = median(am4 - oth), n_am4_higher = sum(am4 > oth),
    p_wilcoxon_paired = wilcox.test(am4, oth, paired = TRUE, exact = .N <= 30)$p.value
  ), by = .(anchor, radius)]
}))
paired[, q_BH_within_scope := p.adjust(p_wilcoxon_paired, 'BH'), by = scope]
fwrite(paired, file.path(out, 'data/m5_proximity_paired_tests.csv'))
s21 <- summ[cohort == 'Post' & main == TRUE]
resp <- rbindlist(lapply(c('am4_pct_mac', 'rel_dist_tumor', 'rel_dist_cd8'), function(m)
  data.table(
    metric = m, pair = c('NR-MPR', 'MPR-PCR', 'NR-PCR'),
    p_stratified = c(strat_p(s21, m, 'NR', 'MPR'), strat_p(s21, m, 'MPR', 'PCR'), strat_p(s21, m, 'NR', 'PCR')),
    median_NR = median(s21[response == 'NR'][[m]]), median_MPR = median(s21[response == 'MPR'][[m]]), median_PCR = median(s21[response == 'PCR'][[m]])
  )))
fwrite(resp, file.path(out, 'data/m5_response_tests.csv'))
print(paired[radius == 100]); print(resp)
log_phase('phase5_m5_stats')
