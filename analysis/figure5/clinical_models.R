# clinical models
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(data.table); library(survival); library(ggplot2)
})
out_dir <- am_data("vibe/186_fig5_clinical_lock_sap"); D <- file.path(out_dir, "data")
dir.create(D, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
log_phase <- function(phase) {
  z <- sprintf("%s\t%s\t%.1f", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), phase, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  cat(z, "\n"); cat(z, "\n", file = file.path(out_dir, "run.log"), append = TRUE)
}
H <- am_data("vibe/141_figure5_s6_clinical_prediction/data/harmonized")
EGA <- am_data("HTAN/EGA")

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
O <- readRDS(file.path(H, "OAK.rds")); P <- readRDS(file.path(H, "POPLAR.rds"))
G <- intersect(rownames(O$expr), rownames(P$expr))
EO <- O$expr[G, ]; EP <- P$expr[G, ]
gs <- readRDS(file.path(H, "gene_sets.rds"))
tw <- intersect(fread(am_data("vibe/145_figure5_clinical_deepdive/data/E0c_signature_genes.csv"))$gene, G)

set.seed(1)
avg <- sort(rowMeans(EO)); cutb <- ggplot2::cut_number(avg + rnorm(length(avg)) / 1e30, n = 24, labels = FALSE, right = FALSE); names(cutb) <- names(avg)
ctrl <- unique(unlist(lapply(tw, function(g) names(sample(cutb[cutb == cutb[g]], 100, replace = FALSE)))))
twas_raw <- function(E) colMeans(E[tw, , drop = FALSE]) - colMeans(E[ctrl, , drop = FALSE])
teff_raw <- function(E, ref) {
  g <- intersect(gs$Teff8, rownames(E)); m <- rowMeans(ref[g, ]); s <- apply(ref[g, ], 1, sd); colMeans((E[g, ] - m) / s)
}
raw <- function(o, E) {
  re <- rank_enrichment(E, gs[c("AM2", "AM4", "AM2_CD4")])
  data.table(
    sample = colnames(E), AM2r = re[, "AM2"], AM4r = re[, "AM4"], AM2CD4r = re[, "AM2_CD4"], TWASr = twas_raw(E),
    Teffr = teff_raw(E, EO), CXr = E["CXCL9", ] - E["SPP1", ]
  )
}
S <- rbind(cbind(trial = "OAK", raw(O, EO)), cbind(trial = "POPLAR", raw(P, EP)))
ref <- S[trial == "OAK"]
zo <- function(x, col) (x - mean(ref[[col]])) / sd(ref[[col]])
S[, `:=`(
  zAM2 = zo(AM2r, "AM2r"), zAM4 = zo(AM4r, "AM4r"), zAM2CD4 = zo(AM2CD4r, "AM2CD4r"), TWASz = zo(TWASr, "TWASr"),
  Teff8 = zo(Teffr, "Teffr"), CX = zo(CXr, "CXr")
)]
S[, AM := zAM2 - zAM4]; sdAM <- sd(S[trial == "OAK", AM]); S[, AMs := AM / sdAM]
cl <- rbind(cbind(trial = "OAK", O$clin), cbind(trial = "POPLAR", P$clin), fill = TRUE)
d <- merge(S, cl[, .(sample, arm, trt = io, pfs_t = pfs_time, pfs_e = pfs_event, os_t = os_time, os_e = os_event, resp, bor, pdl1, tmb)], by = "sample")
d[, arm := fifelse(trt == 1, "Atezolizumab", "Docetaxel")][, grp := factor(fifelse(AM > 0, "AM score-high", "AM score-low"), c("AM score-high", "AM score-low"))]
d[, low := as.integer(AM <= 0)][, twas_hi := as.integer(TWASz > 0)]
stopifnot(nrow(d) == 421, d[trial == "OAK", .N] == 319)
fwrite(d, file.path(D, "locked_scores.csv"))

ob <- fread(am_data("vibe/161_figure5_v5_am_balance/data/patients_balance.csv"))
cc <- merge(d, ob[, .(sample, AM_orig = balance, TWAS_orig = TWASz, Teff_orig = Teff8z)], by = "sample")
conc <- cc[, .(
  n = .N, r_AM = cor(AM, AM_orig), r_TWAS = cor(TWASz, TWAS_orig), r_Teff8 = cor(Teff8, Teff_orig),
  group_agreement = mean((AM > 0) == (AM_orig > 0)), n_high_locked = sum(AM > 0), n_high_orig = sum(AM_orig > 0),
  mean_AM_in_OAK_SD = mean(AMs), mean_TWASz = mean(TWASz), sd_TWASz = sd(TWASz)
), by = trial]
conc[, `:=`(n_twas_genes = length(tw), n_ctrl_genes = length(ctrl), n_common_genes = length(G))]
fwrite(conc, file.path(D, "concordance.csv")); print(conc)
log_phase("lock_scores")

hr <- function(f, term, extra) {
  s <- summary(f); co <- s$coefficients; ci <- s$conf.int
  data.table(term = term, HR = ci[term, 1], lo = ci[term, 3], hi = ci[term, 4], P = co[term, "Pr(>|z|)"], n = f$n, events = f$nevent)[, names(extra) := extra]
}
R <- list(); add <- function(x) R[[length(R) + 1]] <<- x
ep <- list(PFS = c("pfs_t", "pfs_e"), OS = c("os_t", "os_e"))
fsurv <- function(e, rhs) as.formula(sprintf("Surv(%s, %s) ~ %s", ep[[e]][1], ep[[e]][2], rhs))

for (e in names(ep)) {
  f <- coxph(fsurv(e, "trt * AMs + strata(trial)"), d)
  add(hr(f, "trt:AMs", list(analysis = ifelse(e == "PFS", "PRIMARY", "S2"), endpoint = e, trial = "Pooled", model = "trt x AM score (continuous, per OAK SD)")))
  for (g in levels(d$grp)) add(hr(
    coxph(fsurv(e, "trt + strata(trial)"), d[grp == g]), "trt",
    list(analysis = ifelse(e == "PFS", "S1", "S2"), endpoint = e, trial = "Pooled", model = paste0(g, ": atezolizumab vs docetaxel"))
  ))
  add(hr(
    coxph(fsurv(e, "trt * low + strata(trial)"), d), "trt:low",
    list(analysis = ifelse(e == "PFS", "S1", "S2"), endpoint = e, trial = "Pooled", model = "trt x AM group (ratio of HRs, low vs high)")
  ))
  for (tr in c("OAK", "POPLAR")) {
    dd <- d[trial == tr]
    add(hr(coxph(fsurv(e, "trt * AMs"), dd), "trt:AMs", list(analysis = "S6", endpoint = e, trial = tr, model = "trt x AM score (continuous, per OAK SD)")))
    for (g in levels(d$grp)) add(hr(coxph(fsurv(e, "trt"), dd[grp == g]), "trt", list(analysis = "S6", endpoint = e, trial = tr, model = paste0(g, ": atezolizumab vs docetaxel"))))
    add(hr(coxph(fsurv(e, "trt * low"), dd), "trt:low", list(analysis = "S6", endpoint = e, trial = tr, model = "trt x AM group (ratio of HRs, low vs high)")))
  }
  fh <- coxph(fsurv(e, "trt + AMs + trt:AMs + trt:trial + AMs:trial + trt:AMs:trial + strata(trial)"), d)
  add(data.table(
    term = "trt:AMs:trialPOPLAR", P = summary(fh)$coefficients["trt:AMs:trialPOPLAR", "Pr(>|z|)"], n = fh$n,
    analysis = "S6", endpoint = e, trial = "OAK vs POPLAR", model = "heterogeneity of trt x AM score between trials"
  ))
}
log_phase("primary")

dr <- d[!is.na(resp)]
g3 <- glm(resp ~ trt * AMs + trial, family = binomial, data = dr); cg <- summary(g3)$coefficients
add(data.table(
  term = "trt:AMs", HR = exp(cg["trt:AMs", 1]), lo = exp(cg["trt:AMs", 1] - 1.96 * cg["trt:AMs", 2]), hi = exp(cg["trt:AMs", 1] + 1.96 * cg["trt:AMs", 2]),
  P = cg["trt:AMs", 4], n = nrow(dr), analysis = "S3", endpoint = "ORR", trial = "Pooled", model = "trt x AM score, logistic (odds ratio)"
))
orr <- dr[, .(responders = sum(resp), evaluable = .N, ORR = mean(resp)), by = .(grp, arm)][order(grp, arm)]
orr[, fisher_P := {
  g <- grp; sapply(g, function(gg) {
    x <- dr[grp == gg]; fisher.test(table(x$trt, x$resp))$p.value
  })
}]
fwrite(orr, file.path(D, "orr_by_group.csv"))
log_phase("secondary")

ab <- list()
for (e in names(ep)) for (g in levels(d$grp)) {
  x <- d[grp == g]
  for (a in c("Atezolizumab", "Docetaxel")) {
    sf <- survfit(fsurv(e, "1"), x[arm == a]); s12 <- summary(sf, times = 12, extend = TRUE)
    ab[[length(ab) + 1]] <- data.table(
      endpoint = e, grp = g, arm = a, n = sf$n, median = unname(summary(sf)$table["median"]),
      rate12 = s12$surv, rate12_se = s12$std.err
    )
  }
}
ab <- rbindlist(ab)
dif <- ab[, .(
  median_atezo = median[arm == "Atezolizumab"], median_doce = median[arm == "Docetaxel"],
  rate12_diff = rate12[arm == "Atezolizumab"] - rate12[arm == "Docetaxel"],
  rate12_diff_se = sqrt(sum(rate12_se^2))
), by = .(endpoint, grp)]
rm <- list()
for (e in names(ep)) for (g in levels(d$grp)) for (tau in c(12, 18)) {
  per <- rbindlist(lapply(c("OAK", "POPLAR"), function(tr) {
    x <- d[grp == g & trial == tr]
    tb <- summary(survfit(fsurv(e, "arm"), x), rmean = tau)$table
    data.table(
      trial = tr, diff = tb["arm=Atezolizumab", "rmean"] - tb["arm=Docetaxel", "rmean"],
      se = sqrt(tb["arm=Atezolizumab", "se(rmean)"]^2 + tb["arm=Docetaxel", "se(rmean)"]^2)
    )
  }))
  w <- 1 / per$se^2; est <- sum(w * per$diff) / sum(w); se <- sqrt(1 / sum(w))
  rm[[length(rm) + 1]] <- data.table(
    endpoint = e, grp = g, tau = tau, RMST_diff = est, lo = est - 1.96 * se, hi = est + 1.96 * se,
    P = 2 * pnorm(-abs(est / se)), OAK_diff = per$diff[1], POPLAR_diff = per$diff[2]
  )
}
rm <- rbindlist(rm)

rmi <- rm[, .(
  RMST_diff_low_minus_high = RMST_diff[grp == "AM score-low"] - RMST_diff[grp == "AM score-high"],
  se = sqrt((((hi - lo) / 3.92)^2)[1] + (((hi - lo) / 3.92)^2)[2])
), by = .(endpoint, tau)][, P := 2 * pnorm(-abs(RMST_diff_low_minus_high / se))]
fwrite(ab, file.path(D, "km_medians_rates.csv")); fwrite(dif, file.path(D, "absolute_effects.csv")); fwrite(rm, file.path(D, "rmst.csv")); fwrite(rmi, file.path(D, "rmst_interaction.csv"))
ph <- rbindlist(lapply(names(ep), function(e) {
  z1 <- cox.zph(coxph(fsurv(e, "trt * AMs + strata(trial)"), d))$table
  z2 <- cox.zph(coxph(fsurv(e, "trt * low + strata(trial)"), d))$table
  rbind(
    data.table(endpoint = e, model = "trt x AM score (continuous)", term = rownames(z1), chisq = z1[, 1], P = z1[, 3]),
    data.table(endpoint = e, model = "trt x AM group", term = rownames(z2), chisq = z2[, 1], P = z2[, 3])
  )
}))
fwrite(ph, file.path(D, "ph_tests.csv")); print(ph)
log_phase("absolute_ph")

fo <- coxph(Surv(pfs_t, pfs_e) ~ trt * TWASz + trt * AMs, d[trial == "OAK"]); b <- coef(fo)
coefs <- data.table(term = c("trt", "trt:TWASz", "trt:AMs"), locked = b[c("trt", "trt:TWASz", "trt:AMs")], original = c(0.077, 0.198, -0.233))
fwrite(coefs, file.path(D, "classifier_coefficients.csv")); print(coefs)
d[, pred_locked := b["trt"] + b["trt:TWASz"] * TWASz + b["trt:AMs"] * AMs]
d[, pred_origcoef := 0.077 + 0.198 * TWASz - 0.233 * AMs]
for (v in c("pred_locked", "pred_origcoef")) for (tr in c("OAK", "POPLAR")) {
  x <- d[trial == tr][, pg := as.integer(get(v) > 0)]
  for (k in 0:1) add(hr(
    coxph(Surv(pfs_t, pfs_e) ~ trt, x[pg == k]), "trt",
    list(analysis = "S7", endpoint = "PFS", trial = tr, model = paste0(v, ": predicted HR ", ifelse(k == 0, "< 1", "> 1"), ", atezolizumab vs docetaxel"))
  ))
  add(hr(coxph(Surv(pfs_t, pfs_e) ~ trt * pg, x), "trt:pg", list(analysis = "S7", endpoint = "PFS", trial = tr, model = paste0(v, ": trt x predicted group"))))
  x[, pl := get(v)]
  add(hr(coxph(Surv(pfs_t, pfs_e) ~ trt * pl, x), "trt:pl", list(analysis = "S7", endpoint = "PFS", trial = tr, model = paste0(v, ": trt x predicted logHR (continuous)"))))
}
d[, pred_group := fifelse(pred_locked > 0, "Predicted HR > 1", "Predicted HR < 1")]
log_phase("classifier")

tk <- list(); km <- list()
for (tr in c("OAK", "POPLAR")) for (a in c("Atezolizumab", "Docetaxel")) {
  x <- d[trial == tr & arm == a]
  f <- coxph(Surv(pfs_t, pfs_e) ~ twas_hi, x); lr <- survdiff(Surv(pfs_t, pfs_e) ~ twas_hi, x)
  tk[[length(tk) + 1]] <- data.table(
    trial = tr, arm = a, n = nrow(x), n_high = sum(x$twas_hi), HR = exp(coef(f)), lo = summary(f)$conf.int[1, 3],
    hi = summary(f)$conf.int[1, 4], logrank_P = 1 - pchisq(lr$chisq, 1)
  )
  sf <- survfit(Surv(pfs_t, pfs_e) ~ twas_hi, x); ss <- summary(sf, censored = TRUE)
  km[[length(km) + 1]] <- data.table(
    trial = tr, arm = a, twas = fifelse(grepl("=1", ss$strata), "High", "Low"), time = ss$time, surv = ss$surv,
    n.risk = ss$n.risk, n.event = ss$n.event, n.censor = ss$n.censor, lower = ss$lower, upper = ss$upper
  )
}
fwrite(rbindlist(tk), file.path(D, "twas_km.csv")); fwrite(rbindlist(km), file.path(D, "twas_km_steps.csv"))
add(hr(coxph(Surv(pfs_t, pfs_e) ~ trt * TWASz + strata(trial), d), "trt:TWASz", list(analysis = "S8", endpoint = "PFS", trial = "Pooled", model = "trt x TWAS score (continuous)")))
add(hr(coxph(Surv(pfs_t, pfs_e) ~ trt * twas_hi + strata(trial), d), "trt:twas_hi", list(analysis = "S8", endpoint = "PFS", trial = "Pooled", model = "trt x TWAS group (high vs low)")))
for (a in c("Atezolizumab", "Docetaxel")) add(hr(
  coxph(Surv(pfs_t, pfs_e) ~ twas_hi + strata(trial), d[arm == a]), "twas_hi",
  list(analysis = "S8", endpoint = "PFS", trial = "Pooled", model = paste0(a, " arm: TWAS high vs low"))
))
for (a in c("Atezolizumab", "Docetaxel")) add(hr(
  coxph(Surv(pfs_t, pfs_e) ~ TWASz + AMs + strata(trial), d[arm == a]), c("TWASz", "AMs")[1],
  list(analysis = "S8b", endpoint = "PFS", trial = "Pooled", model = paste0(a, " arm: TWAS score per SD, adjusted for AM score"))
))
for (a in c("Atezolizumab", "Docetaxel")) add(hr(
  coxph(Surv(pfs_t, pfs_e) ~ TWASz + AMs + strata(trial), d[arm == a]), "AMs",
  list(analysis = "S8b", endpoint = "PFS", trial = "Pooled", model = paste0(a, " arm: AM score per SD, adjusted for TWAS score"))
))
log_phase("twas")

lrt <- function(dd, base, lab) {
  f0 <- coxph(as.formula(paste("Surv(pfs_t, pfs_e) ~", base)), dd); f1 <- coxph(as.formula(paste("Surv(pfs_t, pfs_e) ~", base, "+ AMs + trt:AMs")), dd)
  stopifnot(f0$n == f1$n)
  a <- anova(f0, f1); s <- summary(f1)
  data.table(
    analysis = "S9", cohort = lab, base_model = base, n = f1$n, events = f1$nevent, LRT_chisq = a$Chisq[2], df = a$Df[2], LRT_P = a$`Pr(>|Chi|)`[2],
    trtxAM_HR = s$conf.int["trt:AMs", 1], trtxAM_lo = s$conf.int["trt:AMs", 3], trtxAM_hi = s$conf.int["trt:AMs", 4], trtxAM_P = s$coefficients["trt:AMs", 5]
  )
}
d[, pdl1_pos := fifelse(is.na(pdl1), NA_integer_, as.integer(pdl1 %in% c("[1,50)", "[50,100]")))][, tmb_hi := fifelse(is.na(tmb), NA_integer_, as.integer(tmb == ">=16"))]
inc <- rbind(
  lrt(d[is.finite(Teff8) & is.finite(CX)], "trt * Teff8 + trt * CX + strata(trial)", "OAK + POPLAR, all"),
  lrt(d[trial == "OAK" & !is.na(pdl1_pos)], "trt * pdl1_pos + trt * Teff8", "OAK, PD-L1 available"),
  lrt(d[trial == "OAK" & !is.na(tmb_hi)], "trt * tmb_hi + trt * Teff8", "OAK, TMB available")
)
fb <- coxph(Surv(pfs_t, pfs_e) ~ trt * pdl1_pos + trt * tmb_hi + trt * Teff8 + trt * AMs, d[trial == "OAK" & !is.na(pdl1_pos) & !is.na(tmb_hi)])
inc <- rbind(inc, data.table(
  analysis = "S9 (descriptive)", cohort = "OAK, PD-L1 and TMB available", base_model = "trt*pdl1_pos + trt*tmb_hi + trt*Teff8 + trt*AMs",
  n = fb$n, events = fb$nevent, trtxAM_HR = summary(fb)$conf.int["trt:AMs", 1], trtxAM_lo = summary(fb)$conf.int["trt:AMs", 3],
  trtxAM_hi = summary(fb)$conf.int["trt:AMs", 4], trtxAM_P = summary(fb)$coefficients["trt:AMs", 5]
), fill = TRUE)
fwrite(inc, file.path(D, "incremental_lrt.csv")); print(inc)
log_phase("incremental")

b699 <- fread(file.path(EGA, "OAK/BYN_n699.csv")); g699 <- fread(file.path(EGA, "OAK/go28915_updated.csv"))
b192 <- fread(file.path(EGA, "POPLAR/BYN_n192.csv"))

bmk <- fread(file.path(EGA, "OAK/biomarker.txt"))[, UNI_ID := sub("^[^-]+-", "", full_anonymized_patient_id)]
lk <- merge(merge(g699, bmk[, .(trunc_UNI_ID = trunc_anonymized_patient_id, UNI_ID)], by = "trunc_UNI_ID"), b699, by = "UNI_ID")
flow <- rbind(
  data.table(trial = "OAK", step = c("Randomised (ITT, trial report)", "Primary analysis population (trial report)"), n = c(1225, 850)),
  data.table(trial = "OAK", step = "Pre-treatment tumour RNA-seq (biomarker-evaluable)", n = nrow(b699)),
  data.table(trial = "OAK", step = "Atezolizumab or docetaxel arm", n = b699[ACTARM %in% c("MPDL3280A", "Docetaxel"), .N]),
  data.table(trial = "OAK", step = "Non-squamous histology", n = b699[ACTARM %in% c("MPDL3280A", "Docetaxel") & HIST == "NON-SQUAMOUS", .N]),
  data.table(trial = "OAK", step = "Specimen from lung (primary site)", n = lk[ACTARM %in% c("MPDL3280A", "Docetaxel") & HIST == "NON-SQUAMOUS" & TISSUE_TYPE.2 == "lung", .N]),
  data.table(trial = "OAK", step = "Analysed (expression and outcome available)", n = d[trial == "OAK", .N]),
  data.table(trial = "POPLAR", step = "Randomised (ITT, trial report)", n = 287),
  data.table(trial = "POPLAR", step = "Pre-treatment tumour RNA-seq (biomarker-evaluable)", n = nrow(b192)),
  data.table(trial = "POPLAR", step = "Non-squamous histology", n = b192[HIST == "NON-SQUAMOUS", .N]),
  data.table(trial = "POPLAR", step = "Analysed (expression and outcome available)", n = d[trial == "POPLAR", .N])
)
fwrite(flow, file.path(D, "flow_counts.csv")); print(flow)
log_phase("flow")

res <- rbindlist(R, fill = TRUE); fwrite(res, file.path(D, "primary_secondary.csv"))
fwrite(d, file.path(D, "locked_scores.csv"))
print(res[, .(analysis, endpoint, trial, model, n, HR = round(HR, 3), lo = round(lo, 3), hi = round(hi, 3), P = signif(P, 3))])
print(rm); print(rmi); print(dif); print(orr)
log_phase("export")
