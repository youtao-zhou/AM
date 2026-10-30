# AM2CD4 killing T 01 screen
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(parallel)

base_dir <- am_data("HTAN/am/cytosig")
n_perm <- 50000
n_cores <- 16
seed <- 125

states <- c("CD8 CTL", "CD8 Tm", "CD8 Tex", "T-IFN", "CD4 Tem")
genes <- c("GZMB", "PRF1", "GNLY", "NKG7", "GZMA", "IFNG")
min_state_cells <- 10

eight <- data.table(
  panel = 1:8,
  response = c("TNFA", "CXCL12", "CD40L", "IL21", "NO", "IL1B", "IFN1", "IL15"),
  axis = c(
    "CD8 CTL | NKG7 expression", "CD4 Tem | IFNG expression", "CD8 CTL | NKG7 expression",
    "CD8 CTL | IFNG expression", "All CD8 | CD8/Treg ratio (log2)", "CD8 CTL | NKG7 expression",
    "T-IFN | fraction", "CD8 Tm | PRF1 expression"
  ),
  ylab = c(
    "CD8 CTL NKG7 expression", "CD4 Tem IFNG expression", "CD8 CTL NKG7 expression",
    "CD8 CTL IFNG expression", "CD8/Treg ratio (log2)", "CD8 CTL NKG7 expression",
    "T-IFN fraction", "CD8 Tm PRF1 expression"
  )
)

f <- function(name) file.path(base_dir, paste0("AM2CD4_killing_T_", name))
am <- fread(f("input_AM2CD4_CytoSig_sample_means.csv"))
tt <- fread(f("input_T_state_sample_values.csv"))
gg <- fread(f("input_T_state_gene_mean_logRNA.csv.gz"))
responses <- setdiff(names(am), c("sample", "n_AM", "patient", "study", "Dataset2", "platform", "n_T", "scope_labels"))
stopifnot(length(responses) == 43)

frac <- tt[outcome_class == "T-cell composition" & T_group %in% c(states, "Treg"), .(sample, T_group, value)]
wide <- dcast(frac, sample ~ T_group, value.var = "value")
readouts <- rbindlist(list(
  frac[T_group %in% states, .(sample, state = T_group, readout = "fraction", y = 100 * value)],
  tt[
    outcome_class == "Within-state expression" & n_state >= min_state_cells & T_group %in% states & module == "Cytotoxic",
    .(sample, state = T_group, readout = "cytotoxic program", y = value)
  ],
  melt(gg[n_cells >= min_state_cells & T_group %in% states, c("sample", "T_group", genes), with = FALSE],
    id.vars = c("sample", "T_group"), variable.name = "gene", value.name = "y"
  )[
    , .(sample, state = T_group, readout = paste(gene, "expression"), y)
  ],
  wide[, .(sample, state = "All CD8", readout = "fraction", y = 100 * (`CD8 CTL` + `CD8 Tm` + `CD8 Tex`))],
  wide[, .(sample,
    state = "All CD8", readout = "CD8/Treg ratio (log2)",
    y = log2((`CD8 CTL` + `CD8 Tm` + `CD8 Tex` + 0.005) / (Treg + 0.005))
  )]
))
readouts[, axis := paste(state, readout, sep = " | ")]

x_long <- melt(am,
  id.vars = c("sample", "patient", "study"), measure.vars = responses,
  variable.name = "response", value.name = "x"
)
x_long[, response := as.character(response)]
long <- merge(x_long, readouts, by = "sample", allow.cartesian = TRUE)

within_study <- function(x, y, study, test_seed) {
  keep <- study %in% names(which(table(study) >= 4))
  x <- x[keep]; y <- y[keep]; study <- study[keep]
  out <- list(n_within = length(x), n_studies = uniqueN(study), rho_within = NA_real_, p_within = NA_real_)
  if (length(x) < 8 || uniqueN(y) < 3) return(out)
  centred_rank <- function(v) (rank(v) - mean(rank(v))) / length(v)
  rx <- ave(x, study, FUN = centred_rank)
  ry <- ave(y, study, FUN = centred_rank)
  if (sd(rx) == 0 || sd(ry) == 0) return(out)
  obs <- cor(rx, ry)
  groups <- split(seq_along(ry), study)
  set.seed(test_seed)
  perm <- vapply(seq_len(n_perm), function(k) {
    yp <- ry
    for (g in groups) yp[g] <- ry[g][sample.int(length(g))]
    cor(rx, yp)
  }, numeric(1))
  out$rho_within <- obs
  out$p_within <- (1 + sum(abs(perm) >= abs(obs) - 1e-12)) / (n_perm + 1)
  out
}

setorder(long, response, axis, sample)
tests <- unique(long[, .(response, axis)])
tests[, test_id := .I]
long <- merge(long, tests, by = c("response", "axis"))
parts <- split(long, by = "test_id", sorted = TRUE)
screen <- rbindlist(mclapply(seq_len(nrow(tests)), function(i) {
  z <- parts[[as.character(i)]]
  pooled <- suppressWarnings(cor.test(z$x, z$y, method = "spearman"))
  c(
    list(
      response = tests$response[i], axis = tests$axis[i], state = z$state[1], readout = z$readout[1],
      n_samples = nrow(z), n_patients = uniqueN(z$patient),
      rho_pooled = unname(pooled$estimate), p_pooled = pooled$p.value
    ),
    within_study(z$x, z$y, z$study, seed * 1000 + tests$test_id[i])
  )
}, mc.cores = n_cores))
screen[, q_within := p.adjust(fifelse(is.na(p_within), 1, p_within), "BH")]
setorder(screen, p_within)
fwrite(screen, f("screen.csv"))

stats <- merge(eight, screen, by = c("response", "axis"), sort = FALSE)[order(panel)]
stopifnot(nrow(stats) == 8)
fwrite(stats, f("eight_panel_statistics.csv"))

points <- rbindlist(lapply(seq_len(nrow(stats)), function(i) {
  s <- stats[i]
  d <- long[response == s$response & axis == s$axis, .(sample, patient, study, x, y)]
  d <- d[study %in% d[, .N, by = study][N >= 4, study]]
  d[, `:=`(x_rank_in_study = (frank(x) - 0.5) / .N, y_rank_in_study = (frank(y) - 0.5) / .N), by = study]
  d[, group := fifelse(x_rank_in_study >= 0.5, "High", "Low")]
  cbind(s[, .(panel, response, axis, ylab)], d)
}))
stopifnot(points[, .N, by = panel]$N == stats$n_within)
fwrite(points, f("eight_panel_points.csv"))

print(stats[, .(panel, response, axis, n_within,
  rho_within = round(rho_within, 3),
  p_within = signif(p_within, 3), q_within = signif(q_within, 3)
)])
message("Wrote: ", f("screen.csv"), ", ", f("eight_panel_statistics.csv"), ", ", f("eight_panel_points.csv"))
