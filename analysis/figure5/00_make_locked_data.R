# 00 make locked data
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)

setwd(am_data("HTAN/am/TWAS_signature"))
dir.create("data", showWarnings = FALSE)
log_phase <- function(phase) cat(sprintf("%s\t%s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), phase),
  file = am_data("vibe/197_fig5_user_scripts_locked/run.log"), append = TRUE
)

s <- fread(am_data("vibe/186_fig5_clinical_lock_sap/data/locked_scores.csv"))
stopifnot(nrow(s) == 421, all(s$twas_hi == (s$TWASz > 0)), all((s$grp == "AM score-high") == (s$AM > 0)))
base <- s[, .(cohort = trial, sample, arm, pfs_months = pfs_t, pfs_event = pfs_e)]

fwrite(cbind(base, s[, .(os_months = os_t, os_event = os_e, TWAS_z = TWASz, TWAS_group = fifelse(TWASz > 0, "High", "Low"))]), "data/Fig5c_patients.csv")

fwrite(cbind(base, s[, .(
  objective_response = resp, AM2_z = zAM2, AM4_z = zAM4, AM_score = AM,
  AM_group = fifelse(AM > 0, "AM score-high", "AM score-low")
)]), "data/Fig5cd_patients.csv")

fwrite(cbind(base, s[, .(AM_z = AMs, TWAS_z = TWASz)]), "data/Fig5f_patients.csv")
cat("high/low:", sum(s$AM > 0), "/", sum(s$AM <= 0), "; TWAS high:", sum(s$TWASz > 0), "\n")
log_phase("data")
