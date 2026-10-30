# 03 AM T CellChat refit
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(CellChat); library(future); library(data.table)
})
base <- am_data("HTAN/am/cellchat")
out <- file.path(base, "data/refit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
z <- readRDS(file.path(base, "data/CellChat_full_AM_full100k_T_prefit.rds"))
stopifnot(length(z@idents) == 170985L, sum(z@idents == "AM-2.CD4") == 9495L)
workers <- as.integer(Sys.getenv("CELLCHAT_WORKERS", "8"))
options(future.globals.maxSize = 128 * 1024^3)
plan(multicore, workers = workers)
z <- computeCommunProb(z,
  type = "truncatedMean", trim = .1, raw.use = TRUE,
  population.size = FALSE, nboot = 2000L, seed.use = 60L, Kh = .5, n = 1
)
z <- filterCommunication(z, min.cells = 20L)
plan(sequential)
saveRDS(z, file.path(out, "CellChat_full_AM_full100k_T.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
message("重算完成。新模型在data/refit；正式图未自动替换。")
