source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "run_figure.R"))
run_figure(2, commandArgs(trailingOnly = TRUE))
