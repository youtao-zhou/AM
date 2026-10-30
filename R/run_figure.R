run_figure <- function(figure, panels = character()) {
  source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
  scripts <- read.csv(am_code("scripts.csv"), stringsAsFactors = FALSE)
  scripts <- scripts[scripts$figure == figure, ]
  available <- unique(unlist(strsplit(scripts$panels, ";", fixed = TRUE)))
  if (length(panels) && any(!panels %in% available)) {
    stop(
      "Available panels: ", paste(available, collapse = ", "),
      ". See docs/figures.md for panels without a plotting script."
    )
  }
  if (length(panels)) {
    keep <- vapply(
      strsplit(scripts$panels, ";", fixed = TRUE),
      function(x) any(x %in% panels), logical(1)
    )
    scripts <- scripts[keep, ]
  }
  inputs <- read.csv(am_code("data/inputs.csv"), stringsAsFactors = FALSE)
  needed <- unique(inputs$path[inputs$script %in% scripts$script])
  missing <- needed[!file.exists(am_data(needed))]
  if (length(missing)) {
    stop("Missing inputs under ", am_data_root, ":\n",
         paste(missing, collapse = "\n"), "\nSee data/README.md.")
  }
  for (i in seq_len(nrow(scripts))) {
    message("Figure ", figure, ": ", scripts$panels[i])
    status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(am_code(scripts$script[i])))
    if (status != 0) stop("Failed: ", scripts$script[i], " (exit ", status, ")")
  }
}
