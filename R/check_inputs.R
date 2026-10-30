source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
inputs <- read.csv(am_code("data/inputs.csv"), stringsAsFactors = FALSE)
missing <- unique(inputs$path[!file.exists(am_data(inputs$path))])
if (length(missing)) {
  cat("Missing files under ", am_data_root, ":\n", sep = "")
  cat(paste0("  ", missing, collapse = "\n"), "\n")
  quit(status = 1)
}
cat(length(unique(inputs$path)), "plotting inputs found.\n")
