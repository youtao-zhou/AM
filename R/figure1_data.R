INPUT_ROOT <- am_data("Figure1/plot_inputs")
SOURCE_ROOT <- am_data("Figure1/source_data")
OUT_ROOT <- am_output("Figure1")
DERIVED_ROOT <- am_output("Figure1", "data")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(DERIVED_ROOT, recursive = TRUE, showWarnings = FALSE)

read_source <- function(key) {
  read.csv(gzfile(file.path(SOURCE_ROOT, paste0(key, ".csv.gz"))), check.names = FALSE)
}
