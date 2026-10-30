am_repo <- normalizePath(Sys.getenv("AM_REPO_DIR", unset = "."), mustWork = TRUE)
Sys.setenv(AM_REPO_DIR = am_repo)
am_data_root <- normalizePath(
  Sys.getenv("AM_DATA_DIR", unset = file.path(am_repo, "data")),
  mustWork = FALSE
)
am_output_root <- normalizePath(
  Sys.getenv("AM_OUTPUT_DIR", unset = file.path(am_repo, "output")),
  mustWork = FALSE
)

am_data <- function(...) file.path(am_data_root, ...)
am_code <- function(...) file.path(am_repo, ...)
am_output <- function(figure, file = "") {
  path <- file.path(am_output_root, figure, file)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}
