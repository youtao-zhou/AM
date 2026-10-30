# AM2 CD8 pCR NR zoom extract
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(RANN)

setwd(am_data("tsa"))

obj_dir <- "20260707AM2/ObjectData"
cases <- data.table(
  patient_id = c("D212640", "2009455"),
  response = c("pCR", "NR"),
  file = c(
    "D212640-PCR-AM2-260704.kfbf_job2365Indica Labs - HighPlex FL v4.2.3 20260707 true AM2_object_results.csv",
    "2009455-NR-AM2-260704.kfbf_job2369Indica Labs - HighPlex FL v4.2.3 20260707 true AM2_object_results.csv"
  ),
  x0 = c(23500, 19000),
  y0 = c(21500, 38000)
)
win <- 1000
radius <- 100
am2_col <- "AM2 marophagee"; cd8_col <- "CD8A"

one_case <- function(k) {
  cs <- cases[k]
  raw <- fread(file.path(obj_dir, cs$file),
    showProgress = FALSE,
    select = c("XMin", "XMax", "YMin", "YMax", "all macrophage", "tumor cells", am2_col, cd8_col)
  )
  raw[, `:=`(
    x = (XMin + XMax) / 2, y = (YMin + YMax) / 2,
    mac = get("all macrophage") > 0, tum = get("tumor cells") > 0,
    am2 = get(am2_col) > 0, cd8 = get(cd8_col) > 0
  )]

  raw[, cell_type := fifelse(
    cd8 & !mac & !tum, "CD8",
    fifelse(
      mac & am2 & !cd8, "AM-2",
      fifelse(
        mac & !am2 & !cd8, "Other macrophages",
        fifelse(tum & !mac & !cd8, "Tumor", "Other cells")
      )
    )
  )]
  z <- raw[, .(x, y, cell_type)]; rm(raw)

  anchors <- as.matrix(z[cell_type == "CD8", .(x, y)])
  is_mac <- z$cell_type %in% c("AM-2", "Other macrophages")
  z[, dist_to_CD8_px := NA_real_]
  z[is_mac, dist_to_CD8_px := RANN::nn2(anchors, as.matrix(z[is_mac, .(x, y)]), k = 1)$nn.dists[, 1]]

  inw <- z[x >= cs$x0 & x < cs$x0 + win & y >= cs$y0 & y < cs$y0 + win]
  inw[, `:=`(
    patient_id = cs$patient_id, response = cs$response,
    x_in_window = x - cs$x0, y_in_window = y - cs$y0
  )]
  near <- function(d, tp) 100 * mean(d[cell_type == tp, dist_to_CD8_px] <= radius)
  wt <- data.table(
    patient_id = cs$patient_id, response = cs$response, x0 = cs$x0, y0 = cs$y0, width_px = win,
    n_AM2 = sum(inw$cell_type == "AM-2"), n_other_macrophages = sum(inw$cell_type == "Other macrophages"),
    n_CD8 = sum(inw$cell_type == "CD8"), n_tumor_label = sum(inw$cell_type == "Tumor"),
    n_other_cells = sum(inw$cell_type == "Other cells"),
    window_AM2_within100_CD8_pct = near(inw, "AM-2"),
    window_other_within100_CD8_pct = near(inw, "Other macrophages"),
    section_n_AM2 = sum(z$cell_type == "AM-2"), section_n_other_macrophages = sum(z$cell_type == "Other macrophages"),
    section_n_CD8 = nrow(anchors),
    section_AM2_within100_CD8_pct = near(z, "AM-2"),
    section_other_within100_CD8_pct = near(z, "Other macrophages"),
    source_file = file.path(obj_dir, cs$file)
  )
  list(cells = inw[, .(patient_id, response, cell_type, x, y, x_in_window, y_in_window, dist_to_CD8_px)], window = wt)
}

res <- lapply(seq_len(nrow(cases)), one_case)
cells <- rbindlist(lapply(res, `[[`, "cells"))
wins <- rbindlist(lapply(res, `[[`, "window"))
wins[, `:=`(
  window_diff_pp = window_AM2_within100_CD8_pct - window_other_within100_CD8_pct,
  section_diff_pp = section_AM2_within100_CD8_pct - section_other_within100_CD8_pct
)]

fwrite(cells, "data/AM2_CD8_pCR_NR_zoom_cells.csv")
fwrite(wins, "data/AM2_CD8_pCR_NR_zoom_windows.csv")
print(wins[, .(patient_id, response, n_AM2, n_other_macrophages, n_CD8, n_tumor_label,
  window_AM2 = round(window_AM2_within100_CD8_pct, 1), window_other = round(window_other_within100_CD8_pct, 1),
  section_diff = round(section_diff_pp, 2)
)])
