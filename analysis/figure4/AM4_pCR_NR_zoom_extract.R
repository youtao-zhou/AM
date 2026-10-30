# AM4 pCR NR zoom extract
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(RANN)

setwd(am_data("tsa"))

obj_dir <- "20260707AM4/ObjectData"
cases <- data.table(
  patient_id = c("2007957", "2009455"),
  response = c("pCR", "NR"),
  file = c(
    "2007957-PCR-AM4-260704.kfbf_job2350Indica Labs - HighPlex FL v4.2.3 20260706 True less strict AM4_object_results.csv",
    "2009455-NR-AM4-260705.kfbf_job2345Indica Labs - HighPlex FL v4.2.3 20260706 True less strict AM4_object_results.csv"
  ),
  x0 = c(20500, 53000),
  y0 = c(35500, 43000)
)
win <- 1000
radius <- 100
k_embed <- 3

need <- c(
  "XMin", "XMax", "YMin", "YMax", "DAPI Positive Classification", "GREEN Positive Classification",
  "ORANGE Positive Classification", "RED Positive Classification", "CY-7 Positive Classification"
)

one_case <- function(k) {
  cs <- cases[k]
  z <- fread(file.path(obj_dir, cs$file), select = need, showProgress = FALSE)
  dapi <- z[["DAPI Positive Classification"]] > 0; green <- z[["GREEN Positive Classification"]] > 0
  red <- z[["RED Positive Classification"]] > 0;  cy7 <- z[["CY-7 Positive Classification"]] > 0
  orange <- z[["ORANGE Positive Classification"]] > 0

  mac <- dapi & cy7; am4 <- mac & orange & !green; oth <- mac & !orange & !green
  tum <- dapi & red & !cy7; cd8 <- dapi & green & !mac & !tum
  typ <- fifelse(cd8, "CD8", fifelse(am4, "AM-4", fifelse(oth, "Other macrophages", fifelse(tum, "Tumor", "Other cells"))))
  d <- data.table(x = (z$XMin + z$XMax) / 2, y = (z$YMin + z$YMax) / 2, cell_type = typ)

  pad <- d[x >= cs$x0 - radius & x < cs$x0 + win + radius & y >= cs$y0 - radius & y < cs$y0 + win + radius]
  inw <- pad[, x >= cs$x0 & x < cs$x0 + win & y >= cs$y0 & y < cs$y0 + win]
  tu <- as.matrix(pad[cell_type == "Tumor", .(x, y)])
  pad[, n_tumor_100px := NA_integer_]
  q <- which(inw & pad$cell_type %in% c("AM-4", "Other macrophages", "CD8"))
  if (length(q) && nrow(tu)) {
    nn <- nn2(tu, as.matrix(pad[q, .(x, y)]), k = min(20, nrow(tu)), searchtype = "radius", radius = radius)$nn.idx
    pad[q, n_tumor_100px := rowSums(nn > 0)]
  }
  cells <- pad[inw][, `:=`(patient_id = cs$patient_id, response = cs$response, xp = x - cs$x0, yp = cs$y0 + win - y)]
  win_stats <- cells[cell_type %in% c("AM-4", "Other macrophages", "CD8"),
    .(n = .N, pct_in_tumour = 100 * mean(n_tumor_100px >= k_embed)),
    by = cell_type
  ]
  win_stats <- dcast(win_stats[, .(cell_type, v = sprintf("%d|%.1f", n, pct_in_tumour))], . ~ cell_type, value.var = "v")
  list(
    cells = cells[, .(patient_id, response, cell_type, x, y, xp, yp, n_tumor_100px)],
    window = data.table(
      patient_id = cs$patient_id, response = cs$response, x0 = cs$x0, y0 = cs$y0, win = win,
      n_tumor = cells[cell_type == "Tumor", .N],
      n_am4 = cells[cell_type == "AM-4", .N], n_other_mac = cells[cell_type == "Other macrophages", .N], n_cd8 = cells[cell_type == "CD8", .N],
      pct_am4_in_tumour = cells[cell_type == "AM-4", 100 * mean(n_tumor_100px >= k_embed)],
      pct_other_in_tumour = cells[cell_type == "Other macrophages", 100 * mean(n_tumor_100px >= k_embed)],
      pct_cd8_in_tumour = cells[cell_type == "CD8", 100 * mean(n_tumor_100px >= k_embed)]
    )
  )
}

res <- lapply(seq_len(nrow(cases)), one_case)
fwrite(rbindlist(lapply(res, `[[`, "cells")), "data/AM4_pCR_NR_zoom_cells.csv")
fwrite(rbindlist(lapply(res, `[[`, "window")), "data/AM4_pCR_NR_zoom_windows.csv")
print(rbindlist(lapply(res, `[[`, "window")))
