# AM2 AM4 tumor proximity untreated extract
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(data.table); library(RANN)
})
tsa <- am_data("tsa")
radii <- c(25, 50, 100, 200)
old <- c("D2407055", "D2408148", "D2411231")
new <- c("D2406926", "D2407910", "D2409692", "D2411944", "D2413695", "D2414554", "D2414722")
code <- c(
  D2406926 = "U01", D2407055 = "U02", D2407910 = "U03", D2408148 = "U04", D2409692 = "U05", D2411231 = "U06",
  D2411944 = "U07", D2413695 = "U08", D2414554 = "U09", D2414722 = "U10"
)
obj <- function(dir, id) {
  f <- list.files(file.path(tsa, dir, "ObjectData"), paste0("^", id, ".*object_results\\.csv$"), full.names = TRUE); stopifnot(length(f) == 1); f
}

cells <- function(id, panel) {
  if (id %in% old) {
    f <- obj(if (panel == "AM-2") "Halo archive 2025-12-08 11-51 - v3.5.3577" else "Halo archive 2025-12-05 16-27 - v3.5.3577", id)
    amcol <- if (panel == "AM-2") "AM2 marophagee" else "AM4 marophagee"
    d <- fread(f, select = c("XMin", "XMax", "YMin", "YMax", "all macrophage", "tumor cells", amcol, "CD8a positive cells"), showProgress = FALSE)
    mac <- d[["all macrophage"]] > 0; am <- d[[amcol]] > 0; cd8 <- d[["CD8a positive cells"]] > 0; tum <- d[["tumor cells"]] > 0 & !mac
  } else {
    f <- obj(if (panel == "AM-2") "Halo_archive_20261009_new2" else "Halo_archive_20261009_new", id)
    hdr <- names(fread(f, nrows = 0))
    d <- fread(f, select = c(
      "XMin", "XMax", "YMin", "YMax", "CY-7 Cytoplasm Intensity", "ORANGE Cytoplasm Intensity",
      grep("Positive Classification$", hdr, value = TRUE)
    ), showProgress = FALSE)
    b <- function(ch) d[[paste(ch, "Positive Classification")]] > 0
    if (panel == "AM-2") {
      cy7h <- d[["CY-7 Cytoplasm Intensity"]] > 20
      mac <- b("DAPI") & cy7h; am <- mac & b("ORANGE") & b("CY-5"); cd8 <- b("DAPI") & b("GREEN") & !cy7h; tum <- b("DAPI") & b("RED") & !cy7h
    } else {
      cy7 <- b("CY-7")
      mac <- b("DAPI") & cy7; am <- mac & d[["ORANGE Cytoplasm Intensity"]] > 20; cd8 <- b("GREEN"); tum <- b("DAPI") & b("RED") & !cy7
    }
  }
  qa <- mac & am & !cd8; qo <- mac & !am & !cd8
  tt <- if (panel == "AM-2") tum & !mac & !cd8 else tum & !mac
  x <- (d$XMin + d$XMax) / 2; y <- (d$YMin + d$YMax) / 2; q <- qa | qo
  dist <- RANN::nn2(cbind(x[tt], y[tt]), cbind(x[q], y[q]), k = 1)$nn.dists[, 1]
  isam <- qa[q]
  rbindlist(lapply(radii, function(r) data.table(
    patient_id = id, patient_code = code[[id]], panel, radius_px = r,
    n_AM = sum(isam), n_other_macrophages = sum(!isam), n_tumour = sum(tt),
    AM_pct = 100 * mean(dist[isam] <= r), other_macrophages_pct = 100 * mean(dist[!isam] <= r)
  )))
}
P <- rbindlist(lapply(c(old, new), function(id) rbind(cells(id, "AM-2"), cells(id, "AM-4"))))
P[, AM_minus_other_pp := AM_pct - other_macrophages_pct]
setorder(P, patient_code, panel, radius_px)
fwrite(P, file.path(tsa, "data/AM2_AM4_tumor_proximity_untreated_patient_data.csv"))

S <- rbindlist(lapply(radii, function(r) {
  g <- P[radius_px == r]
  one <- function(k) {
    h <- g[panel == k]
    data.table(
      radius_px = r, comparison = paste(k, "vs other macrophages (same section)"), n = nrow(h), n_higher = sum(h$AM_pct > h$other_macrophages_pct),
      median_1 = median(h$AM_pct), median_2 = median(h$other_macrophages_pct),
      P = wilcox.test(h$AM_pct, h$other_macrophages_pct, paired = TRUE, exact = TRUE)$p.value
    )
  }
  w <- dcast(g, patient_id ~ panel, value.var = "AM_minus_other_pp")
  rbind(
    one("AM-2"), one("AM-4"),
    data.table(
      radius_px = r, comparison = "AM-4 vs AM-2 (each minus other macrophages of its section)", n = nrow(w), n_higher = sum(w$`AM-4` > w$`AM-2`),
      median_1 = median(w$`AM-4`), median_2 = median(w$`AM-2`), P = wilcox.test(w$`AM-4`, w$`AM-2`, paired = TRUE, exact = TRUE)$p.value
    )
  )
}))
S[, test := "two-sided exact paired Wilcoxon signed-rank"]
fwrite(S, file.path(tsa, "data/AM2_AM4_tumor_proximity_untreated_statistics.csv"))
print(S[, .(radius_px, comparison, n, n_higher, P = signif(P, 3))])
