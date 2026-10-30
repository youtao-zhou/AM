# 03 GRSA by state paired
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(Seurat);library(Matrix);library(data.table);library(ReporterScore)
})
mode <- commandArgs(TRUE)[1]; stopifnot(mode %in% c('AM2', 'OtherAM', 'AM4'))
dd <- am_data("HTAN/am/GRSA/data")
getX <- function(o) tryCatch(LayerData(o, assay = 'RNA', layer = 'data'), error = function(e) GetAssayData(o, assay = 'RNA', slot = 'data'))
key <- function(s) {
  p <- rep(NA_character_, length(s)); tis <- rep(NA_character_, length(s))
  for (i in seq_along(s)) {
    m <- regmatches(s[i], regexec('^(.*)_(tumor_primary|normal_adjacent)$', s[i]))[[1]]
    if (length(m) == 3) {
      p[i] <- m[2]; tis[i] <- ifelse(m[3] == 'normal_adjacent', 'normal', 'tumour'); next
    }
    k <- regmatches(s[i], regexec('^Kim_Lee_2020_LUNG_([NT])([0-9]+)$', s[i]))[[1]]
    if (length(k) == 3) {
      p[i] <- paste0('KimLee_', k[3]); tis[i] <- ifelse(k[2] == 'N', 'normal', 'tumour'); next
    }
    h <- regmatches(s[i], regexec('^HTA3_([0-9]+)_([0-9]+)$', s[i]))[[1]]
    if (length(h) == 3) {
      p[i] <- paste0('HTABU_', h[2]); tis[i] <- ifelse(h[3] == '1001', 'normal', 'tumour')
    }
  }
  data.table(patient = p, tissue = tis)
}
pb <- function(X, cells_by_sample) sapply(cells_by_sample, function(cc) Matrix::rowMeans(X[, cc, drop = FALSE]))

lu <- readRDS(am_data("HTAN/am/object/Rename_am.rds"))
ml <- data.table(cell = colnames(lu), sample = as.character(lu$sample), state = as.character(lu$cell_type))
ml[, cls := sub('^(AM-[0-9]).*$', '\\1', state)]; ml <- cbind(ml, key(ml$sample))
Xl <- getX(lu); rm(lu); invisible(gc())

if (mode %in% c('AM2', 'OtherAM')) {
  tsel <- if (mode == 'AM2') ml[tissue == 'tumour' & cls == 'AM-2'] else ml[tissue == 'tumour' & !cls %in% c('AM-2', 'AM-4')]
  tsel <- tsel[!grepl('^HTABU_', patient)]
  nsel <- ml[tissue == 'normal' & !grepl('^HTABU_', patient)]
  tn <- tsel[, .N, by = .(patient, sample)][N >= 20]; nn <- nsel[, .N, by = .(patient, sample)][N >= 20]
  a2 <- ml[tissue == 'tumour' & cls == 'AM-2' & !grepl('^HTABU_', patient), .N, by = .(patient, sample)][N >= 20]
  pts <- intersect(intersect(unique(a2$patient), unique(nn$patient)), unique(tn$patient))
  T <- pb(Xl, split(tsel[patient %in% pts & sample %in% tn$sample, cell], tsel[patient %in% pts & sample %in% tn$sample, sample]))
  N <- pb(Xl, split(nsel[patient %in% pts & sample %in% nn$sample, cell], nsel[patient %in% pts & sample %in% nn$sample, sample]))
  smp2pt <- unique(ml[, .(sample, patient)]); setkey(smp2pt, sample)
  collapse <- function(M) {
    p <- smp2pt[colnames(M), patient]; sapply(split(seq_len(ncol(M)), p), function(ix) rowMeans(M[, ix, drop = FALSE]))
  }
  T <- collapse(T); N <- collapse(N); genes <- rownames(Xl)
} else {
  no <- readRDS(am_data("HTAN/am_normal/Rename_am_normal.rds"))
  k <- intersect(c('sample', 'Sample', 'orig.ident'), colnames(no@meta.data))[1]
  mn <- data.table(cell = colnames(no), sample = as.character(no@meta.data[[k]])); mn <- cbind(mn, key(mn$sample))
  Xn <- getX(no); rm(no); invisible(gc())
  tsel <- ml[tissue == 'tumour' & cls == 'AM-4' & grepl('^HTABU_', patient)]
  nsel <- mn[tissue == 'normal' & grepl('^HTABU_', patient)]
  tn <- tsel[, .N, by = .(patient, sample)][N >= 20]; nn <- nsel[, .N, by = .(patient, sample)][N >= 20]
  pts <- intersect(unique(tn$patient), unique(nn$patient))
  genes <- intersect(rownames(Xl), rownames(Xn))
  T <- pb(Xl[genes, ], split(tsel[patient %in% pts & sample %in% tn$sample, cell], tsel[patient %in% pts & sample %in% tn$sample, sample]))
  N <- pb(Xn[genes, ], split(nsel[patient %in% pts & sample %in% nn$sample, cell], nsel[patient %in% pts & sample %in% nn$sample, sample]))
  collapse <- function(M) {
    p <- sub('^HTA3_([0-9]+)_.*$', 'HTABU_\\1', colnames(M)); sapply(split(seq_len(ncol(M)), p), function(ix) rowMeans(M[, ix, drop = FALSE]))
  }
  T <- collapse(T); N <- collapse(N)
}
T <- T[genes, , drop = FALSE]; N <- N[genes, , drop = FALSE]
colnames(T) <- paste0('Tumour|', colnames(T)); colnames(N) <- paste0('Normal|', colnames(N))
saveRDS(list(T = T, N = N), file.path(dd, sprintf('Fig1j_pseudobulk_%s.rds', mode)))
M <- cbind(T, N); M <- M[rowSums(M > 0) >= .5 * ncol(M), ]
fwrite(data.table(pseudobulk_id = colnames(M)), file.path(dd, sprintf('Fig1j_pseudobulk_ids_%s.csv', mode)))
cat(sprintf('%s: patients = %d; genes = %d\n', mode, length(pts), nrow(M)))
meta <- data.frame(Group = ifelse(grepl('^Tumour', colnames(M)), 'am_tumor', 'am_normal'), row.names = colnames(M))
set.seed(1)
res <- GRSA(as.data.frame(M), 'Group', meta, mode = 'directed', method = 'wilcox.test', perm = 999, type = 'hsa', feature = 'gene', threads = 1, verbose = FALSE)
rs <- as.data.table(res$reporter_s); rs[, Description := sub(' - Homo sapiens \\(human\\)', '', Description)]
setorder(rs, -ReporterScore)
fwrite(rs, file.path(dd, sprintf('Fig1j_GRSA_%s.csv', mode)))
