# Fig5b 01 score cohorts
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(data.table); library(Matrix)
})
args <- commandArgs(trailingOnly = FALSE); here <- dirname(normalizePath(sub('^--file=', '', args[grep('^--file=', args)])))
D <- normalizePath(file.path(here, '..', 'data'))

HZ <- am_data("vibe/141_figure5_s6_clinical_prediction/data/harmonized")
TCGA <- am_data("HTAN/Scissor/TCGA_LUAD_exp2.RData")
SP <- am_data("vibe/15_tumor_states_spatial_screen/spatial")

rank_enrichment <- function(expr, sets) {
  N <- nrow(expr); r <- seq.int(N, 1); w <- r^.25
  masks <- lapply(sets, function(g) rownames(expr) %in% g)
  ans <- matrix(NA_real_, ncol(expr), length(sets), dimnames = list(colnames(expr), names(sets)))
  for (j in seq_len(ncol(expr))) {
    ord <- order(expr[, j], decreasing = TRUE)
    for (k in seq_along(masks)) {
      I <- masks[[k]][ord]; ans[j, k] <- sum(w[I] * r[I]) / sum(w[I]) - mean(r[!I])
    }
  }
  ans
}
zs <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

gs <- readRDS(file.path(HZ, 'gene_sets.rds'))
blocks <- list(
  'AM-2 activation' = c('CCL3', 'CCL4', 'IL1B', 'CCL20'), 'Cytotoxic T cell' = c('CD8A', 'GZMB', 'PRF1', 'IFNG'),
  'IFN-γ response' = c('CXCL9', 'CXCL10', 'STAT1'), 'Checkpoints' = c('CD274', 'PDCD1', 'LAG3', 'CTLA4')
)
G15 <- unlist(blocks, use.names = FALSE)
fwrite(
  rbind(
    data.table(set = 'AM2_signature', block = NA_character_, gene = gs$AM2), data.table(set = 'AM4_signature', block = NA_character_, gene = gs$AM4),
    rbindlist(lapply(names(blocks), function(b) data.table(set = 'heatmap_gene', block = b, gene = blocks[[b]])))
  ),
  file.path(D, 'Fig5b_gene_sets.csv')
)

cohorts <- list()
for (k in c('OAK', 'POPLAR', 'SU2C_MARK', 'France3', 'France4', 'GSE274975', 'GSE283829', 'GSE207422', 'GSE161537')) {
  o <- readRDS(file.path(HZ, paste0(k, '.rds'))); keep <- o$clin$histology == 'LUAD' & !is.na(o$clin$histology)
  cohorts[[k]] <- list(expr = o$expr[, o$clin$sample[keep], drop = FALSE], label = k, context = unique(o$clin$setting)[1])
}
e <- new.env(); load(TCGA, envir = e)
cohorts[['TCGA-LUAD']] <- list(expr = log2(e$bulk_dataset + 1), label = 'TCGA-LUAD', context = 'untreated_resection'); rm(e)
sp <- list.files(SP); pat <- sub('_LUAD[0-9]*$', '', sp)
pb <- lapply(split(sp, pat), function(ss) Reduce(`+`, lapply(ss, function(s) Matrix::rowSums(readRDS(file.path(SP, s, 'counts.rds'))))))
ug <- Reduce(intersect, lapply(pb, names)); PB <- sapply(pb, function(v) v[ug]); PB <- log2(t(t(PB) / colSums(PB)) * 1e6 + 1)
cohorts[['GSE307534']] <- list(expr = PB, label = 'GSE307534 (Visium pseudobulk)', context = 'untreated_resection_spatial')

inc <- list(); out <- list()
for (k in names(cohorts)) {
  E <- cohorts[[k]]$expr; E <- E[rowSums(is.finite(E)) == ncol(E), , drop = FALSE]
  c2 <- mean(gs$AM2 %in% rownames(E)); c4 <- mean(gs$AM4 %in% rownames(E)); n <- ncol(E)
  ok <- n >= 12 && c2 >= .6 && c4 >= .6
  inc[[k]] <- data.table(
    cohort = k, label = cohorts[[k]]$label, context = cohorts[[k]]$context, n_LUAD = n, AM2_cov = c2, AM4_cov = c4, included = ok,
    reason = if (ok) '' else paste(c(if (n < 12) 'n<12', if (c2 < .6) 'AM2 coverage<60%', if (c4 < .6) 'AM4 coverage<60%'), collapse = '; ')
  )
  if (!ok) next
  r <- rank_enrichment(E, list(AM2 = gs$AM2, AM4 = gs$AM4))
  x <- data.table(cohort = k, sample = colnames(E), AM2_enrichment = r[, 'AM2'], AM4_enrichment = r[, 'AM4'], AM_score = zs(r[, 'AM2']) - zs(r[, 'AM4']))
  for (g in G15) set(x, j = g, value = if (g %in% rownames(E)) unname(E[g, ]) else NA_real_)
  out[[k]] <- x
}
inc <- rbindlist(inc); print(inc)
fwrite(inc, file.path(D, 'Fig5b_cohort_inclusion.csv'))
fwrite(rbindlist(out), file.path(D, 'Fig5b_sample_scores_genes.csv'))
cat('done\n')
