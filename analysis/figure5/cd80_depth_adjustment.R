# cd80 depth adjustment
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(data.table); library(Seurat); library(coin)
})
out_dir <- am_data("vibe/157_cd80_macrophage_icb_mif"); dd <- file.path(out_dir, "data"); dir.create(dd, recursive = TRUE, showWarnings = FALSE)
src <- am_data("vibe/140_figure4_v2_am4_niche/data/gse243013")
t0 <- Sys.time()
log_phase <- function(p, note = "") cat(sprintf("%s\t%s\t%d\t%s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), p, as.integer(difftime(Sys.time(), t0, units = "secs")), note),
  file = file.path(out_dir, "run.log"), append = TRUE
)
q <- readRDS(file.path(src, "gse243013_macrophages_seurat.rds"))
genes <- c("CD80", "CD86", "CD274", "CXCL10", "CXCL9", "IL10", "CCL22", "CD68")
cnt <- GetAssayData(q, assay = "RNA", layer = "counts")[intersect(genes, rownames(q)), , drop = FALSE]
meta <- fread(file.path(src, "macrophage_metadata.csv.gz"))
stopifnot(identical(meta$cellID, colnames(q)))
rm(q); invisible(gc()); log_phase("A_load", sprintf("cells=%d", ncol(cnt)))
cell <- data.table(
  sampleID = meta$sampleID, sub = meta$sub_cell_type, umi = meta$total_counts, CD80 = as.numeric(cnt["CD80", ]),
  CD86 = as.numeric(cnt["CD86", ]), CD274 = as.numeric(cnt["CD274", ]), CXCL10 = as.numeric(cnt["CXCL10", ])
)
cell[, `:=`(cd80pos = CD80 > 0, lumi = log(umi))]
fit <- glm(cd80pos ~ lumi, family = binomial, data = cell); cell[, p_exp := fitted(fit)]
cat("CD80+ overall:", round(100 * mean(cell$cd80pos), 2), "%; depth OR per log-UMI:", round(exp(coef(fit)[2]), 2), "\n")

a3 <- cell[, .(n = .N, CD80pos_pct = 100 * mean(cd80pos), share_of_CD80pos = sum(cd80pos)), by = sub][, share_of_CD80pos := 100 * share_of_CD80pos / sum(share_of_CD80pos)][order(-CD80pos_pct)]
fwrite(a3, file.path(dd, "A3_cd80pos_by_author_subtype.csv")); print(a3)

cl <- fread(file.path(src, "sample_clinical.csv"))[pathological_response != "unknowm"]
cl[, `:=`(
  response = factor(pathological_response, c("pCR", "MPR", "non-MPR"), ordered = TRUE), regression = as.numeric(pathological_response_rate),
  histology = cancer_type, io = `anti-PD1_therapy` != "No"
)]
pt <- cell[, .(
  n_mac = .N, CD80pos_pct = 100 * mean(cd80pos), CD80_OE = sum(cd80pos) / sum(p_exp), med_umi = median(umi),
  CD86pos_pct = 100 * mean(CD86 > 0), CD274pos_pct = 100 * mean(CD274 > 0)
), by = sampleID]
ct <- fread(file.path(src, "sample_celltype_counts.csv.gz"))
tc <- ct[, .(
  n_CD4 = sum(n[grepl("^CD4T_", sub_cell_type)]), n_CD8 = sum(n[grepl("^CD8T_", sub_cell_type)]),
  n_Treg = sum(n[grepl("^CD4T_Treg_", sub_cell_type)]), n_TregCCR8 = sum(n[sub_cell_type == "CD4T_Treg_CCR8"]),
  n_Tex = sum(n[sub_cell_type %in% c("CD8T_Tex_CXCL13", "CD8T_terminal_Tex_LAYN")]),
  n_Tem = sum(n[sub_cell_type %in% c("CD8T_Tem_GZMK+GZMH+", "CD8T_Tem_GZMK+NR4A1+")]),
  n_mac_all = sum(n[grepl("^M", sub_cell_type) & major_cell_type == "Myeloid cell" & sub_cell_type != "Mast cell"]),
  n_immune = sum(n)
), by = sampleID]
tc[, `:=`(
  Treg_of_CD4 = 100 * n_Treg / n_CD4, TregCCR8_of_CD4 = 100 * n_TregCCR8 / n_CD4, Tex_of_CD8 = 100 * n_Tex / n_CD8, Tem_of_CD8 = 100 * n_Tem / n_CD8,
  CD8_of_immune = 100 * n_CD8 / n_immune, Treg_to_CD8 = n_Treg / pmax(n_CD8, 1)
)]
tcr <- fread(am_data("am_xinfuzhu/zhangzemin/GSE243013_T_with_TCR_annotation.csv.gz"), select = c("sampleID", "sub_cell_type", "expansion"))
tcr <- tcr[grepl("^CD8T_", sub_cell_type), .(CD8_expanded_pct = 100 * mean(expansion != "non-expanded")), by = sampleID]
ps <- Reduce(function(a, b) merge(a, b, by = "sampleID", all.x = TRUE), list(merge(pt, cl, by = "sampleID"), tc, tcr))
ps[, CD80pos_per_immune := CD80pos_pct * n_mac / n_immune]
fwrite(ps, file.path(dd, "A_patient_level.csv"))
main <- ps[io == TRUE & n_mac >= 100]; main[, h := factor(histology)]
cat("main n =", nrow(main), "|", paste(names(table(main$response)), table(main$response)), "|", paste(names(table(main$histology)), table(main$histology)), "\n")

a1 <- rbindlist(lapply(c("CD80pos_pct", "CD80_OE", "CD80pos_per_immune", "CD86pos_pct", "CD274pos_pct"), function(v) rbindlist(lapply(c("All (stratified)", "LUAD", "LUSC"), function(sub) {
  d <- if (sub == "All (stratified)") copy(main) else main[histology == sub]; d[, h := factor(histology)]
  tr <- pvalue(independence_test(if (sub == "All (stratified)") as.formula(paste(v, "~ response | h")) else as.formula(paste(v, "~ response")), data = d, scores = list(response = 1:3)))
  d2 <- d[response %in% c("pCR", "non-MPR")][, r2 := factor(as.character(response))]
  pw <- pvalue(wilcox_test(if (sub == "All (stratified)") as.formula(paste(v, "~ r2 | h")) else as.formula(paste(v, "~ r2")), data = d2))
  rs <- if (sub == "All (stratified)") {
    a <- resid(lm(rank(get(v)) ~ h, d)); b <- resid(lm(rank(regression) ~ h, d)); ct <- cor.test(a, b); c(ct$estimate, ct$p.value)
  } else {
    ct <- cor.test(d[[v]], d$regression, method = "spearman", exact = FALSE); c(ct$estimate, ct$p.value)
  }
  data.table(
    variable = v, subset = sub, n = nrow(d), med_pCR = median(d[response == "pCR"][[v]]), med_MPR = median(d[response == "MPR"][[v]]),
    med_nonMPR = median(d[response == "non-MPR"][[v]]), P_trend = as.numeric(tr), P_nonMPR_vs_pCR = as.numeric(pw), rho_regression = rs[1], P_regression = rs[2]
  )
}))))
fwrite(a1, file.path(dd, "A1_response_tests.csv")); print(a1, digits = 3)

tf <- c("Treg_of_CD4", "TregCCR8_of_CD4", "Treg_to_CD8", "Tex_of_CD8", "Tem_of_CD8", "CD8_of_immune", "CD8_expanded_pct")
a2 <- rbindlist(lapply(c("CD80pos_pct", "CD80_OE", "CD80pos_per_immune"), function(v) rbindlist(lapply(tf, function(t) {
  d <- main[is.finite(get(t)) & n_CD4 >= 20 & n_CD8 >= 20]
  a <- resid(lm(rank(d[[v]]) ~ histology + response, d)); b <- resid(lm(rank(d[[t]]) ~ histology + response, d)); c1 <- cor.test(a, b)
  a <- resid(lm(rank(d[[v]]) ~ histology + response + log(med_umi), d)); b <- resid(lm(rank(d[[t]]) ~ histology + response + log(med_umi), d)); c2 <- cor.test(a, b)
  data.table(variable = v, T_feature = t, n = nrow(d), partial_rho = unname(c1$estimate), P = c1$p.value, partial_rho_depth = unname(c2$estimate), P_depth = c2$p.value)
}))))
a2[, q := p.adjust(P, "BH"), by = variable]; fwrite(a2, file.path(dd, "A2_tcell_links.csv")); print(a2, digits = 3)
log_phase("A_done")
