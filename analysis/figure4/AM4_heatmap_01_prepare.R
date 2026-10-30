# AM4 heatmap 01 prepare
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)

oak_dir  <- am_data("HTAN/EGA/OAK")
data_dir <- file.path(oak_dir, "AM4CD4", "data")
n_groups <- 10

scores <- fread(file.path(data_dir, "AM4_heatmap_input_scores.csv"))
genes  <- fread(file.path(data_dir, "AM4_heatmap_genes.csv"))
expr   <- readRDS(file.path(oak_dir, "OAK.rds"))
stopifnot(nrow(scores) == 319, !anyDuplicated(scores$sample),
          all(scores$sample %in% colnames(expr)), all(genes$gene %in% rownames(expr)))

x <- as.matrix(expr[genes$gene, scores$sample, drop = FALSE])
z <- (x - genes$mean) / genes$SD
rm(expr, x); gc()

setorder(scores, AM4)
scores[, group := ceiling(n_groups * frank(AM4, ties.method = "first") / .N)]

tumour_z <- data.table(sample = rep(colnames(z), each = nrow(z)), gene = rep(rownames(z), ncol(z)), z = as.vector(z))
tumour_z <- merge(tumour_z, scores, by = "sample")
tumour_z <- merge(tumour_z, genes[, .(gene, axis, gene_order)], by = "gene")
setorder(tumour_z, AM4, gene_order)
fwrite(tumour_z[, .(sample, AM4, group, axis, gene, z)], file.path(data_dir, "AM4_heatmap_tumour_gene_z.csv"))

plot_data <- tumour_z[, .(mean_z = mean(z), n_tumours = .N), by = .(group, axis, gene, gene_order)]
setorder(plot_data, gene_order, group)
fwrite(plot_data, file.path(data_dir, "AM4_heatmap_plot_data.csv"))

print(scores[, .(n = .N, mean_signature = round(mean(AM4), 2)), by = group])
message("Wrote: AM4_heatmap_tumour_gene_z.csv, AM4_heatmap_plot_data.csv in ", data_dir)
