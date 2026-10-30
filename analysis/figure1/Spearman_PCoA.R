# Spearman PCoA
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(ggrepel)
})

options(stringsAsFactors = FALSE)
set.seed(1)

INPUT_RDS <- am_data("HTAN/am/object/Rename_am.rds")
N_HVG <- 2000L

AM_ORDER <- paste0("AM-", 1:5)
AM2_COLOR <- "#326FAC"
AM4_COLOR <- "#C44D3F"
GREY_COLOR <- "#9CA3AA"

is_technical_gene <- function(g) {
  grepl("^MT-", g) | grepl("^RPL[0-9]", g) | grepl("^RPS[0-9]", g) |
    grepl("^HB[ABDEGQZ][0-9]", g) | g %in% "MALAT1"
}

align_to_reference <- function(coords, reference) {
  a <- scale(coords, center = TRUE, scale = FALSE)
  b <- scale(reference, center = TRUE, scale = FALSE)
  sv <- svd(t(a) %*% b)
  aligned <- a %*% (sv$u %*% t(sv$v))
  aligned * sqrt(sum(b^2) / sum(aligned^2))
}

run_pcoa <- function(x, reference) {
  d <- as.dist(1 - cor(t(x), method = "spearman"))
  fit <- cmdscale(d, k = 2, eig = TRUE, add = TRUE)
  coords <- align_to_reference(fit$points, reference)
  colnames(coords) <- c("Coordinate1", "Coordinate2")
  list(coords = coords, fit = fit)
}

sc <- readRDS(INPUT_RDS)
DefaultAssay(sc) <- "RNA"
md <- sc[[]]
fine_type <- as.character(md$cell_type)
am <- sub("^(AM-[1-5]).*$", "\\1", fine_type)
stopifnot(all(am %in% AM_ORDER))

fine_table <- unique(data.frame(fine_type, AM = am))
fine_table <- fine_table[order(match(fine_table$AM, AM_ORDER), fine_table$fine_type), ]
fine_levels <- fine_table$fine_type

counts <- LayerData(sc, assay = "RNA", layer = "counts")
if (!inherits(counts, "sparseMatrix")) counts <- as(counts, "dgCMatrix")
counts <- counts[, rownames(md), drop = FALSE]

fine_factor <- factor(fine_type, levels = fine_levels)
design <- sparse.model.matrix(~ 0 + fine_factor)
colnames(design) <- fine_levels
pb_counts <- counts %*% design
fine_meta <- data.frame(
  fine_type = fine_levels,
  AM = sub("^(AM-[1-5]).*$", "\\1", fine_levels)
)

lib <- Matrix::colSums(pb_counts)
cpm <- as.matrix(t(t(pb_counts) / lib) * 1e6)
logcpm <- log2(cpm + 1)

gene_var <- apply(logcpm, 1, var)
eligible <- rowSums(cpm >= 1) >= 3 & !is_technical_gene(rownames(cpm)) &
  is.finite(gene_var) & gene_var > 0
hvg <- head(names(sort(gene_var[eligible], decreasing = TRUE)), N_HVG)

x <- scale(t(logcpm[hvg, , drop = FALSE]), center = TRUE, scale = TRUE)

pca_ref <- prcomp(x, center = FALSE, scale. = FALSE)$x[, 1:2, drop = FALSE]
pcoa <- run_pcoa(x, pca_ref)

cdf <- data.frame(AM = fine_meta$AM, as.data.frame(pcoa$coords))
centroids <- aggregate(cbind(Coordinate1, Coordinate2) ~ AM, data = cdf, FUN = mean)
rownames(centroids) <- centroids$AM
centroids <- centroids[AM_ORDER, ]

plot_data <- centroids
plot_data$AM <- factor(plot_data$AM, levels = AM_ORDER)
plot_data$focus <- ifelse(plot_data$AM == "AM-2", "AM-2",
  ifelse(plot_data$AM == "AM-4", "AM-4", "Other")
)

p <- ggplot(plot_data, aes(Coordinate1, Coordinate2)) +
  geom_hline(yintercept = 0, color = "#D8DCE0", linewidth = 0.42) +
  geom_vline(xintercept = 0, color = "#D8DCE0", linewidth = 0.42) +
  geom_point(
    data = subset(plot_data, AM %in% c("AM-2", "AM-4")),
    aes(color = focus), size = 9, alpha = 0.12, stroke = 0, show.legend = FALSE
  ) +
  geom_point(aes(fill = focus, size = focus),
    shape = 21, color = "white",
    stroke = 0.8, show.legend = FALSE
  ) +
  geom_text_repel(
    aes(label = AM),
    size = 4,
    fontface = ifelse(as.character(plot_data$AM) %in% c("AM-2", "AM-4"), "bold", "plain"),
    color = "#252A30",
    box.padding = 0.45, point.padding = 0.35,
    min.segment.length = Inf, seed = 1, max.overlaps = Inf
  ) +
  scale_color_manual(values = c("AM-2" = AM2_COLOR, "AM-4" = AM4_COLOR, "Other" = GREY_COLOR), guide = "none") +
  scale_fill_manual(values = c("AM-2" = AM2_COLOR, "AM-4" = AM4_COLOR, "Other" = GREY_COLOR), guide = "none") +
  scale_size_manual(values = c("AM-2" = 4.7, "AM-4" = 4.7, "Other" = 3.6), guide = "none") +
  labs(x = "Coordinate 1", y = "Coordinate 2") +
  theme_classic(base_size = 11) +
  theme(
    axis.title  = element_text(size = 12, color = "black"),
    axis.text   = element_text(size = 12, color = "black"),
    axis.line   = element_line(color = "black", linewidth = 0.2),
    plot.margin = margin(9, 14, 8, 9)
  )
p
