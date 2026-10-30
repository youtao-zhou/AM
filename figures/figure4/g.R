# TWAS Manhattan plot from the saved 23,455-gene table; threshold: 2e-6.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

root <- am_data("twaslca")
setwd(root)

font <- "Noto Sans"
y_cap <- 25
thr <- 2e-6
col_sig <- "#1D91C0"
col_dir <- c(up = "#C1007A", down = "#436FB8")
chr_cols <- setNames(c(brewer.pal(9, "Set3"), brewer.pal(8, "Set2"), brewer.pal(5, "Set2"))[1:22], 1:22)
data_file <- "Fig4h_TWAS_manhattan_all_genes.csv"

stopifnot(file.exists(data_file))

g <- fread(data_file)
g[, `:=`(sig = OTTERS_p < thr, ly = pmin(-log10(OTTERS_p), y_cap), capped = -log10(OTTERS_p) > y_cap)]
ax <- g[, .(center = mean(range(BPcum))), by = CHROM]

pick <- function(win, drop_noncoding) {
  x <- g[sig == TRUE & !is.na(SYMBOL) & SYMBOL != ""]
  if (drop_noncoding) x <- x[!grepl("^LOC|^LINC|-AS[0-9]*$|-DT$", SYMBOL)]
  x[order(OTTERS_p)][!duplicated(paste(CHROM, GeneStart %/% win))]
}

manhattan <- function(lab, pt, lab_size, base, x_show = 1:22, title = NULL) {
  ggplot(g, aes(BPcum, ly)) +
    geom_point(data = g[sig == FALSE], aes(colour = factor(CHROM)), alpha = 0.8, size = pt) +
    geom_hline(yintercept = -log10(thr), linetype = "dotted", colour = "red", linewidth = base / 20) +
    geom_point(data = g[sig == TRUE], aes(shape = capped), colour = col_sig, size = pt * 2) +
    geom_label_repel(
      data = lab, aes(label = SYMBOL, colour = dir), size = lab_size, family = font,
      box.padding = 0.2, label.padding = 0.08, label.size = 0.15, max.overlaps = Inf,
      segment.size = 0.2, segment.color = "grey50", min.segment.length = 0, seed = 1
    ) +
    scale_colour_manual(values = c(chr_cols, col_dir), guide = "none") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    scale_x_continuous(breaks = ax$center, labels = ifelse(ax$CHROM %in% x_show, ax$CHROM, ""), expand = expansion(mult = 0.005)) +
    scale_y_continuous(limits = c(0, y_cap + 1.5), breaks = seq(0, y_cap, 5), expand = c(0, 0)) +
    labs(x = "Chromosome", y = expression(-log[10](italic(P))), title = title) +
    theme_bw(base_family = font, base_size = base) +
    theme(
      panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90"),
      axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25),
      axis.text = element_text(colour = "black"), plot.title = element_text(hjust = 0.5),
      plot.margin = margin(2, 4, 2, 2)
    )
}

p_main <- manhattan(pick(5e6, TRUE),
  pt = 0.45, lab_size = 2, base = 8, x_show = c(1:12, 14, 16, 18, 20, 22),
  title = sprintf("Lung-cancer TWAS: %d of %d genes with P < 2×10⁻⁶", sum(g$sig), nrow(g))
)
p_full <- manhattan(pick(1e6, FALSE), pt = 1.6, lab_size = 3.2, base = 13)
if (interactive()) {
  print(p_main); print(p_full)
}

save <- function(p, file, w, h, ...) ggsave(file.path(am_output("Figure4", ""), file), p, width = w, height = h, bg = "white", ...)
save(p_main, "Fig4h_TWAS_manhattan_main.pdf", 7.0, 2.0, device = cairo_pdf)
save(p_main, "Fig4h_TWAS_manhattan_main.png", 7.0, 2.0, dpi = 200)
save(p_full, "TWAS_manhattan_full.pdf", 14, 3, device = cairo_pdf)
