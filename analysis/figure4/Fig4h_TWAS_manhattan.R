# Fig4h TWAS manhattan
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

root <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) normalizePath(dirname(sub("^--file=", "", a[1]))) else am_data("twaslca")
})
setwd(root)
dir.create("Figure", showWarnings = FALSE)

font <- "Noto Sans"
y_cap <- 25
thr <- 2e-6
col_sig <- "#1D91C0"
col_dir <- c(up = "#C1007A", down = "#436FB8")
chr_cols <- setNames(c(brewer.pal(9, "Set3"), brewer.pal(8, "Set2"), brewer.pal(5, "Set2"))[1:22], 1:22)
data_file <- "Fig4h_TWAS_manhattan_all_genes.csv"

if (!file.exists(data_file)) {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  acat <- function(p) {
    p <- p[!is.na(p)]; if (!length(p)) return(NA_real_)
    p <- pmin(pmax(p, 1e-300), 1 - 1e-15)
    t <- mean(ifelse(p < 1e-16, 1 / (p * pi), tan((0.5 - p) * pi)))
    if (t > 1e15) 1 / (t * pi) else 0.5 - atan(t) / pi
  }
  meth <- c(lassosum = "lassosum", P0.001 = "PT_0.001", P0.05 = "PT_0.05", PRScs = "PRScs", SDPR = "SDPR")
  f <- CJ(d = list.dirs("stage2", recursive = FALSE), m = names(meth))[, path := file.path(d, paste0(m, ".txt"))][file.exists(path)]
  res <- rbindlist(Map(function(path, m) {
    fread(path, colClasses = list(character = "TargetID"))[, .(CHROM, GeneStart,
      Gene = TargetID, z = as.numeric(FUSION_Z),
      p = as.numeric(FUSION_PVAL), method = meth[[m]]
    )]
  }, f$path, f$m))[!is.na(p)]
  g <- res[, .(CHROM = as.integer(CHROM[1]), GeneStart = GeneStart[1], OTTERS_p = acat(p), z_minp = z[which.min(p)]), by = Gene]

  ref <- fread("TWAS_Result.csv", colClasses = list(character = "Gene"))
  chk <- merge(g, ref[, .(Gene, ref_p = OTTERS_p)], by = "Gene")
  cat(
    "genes tested:", nrow(g), "| matched to TWAS_Result.csv:", nrow(chk),
    "| max |log10 P diff|:", signif(max(abs(log10(chk$OTTERS_p / chk$ref_p))), 3), "\n"
  )

  g[, SYMBOL := unname(suppressMessages(mapIds(org.Hs.eg.db, Gene, "SYMBOL", "ENSEMBL", multiVals = "first")))]
  setorder(g, CHROM, GeneStart)
  g[, BPcum := GeneStart + c(0, cumsum(as.numeric(tapply(GeneStart, CHROM, max))))[CHROM]]
  g[, dir := fifelse(z_minp > 0, "up", "down")]
  fwrite(g, data_file)
}

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

save <- function(p, file, w, h, ...) ggsave(file.path("Figure", file), p, width = w, height = h, bg = "white", ...)
save(p_main, "Fig4h_TWAS_manhattan_main.pdf", 7.0, 2.0, device = cairo_pdf)
save(p_main, "Fig4h_TWAS_manhattan_main.png", 7.0, 2.0, dpi = 200)
save(p_full, "TWAS_manhattan_full.pdf", 14, 3, device = cairo_pdf)
