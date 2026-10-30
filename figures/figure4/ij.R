# CD80-positive macrophages and T-cell markers in 52 LUAD patients.
# Positivity is based on RNA expression, not CD68 protein gating.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(ggplot2)
library(patchwork)
library(coin)

root <- am_data("HTAN/am/GSE243013")
setwd(root)

font <- "Noto Sans"
ink <- "#20292F"
cols <- c(`non-MPR` = "#5BB5AC", MPR = "#DAC6CF", pCR = "#DE526C")
resp_lab <- c(`non-MPR` = "no MPR/pCR", MPR = "MPR", pCR = "pCR")
n_bins <- 5
markers <- c("CTLA4", "CD28", "CCR8", "TNFRSF4", "TNFRSF18", "TNFRSF9", "ICOS", "IL10", "LAYN", "PDCD1", "HAVCR2", "LAG3", "TIGIT",
             "TOX", "ENTPD1", "CXCL13", "GZMB", "PRF1", "IFNG", "GZMK", "IL7R", "TCF7")
shown <- list(Treg = c("CTLA4", "CCR8", "TNFRSF18", "TNFRSF4", "TNFRSF9", "ICOS", "LAYN", "TCF7"),
              CD8  = c("CTLA4", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "ENTPD1", "LAYN", "GZMB"))
alias <- c(CTLA4 = "CTLA4", CCR8 = "CCR8", TNFRSF18 = "GITR", TNFRSF4 = "OX40", TNFRSF9 = "4-1BB", ICOS = "ICOS", LAYN = "LAYN", TCF7 = "TCF7",
           PDCD1 = "PD-1", HAVCR2 = "TIM-3", LAG3 = "LAG3", TIGIT = "TIGIT", ENTPD1 = "CD39", GZMB = "GZMB")
f_mac <- "data/Fig4g_patient_macrophage.csv"
f_tm  <- "data/Fig4g_patient_T_markers.csv"
f_rat <- "data/Fig4g_patient_CTLA4_CD28_ratio.csv"

stopifnot(all(file.exists(c(f_mac, f_tm, f_rat))))

pm <- fread(f_mac); tm <- fread(f_tm); rt <- fread(f_rat)
luad <- pm[io == TRUE & pathological_response != "unknowm" & n_mac >= 100 & histology == "LUAD"]
stopifnot(nrow(luad) == 52)

part_rho <- function(d, w) {
  ct <- cor.test(resid(lm(rank(CD80pos_pct) ~ regression + log(med_umi), d)), resid(lm(rank(d[[w]]) ~ regression + log(med_umi), d)))
  c(unname(ct$estimate), ct$p.value)
}

links <- rbindlist(lapply(split(tm, by = c("cls", "gene")), function(x) {
  d <- merge(luad[, .(sampleID, CD80pos_pct, regression, med_umi)], x[n >= 20], by = "sampleID")
  if (nrow(d) < 30) return(NULL)
  r <- part_rho(d, "OE")
  data.table(cls = x$cls[1], gene = x$gene[1], n = nrow(d), rho = r[1], P = r[2])
}))
links[, q := p.adjust(P, "BH")]
fwrite(links[order(P)], am_output("Figure4", "Fig4g_stats_marker_correlations_LUAD.csv"))

pz <- tm[cls %in% names(shown) & n >= 20][, z := scale(log2(OE + 0.05))[, 1], by = .(cls, gene)]
luad[, bin := cut(rank(CD80pos_pct, ties.method = "first"), n_bins, labels = FALSE)]
hm <- merge(pz, luad[, .(sampleID, bin)], by = "sampleID")[, .(z = mean(z)), by = .(cls, gene, bin)]
sel <- rbindlist(lapply(names(shown), function(k) data.table(cls = k, gene = shown[[k]], ord = seq_along(shown[[k]]))))
hm <- merge(hm, merge(sel, links, by = c("cls", "gene")), by = c("cls", "gene"))
hm[, `:=`(lab = alias[gene], cls = factor(cls, names(shown)))]
setorder(hm, cls, ord, bin)
hm[, row := factor(paste(cls, lab), rev(unique(paste(cls, lab))))]
top <- luad[, .(nonMPR = 100 * mean(pathological_response == "non-MPR"), median_CD80pos_pct = median(CD80pos_pct), n = .N), by = bin][order(bin)]
fwrite(hm[, .(cls, gene, label = lab, bin, z, partial_rho = rho, P, q)][order(cls, gene, bin)], am_output("Figure4", "Fig4g_plot_heatmap.csv"))
fwrite(top, am_output("Figure4", "Fig4g_plot_heatmap_top.csv"))

rl <- merge(luad[, .(sampleID, CD80pos_pct, regression, med_umi, pathological_response)], rt[n_T >= 200], by = "sampleID")
rl[, `:=`(ro = factor(pathological_response, c("pCR", "MPR", "non-MPR"), ordered = TRUE),
          response = factor(pathological_response, names(cols)), Treg_pct = 100 * n_Treg / n_T)]
r <- part_rho(rl, "ratio")
wp <- function(g) wilcox.test(rl[response == "non-MPR", ratio], rl[response == g, ratio])$p.value
ratio_test <- data.table(n = nrow(rl), rho = r[1], P = r[2],
                         P_trend = as.numeric(pvalue(independence_test(ratio ~ ro, data = rl, scores = list(ro = 1:3)))),
                         P_nonMPR_vs_MPR = wp("MPR"), P_nonMPR_vs_pCR = wp("pCR"),
                         rho_ratio_vs_Treg_pct = cor(rl$ratio, rl$Treg_pct, method = "spearman"))
fwrite(ratio_test, am_output("Figure4", "Fig4g_stats_ratio_LUAD.csv"))
rg <- rl[, .(n = .N, mean = mean(ratio), se = sd(ratio) / sqrt(.N), median = median(ratio)), by = response][order(response)]
fwrite(rg, am_output("Figure4", "Fig4g_plot_ratio_groups.csv"))
fwrite(rl[, .(sampleID, pathological_response, CD80pos_pct, ctla4_treg, cd28_t, n_Treg, n_T, Treg_pct, ratio)], am_output("Figure4", "Fig4g_plot_ratio.csv"))

th <- theme_classic(base_size = 8, base_family = font) +
  theme(axis.text = element_text(colour = "black", size = 7.5), axis.line = element_line(linewidth = 0.15),
        axis.ticks = element_line(linewidth = 0.15))

st <- unique(hm[, .(row, P)])[, txt := as.character(cut(P, c(0, 0.001, 0.01, 0.05, 1), c("***", "**", "*", "")))]
st <- merge(st, unique(hm[, .(row, cls)]), by = "row")
lim <- max(abs(hm$z))
ph <- ggplot(hm, aes(factor(bin), row)) +
  geom_tile(aes(fill = z), colour = "white", linewidth = 0.4) +
  geom_text(data = st, aes(x = n_bins + 0.6, y = row, label = txt), hjust = 0, vjust = 0.75, size = 2.6, family = font, inherit.aes = FALSE) +
  facet_grid(cls ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(labels = function(x) sub("^\\S+ ", "", x)) +
  scale_x_discrete(labels = c("Low", rep("", n_bins - 2), "High")) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C8102E", limits = c(-lim, lim), name = "Z-score") +
  coord_cartesian(clip = "off") + labs(x = "CD68+ CD80+ Macrophage", y = NULL) + th +
  theme(strip.background = element_blank(), strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold",size=8),
        axis.line = element_blank(), axis.ticks = element_blank(), plot.margin = margin(2, 16, 2, 2),
        legend.position = "bottom", legend.key.width = unit(16, "pt"), legend.key.height = unit(5, "pt"),
        legend.title = element_text(size = 7, vjust = 0.9), legend.text = element_text(size = 7))
pt <- ggplot(top, aes(factor(bin), 1)) +
  geom_tile(aes(fill = nonMPR), colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.0f%%", nonMPR)), size = 2.1, family = font) +
  annotate("text", x = 0.42, y = 1, label = resp_lab[["non-MPR"]], hjust = 1, size = 2.6, family = font,face="bold") +
  scale_fill_gradient(low = "#BFE3DF", high = "#2F8C83", limits = c(0, 100), guide = "none") +
  coord_cartesian(clip = "off") + theme_void(base_family = font) + theme(plot.margin = margin(2, 16, 0, 2))
p_heat <- pt / ph + plot_layout(heights = c(1, 17))
if (interactive()) print(p_heat)
ggsave(am_output("Figure4", "CD80_Treg_heatmap.pdf"), p_heat, width = 3.32, height = 3.08, device = cairo_pdf, bg = "white")

fs <- 9
ytop <- max(rl$ratio)
br <- data.table(x2 = c(2, 3), P = c(ratio_test$P_nonMPR_vs_MPR, ratio_test$P_nonMPR_vs_pCR), y = ytop * c(1.07, 1.21))
br[, lab := fifelse(P < 0.05, sprintf("P = %.3f", P), "ns")]
brk <- br[, .(x = c(1, 1, x2, x2), y = c(y - 0.02 * ytop, y, y, y - 0.02 * ytop)), by = .(id = x2)]
p_ratio <- ggplot(rg, aes(response, mean)) +
  geom_col(aes(fill = response), width = 0.62, colour = ink, linewidth = 0.25, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.22, linewidth = 0.35) +
  geom_point(data = rl, aes(response, ratio), position = position_jitter(width = 0.13, seed = 158), size = 0.9, colour = ink, alpha = 0.6) +
  geom_path(data = brk, aes(x, y, group = id), linewidth = 0.3, inherit.aes = FALSE) +
  geom_text(data = br, aes((1 + x2) / 2, y, label = lab), vjust = -0.35, size = fs * 0.9 / .pt, family = font, inherit.aes = FALSE) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_discrete(labels = sprintf("%s\n(n = %d)", resp_lab[as.character(rg$response)], rg$n)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "CTLA4+ Treg / CD28+ non-Treg T") + th +
  theme(axis.text = element_text(size = fs), axis.title.y = element_text(size = fs + 1),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(2, 8, 2, 2))
if (interactive()) print(p_ratio)
ggsave(am_output("Figure4", "CTLA_CD28_Treg_ratio.pdf"), p_ratio, width = 3.36, height = 3.08, device = cairo_pdf, bg = "white")
