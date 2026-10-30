# AM score/gene Spearman correlations; DerSimonian-Laird pooling on Fisher z.
# Sampling variance is 1/(n - 3).
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

root <- am_data("HTAN/am/TWAS_signature")
setwd(root)

font <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "Liberation Sans"
fs <- 9
fs_l <- 10
short <- c(
  `TCGA-LUAD` = "TCGA", OAK = "OAK", POPLAR = "POPLAR", SU2C_MARK = "SU2C-MARK", GSE274975 = "GSE274975",
  GSE283829 = "GSE283829", France3 = "France3", GSE307534 = "GSE307534 (ST)"
)
blab <- c(
  `AM-2 activation` = "AM-2\nactivation", `Cytotoxic T cell` = "Cytotoxic\nT cell",
  `IFN-γ response` = "IFN-γ\nresponse", Checkpoints = "Checkpoints"
)
pal <- c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B")

x <- fread("data/Fig5b_sample_scores_genes.csv")
inc <- fread("data/Fig5b_cohort_inclusion.csv")[included == TRUE][order(-n_LUAD)]
gsets <- fread("data/Fig5b_gene_sets.csv")[set == "heatmap_gene", .(gene, block)]
G15 <- gsets$gene

long <- melt(x, id.vars = c("cohort", "AM_score"), measure.vars = G15, variable.name = "gene", variable.factor = FALSE)
res <- long[, if (isTRUE(sd(value, na.rm = TRUE) > 0)) {
  ct <- suppressWarnings(cor.test(AM_score, value, method = "spearman", exact = FALSE))
  .(n = .N, rho = unname(ct$estimate), P = ct$p.value)
}, by = .(cohort, gene)]

dl <- function(z, v) {
  w <- 1 / v; fe <- sum(w * z) / sum(w); Q <- sum(w * (z - fe)^2); df <- length(z) - 1
  tau2 <- max(0, (Q - df) / (sum(w) - sum(w^2) / sum(w)))
  ws <- 1 / (v + tau2)
  list(est = sum(ws * z) / sum(ws), se = sqrt(1 / sum(ws)), I2 = if (Q > 0) max(0, (Q - df) / Q) else 0)
}
pool <- res[,
  {
    m <- dl(atanh(rho), 1 / (n - 3))
    .(
      rho_pooled = tanh(m$est), lo = tanh(m$est - 1.96 * m$se), hi = tanh(m$est + 1.96 * m$se),
      P = 2 * pnorm(-abs(m$est / m$se)), I2 = m$I2, k = .N, n_pos = sum(rho > 0), n_sig_pos = sum(rho > 0 & P < 0.05)
    )
  },
  by = gene
]
pool <- merge(pool, gsets, by = "gene")[order(match(block, names(blab)), -rho_pooled)]
fwrite(res, am_output("Figure5", "Fig5b_per_cohort_gene_rho.csv"))
fwrite(pool, am_output("Figure5", "Fig5b_pooled_rho.csv"))
print(pool, digits = 3)
res[, cat(sprintf(
  "cohort-gene pairs: %d; rho > 0: %d; rho > 0 & P < .05: %d; rho < 0 & P < .05: %d\n",
  .N, sum(rho > 0), sum(rho > 0 & P < 0.05), sum(rho < 0 & P < 0.05)
))]

colab <- c(setNames(sprintf("%s (%d)", short[inc$cohort], inc$n_LUAD), inc$cohort), Pooled = "Meta")
hm <- rbind(res[, .(cohort, gene, rho, P)], pool[, .(cohort = "Pooled", gene, rho = rho_pooled, P)])
hm <- merge(hm, gsets, by = "gene")
hm[, `:=`(
  col = factor(colab[cohort], colab), gene = factor(gene, rev(G15)),
  block = factor(blab[block], blab), pooled = factor(cohort == "Pooled", c(FALSE, TRUE))
)]

xlab_bottom <- labeller(pooled = c(`FALSE` = "AM score", `TRUE` = ""))

p <- ggplot(hm, aes(col, gene)) +
  geom_tile(aes(fill = rho), colour = "white", linewidth = 0.4) +
  geom_point(data = hm[P < 0.05], aes(shape = "P < 0.05"), size = 0.9, colour = "black") +
  facet_grid(block ~ pooled, scales = "free", space = "free", switch = "both", labeller = xlab_bottom) +
  scale_fill_gradientn(
    colours = pal, limits = c(-0.8, 0.8), breaks = c(-0.8, -0.4, 0, 0.4, 0.8),
    labels = c("−0.8", "−0.4", "0", "0.4", "0.8"), oob = scales::squish,
    name = "Spearman ρ",
    guide = guide_colourbar(order = 1, ticks.colour = "black", title.vjust = 1)
  ) +
  scale_shape_manual(
    values = 16, name = NULL,
    guide = guide_legend(order = 2, override.aes = list(size = 2))
  ) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = fs, base_family = font) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(size = fs, colour = "black"),
    axis.text.y = element_text(margin = margin(r = 2)),
    axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, margin = margin(b = 2)),
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x.bottom = element_text(size = fs_l, colour = "black", margin = margin(t = 5)),
    strip.text.y.left = element_text(
      angle = 0, hjust = 1, vjust = 0.5, size = fs_l,
      lineheight = 0.95, colour = "black", margin = margin(r = 5)
    ),
    strip.switch.pad.grid = unit(0, "pt"),
    panel.spacing.y = unit(4, "pt"),
    panel.spacing.x = unit(6, "pt"),
    legend.position = "right",
    legend.justification = c(0, 0.5),
    legend.title = element_text(size = fs_l, margin = margin(b = 4)),
    legend.text = element_text(size = fs),
    legend.key.width = unit(9, "pt"),
    legend.key.height = unit(18, "pt"),
    legend.spacing.y = unit(8, "pt"),
    legend.margin = margin(0, 0, 0, 4),
    plot.margin = margin(4, 4, 4, 4)
  )
if (interactive()) print(p)

ggsave(am_output("Figure5", "AMScore_multicohort_heatmap.pdf"), p, width = 4.82, height = 5.01, device = cairo_pdf, bg = "white")
