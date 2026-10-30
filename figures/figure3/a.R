# SCENIC regulon activity, averaged within sample and AM state.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(ggplot2)

setwd(am_data("HTAN/am/SCENIC"))

tfs <- c(
  "ELF4(+)", "EP300(+)", "RXRA(+)", "SMARCA1(+)",
  "ZBTB7A(+)", "FOXO1(+)", "FOXO3(+)", "NFIL3(+)"
)
ams <- c("AM-1", "AM-2", "AM-3", "AM-4")
tissues <- c("Normal", "Tumor")

d <- fread("SCENIC_S4A_Tumor_Normal_original_sample_mean_auc.csv")
d <- d[n_cells >= 20 & AM %in% ams & tissue %in% tissues]

long <- melt(d,
  id.vars = c("tissue", "AM"), measure.vars = tfs,
  variable.name = "TF", value.name = "auc", variable.factor = FALSE
)
m <- long[, .(auc = if (all(is.na(auc))) NA_real_ else mean(auc, na.rm = TRUE)),
  by = .(tissue, TF, AM)
]
m <- m[CJ(tissue = tissues, TF = tfs, AM = ams), on = .(tissue, TF, AM)]

m[, z := {
  mu <- mean(auc, na.rm = TRUE)
  s <- sqrt(mean((auc - mu)^2, na.rm = TRUE))
  if (is.finite(s) && s > 0) (auc - mu) / s else NA_real_
}, by = .(tissue, TF)]

m[, TF := factor(TF, rev(tfs))]
m[, AM := factor(AM, ams)]

rdbu <- c(
  "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"
)
j <- match("AM-2", ams)

p <- ggplot(m, aes(AM, TF, fill = z)) +
  geom_tile(color = "white", linewidth = 0.25) +
  annotate("rect",
    xmin = j - 0.5, xmax = j + 0.5, ymin = 0.5, ymax = length(tfs) + 0.5,
    fill = NA, color = "#258A91", linewidth = 0.65
  ) +
  facet_wrap(~tissue, nrow = 1) +
  scale_fill_gradientn(
    colors = rdbu, limits = c(-2, 2), breaks = -2:2,
    oob = scales::squish, na.value = "#E4E4E4", name = "Zscore"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  guides(fill = guide_colorbar(
    title.position = "top", title.hjust = 0.5,
    barwidth = unit(120, "pt"), barheight = unit(7, "pt")
  )) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(color = "black", size = 12),
    strip.text = element_text(size = 12, color = "black", face = "bold"),
    panel.spacing = unit(12, "pt"),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

if (interactive()) print(p)
ggsave(am_output("Figure3", "Figure3a.pdf"), p, width = 6.5, height = 4, device = cairo_pdf)
fwrite(m, am_output("Figure3", "Figure3a_values.csv"))
