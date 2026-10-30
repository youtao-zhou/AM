# cNMF usage: cell-count-weighted means, then population z-scores within each program.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
source(am_code("R/figure1_data.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

data_dir <- file.path(INPUT_ROOT, "g")

MACRO <- c("AM-1", "AM-2", "AM-3", "AM-4", "AM-5")
COL_AM2 <- "#2C6A9F"
COL_AM4 <- "#B94E43"

usage <- fread(file.path(data_dir, "program_usage_by_state.csv"))
usage[, macro := sub("\\..*$", "", cell_type)]
stopifnot(all(usage$macro %in% MACRO))

agg <- usage[, .(mean_usage = sum(mean_usage * cells) / sum(cells)),
  by = .(program, macro)
]

agg[, z := (mean_usage - mean(mean_usage)) /
  (sqrt(mean((mean_usage - mean(mean_usage))^2)) + 1e-12),
by = program
]

own_dt <- agg[
  agg[, .I[which.max(mean_usage)], by = program]$V1,
  .(program, owner = macro)
]
agg <- merge(agg, own_dt, by = "program")

ord <- copy(own_dt)
ord[, own_z := agg[.SD, on = .(program, macro = owner), x.z]]
ord <- ord[order(match(owner, MACRO), -own_z)]

agg[, program_lab := factor(paste0("P", program),
  levels = rev(paste0("P", ord$program))
)]
agg[, macro := factor(macro, levels = MACRO)]

hl <- agg[owner %in% c("AM-2", "AM-4") & macro == owner]
hl[, col := ifelse(owner == "AM-2", COL_AM2, COL_AM4)]

lab_cols <- setNames(
  rep("grey20", nlevels(agg$program_lab)),
  levels(agg$program_lab)
)
for (i in seq_len(nrow(ord))) {
  key <- paste0("P", ord$program[i])
  if (ord$owner[i] == "AM-2") lab_cols[key] <- COL_AM2
  if (ord$owner[i] == "AM-4") lab_cols[key] <- COL_AM4
}
lab_face <- ifelse(lab_cols == "grey20", "plain", "bold")

x_cols <- ifelse(MACRO == "AM-2", COL_AM2, ifelse(MACRO == "AM-4", COL_AM4, "grey20"))
x_face <- ifelse(MACRO %in% c("AM-2", "AM-4"), "bold", "plain")

p <- ggplot(agg, aes(x = macro, y = program_lab, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_tile(
    data = hl, aes(x = macro, y = program_lab),
    fill = NA, colour = hl$col, linewidth = 1.1
  ) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-2, 2), oob = scales::squish,
    name = "Zscore"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = "cNMF Program") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.6),
    axis.text.x = element_text(colour = x_cols, face = x_face, size = 12),
    axis.text.y = element_text(
      colour = lab_cols[levels(agg$program_lab)],
      face = lab_face[levels(agg$program_lab)], size = 12
    ),
    axis.ticks = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "bottom",
    legend.title = element_text(size = 12, vjust = 0.9, color = "black"),
    legend.key.width = unit(2.2, "cm"),
    legend.key.height = unit(0.42, "cm")
  ) +
  guides(fill = guide_colourbar(title.position = "left"))

if (interactive()) print(p)
ggplot2::ggsave(file.path(OUT_ROOT, "Figure1g.pdf"), p, width = 5.75, height = 5.40277778, device = grDevices::cairo_pdf, bg = "white")
ggplot2::ggsave(file.path(OUT_ROOT, "Figure1g.png"), p, width = 5.75, height = 5.40277778, dpi = 180, bg = "white")
fwrite(agg[, .(program, macro, mean_usage, z, owner)], file.path(DERIVED_ROOT, "g_plotted_values.csv"))
