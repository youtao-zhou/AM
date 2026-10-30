# AM abundance along the tumor boundary; mean and SEM within distance bins.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(RANN); library(scales); library(grid)
})
setwd(am_data("sp/BANKSY"))
base_family <- "sans"
if (requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE)) {
  tryCatch(
    {
      ff <- sysfonts::font_files()
      sysfonts::font_add("Helvetica", ff$path[grepl("Helvetica|Arial", ff$family)][1])
      showtext::showtext_auto(); base_family <- "Helvetica"
    },
    error = function(e) NULL
  )
}
FEAT <- c(AM4 = "AM-4", AM2 = "AM-2", CD8 = "CD8 T")
COL <- c("AM-4" = "#EF767B", "AM-2" = "#43A3EF", "CD8 T" = "#E8C24A")

D <- read.csv(gzfile("data_domain_am_cd8_26.csv.gz"), stringsAsFactors = FALSE)

rows <- list()
for (s in unique(D$sample)) {
  d <- D[D$sample == s, ]
  tum <- d$compartment == "Tumor"
  if (sum(tum) < 30 || sum(!tum) < 30) next
  XY <- cbind(d$array_col * 0.5, as.numeric(d$array_row))
  dist <- numeric(nrow(d))
  dist[!tum] <- RANN::nn2(XY[tum, , drop = FALSE], XY[!tum, , drop = FALSE], k = 1)$nn.dists[, 1]
  dist[tum] <- -RANN::nn2(XY[!tum, , drop = FALSE], XY[tum, , drop = FALSE], k = 1)$nn.dists[, 1]
  dd <- data.frame(sample = s, dist = dist)
  for (fc in names(FEAT)) dd[[fc]] <- d[[fc]]
  rows[[length(rows) + 1]] <- dd
}
ALL <- bind_rows(rows)
nL <- length(unique(ALL$sample)); message("LUAD used = ", nL, " spots = ", nrow(ALL))

core_min <- -8;  peri_max <- 8
FEATcols <- names(FEAT)

PF <- ALL %>%
  mutate(bin = pmin(pmax(round(dist), core_min), peri_max)) %>%
  group_by(sample) %>%
  mutate(across(all_of(FEATcols), ~ .x / mean(.x, na.rm = TRUE), .names = "fold_{.col}")) %>%
  ungroup() %>%
  pivot_longer(starts_with("fold_"),
    names_to = "feat",
    values_to = "fold", names_prefix = "fold_"
  ) %>%
  group_by(feat, bin) %>%
  summarise(
    n = sum(is.finite(fold)),
    m = mean(fold, na.rm = TRUE),
    sem = sd(fold, na.rm = TRUE) / sqrt(n), .groups = "drop"
  ) %>%
  mutate(feat = factor(dplyr::recode(feat, !!!FEAT),
    levels = c("AM-4", "AM-2", "CD8 T")
  ))

xr <- range(PF$bin)
ytop <- max(PF$m + PF$sem); ymin <- min(PF$m - PF$sem)
y_arr <- ytop + 0.07
y_lab <- ytop + 0.15

th <- theme_minimal(base_size = 15, base_family = base_family) +
  theme(
    text               = element_text(color = "black"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = .2),
    axis.line          = element_line(color = "black", linewidth = .2),
    axis.ticks         = element_line(color = "black", linewidth = .2),
    axis.ticks.length  = unit(3.5, "pt"),
    axis.title.x       = element_text(size = 16, color = "black", margin = margin(t = 10)),
    axis.title.y       = element_text(size = 16, color = "black", margin = margin(r = 10)),
    axis.text          = element_text(size = 14, color = "black"),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 14, color = "black"),
    legend.key.width   = unit(30, "pt"),
    legend.margin      = margin(t = 6),
    plot.margin        = margin(14, 18, 10, 14)
  )

pA1 <- ggplot(PF, aes(bin, m, color = feat, fill = feat)) +
  geom_vline(xintercept = 0, linetype = "22", color = "grey55", linewidth = .5) +
  geom_ribbon(aes(ymin = m - sem, ymax = m + sem), alpha = .1, color = NA) +
  geom_line(linewidth = 1.2, lineend = "round") +
  geom_point(shape = 21, size = 2.2, color = "grey25", stroke = .5) +
  annotate("segment",
    x = -0.4, xend = xr[1] + 0.3, y = y_arr, yend = y_arr,
    arrow = arrow(length = unit(8, "pt"), type = "closed"),
    linewidth = .35, color = "grey30"
  ) +
  annotate("segment",
    x = 0.4, xend = xr[2] - 0.3, y = y_arr, yend = y_arr,
    arrow = arrow(length = unit(8, "pt"), type = "closed"),
    linewidth = .35, color = "grey30"
  ) +
  annotate("text",
    x = xr[1], y = y_lab, label = "Tumor core", hjust = 0,
    size = 5, color = "black", family = base_family
  ) +
  annotate("text",
    x = xr[2], y = y_lab, label = "Periphery", hjust = 1,
    size = 5, color = "black", family = base_family
  ) +
  scale_color_manual(values = COL) +
  scale_fill_manual(values = COL) +
  scale_x_continuous(breaks = seq(-8, 8, 2), expand = expansion(mult = .012)) +
  scale_y_continuous(breaks = seq(0.5, 1.5, 0.5)) +
  coord_cartesian(ylim = c(ymin - 0.03, y_lab + 0.05), clip = "off") +
  labs(x = "Distance to tumor boundary", y = "Abundance") +
  th

if (interactive()) print(pA1)
ggsave(am_output("Figure2", "Figure2d.pdf"), pA1, width = 6.5, height = 4.5, device = cairo_pdf)
write.csv(PF, am_output("Figure2", "Figure2d_values.csv"), row.names = FALSE)
