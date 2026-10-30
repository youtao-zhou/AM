# AM abundance and T-cell programs in 25 patients, using the atlas reference.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(ggplot2)

setwd(am_data("sp"))

REFERENCE <- "atlas"
font <- "Noto Sans"
ink <- "#20292F"
axis_lwd <- 0.15
cols <- c("AM-2-rich" = "#E25563", "AM-4-rich" = "#2E44B8")
programmes <- c(CD8_CTL = "CD8 CTL", CD8_Tm = "CD8 Tm", Treg = "Treg")

fig_width <- 7.17
fig_height <- 3.20

sfx <- if (REFERENCE == "atlas") "_atlas" else ""

d <- fread(sprintf("AM4_AM2_visium_T_programme_scores_patient_data%s.csv", sfx))[
  program %in% names(programmes)
]
d[, group := factor(group, levels = c("AM-2-rich", "AM-4-rich"))]

w <- dcast(d, patient + program ~ group, value.var = "score")

st <- w[, .(
  n = .N,
  median_AM2_rich = median(`AM-2-rich`),
  median_AM4_rich = median(`AM-4-rich`),
  n_AM4_higher = sum(`AM-4-rich` > `AM-2-rich`),
  P = wilcox.test(
    `AM-4-rich`, `AM-2-rich`,
    paired = TRUE, exact = TRUE
  )$p.value
), by = program]

st[, label := ifelse(
  P < 0.001,
  sprintf("italic(P) == %.1f %%*%% 10^%d", P / 10^floor(log10(P)), as.integer(floor(log10(P)))),
  sprintf("italic(P) == \"%.3f\"", P)
)]

fwrite(st, am_output("Figure4", "Figure4f_statistics.csv"))

d[, program := factor(programmes[program], levels = unname(programmes))]
st[, program := factor(programmes[program], levels = unname(programmes))]

p <- ggplot(d, aes(group, score)) +
  geom_line(
    aes(group = patient),
    colour = "grey82", linewidth = 0.25
  ) +
  geom_boxplot(
    aes(fill = group),
    width = 0.55, outlier.shape = NA,
    alpha = 0.5, linewidth = 0.35
  ) +
  geom_point(
    aes(fill = group),
    shape = 21, size = 1.1, stroke = 0.2, colour = ink
  ) +
  facet_wrap(~program, nrow = 1, scales = "free_y") +
  geom_text(
    data = st,
    aes(x = 1.5, y = Inf, label = label),
    vjust = 1.3, size = 2.8, family = font, parse = TRUE,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_discrete(labels = c(
    "AM-2-rich" = "AM-2",
    "AM-4-rich" = "AM-4"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = "Programme score") +
  theme_classic(base_size = 8, base_family = font) +
  theme(
    axis.text = element_text(colour = "black", size = 9),
    axis.title = element_text(colour = "black", size = 9.5),
    axis.line = element_line(linewidth = axis_lwd),
    axis.ticks = element_line(linewidth = axis_lwd),
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "bold")
  )

if (interactive()) print(p)

stem <- am_output("Figure4", "Figure4f")
ggsave(paste0(stem, ".pdf"), p, width = fig_width, height = fig_height, device = cairo_pdf, bg = "white")
ggsave(paste0(stem, ".png"), p, width = fig_width, height = fig_height, dpi = 300, bg = "white")
print(st[, .(program, n, n_AM4_higher, P)])
