# Paired pathway scores for AM-2, AM-4 and other AM states.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
source(am_code("R/figure1_data.R"))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

DEVICE <- "pdf"
font <- if (DEVICE == "pdf") "Helvetica" else if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "Liberation Sans"
fs <- 12
states <- c(AM2 = "AM-2", OtherAM = "Other AMs", AM4 = "AM-4")
if (DEVICE == "pdf") states <- gsub("-", "\u00ad", states, fixed = TRUE)
pathways <- c(
  "Th1 and Th2 cell differentiation", "Th17 cell differentiation", "Antigen processing and presentation",
  "Phagosome", "Cell adhesion molecules", "Endocytosis", "Natural killer cell mediated cytotoxicity"
)

rs <- rbindlist(lapply(names(states), function(s) fread(file.path(INPUT_ROOT, "j", sprintf("Fig1j_GRSA_%s.csv", s)))[, set := s]), fill = TRUE)
df <- rs[Description %in% pathways]
df[, Description := factor(Description, rev(pathways))]
df[, state := factor(states[set], states)]
df[, dir := fcase(ReporterScore >= 1.64, "Up in Tumor", ReporterScore <= -1.64, "Down in Tumor", default = "ns")]
fwrite(df[, .(state, Description, ReporterScore, p.value, dir)], file.path(DERIVED_ROOT, "j_plotted_values.csv"))
p <- ggplot(df, aes(ReporterScore, Description)) +
  geom_vline(xintercept = c(-1.64, 1.64), linetype = 2, linewidth = 0.18, colour = "grey50") +
  geom_vline(xintercept = 0, linewidth = 0.18) +
  geom_segment(aes(x = 0, xend = ReporterScore, yend = Description), colour = "grey70", linewidth = 0.36) +
  geom_point(aes(fill = dir), shape = 21, size = 3.1, stroke = 0.45) +
  facet_wrap(~state, nrow = 1) +
  scale_fill_manual(
    values = c("Up in Tumor" = "#C0392B", "Down in Tumor" = "#2E6DB4", ns = "white"),
    breaks = c("Up in Tumor", "Down in Tumor", "ns"), name = NULL
  ) +
  guides(fill = guide_legend(override.aes = list(size = 4.5, stroke = 0.6))) +
  scale_x_continuous(breaks = c(-6, -3, 0, 3)) +
  labs(x = "Reporter score", y = NULL) +
  theme_classic(base_size = fs, base_family = font) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(size = fs, colour = "black"),
    axis.line = element_line(linewidth = 0.36),
    axis.ticks = element_line(linewidth = 0.36),
    strip.background = element_blank(),
    strip.text = element_text(size = fs),
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.text = element_text(size = fs, margin = margin(l = 3, r = 12)),
    legend.key.size = unit(16, "pt")
  )
if (interactive()) print(p)
if (DEVICE == "pdf") {
  ggsave(file.path(OUT_ROOT, "Figure1j.pdf"), p, width = 8.4, height = 3.35, device = pdf, useDingbats = FALSE, version = "1.4", bg = "white")
} else {
  ggsave(file.path(OUT_ROOT, "Figure1j.pdf"), p, width = 8.4, height = 3.35, device = cairo_pdf, bg = "white")
}

df_png <- copy(df)[, state := factor(gsub("\u00ad", "-", as.character(state), fixed = TRUE), gsub("\u00ad", "-", levels(state), fixed = TRUE))]
ggsave(file.path(OUT_ROOT, "Figure1j.png"), p + df_png, width = 8.4, height = 3.35, dpi = 300, bg = "white")
