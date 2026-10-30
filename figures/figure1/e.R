# Paired tumor/normal pathway scores.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
source(am_code("R/figure1_data.R"))
library(data.table)
library(ggplot2)
library(stringr)

rs <- fread(file.path(INPUT_ROOT, "e", "Fig1e_GRSA_pseudobulk_all_pathways.csv"))
luad_set <- c(
  "Antigen processing and presentation", "Phagosome", "Cell adhesion molecules",
  "Th1 and Th2 cell differentiation", "Th17 cell differentiation"
)
luad <- rs[Description %in% luad_set]
luad$side <- "LUAD"
normal <- rs[ReporterScore < 0][order(ReporterScore)][1:5]
normal$side <- "Normal"
d <- rbind(luad, normal)
stopifnot(nrow(d) == 10, all(d[side == "LUAD"]$ReporterScore > 0))
d$Group <- factor(ifelse(d$ReporterScore > 0, "Tumor AM", "Normal AM"), levels = c("Tumor AM", "Normal AM"))

setorder(d, ReporterScore)
d$label <- factor(str_wrap(d$Description, 34), levels = str_wrap(d$Description, 34))
fwrite(
  d[, .(Description, side, ReporterScore, Z_score, p.value, p.adjust, Exist_K_num, Significant_K_num)],
  file.path(DERIVED_ROOT, "e_plotted_values.csv")
)

font <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "Liberation Sans"
fs <- 11
axis_lw <- 0.2
cols <- c("Normal AM" = "#4DBBD5", "Tumor AM" = "#E64B35")
xmax <- ceiling(max(abs(d$ReporterScore)))

p1e <- ggplot(d, aes(x = ReporterScore, y = label, colour = Group)) +
  geom_vline(xintercept = c(-1.64, 1.64), linetype = "dashed", colour = "grey70", linewidth = axis_lw) +
  geom_vline(xintercept = 0, colour = "black", linewidth = axis_lw) +
  geom_segment(aes(x = 0, xend = ReporterScore, yend = label), colour = "grey82", linewidth = 1.2, lineend = "round") +
  geom_point(size = 3) +
  scale_colour_manual(values = cols, name = "Enriched in") +
  scale_x_continuous(limits = c(-6, xmax), breaks = seq(-4, xmax, 4), expand = expansion(mult = c(0.02, 0.03))) +
  labs(x = "ReporterScore", y = NULL) +
  theme_classic(base_size = fs, base_family = font) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black", size = fs),
    axis.title = element_text(size = fs),
    axis.line.x = element_line(linewidth = axis_lw),
    axis.ticks.x = element_line(linewidth = axis_lw),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    legend.title = element_text(size = fs),
    legend.text = element_text(size = fs),
    legend.position = "inside", legend.position.inside = c(1, 0.02), legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.size = unit(11, "pt"), legend.key.spacing.y = unit(0, "pt"),
    plot.margin = margin(4, 8, 2, 2)
  )

print(d[, .(Description, side, RS = round(ReporterScore, 2), p.value)])
if (interactive()) print(p1e)
ggsave(file.path(OUT_ROOT, "Figure1e.pdf"), p1e, width = 4.66, height = 3.69, device = cairo_pdf, bg = "white")
