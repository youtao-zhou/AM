# Tumor contact and expression contrasts from the saved single-cell tables.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(ggplot2)

setwd(am_data("sp"))

data_dir <- "Figure2_AM24_ABC/data"
out_dir <- am_output("Figure2", "jk")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pairs <- read.csv(
  file.path(data_dir, "fixed_pairs.csv"),
  stringsAsFactors = FALSE
)

dot <- read.csv(
  file.path(data_dir, "A_pair_expression.csv"),
  stringsAsFactors = FALSE
)

dot$display <- factor(
  dot$display,
  levels = pairs$display
)

dot$role <- factor(
  dot$role,
  levels = c("AM-4", "AM-2", "Tumor")
)

pair_labels <- setNames(
  paste0(
    pairs$display,
    ifelse(pairs$scope == "extra4", " \u2020", "")
  ),
  pairs$display
)

pA <- ggplot(dot, aes(x = display, y = role)) +
  geom_tile(
    fill = "white",
    color = "#D5D5D5",
    linewidth = 0.3,
    width = 1,
    height = 1
  ) +
  geom_point(
    data = dot[!is.na(dot$detection) & dot$detection > 0, ],
    aes(
      size = detection * 100,
      fill = expression
    ),
    shape = 21,
    color = "#454545",
    stroke = 0.3
  ) +
  scale_size_area(
    name = "Cells expressing (%)",
    max_size = 8,
    limits = c(0, 100),
    breaks = c(10, 30, 60)
  ) +
  scale_fill_gradientn(
    name = "Expression",
    colours = c("#F5EDE7", "#DE8F80", "#A52A42"),
    space = "Lab",
    na.value = "#BDBDBD"
  ) +
  scale_x_discrete(
    labels = pair_labels,
    drop = FALSE,
    expand = expansion(mult = 0, add = 0)
  ) +
  scale_y_discrete(
    drop = FALSE,
    expand = expansion(mult = 0, add = 0)
  ) +
  guides(
    size = guide_legend(
      order = 1,
      direction = "horizontal",
      nrow = 1,
      title.position = "top",
      title.hjust = 0,
      override.aes = list(
        shape = 21,
        fill = "#DE8F80",
        color = "#454545",
        stroke = 0.3
      )
    ),
    fill = guide_colorbar(
      order = 2,
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0,
      barwidth = grid::unit(30, "mm"),
      barheight = grid::unit(3, "mm")
    )
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_family = "Arial",
    base_size = 11
  ) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9,
      color = "#222222",
      margin = margin(t = 4)
    ),
    axis.text.y = element_text(
      size = 10,
      color = "#222222",
      margin = margin(r = 4)
    ),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(size = 8),
    legend.key.width = grid::unit(8, "mm"),
    legend.key.height = grid::unit(8, "mm"),
    legend.spacing.x = grid::unit(3, "mm"),
    legend.box.spacing = grid::unit(1, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(4, 4, 4, 8)
  )

if (interactive()) print(pA)

ggsave(
  filename = am_output("Figure2", "scRNA_AM24_Tumor_Contact.pdf"),
  plot = pA,
  width = 7.36,
  height = 3.24,
  units = "in",
  device = cairo_pdf,
  family = "Arial",
  bg = "white"
)

library(ggplot2)

setwd(am_data("sp"))

data_dir <- "Figure2_AM24_ABC/data"
out_dir <- am_output("Figure2", "jk")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr <- read.csv(
  file.path(data_dir, "all_patient_expression.csv"),
  stringsAsFactors = FALSE
)

stat <- read.csv(
  file.path(data_dir, "all20_unadjusted_between_group_tests.csv"),
  stringsAsFactors = FALSE
)

expr$AM <- factor(
  expr$AM,
  levels = c("AM-2", "AM-4")
)

contact_genes <- c("PTPRM", "ADGRE5", "NOTCH2")
immune_genes <- c("HLA-DRA", "HLA-DPA1", "CD74")

cols <- c(
  "AM-2" = "#258CA5",
  "AM-4" = "#C9503F"
)

gene_labels <- c(
  "PTPRM" = "PTPRM",
  "ADGRE5" = "CD97",
  "NOTCH2" = "NOTCH2",
  "HLA-DRA" = "HLA-DRA",
  "HLA-DPA1" = "HLA-DPA1",
  "CD74" = "CD74"
)

expression_plot <- function(genes, from_zero = TRUE) {
  df <- expr[expr$entity %in% genes, ]
  df$entity <- factor(df$entity, levels = genes)

  ann <- data.frame(
    entity = genes,
    p = stat$p[match(genes, stat$entity)]
  )

  ymax <- sapply(genes, function(g) {
    max(df$AM_expression[df$entity == g], na.rm = TRUE)
  })

  ymin <- sapply(genes, function(g) {
    min(df$AM_expression[df$entity == g], na.rm = TRUE)
  })

  span <- if (from_zero) ymax else ymax - ymin
  span[span == 0] <- pmax(abs(ymax[span == 0]) * 0.1, 0.1)

  ann$bracket <- ymax + 0.08 * span
  ann$tip <- ymax + 0.04 * span
  ann$text_y <- ymax + 0.15 * span
  ann$entity <- factor(ann$entity, levels = genes)

  ann$label <- paste0(
    "p = ",
    ifelse(
      ann$p < 0.001,
      formatC(ann$p, format = "e", digits = 1),
      formatC(ann$p, format = "f", digits = 3)
    )
  )

  p <- ggplot(
    df,
    aes(x = AM, y = AM_expression, color = AM, fill = AM)
  ) +
    geom_boxplot(
      width = 0.5,
      outlier.shape = NA,
      alpha = 0.18,
      linewidth = 0.5
    ) +
    geom_point(
      position = position_jitter(
        width = 0.12,
        height = 0,
        seed = 38
      ),
      size = 1.5,
      alpha = 0.65
    ) +
    geom_segment(
      data = ann,
      aes(x = 1, xend = 2, y = bracket, yend = bracket),
      inherit.aes = FALSE,
      linewidth = 0.25
    ) +
    geom_segment(
      data = ann,
      aes(x = 1, xend = 1, y = tip, yend = bracket),
      inherit.aes = FALSE,
      linewidth = 0.25
    ) +
    geom_segment(
      data = ann,
      aes(x = 2, xend = 2, y = tip, yend = bracket),
      inherit.aes = FALSE,
      linewidth = 0.25
    ) +
    geom_text(
      data = ann,
      aes(x = 1.5, y = text_y, label = label),
      inherit.aes = FALSE,
      family = "Arial",
      size = 3.3,
      vjust = 0
    ) +
    facet_wrap(
      ~entity,
      ncol = 3,
      scales = "free_y",
      labeller = as_labeller(gene_labels)
    ) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    scale_y_continuous(
      expand = expansion(
        mult = c(if (from_zero) 0 else 0.06, 0.15)
      )
    ) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = NULL,
      y = "Mean Expression"
    ) +
    theme_classic(
      base_family = "Arial",
      base_size = 12
    ) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(
        size = 11,
        face = "bold",
        color = "#222222"
      ),
      axis.text.x = element_text(
        size = 11,
        color = "#222222",
        margin = margin(t = 5)
      ),
      axis.text.y = element_text(
        size = 11,
        color = "#222222"
      ),
      axis.title.y = element_text(
        size = 12,
        margin = margin(r = 8)
      ),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      panel.border = element_blank(),
      panel.spacing.x = grid::unit(6, "mm"),
      plot.margin = margin(8, 10, 8, 8)
    )

  if (from_zero) {
    p <- p + geom_blank(aes(y = 0))
  }

  return(p)
}

pB <- expression_plot(
  genes = contact_genes,
  from_zero = TRUE
)

pC <- expression_plot(
  genes = immune_genes,
  from_zero = FALSE
)

if (interactive()) print(pB)
if (interactive()) print(pC)

ggsave(
  filename = am_output("Figure2", "scRNA_AM24_contact_boxplot1.pdf"),
  plot = pB,
  width = 8.76,
  height = 3.74,
  units = "in",
  device = cairo_pdf,
  family = "Arial",
  bg = "white"
)

ggsave(
  filename = am_output("Figure2", "scRNA_AM24_contact_boxplot2.pdf"),
  plot = pC,
  width = 8.76,
  height = 3.74,
  units = "in",
  device = cairo_pdf,
  family = "Arial",
  bg = "white"
)
