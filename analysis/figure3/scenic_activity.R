# scenic activity
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork)
})

data_dir <- am_data("HTAN/am/SCENIC")
out_dir <- data_dir
prefix <- "SCENIC_S4A_Tumor_Normal_original"
tf_order <- c(
  "ELF4(+)", "EP300(+)", "RXRA(+)", "SMARCA1(+)",
  "ZBTB7A(+)", "FOXO1(+)", "FOXO3(+)", "NFIL3(+)"
)
am_order <- c("AM-1", "AM-2", "AM-3", "AM-4")
environment_order <- c("Normal", "Tumor")
min_cells <- 20
font_size <- 11
font_family <- "sans"
color_limits <- c(-2, 2)

heat_colors <- c(
  "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
  "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"
)
missing_color <- "#E4E4E4"
highlight_color <- "#258A91"
highlight_AM <- "AM-2"
figure_width <- 9
figure_height <- 4.7
png_dpi <- 300
save_files <- TRUE
show_note <- TRUE

auc <- read.csv(file.path(data_dir, paste0(prefix, "_sample_mean_auc.csv")),
  check.names = FALSE, stringsAsFactors = FALSE
)
stopifnot(all(c("sample", "AM", "n_cells", "tissue", tf_order) %in% names(auc)))
selected <- auc[auc$n_cells >= min_cells & auc$AM %in% am_order &
  auc$tissue %in% environment_order, ]

long_list <- list()
raw_matrices <- list()
z_matrices <- list()
for (env in environment_order) {
  x <- selected[selected$tissue == env, ]
  m <- matrix(NA_real_,
    nrow = length(tf_order), ncol = length(am_order),
    dimnames = list(tf_order, am_order)
  )
  for (g in tf_order) for (am in am_order) {
    v <- x[x$AM == am, g]
    if (any(is.finite(v))) m[g, am] <- mean(v, na.rm = TRUE)
  }
  z <- m
  for (g in tf_order) {
    v <- m[g, ]; center <- mean(v, na.rm = TRUE)
    pop_sd <- sqrt(mean((v - center)^2, na.rm = TRUE))
    z[g, ] <- if (is.finite(pop_sd) && pop_sd > 0) (v - center) / pop_sd else NA_real_
  }
  raw_matrices[[env]] <- m
  z_matrices[[env]] <- z
  long_list[[env]] <- data.frame(
    environment = env,
    TF = rep(tf_order, times = length(am_order)),
    AM = rep(am_order, each = length(tf_order)),
    mean_auc = as.vector(m), within_environment_z = as.vector(z)
  )
}
plot_data <- do.call(rbind, long_list)
rownames(plot_data) <- NULL
plot_data$TF <- factor(plot_data$TF, levels = rev(tf_order))
plot_data$AM <- factor(plot_data$AM, levels = am_order)

make_heatmap <- function(env) {
  d <- plot_data[plot_data$environment == env, ]
  p <- ggplot(d, aes(AM, TF, fill = within_environment_z)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradientn(
      colors = heat_colors, limits = color_limits,
      breaks = seq(-2, 2, 1), oob = scales::squish,
      na.value = missing_color, name = "Within-tissue Z"
    ) +
    scale_x_discrete(limits = am_order, expand = c(0, 0)) +
    scale_y_discrete(limits = rev(tf_order), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = env) +
    theme_minimal(base_size = font_size, base_family = font_family) +
    theme(
      panel.grid = element_blank(), axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(color = "black"),
      plot.margin = margin(5, 12, 5, 5), legend.position = "right"
    )
  if (!is.null(highlight_AM) && highlight_AM %in% am_order) {
    j <- match(highlight_AM, am_order)
    p <- p + annotate("rect",
      xmin = j - 0.5, xmax = j + 0.5,
      ymin = 0.5, ymax = length(tf_order) + 0.5,
      fill = NA, color = highlight_color, linewidth = 0.65
    )
  }
  p
}
p_Normal <- make_heatmap("Normal")
p_Tumor <- make_heatmap("Tumor")
panel_list <- lapply(environment_order, make_heatmap)
fig_S4A <- wrap_plots(panel_list, nrow = 1, guides = "collect")
if (show_note) fig_S4A <- fig_S4A + plot_annotation(
  caption = paste(
    "Separate regulon networks; target sets differ.",
    "Colours are not direct Normal-Tumor AUC differences."
  ),
  theme = theme(plot.caption = element_text(size = 8, color = "grey40", hjust = 0.5))
)
print(fig_S4A)

if (save_files) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- function(s) file.path(out_dir, paste0(prefix, s))
  ggsave(out("_R.pdf"), fig_S4A,
    width = figure_width, height = figure_height,
    device = grDevices::cairo_pdf, bg = "white"
  )
  ggsave(out("_R.png"), fig_S4A,
    width = figure_width, height = figure_height,
    dpi = png_dpi, bg = "white"
  )
  write.csv(plot_data, out("_plot_values.csv"), row.names = FALSE, na = "")
  for (env in environment_order) {
    write.csv(raw_matrices[[env]], out(paste0("_", env, "_mean_auc.csv")), na = "")
    write.csv(z_matrices[[env]], out(paste0("_", env, "_z_matrix.csv")), na = "")
  }
  cov <- aggregate(sample ~ tissue + AM, selected, length)
  names(cov)[3] <- "n_sample_AM_groups"
  write.csv(cov, out("_used_sample_counts.csv"), row.names = FALSE)
  message("完成：最初S4A，Normal/Tumor并列；不是25患者配对版。")
}
