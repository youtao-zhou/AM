# Spatial TWAS maps. Smoothing is used for display only; statistics use raw scores.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(ggplot2)
library(patchwork)

setwd(am_data("sp"))

main_section <- "P15_LUAD"
supp_sections <- c("P5_LUAD", "P11_LUAD", "P4_LUAD", "P8_LUAD", "P21_LUAD")
spot_file <- "AM4_TWAS_visium_maps_spot_data.csv.gz"
font <- "Noto Sans"

stopifnot(file.exists(spot_file))

spot <- fread(spot_file)

st <- spot[, .(
  n_spots = .N,
  n_twas_genes = n_twas_genes[1],
  rho_AM4_TWAS = cor(AM4, twas, method = "spearman", use = "complete.obs"),
  rho_AM2_TWAS = cor(AM2, twas, method = "spearman", use = "complete.obs"),
  twas_tumour_minus_rest =
    median(twas[tumour], na.rm = TRUE) -
      median(twas[!tumour], na.rm = TRUE)
), by = section]

sec_rho <- fread("AM4_TWAS_visium_maps_section_rho.csv")
st <- merge(st, sec_rho[, .(section, partial_rho_TWAS_AM4)], by = "section")
fwrite(st, am_output("Figure4", "Figure4h_statistics.csv"))

pal <- list(
  AM4 = c("#F2F4FB", "#6888F5", "#1F2F8C"),
  AM2 = c("#FBF2F3", "#E25563", "#8C1F2B"),
  twas_smooth = c("#FBF7EE", "#D9A43A", "#7A5410")
)

map_theme <- theme_void(base_family = font, base_size = 10) +
  theme(
    text = element_text(colour = "black"),
    plot.title = element_text(size = 11, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = grid::unit(2, "pt"),
    plot.margin = margin(4, 4, 4, 4)
  )

cont_map <- function(d, v, title, legend_title) {
  lim <- quantile(d[[v]], c(0.01, 0.99), na.rm = TRUE)
  pad <- max(diff(lim) * 1e-6, 1e-8)

  ggplot(d, aes(x, y, colour = .data[[v]])) +
    geom_point(size = 0.7, shape = 15) +
    coord_equal(expand = FALSE) +
    scale_colour_gradientn(
      colours = pal[[v]],
      limits = lim + c(-pad, pad),
      oob = scales::squish,
      breaks = lim,
      labels = sprintf("%.2f", lim),
      name = legend_title,
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        direction = "horizontal",
        barwidth = grid::unit(2, "cm"),
        barheight = grid::unit(0.22, "cm"),
        ticks = FALSE
      )
    ) +
    labs(title = title) +
    map_theme
}

section_maps <- function(s) {
  d <- spot[section == s]

  tumor_map <- ggplot(d, aes(x, y, colour = tumour)) +
    geom_point(size = 0.7, shape = 15) +
    coord_equal(expand = FALSE) +
    scale_colour_manual(
      values = c(`TRUE` = "#7A6A55", `FALSE` = "#EDE9E3"),
      breaks = c("TRUE", "FALSE"),
      labels = c("Tumor", "Non-tumor"),
      name = NULL,
      guide = guide_legend(
        direction = "horizontal",
        override.aes = list(size = 2.5)
      )
    ) +
    labs(title = "Tumor") +
    map_theme

  wrap_plots(
    cont_map(d, "AM4", "AM-4", "RCTD weight"),
    cont_map(d, "AM2", "AM-2", "RCTD weight"),
    tumor_map,
    cont_map(d, "twas_smooth", "TWAS Program", "Programme score"),
    ncol = 2
  ) +
    plot_annotation(
      title = sub("_LUAD$", "", s),
      theme = theme(
        plot.title = element_text(
          family = font, colour = "black",
          size = 11, hjust = 0.5, face = "bold"
        )
      )
    )
}

p_main <- section_maps(main_section)
if (interactive()) print(p_main)
ggsave(
  am_output("Figure4", "AM4_TWAS_visium_main.pdf"),
  p_main,
  width = 5.63, height = 4.76,
  device = cairo_pdf, bg = "white"
)

for (s in supp_sections) {
  p <- section_maps(s)

  ggsave(
    file.path(am_output("Figure4", ""), paste0("AM4_TWAS_visium_", s, ".pdf")),
    p,
    width = 5.63, height = 4.76,
    device = cairo_pdf, bg = "white"
  )
}
