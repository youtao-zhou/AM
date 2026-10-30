# Predicted treatment HR from the OAK model with both score interactions.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(survival)
library(ggplot2)

setwd(am_data("HTAN/am/TWAS_signature"))

d <- fread("data/Fig5f_patients.csv")

d[, TWz := TWAS_z]
d[, AMz := AM_z]
d$io <- ifelse(d$arm == "Atezolizumab", 1, 0)

fit <- coxph(Surv(pfs_months, pfs_event) ~ io * TWz + io * AMz, data = d[cohort == "OAK" & pfs_months > 0])
b <- coef(fit)
b0 <- b[["io"]]
bT <- b[["io:TWz"]]
bA <- b[["io:AMz"]]

stopifnot(abs(c(b0, bT, bA) - c(0.0772351988963016, 0.199310968241257, -0.232134428875549)) < 1e-6)
coef_tab <- data.table(term = names(b), coef = b, HR = exp(b), P = summary(fit)$coefficients[, 5], n = fit$n)
fwrite(coef_tab, am_output("Figure5", "Fig5e_oak_model_coefficients.csv"))
print(coef_tab)

d$pred_HR <- exp(b0 + bT * d$TWz + bA * d$AMz)
fwrite(d[, .(cohort, sample, arm, AM_z = AMz, TWAS_z = TWz, pred_HR)], am_output("Figure5", "Fig5e_patient_predicted_hr.csv"))
rho <- cor(d$AMz, d$TWz, method = "spearman")
cat(sprintf("OAK model: log HR = %.3f + %.3f * TWAS_z %+.3f * AM_z\n", b0, bT, bA))
cat(sprintf("Spearman rho (AM vs TWAS) = %.2f\n", rho))

font <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "Liberation Sans"
fs <- 11
ts <- fs / .pt
axis_lw <- 0.3
col_icb <- "#E25563"
col_chemo <- "#2E44B8"

tint <- function(col, w) rgb(t(col2rgb(col) * (1 - w) + 255 * w), maxColorValue = 255)

xr <- c(-3.6, 3.6)
yr <- c(-5.4, 3.6)
grd <- CJ(x = seq(xr[1], xr[2], length.out = 220), y = seq(yr[1], yr[2], length.out = 270))
cat("Points outside the plotted range:", d[AMz < xr[1] | AMz > xr[2] | TWz < yr[1] | TWz > yr[2], .N], "\n")
stopifnot(d[AMz < xr[1] | AMz > xr[2] | TWz < yr[1] | TWz > yr[2], .N] == 0)
grd$lhr <- b0 + bT * grd$y + bA * grd$x

lev <- c(0.5, 0.7, 1, 1.4, 2)
iso <- data.table()
for (h in lev) {
  yt <- 2.55
  xt <- (log(h) - b0 - bT * yt) / bA
  if (abs(xt) > 3.1) {
    xt <- sign(xt) * 3.1
    yt <- (log(h) - b0 - bA * xt) / bT
  }
  iso <- rbind(iso, data.table(h = h, x = xt, y = yt, lab = format(h)))
}
iso <- iso[y <= 3.4 & y >= yr[1] + 0.2]

corner_label <- function(x, y, lab, col, hj, vj) {
  annotate("label",
    x = x, y = y, label = lab, hjust = hj, vjust = vj, colour = col,
    fontface = "bold", size = ts, family = font, fill = scales::alpha("white", 0.85),
    border.colour = NA, label.padding = unit(2, "pt"), label.r = unit(1.5, "pt")
  )
}
p5e <- ggplot(grd, aes(x, y)) +
  geom_raster(aes(fill = lhr), interpolate = TRUE) +
  geom_contour(aes(z = lhr), breaks = log(c(0.5, 0.7, 1.4, 2)), colour = "grey35", linetype = "22", linewidth = 0.2) +
  geom_contour(aes(z = lhr), breaks = 0, colour = "grey15", linewidth = 0.35) +
  geom_point(
    data = d, aes(AMz, TWz, shape = cohort), inherit.aes = FALSE,
    fill = "grey20", colour = "white", stroke = 0.15, size = 1.3, alpha = 0.85
  ) +
  geom_label(
    data = iso, aes(x, y, label = lab), size = ts, family = font,
    fontface = ifelse(iso$h == 1, "bold", "plain"), fill = scales::alpha("white", 0.85),
    border.colour = NA, label.padding = unit(1, "pt"), label.r = unit(1, "pt")
  ) +
  annotate("rect", xmin = -3.55, xmax = -1.7, ymin = 2.45, ymax = 3.55, fill = scales::alpha("white", 0.85)) +
  annotate("point", x = -3.3, y = c(3.25, 2.75), shape = c(21, 24), fill = "grey20", colour = "white", stroke = 0.15, size = 1.6) +
  annotate("text", x = -3.05, y = c(3.25, 2.75), label = c("OAK", "POPLAR"), hjust = 0, size = ts, family = font) +
  scale_shape_manual(values = c(OAK = 21, POPLAR = 24), guide = "none") +
  scale_fill_gradient2(
    low = tint(col_icb, 0.45), mid = "white", high = tint(col_chemo, 0.45),
    midpoint = 0, limits = log(c(0.5, 2)), oob = scales::squish,
    breaks = log(c(0.5, 1, 2)), labels = c("0.5", "1", "2"),
    name = "Predicted PFS HR"
  ) +
  scale_x_continuous(breaks = c(-3, 0, 3)) +
  scale_y_continuous(breaks = c(-3, 0, 3)) +
  coord_cartesian(xlim = xr, ylim = yr, expand = FALSE) +
  labs(x = "AM Score", y = "TWAS Score") +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, ticks.colour = "black")) +
  theme_classic(base_size = fs, base_family = font) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black", size = fs),
    axis.title = element_text(size = fs),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = axis_lw),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = axis_lw),
    legend.position = "bottom",
    legend.title = element_text(size = fs),
    legend.text = element_text(size = fs),
    legend.key.width = unit(14, "pt"),
    legend.key.height = unit(7, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 6, 2, 2)
  )

if (interactive()) print(p5e)
ggsave(am_output("Figure5", "predicted_hr_map.pdf"), p5e, width = 3.92, height = 4.36, device = cairo_pdf, bg = "white")
ggsave(am_output("Figure5", "predicted_hr_map.png"), p5e, width = 3.92, height = 4.36, dpi = 300, bg = "white")

fwrite(grd, am_output("Figure5", "Figure5f_grid.csv"))
