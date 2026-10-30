# AM-4 abundance by response. Exact trend P = 1940153/37837800 (> 0.05).
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(ggplot2)

setwd(am_data("tsa"))

font <- "Noto Sans"
ink <- "#20292F"
axis_lwd <- 0.25
cols <- c(NR = "#5BB5AC", MPR = "#DAC6CF", PCR = "#DE526C")
p_digits <- 2
fig_width <- 2.6
fig_height <- 2.8

d <- fread("data/AM4_post21_NR_MPR_pCR_patient_data.csv", colClasses = c(patient_id = "character"))
d[, response := factor(response, c("NR", "MPR", "PCR"))]

rho <- cor(d$am4_pct_mac, as.integer(d$response), method = "spearman")
p_exact <- 1940153 / 37837800
stopifnot(abs(rho - (-0.4318)) < 1e-3)
st <- data.table(
  comparison = "Trend across NR, MPR, pCR (coded 1, 2, 3)",
  test = "Spearman correlation, exact permutation of response within delivery batch, two-sided",
  n_NR = sum(d$response == "NR"), n_MPR = sum(d$response == "MPR"), n_pCR = sum(d$response == "PCR"),
  median_NR = median(d[response == "NR", am4_pct_mac]), median_MPR = median(d[response == "MPR", am4_pct_mac]),
  median_pCR = median(d[response == "PCR", am4_pct_mac]), rho = rho, P = p_exact
)
fwrite(st, am_output("Figure4", "Figure4d_statistics.csv"))
lab <- sprintf(paste0("P = %.", p_digits, "f"), st$P)
n_by <- d[, .N, by = response]
xlab <- setNames(sprintf("%s\n(n = %d)", c(NR = "NR", MPR = "MPR", PCR = "pCR")[as.character(n_by$response)], n_by$N), as.character(n_by$response))

p <- ggplot(d, aes(response, am4_pct_mac, fill = response)) +
  geom_boxplot(
    width = 0.55, outlier.shape = NA,
    alpha = 0.55, linewidth = 0.35
  ) +
  stat_summary(
    aes(group = 1),
    fun = median, geom = "line",
    linewidth = 0.4, linetype = "22", colour = ink
  ) +
  geom_point(
    shape = 21, size = 1.3, stroke = 0.25, colour = ink,
    position = position_jitter(width = 0.09, seed = 140)
  ) +
  annotate(
    "text",
    x = 2, y = 109, label = lab,
    size = 3.1, family = font, colour = "black"
  ) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_discrete(labels = xlab) +
  scale_y_continuous(limits = c(0, 113), breaks = seq(0, 100, 25)) +
  labs(
    x = NULL,
    y = "AM-4 / macrophages (%)",
    subtitle = "Neoadjuvant mIF"
  ) +
  theme_classic(base_size = 9, base_family = font) +
  theme(
    axis.text = element_text(colour = "black", size = 9),
    axis.title = element_text(colour = "black", size = 10),
    axis.line = element_line(linewidth = 0.15),
    axis.ticks = element_line(linewidth = 0.15),
    plot.subtitle = element_text(
      colour = "black", size = 9, hjust = 0.5
    )
  )

if (interactive()) print(p)

ggsave(
  am_output("Figure4", "AM4_post21_NR_MPR_pCR_boxplot.pdf"),
  p,
  width = 4.67, height = 3.20,
  device = cairo_pdf, bg = "white"
)
ggsave(
  am_output("Figure4", "AM4_post21_NR_MPR_pCR_boxplot.png"),
  p,
  width = 4.67, height = 3.20,
  dpi = 300, bg = "white", type = "cairo"
)
