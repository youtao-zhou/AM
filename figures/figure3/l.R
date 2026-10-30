# AM-2/CD8 distances by pathological response; exact stratified comparisons.
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(ggplot2)

setwd(am_data("tsa"))

d <- fread("data/AM2_CD8_post21_relative_distance_patient_data.csv",
  colClasses = list(character = "patient_id")
)
groups <- c("NR", "MPR", "PCR")
d[, response := factor(response_group, groups)]
d[, y := gm1p_percent]

strat_p <- function(g1, g2) {
  x <- d[response_group %in% c(g1, g2)]
  sc <- numeric(nrow(x))
  alloc <- list()
  for (ids in split(seq_len(nrow(x)), x$delivery_source)) {
    n <- length(ids); k <- sum(x$response_group[ids] == g1)
    sc[ids] <- (rank(x$y[ids]) - (n + 1) / 2) / (n + 1)
    alloc[[length(alloc) + 1]] <- if (k == 0) 0 else if (k == n) sum(sc[ids]) else
      combn(ids, k, function(j) sum(sc[j]))
  }
  null <- Reduce(function(a, b) as.vector(outer(a, b, "+")), alloc)
  obs <- sum(sc[x$response_group == g1])
  mean(abs(null) >= abs(obs) - 1e-12)
}

cmp <- data.table(g1 = c("NR", "MPR", "NR"), g2 = c("MPR", "PCR", "PCR"))
cmp[, p := mapply(strat_p, g1, g2)]
print(cmp)
cmp <- cmp[p < 0.05]

sm <- d[, .(mean = mean(y), se = sd(y) / sqrt(.N)), by = response]

top <- max(d$y, sm$mean + sm$se)
tick <- 0.02 * top
cmp[, `:=`(x1 = match(g1, groups), x2 = match(g2, groups))]
cmp <- cmp[order(x2 - x1)]
cmp[, yb := top * (1.08 + 0.1 * (seq_len(.N) - 1))]
ytop <- if (nrow(cmp)) max(cmp$yb) + 0.1 * top else top * 1.05

fmt_p <- function(p) ifelse(p < 0.01, "P < 0.01", sprintf("P = %.2f", p))

cols <- c(NR = "#5BB5AC", MPR = "#DAC6CF", PCR = "#DE526C")
ink <- "#20292F"
fnt <- "Noto Sans"
n <- table(d$response)
labs_x <- setNames(sprintf("%s\n(n = %d)", c("NR", "MPR", "pCR"), n[groups]), groups)

p <- ggplot(d, aes(response, y)) +
  geom_col(
    data = sm, aes(response, mean, fill = response), width = 0.55,
    colour = ink, linewidth = 0.4, alpha = 0.85
  ) +
  geom_errorbar(
    data = sm, aes(response, ymin = mean - se, ymax = mean + se), inherit.aes = FALSE,
    width = 0.18, linewidth = 0.4, colour = ink
  ) +
  geom_hline(yintercept = 100, colour = "grey55", linetype = "22", linewidth = 0.4) +
  geom_point(aes(fill = response),
    shape = 21, size = 2.3, stroke = 0.35, colour = ink,
    position = position_jitter(width = 0.1, height = 0, seed = 109)
  ) +
  geom_segment(
    data = cmp, aes(x = x1, xend = x2, y = yb, yend = yb), inherit.aes = FALSE,
    linewidth = 0.35, colour = ink
  ) +
  geom_segment(
    data = cmp, aes(x = x1, xend = x1, y = yb, yend = yb - tick), inherit.aes = FALSE,
    linewidth = 0.35, colour = ink
  ) +
  geom_segment(
    data = cmp, aes(x = x2, xend = x2, y = yb, yend = yb - tick), inherit.aes = FALSE,
    linewidth = 0.35, colour = ink
  ) +
  geom_text(
    data = cmp, aes(x = (x1 + x2) / 2, y = yb + 2.5 * tick, label = fmt_p(p)),
    inherit.aes = FALSE, size = 3.4, family = fnt, colour = ink
  ) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_discrete(labels = labs_x) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, ytop), clip = "off") +
  labs(x = NULL, y = "AM2-CD8 distance (%)") +
  theme_classic(base_size = 11, base_family = fnt) +
  theme(
    axis.text = element_text(colour = ink, size = 10),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.line = element_line(colour = ink, linewidth = 0.25),
    axis.ticks = element_line(colour = ink, linewidth = 0.25),
    plot.margin = margin(10, 14, 10, 10)
  )
if (interactive()) print(p)
ggsave(am_output("Figure3", "AM2_CD8_post21_distance_bar.pdf"), p, width = 5.44, height = 3.30, device = cairo_pdf, bg = "white")

fwrite(sm, am_output("Figure3", "Figure3l_summary.csv"))
fwrite(cmp, am_output("Figure3", "Figure3l_tests.csv"))
