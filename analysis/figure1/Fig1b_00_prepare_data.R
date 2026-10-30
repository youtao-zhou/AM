# Fig1b 00 prepare data
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(dplyr)

setwd(am_data("HTAN/scPagwas/scPagwas_lc"))
dir.create("data", showWarnings = FALSE)
src <- am_data("vibe/191_scpagwas_hg19")

rd <- function(f) read.csv(f, check.names = FALSE)[, -1]
tumor <- rd(file.path(src, "tumor_allcelltype/subset1_celltypes_bootstrap_results_hg19.csv"))
normal <- rd(file.path(src, "normal_allcelltype/subset1_celltypes_bootstrap_results_hg19.csv"))
stopifnot(all(c("annotation", "bp_value", "bias_corrected_estimate", "CI_lo", "CI_hi") %in% names(tumor)))
write.csv(tumor, "data/Tumor_celltypes_bootstrap_results.csv", row.names = FALSE)
write.csv(normal, "data/Normal_celltypes_bootstrap_results.csv", row.names = FALSE)

alpha <- 0.05
sig_t <- tumor$annotation[tumor$bp_value < alpha & tumor$annotation != "Intercept"]
sig_n <- normal$annotation[normal$bp_value < alpha & normal$annotation != "Intercept"]
plot_data <- bind_rows(
  tumor %>% filter(annotation %in% setdiff(sig_t, sig_n)) %>% mutate(side = "Tumor"),
  normal %>% filter(annotation %in% setdiff(sig_n, sig_t)) %>% mutate(side = "Normal")
) %>%
  transmute(side, annotation, bias_corrected_estimate, CI_lo, CI_hi, bp_value, neglog10_p = -log10(pmax(bp_value, 1e-300)))
write.csv(plot_data, "data/Fig1b_plot_data.csv", row.names = FALSE)
cat("Significant in both (not plotted):", paste(intersect(sig_t, sig_n), collapse = "; "), "\n")
print(plot_data)
