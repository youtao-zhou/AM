## Figure5

library(data.table)
library(dplyr)
library(stringr)
library(spatialTIME)
library(tidyr)
library(ggplot2)
library(patchwork)
library(cowplot)
setwd("~/am/tsa/Halo archive 2025-12-08 11-51 - v3.5.3577/ObjectData/")
select <- dplyr::select

files <- list.files(pattern = ".csv$")
dt_list <- lapply(files, function(f) {
  x <- fread(f)
  x[, Source_File := f]  
  x
})

AM2_Total <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
fwrite(AM2_Total, "AM2_Total_Cells_Result.csv")

data <- fread("AM2_Total_Cells_Result.csv")
dt <- as.data.table(data)

dt <- janitor::clean_names(dt)
setnames(dt,
         old = c("analysis_region","object_id","x_min", "x_max", "y_min", "y_max"),
         new = c("Analysis_Region","Object_Id","XMin", "XMax", "YMin", "YMax"),
         skip_absent = TRUE)

dt[, file_name := basename(image_location)]
dt[, deidentified_sample := paste(file_name, Analysis_Region, sep = "_")]
dt[, deidentified_id := str_extract(file_name, "^D\\d+")]
dt[is.na(deidentified_id) | deidentified_id == "", deidentified_id := file_name]

dt[, `:=`(
  x = (XMin + XMax)/2,
  y = (YMin + YMax)/2
)]

dt[, am2_macrophage      := (am2_marophagee == 1)]
dt[, tumor               := (tumor_cells == 1)]
dt[, stromal             := (stromal_cells == 1)]
dt[, cd8a_positive_cells := (cd8a_positive_cells == 1)]  

am2_loc <- dt[am2_macrophage == TRUE, .(
  n_am2       = .N,
  n_in_tumor  = sum(tumor,   na.rm = TRUE),
  n_in_stroma = sum(stromal, na.rm = TRUE)
), by = deidentified_sample][
  ,
  `:=`(
    frac_in_tumor  = n_in_tumor  / pmax(n_am2, 1),
    frac_in_stroma = n_in_stroma / pmax(n_am2, 1)
  )
]

summary(am2_loc$frac_in_tumor)
summary(am2_loc$frac_in_stroma)

r_seq <- seq(0, 200, by = 50)  

spatial_list <- split(
  dt[, .(deidentified_sample, x, y,
         am2_macrophage, tumor, stromal, cd8a_positive_cells,
         XMin, XMax, YMin, YMax, Object_Id)],
  by = "deidentified_sample", drop = TRUE
)

sample_data <- dt[, .(
  `Total Cells`  = .N,
  n_am2          = sum(am2_macrophage, na.rm = TRUE),
  n_tumor        = sum(tumor, na.rm = TRUE),
  n_stromal      = sum(stromal, na.rm = TRUE)
), by = .(deidentified_id, deidentified_sample)]

clinical_data <- unique(dt[, .(deidentified_id, deidentified_sample)])

mif <- create_mif(
  clinical_data = clinical_data,
  sample_data   = sample_data,
  spatial_list  = spatial_list,
  patient_id    = "deidentified_id",
  sample_id     = "deidentified_sample"
)

mif <- spatialTIME::bi_ripleys_k(
  mif,
  mnames = c("tumor","am2_macrophage"),
  r_range = r_seq,
  edge_correction = "translation",
  num_permutations = 0,
  permute = FALSE,
  keep_permutation_distribution = FALSE,
  overwrite = TRUE,
  workers = 30,
  xloc = "x",
  yloc = "y",
  force = TRUE
)

bk <- mif$derived$bivariate_Count

bk_tumor_am2 <- bk %>%
  dplyr::filter(Anchor == "tumor", Counted == "am2_macrophage") %>%
  dplyr::rename(DC = `Degree of Clustering Exact`)

bk_sum <- bk_tumor_am2 %>%
  dplyr::group_by(r) %>%
  dplyr::summarise(
    median_DC = median(DC, na.rm = TRUE),
    q25       = quantile(DC, 0.25, na.rm = TRUE),
    q75       = quantile(DC, 0.75, na.rm = TRUE),
    .groups   = "drop"
  )

ggplot(bk_sum, aes(x = r, y = median_DC)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "red") +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3) +
  geom_line(size = 1) +
  labs(x = "r (µm)", y = "Degree of clustering (AM2 around tumor)") +
  theme_bw()

mif <- spatialTIME::bi_ripleys_k(
  mif,
  mnames = c("am2_macrophage", "cd8a_positive_cells"),
  r_range = r_seq,
  permute = TRUE,
  edge_correction = "translation",
  xloc = "x", yloc = "y",
  workers = 30,
  force = TRUE
)

bk2 <- mif$derived$bivariate_Count

bk_eff_anchor <- bk2 %>%
  dplyr::filter(Anchor == "cd8a_positive_cells",
                Counted == "am2_macrophage") %>%
  dplyr::rename(DC = `Degree of Clustering Permutation`)

bk_eff_sum <- bk_eff_anchor %>%
  dplyr::group_by(r) %>%
  dplyr::summarise(
    median_DC = median(DC, na.rm = TRUE),
    q25       = quantile(DC, 0.25, na.rm = TRUE),
    q75       = quantile(DC, 0.75, na.rm = TRUE),
    .groups   = "drop"
  )

ggplot(bk_eff_sum, aes(x = r, y = median_DC)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "red") +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3) +
  geom_line(size = 1) +
  labs(x = "r (µm)", y = "Degree of clustering (AM2 around CD8A+ T, permutation)") +
  theme_bw()

summary <- as.data.frame(fread("../Total_Summary_Results.csv", sep = ","))

bk_sum_am2_tumor <- bk_sum
bk_sum_am2_effT  <- bk_eff_sum
base_size <- 16

sum_df <- summary

am2_frac_am    <- sum(sum_df$`AM2 marophagee Cells`)   / sum(sum_df$`all macrophage Cells`)
tumor_frac_tot <- sum(sum_df$`tumor cells Cells`)      / sum(sum_df$`Total Cells`)
effT_frac_mean <- mean(sum_df$`% CD8a positive cells Positive Cells`, na.rm = TRUE)
effT_frac_mean <- max(min(effT_frac_mean, 1), 0)

right_bar_df <- tibble::tibble(
  bar   = factor(rep("Eff T / total", 2)),
  group = c("Eff T", "No Eff T"),
  value = c(effT_frac_mean, 1 - effT_frac_mean)
)


left_bar_df <- tibble::tibble(
  bar   = factor(rep(c("AM2 / all AM", "Tumor / total"), each = 2),
                 levels = c("AM2 / all AM", "Tumor / total")),
  group = c("AM2", "No AM2",
            "Tumor", "No Tumor"),
  value = c(am2_frac_am, 1 - am2_frac_am,
            tumor_frac_tot, 1 - tumor_frac_tot)
)



fill_cols <- c(
  "AM2"       = "#D6E6F2",
  "No AM2"    = "#D3D2D2",
  "Tumor"     = "#D7F0D0",
  "No Tumor"  = "#FAC4CC",
  "Eff T"     = "#8E7FB8",
  "No Eff T"  = "#A2C9AE"
)

label_scientific_math <- function(digits = 2) {
  function(x) {
    sapply(x, function(val) {
      if (is.na(val) || is.nan(val)) return(NA)
      if (!is.finite(val)) return(as.character(val))
      if (val == 0) return("0")
      sci   <- format(val, scientific = TRUE, digits = digits)
      parts <- strsplit(sci, "e")[[1]]
      coef  <- parts[1]
      expo  <- as.integer(parts[2])
      coef  <- sub("\\.?0+$", "", coef)
      paste0(coef, "×10^", expo)
    })
  }
}

p_median_am2_tumor <- ggplot(bk_sum_am2_tumor,
                             aes(x = r, y = median_DC)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_ribbon(aes(ymin = q25, ymax = q75),
              fill  = "#E48080",
              alpha = 0.6) +
  geom_line(linewidth = 1.3,
            colour   = "#D20000") +
  scale_y_continuous(labels = label_scientific_math(digits = 2))+
  labs(
    title = "AM2 around Tumor",
    x = "r (µm)",
    y = "Degree of Clustering"
  ) +
  theme_bw(base_size = base_size) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      face  = "plain",
      size  = base_size - 2
    ),
    axis.title = element_text(size = base_size),
    axis.text  = element_text(size = base_size - 2)
  )

p_median_am2_effT <- ggplot(bk_sum_am2_effT,
                            aes(x = r, y = median_DC)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_ribbon(aes(ymin = q25, ymax = q75),
              fill  = "#9A99B7",
              alpha = 0.6) +
  scale_y_continuous(labels = label_scientific_math(digits = 2))+
  geom_line(linewidth = 1.3,
            colour   = "#47373F") +
  labs(
    title = "AM2 around CD8A⁺ T cells",
    x = "r (µm)",
    y = NULL
  ) +
  theme_bw(base_size = base_size) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      face  = "plain",
      size  = base_size - 2
    ),
    axis.title.x = element_text(size = base_size),
    axis.text    = element_text(size = base_size - 2)
  )

p_bar_left_full <- ggplot(left_bar_df,
                          aes(x = bar, y = value, fill = group)) +
  geom_col(width = 0.35) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = fill_cols, name = "") +
  theme_void(base_size = base_size - 4) +
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.key.size   = unit(0.35, "cm"),
    legend.text       = element_text(size = base_size - 6),
    legend.margin     = margin(t = 0, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = -4, r = 0, b = 0, l = 0)
  ) +
  guides(fill = guide_legend(
    nrow = 2,
    byrow = TRUE
  ))

left_legend <- cowplot::get_legend(p_bar_left_full)
p_bar_left_small <- p_bar_left_full + theme(legend.position = "none")

left_inset <- cowplot::plot_grid(
  p_bar_left_small,
  left_legend,
  ncol = 1,
  rel_heights = c(1, 0.3)
)

p_bar_right_full <- ggplot(right_bar_df,
                           aes(x = bar, y = value, fill = group)) +
  geom_col(width = 0.35) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = fill_cols, name = "") +
  theme_void(base_size = base_size - 4) +
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.key.size   = unit(0.35, "cm"),
    legend.text       = element_text(size = base_size - 6),
    legend.margin     = margin(t = 0, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = -4, r = 0, b = 0, l = 0)
  )

right_legend <- cowplot::get_legend(p_bar_right_full)
p_bar_right_small <- p_bar_right_full + theme(legend.position = "none")

right_inset <- cowplot::plot_grid(
  p_bar_right_small,
  right_legend,
  ncol = 1,
  rel_heights = c(1, 0.3)
)

p_left_combined <- p_median_am2_tumor +
  patchwork::inset_element(
    left_inset,
    left   = 0.05,
    bottom = 0.63,
    right  = 0.60,
    top    = 0.98
  )

p_right_combined <- p_median_am2_effT +
  patchwork::inset_element(
    right_inset,
    left   = 0.05,
    bottom = 0.63,
    right  = 0.60,
    top    = 0.98
  )

final_fig <- p_left_combined | p_median_am2_effT

final_fig

x <- plot_immunoflo(
  mif,
  plot_title = "deidentified_sample",
  mnames = c("am2_macrophage", "tumor")
)
find_xy_tables <- function(obj, path="mif", max_depth=6) {
  out <- list()
  walk <- function(x, p, depth) {
    if (depth > max_depth) return()
    if (is.data.frame(x) && all(c("x","y") %in% colnames(x))) {
      out[[p]] <<- x
      return()
    }
    if (is.list(x)) {
      nms <- names(x)
      for (i in seq_along(x)) {
        nm <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else paste0("[[", i, "]]")
        walk(x[[i]], paste0(p, "$", nm), depth + 1)
      }
    }
  }
  walk(obj, path, 1)
  out
}

xy_candidates <- find_xy_tables(mif, "mif", max_depth = 6)
names(xy_candidates)

df <- xy_candidates[[1]]
df <- as.data.frame(df)   

plot_AM_tumor_fixed <- function(df,
                                am_col = "am2_macrophage",   
                                tumor_col = "tumor",
                                stromal_col = "stromal",
                                am_n = Inf,
                                tumor_n = 80000,
                                stromal_n = 20000,
                                seed = 1,
                                title = NULL,
                                alpha_tumor = 0.45,
                                alpha_am = 0.90,
                                alpha_stromal = 0.15,
                                size_tumor = 0.22,
                                size_am = 0.55,
                                size_stromal = 0.14) {
  
  set.seed(seed)
  df <- as.data.frame(df)
  safe_true <- function(v) !is.na(v) & (v == TRUE)

  idx_am <- which(safe_true(df[[am_col]]))
  if (is.finite(am_n) && length(idx_am) > am_n) idx_am <- sample(idx_am, am_n)
  df_am <- df[idx_am, c("x","y"), drop = FALSE]
  df_am$label <- "AM"

  idx_t <- which(safe_true(df[[tumor_col]]))
  if (is.finite(tumor_n) && length(idx_t) > tumor_n) idx_t <- sample(idx_t, tumor_n)
  df_tumor <- df[idx_t, c("x","y"), drop = FALSE]
  df_tumor$label <- "tumor"

  df_stromal <- NULL
  if (!is.null(stromal_col) && stromal_col %in% colnames(df)) {
    idx_s <- which(safe_true(df[[stromal_col]]))
    if (length(idx_s) > 0) {
      if (is.finite(stromal_n) && length(idx_s) > stromal_n) idx_s <- sample(idx_s, stromal_n)
      df_stromal <- df[idx_s, c("x","y"), drop = FALSE]
      df_stromal$label <- "stromal"
    }
  }
  x_min <- min(df$x, na.rm = TRUE); x_max <- max(df$x, na.rm = TRUE)
  y_min <- min(df$y, na.rm = TRUE); y_max <- max(df$y, na.rm = TRUE)
  p <- ggplot() +
    coord_equal(xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand = FALSE) +
    theme_classic(base_size = 12) +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(color = "black")
    ) +
    labs(
      title = "AM-2 Macrophage",
      x = "", y = ""
    )
  
  if (!is.null(df_stromal) && nrow(df_stromal) > 0) {
    p <- p + geom_point(data = df_stromal, aes(x, y, color = label),
                        size = size_stromal, alpha = alpha_stromal)
  }
  p <- p + geom_point(data = df_tumor, aes(x, y, color = label),
                      size = size_tumor, alpha = alpha_tumor)
  p <- p + geom_point(data = df_am, aes(x, y, color = label),
                      size = size_am, alpha = alpha_am)
  p <- p + guides(
    color = guide_legend(
      override.aes = list(size = 1, alpha = 1)
    )
  )
  
  p
}
p_am2 <- plot_AM_tumor_fixed(
  df,
  am_col = "am2_macrophage",
  am_n = 50000,
  tumor_n = 50000,
  stromal_n = 20000,
  title = unique(df$deidentified_sample)[2]
)
p_am2

#########
genes_am2<-c('CD4',"C5AR1","CD8A","IFNG")


genes_use <- intersect(genes_am2, rownames(bulk_dataset))
bulk_am2 <- bulk_dataset[genes_use, icb3$trunc_anonymized_sample_ids, drop = FALSE]
bulk_am2_df <- as.data.frame(t(bulk_am2))
bulk_am2_df$trunc_anonymized_sample_ids <- rownames(bulk_am2_df)

icb3_am2 <- icb3 %>%
  dplyr::left_join(
    bulk_am2_df,
    by = "trunc_anonymized_sample_ids"
  )

genes_use <- intersect(genes_use, colnames(icb3_am2))

d <- icb3_am2 %>%
  dplyr::select(pfs_time, pfs_event, TWAS_z_icb, AM2_z, all_of(genes_use)) %>%
  mutate(
    across(
      where(is.numeric),
      ~ ifelse(is.na(.), median(., na.rm = TRUE), .)
    )
  )

d[genes_use] <- scale(d[genes_use])

X <- model.matrix(~ . - pfs_time - pfs_event - 1, data = d)
y <- Surv(d$pfs_time, d$pfs_event)

set.seed(1)
pf <- rep(1, ncol(X))
pf[match("TWAS_z_icb", colnames(X))] <- 0
pf[match("AM2_z",      colnames(X))] <- 0

cvfit <- cv.glmnet(
  X, y,
  family         = "cox",
  alpha          = 0.5,
  nfolds         = 10,
  penalty.factor = pf
)

sel_min <- colnames(X)[which(as.numeric(coef(cvfit, s = cvfit$lambda.min)) != 0)]
cat("Selected @ lambda.min:\n")
print(sel_min)

extra_genes <- setdiff(sel_min, c("TWAS_z_icb", "AM2_z"))
print(extra_genes)
extra_genes<-c('CD4',"C5AR1","CD8A","IFNG")
form_final <- as.formula(
  paste(
    "Surv(pfs_time, pfs_event) ~ TWAS_z_icb + AM2_z +",
    paste(extra_genes, collapse = " + ")
  )
)
form_final

fit_final <- coxph(form_final, data = icb3_am2)
print(broom::tidy(fit_final, exponentiate = TRUE, conf.int = TRUE))

lp_final <- predict(fit_final, newdata = icb3_am2, type = "lp")
Y_final  <- Surv(icb3_am2$pfs_time, icb3_am2$pfs_event)

times_eval <- c(12, 24)

auc_uno_final <- AUC.uno(
  Surv.rsp     = Y_final,
  Surv.rsp.new = Y_final,
  lpnew        = lp_final,
  times        = times_eval
)

td_auc_final <- data.frame(
  time = times_eval,
  AUC  = round(as.numeric(auc_uno_final$auc), 3)
)
cat("TWAS + AM2 + AM2-signature 的 time-dependent AUC:\n")
print(td_auc_final)

tr_final <- timeROC(
  T      = icb3_am2$pfs_time,
  delta  = icb3_am2$pfs_event,
  marker = lp_final,
  cause  = 1,
  times  = times_eval,
  iid    = TRUE
)

print(data.frame(time = times_eval,
                 AUC = round(tr_final$AUC, 3)))

df_roc_final <- do.call(rbind, lapply(seq_along(times_eval), function(i){
  data.frame(
    FPR   = tr_final$FP[, i],
    TPR   = tr_final$TP[, i],
    time  = times_eval[i],
    label = sprintf("%d month (AUC = %.3f)", times_eval[i], tr_final$AUC[i])
  )
}))

p_roc_final <- ggplot(df_roc_final, aes(x = FPR, y = TPR, color = label)) +
  geom_path(linewidth = 1.2, lineend = "round") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.6, color = "grey60") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous("False Positive Rate", labels = percent_format(accuracy = 1)) +
  scale_y_continuous("True Positive Rate",  labels = percent_format(accuracy = 1)) +
  scale_color_manual(
    values = c("#8fbca5", "#7bbdd3"),
    guide  = guide_legend(title = "",
                          override.aes = list(linewidth = 2))
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position   = "bottom",
    legend.box.margin = margin(t = -10),
    legend.title      = element_text(size = 12),
    legend.text       = element_text(size = 12),
    axis.title        = element_text(size = 12, face = "bold"),
    axis.text         = element_text(size = 12, color = "black"),
    axis.ticks.length = unit(3, "pt"),
    plot.margin       = margin(t = 0, r = 50, l = 0, b = 0)
  ) +
  labs(title = "") +
  annotate(
    "text",
    x = 0.68, y = 0.10, 
    label = paste0(
      "Surv(PFS Time, PFS Event) ~ \n",
      "TWAS Score + AM2 Score\n",
      "+ CD4 + C5AR1 + CD8A + IFNG"
    ),
    hjust = 0.6, vjust = 0.4,
    size  =3
  )
print(p_roc_final)
