#### Figure4
library(BisqueRNA)
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(SingleCellExperiment)
library(ggbeeswarm)
library(stringr)
library(survival)
library(broom)
library(GSVA)
library(BiocParallel)
library(tidyr)
library(rms)
library(survminer)
library(survAUC)
library(timeROC)
library(scales)
library(glmnet)
library(gridExtra)
library(grid)
library(ggsurvfit)
library(progeny)
library(pheatmap)
library(future)
library(CellChat)
library(dior)
library(irlba)
library(ggrepel)

am <- readRDS("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am/Rename_am.rds")
meta <- am@meta.data     
expr_sc <- GetAssayData(am, assay = "RNA", slot = "counts")
meta_sc <- am@meta.data[, c("cell_type", "sample")]
colnames(meta_sc) <- c("cellType_raw", "SubjectName")
meta_sc <- meta_sc[colnames(expr_sc), , drop = FALSE]
cellType_major <- str_extract(meta_sc$cellType_raw, "AM-\\d+")
cellType_major[!cellType_major %in% paste0("AM-", 1:5)] <- NA
keep <- !is.na(cellType_major)
expr_sc_use <- expr_sc[, keep, drop = FALSE]
meta_sc_use <- meta_sc[keep, , drop = FALSE]
cellType_major <- cellType_major[keep]


sce <- SingleCellExperiment(assays = list(counts = as.matrix(expr_sc_use)))
colData(sce)$cellType <- cellType_major                
colData(sce)$SubjectName <- meta_sc_use$SubjectName

marker <- readRDS("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am/marker12.rds")  

cluster_celltype_map <- meta %>%
  dplyr::select(seurat_clusters, cell_type) %>%
  distinct() %>%
  group_by(seurat_clusters) %>%
  summarise(cell_type = dplyr::first(na.omit(cell_type)), .groups = "drop")

cluster_celltype_map$seurat_clusters <- as.character(cluster_celltype_map$seurat_clusters)

marker_with_ct <- marker %>%
  mutate(cluster = as.character(cluster)) %>%
  left_join(cluster_celltype_map,
            by = c("cluster" = "seurat_clusters")) %>%
  mutate(cell_type_major = str_extract(cell_type, "AM-\\d+")) %>%
  filter(!is.na(cell_type_major),
         cell_type_major %in% paste0("AM-", 1:5))

marker_list <- marker_with_ct %>%
  group_by(cell_type_major) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 50) %>%
  summarise(genes = list(unique(gene)), .groups = "drop") %>%
  deframe()

bulk <- readRDS("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/EGA/OAK/OAK.rds")
bulk <- as.matrix(bulk)

common_genes <- intersect(rownames(bulk), rownames(sce))
if (length(common_genes) < 2000) {
  warning("< 2000")
}
bulk_use <- bulk[common_genes, , drop = FALSE]
sce_use  <- sce [common_genes, , drop = FALSE]

Bulk_obj <- ExpressionSet(as.matrix(bulk_use))
sc_obj   <- ExpressionSet(as.matrix(counts(sce_use)))

pd <- data.frame(
  cellType     = as.character(colData(sce_use)$cellType),
  SubjectName  = as.character(colData(sce_use)$SubjectName),
  row.names    = colnames(sce_use)
)
phenoData(sc_obj) <- AnnotatedDataFrame(pd)

res_bisque <- ReferenceBasedDecomposition(
  bulk.eset   = Bulk_obj,
  sc.eset     = sc_obj,
  use.overlap = FALSE,
)

res_bisque<-readRDS("res_bisque_AM5.rds")
prop_df<-as.data.frame(t(res_bisque$bulk.props))
head(prop_df)

data1<-read.csv("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/EGA/OAK/go28915_updated.csv")
data2<-read.csv("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/EGA/OAK/BYN_n699.csv")
data3<-read.csv("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/EGA/OAK/ss_ss.csv")

biomarker<-read.table("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/EGA/OAK/biomarker.txt",header=T,sep="\t")

data4<-left_join(data1,data3,by=c("trunc_anonymized_sample_ids"="alias"))

phenotype<-left_join(data4,biomarker,by=c("trunc_UNI_ID"="trunc_anonymized_patient_id"))

phenotype<-separate(phenotype,col = "full_anonymized_patient_id",into = c("PAT","UNI_ID"))

phenotype<-left_join(phenotype,data2,by="UNI_ID")

phenotype <- phenotype[!is.na(phenotype$ACTARM), ]
phenotype<-phenotype[phenotype$HIST%in%"NON-SQUAMOUS",]
phenotype<-phenotype[phenotype$TISSUE_TYPE.2%in%"lung",]
phenotype<-phenotype[phenotype$BCOR%in%c("CR","PR","SD"),]

prop_df <- as.data.frame(prop_df)
prop_df$sample <- rownames(prop_df)

pheno_sub <- phenotype %>%
  dplyr::select(sample =trunc_anonymized_sample_ids, Group=ACTARM)

data_merged <- left_join(pheno_sub,prop_df, by = "sample")

data_long <- melt(
  data_merged,
  id.vars = c("sample", "Group"),
  variable.name = "Cell_Type",
  value.name = "Proportion"
)

data_long <- data_long %>%
  group_by(Cell_Type) %>%
  filter(mean(Proportion) > 0.001) %>%
  ungroup()

data_long<-subset(data_long,data_long$Cell_Type%in%"AM-4")
data_long2<-readRDS("data_long_AM2.rds")
data_long<-rbind(data_long,data_long2)

data_long <- data_long %>%
  mutate(Cell_Type = case_when(
    Cell_Type %in% c("AM-2.CK", "AM-2.CD4") ~ "AM-2",
    TRUE ~ Cell_Type
  ))

pal_fill <- c("Docetaxel"="#84A27C", "MPDL3280A"="#E17357") 
pal_line <- c("Docetaxel"="#84A27C", "MPDL3280A"="#E17357") 
dodge_w <- 0.55

p <- ggplot(data_long, aes(x = Cell_Type, y = Proportion, fill = Group)) +
  ggbeeswarm::geom_quasirandom(
    aes(color = Group),
    dodge.width = dodge_w, width = 0.16, varwidth = FALSE,
    size = 1.4, alpha = 0.55, stroke = 0
  ) +
  geom_boxplot(
    width = 0.52, outlier.shape = NA, linewidth = 0.6,
    color = "black"
  ) +
  stat_summary(
    fun = median, geom = "point",
    shape = 21, size = 2.6, stroke = 0.7, color = "black",
    position = position_dodge(width = dodge_w), fill = "white"
  ) +
  scale_fill_manual(values = pal_fill) +
  scale_color_manual(values = pal_line, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.07))) +
  labs(
    title = NULL,
    x = NULL,
    y = "Estimated by BisqueRNA"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y  = element_text(face = "bold", size = 18, margin = margin(r = 10)),
    axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 15),
    axis.text.y   = element_text(face = "bold", size = 15),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#E6E6E6", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "#F3F3F3", linewidth = 0.4),
    axis.line = element_line(colour = "black", linewidth = 0.6),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"),
    legend.key.size = unit(16, "pt"),   
    legend.key.width  = unit(18, "pt"),
    legend.key.height = unit(14, "pt"),
    plot.margin = margin(6, 6, 6, 6)
  )+
  stat_compare_means(
    aes(x = Cell_Type, y = Proportion, group = Group),
    method = "t.test",
    label = "p.signif",
    hide.ns = TRUE,
    label.y = 0.65,   
    size = 6
  )
p

############
TWAS_signature_score_OAK <- readRDS("FTWAS_signature_score_OAK.rds")
res_bisque               <- readRDS("res_bisque_AM5.rds") 
data1     <- read.csv("go28915_updated.csv")
data2     <- read.csv("BYN_n699.csv")
data3     <- read.csv("ss_ss.csv")
biomarker <- read.table("biomarker.txt",header = TRUE, sep = "\t")
data4 <- dplyr::left_join(data1, data3,
                          by = c("trunc_anonymized_sample_ids" = "alias"))
phenotype <- dplyr::left_join(data4, biomarker, by = c("trunc_UNI_ID" = "trunc_anonymized_patient_id"))
phenotype <- tidyr::separate(
  phenotype,
  col    = "full_anonymized_patient_id",
  into   = c("PAT","UNI_ID"),
  remove = FALSE
)
phenotype <- dplyr::left_join(phenotype, data2, by = "UNI_ID")
phenotype <- phenotype[phenotype$HIST %in% "NON-SQUAMOUS", ]
phenotype <- phenotype[phenotype$TISSUE_TYPE.2 %in% "lung", ]
phenotype <- phenotype[phenotype$BCOR %in% c("CR","PR","SD"), ]
bulk_dataset <- readRDS("OAK.rds")
colnames(bulk_dataset) <- gsub("\\.", "-", colnames(bulk_dataset)) 
common_samps <- intersect(colnames(bulk_dataset),
                          phenotype$trunc_anonymized_sample_ids)
bulk_use  <- bulk_dataset[, common_samps, drop = FALSE]
pheno_use <- phenotype[match(common_samps,
                             phenotype$trunc_anonymized_sample_ids),
                       , drop = FALSE]
stopifnot(identical(colnames(bulk_use),
                    pheno_use$trunc_anonymized_sample_ids))
bulk_use_log <- if (max(bulk_use, na.rm = TRUE) > 100) log2(bulk_use + 1) else as.matrix(bulk_use)
bulk_use_log <- as.matrix(bulk_use_log)

prop_df <- as.data.frame(t(res_bisque$bulk.props))
prop_df$trunc_anonymized_sample_ids <- rownames(prop_df)

bisque_am2 <- prop_df %>%
  transmute(
    trunc_anonymized_sample_ids,
    AM2_bisque = `AM-2`
  )

dat <- pheno_use %>%
  left_join(
    TWAS_signature_score_OAK %>%
      select(trunc_anonymized_sample_ids, TWAS_signature1),
    by = "trunc_anonymized_sample_ids"
  ) %>%
  left_join(
    bisque_am2,
    by = "trunc_anonymized_sample_ids"
  ) %>%
  mutate(
    treatment = case_when(
      str_detect(toupper(as.character(ACTARM)),
                 "MPDL|ATEZ|PD-?L1|IMMUNO|ICB") ~ "ICB",
      TRUE ~ "Chemo"
    ),
    treatment = factor(treatment, levels = c("Chemo","ICB")),
    os_time  = OS_MONTHS,   os_event  = OS_CENSOR,
    pfs_time = PFS_MONTHS,  pfs_event = PFS_CENSOR,
    age      = AGE
  )

icb <- dat %>%
  filter(treatment == "ICB") %>%
  mutate(
    TWAS_z_icb   = as.numeric(scale(TWAS_signature1)),
    AM2_bisque_z = as.numeric(scale(AM2_bisque)),
    AM2_z        = AM2_bisque_z        
  )

fit_bisque_am2 <- coxph(
  Surv(pfs_time, pfs_event) ~ TWAS_z_icb + AM2_bisque_z,
  data = icb
)
print(broom::tidy(fit_bisque_am2, exponentiate = TRUE, conf.int = TRUE))

fit_bisque_am2_age <- coxph(
  Surv(pfs_time, pfs_event) ~ TWAS_z_icb + AM2_bisque_z + age,
  data = icb
)
print(broom::tidy(fit_bisque_am2_age, exponentiate = TRUE, conf.int = TRUE))

icb_dd <- icb %>% select(TWAS_z_icb, age, pfs_time, pfs_event)
dd <- datadist(icb_dd); options(datadist="dd")

fit_rcs <- cph(
  Surv(pfs_time, pfs_event) ~ rcs(TWAS_z_icb, 4) + age,
  data = icb, x = TRUE, y = TRUE, surv = TRUE
)

z_grid <- seq(min(icb$TWAS_z_icb, na.rm=TRUE),
              max(icb$TWAS_z_icb, na.rm=TRUE),
              length.out = 801)

newd   <- data.frame(
  TWAS_z_icb = z_grid,
  age        = median(icb$age, na.rm = TRUE)
)

lp_rcs <- as.numeric(predict(fit_rcs, newdata = newd, type = "lp"))

i_min  <- which.min(lp_rcs)
hr_rel <- exp(lp_rcs - lp_rcs[i_min])
eps    <- 0.02
idx    <- i_min + which(hr_rel[(i_min+1):length(hr_rel)] >= (1 + eps))[1]

ct_twas_icb <- uniroot(
  approxfun(z_grid, hr_rel - (1 + eps)),
  lower = z_grid[idx-1], upper = z_grid[idx]
)$root

message(sprintf("TWAS切点(ICB, RCS) = %.4f", ct_twas_icb))

icb2 <- icb %>%
  mutate(TWAS_bin = ifelse(TWAS_z_icb >= ct_twas_icb, "High","Low"))

am2_med_hi <- median(icb2$AM2_z[icb2$TWAS_bin == "High"], na.rm = TRUE)

icb3 <- icb2 %>%
  mutate(
    AM2_bin = ifelse(
      AM2_z >= am2_med_hi, "High",
      ifelse(!is.na(AM2_z), "Low", NA_character_)
    ),
    risk4 = dplyr::case_when(
      TWAS_bin == "Low"  & AM2_bin == "High" ~ "Best",
      TWAS_bin == "High" & AM2_bin == "Low"  ~ "Worst",
      TWAS_bin == "High" & AM2_bin == "High" ~ "TWAShigh_AM2high",
      TWAS_bin == "Low"  & AM2_bin == "Low"  ~ "TWASlow_AM2low",
      TRUE                                   ~ NA_character_
    ),
    risk4 = factor(
      risk4,
      levels = c("Best","TWAShigh_AM2high","TWASlow_AM2low","Worst")
    )
  )

table(icb3$TWAS_bin, icb3$AM2_bin, useNA = "ifany")
print(table(icb3$risk4, useNA = "ifany"))

icb3$risk3 <- case_when(
  icb3$risk4 == "Both High"        ~ "Intermediate",
  icb3$risk4 == "TWAShigh_AM2high" ~ "Intermediate",
  icb3$risk4 == "TWASlow_AM2low"   ~ "Intermediate",
  icb3$risk4 == "Best"             ~ "Best",
  icb3$risk4 == "Worst"            ~ "Worst"
)

icb3$risk3 <- factor(icb3$risk3, levels = c("Best","Intermediate","Worst"))

fit3 <- survfit(Surv(pfs_time, pfs_event) ~ risk3, data = icb3)
sf3  <- survfit2(Surv(pfs_time, pfs_event) ~ risk3, data = icb3)

max_time <- max(icb3$pfs_time, na.rm = TRUE)

label_df <- data.frame(
  risk3  = levels(icb3$risk3),
  y      = c(0.85, 0.50, 0.20),
  x      = max_time * 0.97
)

pairwise_p <- function(groups) {
  dsub <- icb3[icb3$risk3 %in% groups, ]
  dsub$risk3_pair <- droplevels(dsub$risk3)
  sdiff <- survdiff(Surv(pfs_time, pfs_event) ~ risk3_pair, data = dsub)
  1 - pchisq(sdiff$chisq, df = length(sdiff$n) - 1)
}

p_worst_best   <- pairwise_p(c("Worst","Best"))
p_worst_int    <- pairwise_p(c("Worst","Intermediate"))
p_best_int     <- pairwise_p(c("Best","Intermediate"))

p_tbl <- data.frame(
  Comparison = c(
    "Worst vs Best",
    "Worst vs Intermediate",
    "Best vs Intermediate"
  ),
  `P value` = sprintf("%.3f", c(p_worst_best,p_worst_int,p_best_int))
)

tbl_grob <- tableGrob(
  p_tbl,
  rows = NULL,
  theme = ttheme_minimal(base_size = 9)
)

p3 <- sf3 %>%
  ggsurvfit(linewidth = 0.9) +
  add_censor_mark(size = 3, shape = 73) +
  add_pvalue(
    location = "annotation",
    x = max_time * 0.01,   
    y = 0.10,              
    hjust = 0,             
    size = 4.5,
    caption = "Log-rank {p.value}"
  ) +
  labs(
    title = NULL,
    x = "Months",
    y = "Progression-free survival"
  ) +
  scale_x_continuous(breaks = seq(0, max_time, by = 6), expand = c(0.03,0)) +
  scale_y_continuous(breaks = seq(0,1,0.25), labels = seq(0,100,25)) +
  scale_color_manual(values = c("#5BB5AC", "#D8B365", "#DE526C")) +
  scale_fill_manual(values  =c("#5BB5AC", "#D8B365", "#DE526C")) +
  theme_classic() +
  theme(
    axis.text          = element_text(size=14,color="black"),
    axis.title.x       = element_text(size=14),
    axis.title.y       = element_text(size=14),
    panel.grid.major.y = element_line(color="#DCDDDF"),
    legend.position    = "none"
  ) +
  geom_text(
    data = label_df,
    aes(x=x, y=y, label=risk3, color=risk3),
    size = 5, hjust=1.2, vjust=4.5, show.legend=FALSE
  ) +
  annotation_custom(
    grob = tbl_grob,
    xmin = max_time * 0.75,
    xmax = max_time * 0.98,
    ymin = 0.40,
    ymax = 1.4
  )

print(p3)

fit_g3 <- coxph(
  Surv(pfs_time, pfs_event) ~ TWAS_z_icb + AM2_z,
  data = icb3
)
print(broom::tidy(fit_g3, exponentiate = TRUE, conf.int = TRUE))
df <- broom::tidy(fit_g3, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%  
  mutate(
    Variable = dplyr::recode(term,
                             "TWAS_z_icb" = "TWAS score",
                             "AM2_z"      = "AM2 score"),
    HR       = estimate,
    lower_CI = conf.low,
    upper_CI = conf.high,
    HR_label = sprintf("%.2f (%.2f–%.2f)", HR, lower_CI, upper_CI),
    P_label  = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  )

df$Variable <- factor(df$Variable, levels = rev(unique(df$Variable)))
df <- df %>% arrange(Variable) %>%
  mutate(row_id = as.numeric(Variable))
pal <- c("#5BB5AC", "#DE526C")
col_map <- setNames(pal[seq_along(levels(df$Variable))], levels(df$Variable))

x_var     <- -1
x_hr_text <- -0.3   
x_p_text  <- 2.3  
x_forestL <- 1.0   
x_forestR <- 2.1   

bg <- df %>%
  filter(row_id %% 2 == 0) %>%
  transmute(ymin = row_id - 0.5, ymax = row_id + 0.5)

p_forest <- ggplot(df, aes(y = Variable)) +
  geom_rect(data = bg, inherit.aes = FALSE,
            aes(ymin = ymin, ymax = ymax),
            xmin = x_var - 0.05, xmax = x_forestR + 0.2, fill = "grey95") +
  geom_text(aes(x = x_var, label = Variable),
            hjust = 0, size = 4.5) +
  geom_text(aes(x = x_hr_text, label = HR_label),
            hjust = 0, size = 4.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI, color = Variable),
                 height = 0.20, linewidth = 1) +
  geom_point(aes(x = HR, color = Variable), size = 3.0) +
  geom_text(aes(x = x_p_text, label = P_label),
            hjust = 0, size = 4.5, fontface = "italic") +
  
  scale_color_manual(values = col_map, guide = "none") +
  scale_x_continuous(
    limits = c(-1, 2.5),
    breaks = c(-1, 0, 1, 2.0),
    labels = c("-1", "0", "1.0", "2.0")
  ) +
  
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.text.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    axis.text.x   = element_text(color = "black", face = "bold", size = 13),
    panel.grid    = element_blank(),
    panel.border  = element_blank(),
    axis.line.x   = element_line(color = "black"),
    plot.margin   = margin(10, 10, 10, 10)
  ) +
  annotate("text", x = x_var,     y = max(df$row_id) + 1,
           label = "Variable", hjust = 0, size = 5, fontface = "bold") +
  annotate("text", x = x_hr_text, y = max(df$row_id) + 1,
           label = "HR (95% CI)", hjust = 0, size = 5, fontface = "bold") +
  annotate("text", x = x_p_text,  y = max(df$row_id) + 1,
           label = "P value", hjust = 0, size = 5, fontface = "bold")

p_forest



###############
seurat_object<-readRDS("Rename_am.rds")
metadata<-seurat_object@meta.data

model <- progeny::model_human_full
head(model)

DefaultAssay(seurat_object) <- "RNA"
expr <- GetAssayData(seurat_object, slot = "data") 

genes <- rownames(expr)
dup <- duplicated(genes)
if (any(dup)) {
  expr <- expr[!dup, ]
}
expr_dense <- as.matrix(expr)  
pathway_scores <- progeny(
  expr      = expr_dense,      
  scale     = TRUE,
  organism  = "Human",
  top       = 500,
  z_scores  = TRUE,
  verbose   = TRUE,
  model     = model_use
)


dim(pathway_scores); pathway_scores[1:3,1:3]

all.equal(rownames(metadata), rownames(pathway_scores)) 

pathway_meta <- cbind(metadata, pathway_scores)

pathway_meta <- pathway_meta %>%
  mutate(AM_group = str_extract(cell_type, "^AM-[0-9]+"))
head(pathway_meta[, c("cell_type", "AM_group")])

pathway_cluster_mean <- pathway_meta %>%
  group_by(AM_group) %>%
  summarise(across(Androgen:WNT, mean, na.rm = TRUE)) %>%   
  as.data.frame()

rownames(pathway_cluster_mean) <- pathway_cluster_mean$AM_group
pathway_cluster_mean <- pathway_cluster_mean[ , -1]  
head(pathway_cluster_mean)


pw_am2_immune <- c("NFkB","JAK-STAT","TNFa")
pw_am35_rtk <- c("EGFR","MAPK","PI3K","p53")
pw_am4_repair <- c("TGFb","VEGF","Hypoxia","WNT")
pw_am1_apoptosis <- c("Trail","Androgen","Estrogen")

pw_order <- c(pw_am1_apoptosis, pw_am2_immune, pw_am35_rtk, pw_am4_repair)
pw_order <- pw_order[pw_order %in% colnames(pathway_cluster_mean)]


ppathway_cluster_mean <- pathway_cluster_mean[, pw_order, drop=FALSE]
gaps <- c(length(pw_am1_apoptosis),
          length(pw_am1_apoptosis)+length(pw_am2_immune),
          length(pw_am1_apoptosis)+length(pw_am2_immune)+length(pw_am35_rtk))

library(pheatmap)
pheatmap(pathway_cluster_mean,
         cluster_rows = FALSE, cluster_cols = FALSE,
         fontsize = 15,
         scale = "row",
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
         gaps_col = gaps,
         border_color = NA,
         main = "PROGENy")


##################
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 2 * 1024^3)

groupSize <- readRDS("F:/胸外/twaslc/am/cellchat/groupSize.rds")
cellchat  <- readRDS("F:/胸外/twaslc/am/cellchat/cellchat_am_tumor_computeCommunProb.rds")

if (is.null(cellchat@netP$pathways)) cellchat <- computeCommunProbPathway(cellchat)
if (is.null(cellchat@net$weight))    cellchat <- aggregateNet(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

FOCUS <- "AM-2.CK"
TopK_ck <- 30   
EPS <- 1e-6
topN  <- 15
COMPARE_MODE <- "other_q90"   
Q <- 0.9                     
POS_THRESH <- 0

palette_all <- c(
  "Cytokine_Chemokine"   = "#E684AE",
  "Antigen_presentation" = "#92C5DE",
  "Immune_CoSignaling"   = "#7FBF7B",
  "Inflammatory_DAMP"    = "#C9C7E3",
  "Phagocytosis_Signals" = "#FDB863",
  "Other"                = "grey80"
)

.norm_pair <- function(x) toupper(gsub("[^A-Z0-9]+", "", x))

module_map <- function(x) {
  n <- .norm_pair(x)
  dplyr::case_when(
    grepl("(IL[0-9]+|TNF|IFN|CXCL[0-9]+|CCL[0-9]+|CX3CL1|CSF[0-9]*|TGFB[0-9]*)", n) ~ "Cytokine_Chemokine",
    grepl("(HLA|B2M|CD74)", n) ~ "Antigen_presentation",
    grepl(paste0(
      "(PDCD1|CD274|PDCD1LG2|CTLA4|CD80|CD86|ICOS|ICOSLG|CD40|CD40LG|",
      "CD27|CD70|TNFRSF[0-9]+|TNFSF[0-9]+|TNFRSF4|TNFSF4|TNFRSF9|TNFSF9|",
      "TNFRSF18|TNFSF18|TNFRSF14|TNFSF14)"
    ), n) ~ "Immune_CoSignaling",
    grepl("(MIF|ANXA1|S100A8|S100A9|RETN|HMGB1)", n) ~ "Inflammatory_DAMP",
    grepl("(CD47|SIRPA|CALR|LRP1|MERTK|AXL|TYRO3|GAS|PROS1|MARCO)", n) ~ "Phagocytosis_Signals",
    TRUE ~ "Other"
  )
}

comm_out <- subsetCommunication(cellchat, sources.use = FOCUS, targets.use = NULL) %>%
  mutate(direction = "Outgoing", partner = target)
comm_in  <- subsetCommunication(cellchat, sources.use = NULL, targets.use = FOCUS) %>%
  mutate(direction = "Incoming", partner = source)
comm_all <- bind_rows(comm_out, comm_in) %>% filter(!is.na(partner))

pair_summ <- comm_all %>%
  group_by(direction, partner, interaction_name_2) %>%
  summarise(mean_prob = mean(prob, na.rm = TRUE), .groups = "drop") %>%
  mutate(module = module_map(interaction_name_2))

keep_ints_ck <- pair_summ %>%
  filter(module == "Cytokine_Chemokine") %>%
  group_by(interaction_name_2) %>%
  summarise(mp_overall = mean(mean_prob), .groups = "drop") %>%
  arrange(desc(mp_overall)) %>%
  slice_head(n = TopK_ck) %>%
  pull(interaction_name_2)

pair_ck <- pair_summ %>%
  filter(module == "Cytokine_Chemokine",
         interaction_name_2 %in% keep_ints_ck)

partner_order_ck <- pair_ck %>%
  group_by(direction, partner) %>%
  summarise(tot = sum(mean_prob), .groups = "drop") %>%
  group_by(direction) %>%
  arrange(desc(tot), .by_group = TRUE) %>%
  mutate(order = row_number()) %>%
  ungroup()

pair_ck <- pair_ck %>%
  left_join(partner_order_ck, by = c("direction","partner")) %>%
  mutate(
    partner = reorder(partner, order),
    interaction_name_2 = reorder(interaction_name_2, mean_prob),
    module = factor("Cytokine_Chemokine", levels = "Cytokine_Chemokine") 
  )

p_ck <- ggplot(pair_ck, aes(x = partner, y = interaction_name_2)) +
  geom_point(aes(size = mean_prob, fill = module),
             shape = 21, color = "black", stroke = 0.25) +
  scale_size_continuous(name = "Mean prob", range = c(2.2, 9)) +
  scale_fill_manual(values = c(Cytokine_Chemokine = "#E684AE")) +
  facet_grid(rows = vars(direction), switch = "y", scales = "free_y") +
  labs(
    x = NULL, y = NULL,
    title = paste0(FOCUS, " Communications")
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 40, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

p_ck


out_focus <- subsetCommunication(cellchat, sources.use = FOCUS, targets.use = NULL) %>%
  group_by(target, interaction_name_2) %>%
  summarise(focus_mean = mean(prob, na.rm = TRUE), .groups = "drop") %>%
  rename(partner = target)

out_others_all <- subsetCommunication(cellchat, sources.use = NULL, targets.use = NULL) %>%
  filter(source != FOCUS)

out_others_full <- out_others_all %>%
  group_by(target, interaction_name_2) %>%
  summarise(other_mean = mean(prob, na.rm = TRUE),
            other_max  = max(prob,  na.rm = TRUE),
            other_q90  = quantile(prob, probs = Q, na.rm = TRUE, names = FALSE),
            other_sd   = sd(prob, na.rm = TRUE), .groups = "drop") %>%
  rename(partner = target) %>%
  mutate(across(c(other_mean, other_max, other_q90, other_sd), ~ tidyr::replace_na(., 0)))

uniq_out <- out_focus %>%
  left_join(out_others_full, by = c("partner","interaction_name_2")) %>%
  mutate(delta_max  = focus_mean - other_max,
         delta_mean = focus_mean - other_mean,
         delta_q90  = focus_mean - other_q90,
         log2FC     = log2((focus_mean + EPS)/(other_mean + EPS)),
         zscore     = (focus_mean - other_mean)/pmax(other_sd, EPS),
         spec       = pmax(0, focus_mean - other_mean)/(focus_mean + other_mean + EPS) * log2FC,
         module     = module_map(interaction_name_2)) %>%
  mutate(metric = dplyr::case_when(
    COMPARE_MODE == "other_max"  ~ delta_max,
    COMPARE_MODE == "other_mean" ~ delta_mean,
    COMPARE_MODE == "other_q90"  ~ delta_q90,
    TRUE ~ delta_q90
  ))

sel_out <- uniq_out %>% filter(is.finite(metric), metric > POS_THRESH)
if (nrow(sel_out) == 0) sel_out <- uniq_out %>% filter(is.finite(log2FC), log2FC > 0)
if (nrow(sel_out) == 0) sel_out <- uniq_out

uniq_out_top2 <- sel_out %>%
  arrange(desc(metric), desc(log2FC), desc(focus_mean)) %>%
  slice_head(n = topN) %>%
  mutate(key = paste(partner, interaction_name_2, sep = "\n"),
         key = factor(key, levels = rev(unique(key))),
         module = factor(module))

mods_out <- levels(droplevels(uniq_out_top2$module))
pal_out  <- palette_all[mods_out]

ylab_out <- dplyr::case_when(
  COMPARE_MODE == "other_max"  ~ "Δ vs best other",
  COMPARE_MODE == "other_mean" ~ "Δ vs mean other",
  TRUE                         ~ paste0("Δ vs ", Q*100, "th percentile")
)

p_out <- ggplot(uniq_out_top2, aes(x = key, y = metric, fill = module)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = pal_out, breaks = mods_out, limits = mods_out, drop = TRUE) +
  labs(x = NULL, y = "Signal difference",
       title = paste0(FOCUS, " unique outgoing pairs")) +
  theme_classic(base_size = 13) +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        legend.position = "bottom")
p_out

in_focus <- subsetCommunication(cellchat, sources.use = NULL, targets.use = FOCUS) %>%
  group_by(source, interaction_name_2) %>%
  summarise(focus_mean = mean(prob, na.rm = TRUE), .groups = "drop") %>%
  rename(partner = source)

in_others_all <- subsetCommunication(cellchat, sources.use = NULL, targets.use = NULL) %>%
  filter(target != FOCUS)

in_others_full <- in_others_all %>%
  group_by(source, interaction_name_2) %>%
  summarise(other_mean = mean(prob, na.rm = TRUE),
            other_max  = max(prob,  na.rm = TRUE),
            other_q90  = quantile(prob, probs = Q, na.rm = TRUE, names = FALSE),
            other_sd   = sd(prob, na.rm = TRUE), .groups = "drop") %>%
  rename(partner = source) %>%
  mutate(across(c(other_mean, other_max, other_q90, other_sd), ~ tidyr::replace_na(., 0)))

uniq_in <- in_focus %>%
  left_join(in_others_full, by = c("partner","interaction_name_2")) %>%
  mutate(delta_max  = focus_mean - other_max,
         delta_mean = focus_mean - other_mean,
         delta_q90  = focus_mean - other_q90,
         log2FC     = log2((focus_mean + EPS)/(other_mean + EPS)),
         zscore     = (focus_mean - other_mean)/pmax(other_sd, EPS),
         spec       = pmax(0, focus_mean - other_mean)/(focus_mean + other_mean + EPS) * log2FC,
         module     = module_map(interaction_name_2)) %>%
  mutate(metric = dplyr::case_when(
    COMPARE_MODE == "other_max"  ~ delta_max,
    COMPARE_MODE == "other_mean" ~ delta_mean,
    COMPARE_MODE == "other_q90"  ~ delta_q90,
    TRUE ~ delta_q90
  ))

sel_in <- uniq_in %>% filter(is.finite(metric), metric > POS_THRESH)
if (nrow(sel_in) == 0) sel_in <- uniq_in %>% filter(is.finite(log2FC), log2FC > 0)
if (nrow(sel_in) == 0) sel_in <- uniq_in

uniq_in_top2 <- sel_in %>%
  arrange(desc(metric), desc(log2FC), desc(focus_mean)) %>%
  slice_head(n = topN) %>%
  mutate(key = paste(partner, interaction_name_2, sep = "\n"),
         key = factor(key, levels = rev(unique(key))),
         module = factor(module))

mods_in <- levels(droplevels(uniq_in_top2$module))
pal_in  <- palette_all[mods_in]

p_in <- ggplot(uniq_in_top2, aes(x = key, y = metric, fill = module)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = pal_in, breaks = mods_in, limits = mods_in, drop = TRUE) +
  labs(x = NULL, y = "Signal difference",
       title = paste0(FOCUS, " unique incoming pairs")) +
  theme_classic(base_size = 13) +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        legend.position = "bottom")
p_in

library(patchwork)

combined_plot <- (wrap_elements(p_out) / wrap_elements(p_in)) +
  plot_annotation(title = "") +
  plot_layout(heights = c(1, 1))
combined_plot

#############
core_tfs <- c("KLF10(+)","FOXO3(+)","NFIL3(+)","FOXO1(+)")
cols_AM2 <- c("AM-2_Normal","AM-2_Tumor")

core_tfs <- intersect(core_tfs, rownames(pair_scaled))
cols_AM2 <- intersect(cols_AM2, colnames(pair_scaled))

mat_core <- pair_scaled[core_tfs, cols_AM2, drop = FALSE]
mat_core <- as.matrix(mat_core)
storage.mode(mat_core) <- "numeric"

lim  <- max(abs(mat_core), na.rm = TRUE)
col_fun <- circlize::colorRamp2(
  c(0, 1, 2.195),
  c("#7E9BB7","white","#D6AFB9")
)
colnames(mat_core) <- sub("_", " ", colnames(mat_core))

ht1 <- Heatmap(
  mat_core,
  name = "AM2_core",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "",
  row_title = "",
  column_names_gp = gpar(fontsize = 12),
  row_names_gp = gpar(fontsize = 12),
  column_names_rot = 0,
  heatmap_legend_param = list(
    title  = "Activity",
    at     = c(-lim, 0, lim),
    labels = c("low", "mid", "high")
  ),
  height = unit(3, "cm")
)

ht1


regulons_AM2_CD4 <- c("ELF4(+)","EP300(+)","RXRA(+)","SMARCA1(+)","ZBTB7A(+)")
regulons_AM2_CD4 <- intersect(regulons_AM2_CD4, rownames(auc_tumor))

auc_tumor_mat <- getAUC(auc_tumor)

group_def <- list(
  AM1    = "^AM-1",
  AM2.CK = "^AM-2.*CK",
  AM2.CD4= "^AM-2.*CD4",
  AM3    = "^AM-3",
  AM4    = "^AM-4",
  AM5    = "^AM-5"
)

mat_list <- list()
for(g in names(group_def)){
  pat <- group_def[[g]]
  cells_g <- rownames(anno_tumor)[grepl(pat, anno_tumor$celltype)]
  cells_g <- intersect(cells_g, colnames(auc_tumor_mat))
  if(length(cells_g) == 0) next 
  mat_list[[g]] <- rowMeans(
    auc_tumor_mat[regulons_AM2_CD4, cells_g, drop = FALSE],
    na.rm = TRUE
  )
}

mat_comb_raw <- do.call(cbind, mat_list)
mat_comb_raw <- as.matrix(mat_comb_raw)
storage.mode(mat_comb_raw) <- "numeric"

wanted_order <- c("AM1","AM2.CK","AM2.CD4","AM3","AM4","AM5")
wanted_order <- intersect(wanted_order, colnames(mat_comb))
mat_comb <- mat_comb[, wanted_order, drop = FALSE]

ha_bottom <- HeatmapAnnotation(
  AM = factor(colnames(mat_comb), levels = wanted_order),
  col = list(
    AM = setNames(
      c("#adc9ca","#a7a5d7","#9ED1E2","#ffac5f","#7aa5d2","#e49898"),
      wanted_order
    )
  )
)

col_fun2 <- circlize::colorRamp2(
  c(-0.6, 0, 1),
  c("#7E9BB7","white","#D6AFB9")
)
ha_bottom <- HeatmapAnnotation(
  AM = factor(colnames(mat_comb), levels = wanted_order),
  col = list(
    AM = setNames(
      c("#adc9ca","#a7a5d7","#9ED1E2","#ffac5f","#7aa5d2","#e49898"),
      wanted_order
    )
  ),
  annotation_legend_param = list(
    AM = list(
      title = "AM",
      direction = "vertical" 
    )
  )
)
ht2 <- Heatmap(
  mat_comb,
  name = "AM2CD4_prog",
  col = col_fun2,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  top_annotation = ha_bottom,
  column_title = "",
  row_title = "",
  column_names_gp = gpar(fontsize = 12),
  row_names_gp = gpar(fontsize = 12),
  column_names_rot = 0,
  heatmap_legend_param = list(
    title  = "Activity",
    at     = c(-lim2, 0, lim2),
    labels = c("low", "mid", "high"),
    legend_title_gp  = gpar(fontsize = 14),
    legend_labels_gp = gpar(fontsize = 14) 
  )
)
draw(
  ht2,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE  
)


############
tcell=read_h5("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/cellhint/tcell_100k_symbol.h5")

DefaultAssay(tcell) <- "RNA"
Idents(tcell) <- tcell$cell_type

p1 <- DimPlot(
  tcell,
  reduction = "umap",
  group.by  = "cell_type",
  label = TRUE,
  label.size = 4,
  pt.size = 0.2
) + theme_classic()

p1

tcell <- FindVariableFeatures(
  tcell,
  selection.method = "vst",
  nfeatures = 3000
)

tcell <- ScaleData(
  tcell,
  features = VariableFeatures(tcell),
  verbose = FALSE
)

tcell <- RunPCA(
  tcell,
  features = VariableFeatures(tcell),
  npcs = 50,
  verbose = FALSE
)
tcell <- RunUMAP(tcell, dims = 1:30)
tcell <- FindNeighbors(tcell, dims = 1:30)
tcell <- FindClusters(tcell, resolution = 0.2)
p_new_ct <- DimPlot(
  tcell,
  reduction = "umap",
  label = TRUE,
  label.size = 4,
  pt.size = 0.2
) + theme_classic()

p_new_ct

all_markers<-readRDS("tcell_markers13.rds")

top10 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

top10 %>% arrange(cluster, desc(avg_log2FC)) %>% head(60)

cluster_lists <- top10 %>%
  group_by(cluster) %>%
  summarise(
    genes = paste(gene, collapse = ", ")
  )


tcell$sub_type <- "Other"
cluster2type <- c(
  "0" = "CD4_Tfh",
  "1" = "CD4_Tfh",
  "2" = "CD8_Eff",
  "5" = "CD4_Tcm",
  "7" = "Treg"
)
current_clusters <- as.character(Idents(tcell))
hit <- current_clusters %in% names(cluster2type)
tcell$sub_type[hit] <- cluster2type[current_clusters[hit]]

tcell$sub_type <- factor(
  tcell$sub_type,
  levels = c("CD4_Tfh", "CD4_Tcm", "Treg", "CD8_Eff", "Other")
)
tcell$cell_type <- tcell$sub_type

tcell_palette <- c(
  "CD4_Tfh" = "#EE706F",
  "CD4_Tcm" = "#FEE09E",
  "Treg"    = "#4C938C",
  "CD8_Eff" = "#9ED1E2",
  "Other"   = "#D4E3D0"
)

emb_t <- as.data.frame(Embeddings(tcell, "umap"))
colnames(emb_t)[1:2] <- c("UMAP_1","UMAP_2")

df_t <- emb_t %>%
  mutate(
    ident    = as.character(Idents(tcell)),
    sub_type = as.character(tcell$sub_type)
  )

ellipse_groups <- c("CD4_Tfh","CD4_Tcm","Treg","CD8_Eff")
df_ell <- df_t %>% filter(sub_type %in% ellipse_groups)

label_pos <- df_ell %>%
  group_by(sub_type) %>%
  summarise(
    cx = mean(UMAP_1),
    cy = mean(UMAP_2),
    .groups = "drop"
  )
xr <- range(df_t$UMAP_1, na.rm = TRUE)
yr <- range(df_t$UMAP_2, na.rm = TRUE)
dx <- diff(xr); dy <- diff(yr)
x0 <- xr[1]; y0 <- yr[1]
ax_len <- 0.12 * min(dx, dy)
hx <- ax_len; hy <- ax_len

p_t <- ggplot(df_t, aes(UMAP_1, UMAP_2)) +
  geom_point(
    aes(color = sub_type),
    size = 0.6, alpha = 0.9, show.legend = FALSE
  ) +
  scale_color_manual(values = tcell_palette, drop = FALSE) +
  stat_ellipse(
    data = df_ell,
    aes(group = sub_type),
    type = "norm", level = 0.9,
    linetype = "dashed", linewidth = 0.35, color = "black",
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = label_pos,
    aes(x = cx, y = cy, label = sub_type),
    inherit.aes = FALSE,
    size = 5, fontface = "bold",
    box.padding = 0.15, point.padding = 0.2,
    min.segment.length = 0, segment.size = 0.3
  ) +
  geom_segment(
    inherit.aes = FALSE,
    aes(x = x0, y = y0, xend = x0 + hx, yend = y0),
    linewidth = 0.8,
    arrow = arrow(length = unit(0.2, "cm"), type = "closed")
  ) +
  geom_segment(
    inherit.aes = FALSE,
    aes(x = x0, y = y0, xend = x0, yend = y0 + hy),
    linewidth = 0.8,
    arrow = arrow(length = unit(0.2, "cm"), type = "closed")
  ) +
  annotate("text", x = x0 + hx*0.55, y = y0 - dy*0.03, label = "UMAP1", size = 5) +
  annotate("text", x = x0 - dx*0.03, y = y0 + hy*0.55, label = "UMAP2", size = 5, angle = 90) +
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(legend.position = "none")
p_t

pdf("tcell_umap",w=5.5,h=4.5)
p_t
dev.off

#######
.libPaths("/GPUFS/gyfyy_jxhe_1/mambaforge/envs/sc/lib/R/library")
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(Seurat)
library(ggplot2)
library(ggExtra)
library(ggpubr)

am_path      <- "/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am/Rename_am.rds"
coef_path    <- "/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am/cytosig/cytosig_scores.Coef"
tcell_rdata  <- "tcell.RData"  
sample_col   <- "sample"        

am <- readRDS(am_path)
load(tcell_rdata)   

coef <- fread(coef_path, header = FALSE, sep = "\t")
coef <- as.data.frame(coef)
rownames(coef) <- coef[[1]]
coef <- coef[, -1, drop = FALSE]

stopifnot(ncol(coef) == ncol(am))
colnames(coef) <- colnames(am)


module_list <- list(
  CD4_Tfh_module    = c("CD40LG", "TNFRSF25"),
  CD8_Eff_module    = c("PRF1", "GZMB", "CCL5", "NKG7", "IFNG"),
  CD4_Memory_module = c("IL7R", "CCR7", "KLF2"),
  Treg_module       = c("FOXP3", "CTLA4", "IL2RA", "TNFRSF4")
)

tcell <- AddModuleScore(
  object  = tcell,
  features = module_list,
  name    = names(module_list)
)

meta_t <- tcell@meta.data
meta_t$sample <- meta_t[[sample_col]]

module_cols <- c("CD4_Tfh_module1",
                 "CD8_Eff_module2",
                 "CD4_Memory_module3",
                 "Treg_module4")

tcell_sample_mean <- meta_t %>%
  group_by(sample) %>%
  summarise(
    CD4_Tfh = mean(CD4_Tfh_module1,      na.rm = TRUE),
    CD8_Eff = mean(CD8_Eff_module2,      na.rm = TRUE),
    Memory  = mean(CD4_Memory_module3,   na.rm = TRUE),
    Treg    = mean(Treg_module4,         na.rm = TRUE),
    .groups = "drop"
  )

meta_am <- am@meta.data

am2_cells <- rownames(meta_am)[meta_am$cell_type == "AM-2.CD4"]
meta_am2  <- meta_am[am2_cells, , drop = FALSE]

cyto_am2_mat <- t(coef[, am2_cells, drop = FALSE]) 

gene_var   <- apply(cyto_am2_mat, 2, var, na.rm = TRUE)
cytos_use  <- names(gene_var[gene_var > quantile(gene_var, 0.3)])  
cyto_am2_mat <- cyto_am2_mat[, cytos_use, drop = FALSE]

cyto_am2_df <- as.data.frame(cyto_am2_mat)
cyto_am2_df$sample <- meta_am2[[sample_col]]

am2_sample_mean <- cyto_am2_df %>%
  group_by(sample) %>%
  summarise(across(all_of(cytos_use), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup()



df_merge <- inner_join(am2_sample_mean, tcell_sample_mean, by = "sample")

modules <- c("CD4_Tfh", "CD8_Eff", "Memory", "Treg")


cor_result <- purrr::map_dfr(
  cytos_use,
  function(g) {
    purrr::map_dfr(
      modules,
      function(m) {
        df_test <- df_merge[, c(g, m)]
        colnames(df_test) <- c("cyto", "module")
        if (sd(df_test$cyto, na.rm = TRUE) == 0) return(NULL)
        
        ct <- suppressWarnings(
          cor.test(df_test$cyto, df_test$module, method = "spearman")
        )
        
        tibble(
          cytokine = g,
          module   = m,
          rho      = unname(ct$estimate),
          pval     = ct$p.value
        )
      }
    )
  }
)

p_cut  <- 0.05
top_n    <- 100

top_hits <- cor_result %>%
  filter(pval < p_cut) %>%
  mutate(direction = ifelse(rho > 0, "positive", "negative")) %>%
  group_by(module) %>%
  slice_max(order_by = abs(rho), n = top_n) %>%
  arrange(module, desc(abs(rho))) %>%
  ungroup()

top_hits

plot_cyto_module_median <- function(df_merge,
                                    gene_use,
                                    t_module,
                                    x_label = NULL,
                                    y_label = NULL) {
  df <- df_merge[, c("sample", gene_use, t_module)]
  colnames(df) <- c("sample", "X", "Y")
  df <- df[complete.cases(df), ]
  
  if (nrow(df) < 3) {
    stop(paste("Too few samples for", gene_use, "x", t_module))
  }
  med <- median(df$X, na.rm = TRUE)
  df$Group <- ifelse(df$X >= med, "High", "Low")
  df$Group <- factor(df$Group, levels = c("Low", "High"))
  
  if (is.null(x_label)) {
    x_label_use <- paste0("AM-2.CD4 ", gene_use, " score")
  } else {
    x_label_use <- x_label
  }
  
  if (is.null(y_label)) {
    base_mod <- sub("\\d+$", "", t_module)
    y_label_use <- case_when(
      grepl("CD8_Eff", base_mod)    ~ "CD8 Effector T module",
      grepl("CD4_Tfh", base_mod)    ~ "CD4 Tfh module",
      grepl("Memory",   base_mod)   ~ "CD4 Memory module",
      grepl("Treg",     base_mod)   ~ "Treg module",
      TRUE                          ~ base_mod
    )
  } else {
    y_label_use <- y_label
  }
  
  cor_res <- suppressWarnings(
    cor.test(df$X, df$Y, method = "spearman")
  )
  r_val <- as.numeric(cor_res$estimate)
  p_val <- cor_res$p.value
  
  r_lab <- paste0("r = ", formatC(r_val, digits = 3, format = "f"))
  p_lab <- if (p_val < 0.001) {
    "p < 0.001"
  } else {
    paste0("p = ", formatC(p_val, digits = 3, format = "f"))
  }
  label_text <- paste(r_lab, p_lab, sep = "\n")
  
  x_range <- range(df$X, na.rm = TRUE)
  y_range <- range(df$Y, na.rm = TRUE)
  y_top   <- y_range[1] + 0.80 * diff(y_range)  
  x_left  <- x_range[1] + 0.03 * diff(x_range)
  
  y_top2   <- y_range[1] + 0.95 * diff(y_range)
  x_left2  <- x_range[1] + 0.05 * diff(x_range)
  x_right2 <- x_range[1] + 0.95 * diff(x_range)
  
  p <- ggplot(df, aes(x = X, y = Y)) +
    geom_point(aes(color = Group), size = 2) +
    geom_vline(xintercept = med, linetype = 2) +
    geom_smooth(method = lm, formula = y ~ x, color = "grey40", se = FALSE) +
    annotate("text",
             x = x_left, y = y_top,
             label = label_text,
             hjust = 0, vjust = 1,
             size = 3.5) +
    annotate("text",
             x = x_left2, y = y_top2,
             label = "Low",
             color = "#70CDBE",
             size = 4, hjust = 0) +
    annotate("text",
             x = x_right2, y = y_top2,
             label = "High",
             color = "#AC99D2",
             size = 4, hjust = 1) +
    scale_color_manual(values = c("Low"="#70CDBE","High"="#AC99D2")) +
    theme_test() +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x  = element_text(size = 11),
      axis.text.y  = element_text(size = 11),
      plot.title   = element_blank()
    ) +
    labs(
      x = x_label_use,
      y = y_label_use
    )
  
  p_marginal <- ggMarginal(
    p,
    type        = "density",
    margins     = "y",
    groupColour = TRUE,
    groupFill   = TRUE
  )
  
  return(p_marginal)
}


hits_eff_tfh <- top_hits %>%
  filter(module %in% c("CD4_Tfh","CD8_Eff")) %>%
  arrange(factor(module, levels = c("CD4_Tfh","CD8_Eff")),
          desc(abs(rho)))

hits_eff_tfh
plots_eff_tfh <- purrr::map2(
  hits_eff_tfh$cytokine,
  hits_eff_tfh$module,
  ~ plot_cyto_module_median(
    df_merge = df_merge,
    gene_use = .x,
    t_module = .y  
  )
)
p_eff_tfh_grid <- ggarrange(
  plotlist = plots_eff_tfh,
  nrow = 2,
  ncol = 3,
  align = "hv"
)

p_eff_tfh_grid

