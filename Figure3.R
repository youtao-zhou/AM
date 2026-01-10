### Figure3

library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(tidyr)
library(survival)


select<-dplyr::select
source("get_otters.R")
pbmc<-readRDS("Rename_am.rds")

result_lung<-get_otters("stage2")
result_lung2<-subset(result_lung,result_lung$OTTERS_p<0.05)
twasgene<-result_lung2$Gene
twasgene<-bitr(twasgene,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
twasgene_list <- list(twasgene)


pbmc <- AddModuleScore(
  object = pbmc,
  features = twasgene_list,
  ctrl = 100, 
  name = 'TWAS_Signature'
)

colnames(pbmc@meta.data)

######## umap TWAS Signature
umap_data <- as.data.frame(Embeddings(pbmc, "umap"))
umap_data$TWAS_Signature1 <- pbmc@meta.data$TWAS_Signature1

mydata<- FetchData(pbmc,vars = c("UMAP_1","UMAP_2","TWAS_Signature1"))
colnames(mydata)[3]<-"TWAS_Signature"

seurat_object<-pbmc


data<- FetchData(seurat_object,vars = c("cell_type","TWAS_Signature1"))
exprSet<-data
exprSet2 <- na.omit(exprSet[,c("cell_type","TWAS_Signature1")])
exprSet2$TWAS_Signature1 <- scale(exprSet2$TWAS_Signature1)
exprSet2$AM_group <- factor(ifelse(grepl("AM-1", exprSet2$cell_type), "AM-1", 
                                   ifelse(grepl("AM-2", exprSet2$cell_type), "AM-2", 
                                          ifelse(grepl("AM-3", exprSet2$cell_type), "AM-3",
                                                 ifelse(grepl("AM-4", exprSet2$cell_type), "AM-4", "AM-5")))),
                            levels = c("AM-1", "AM-2", "AM-3", "AM-4", "AM-5"))
exprSet2$marker <- gsub(".*\\.", "", exprSet2$cell_type)
exprSet2_filtered <- exprSet2[!is.na(exprSet2$TWAS_Signature1), ]
colors <- c("#e97257", "#72c5d9", "#58ae9a", "#646e9a", "#efa78e", 
            "#f3a9c4", "#ffbf69", "#b4d33e", "#ed7d31", "#6b83b3",   
            "#b7d3ff", "#ff72cc", "#d4a4ff")  
theme <- theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 315, size = 16, hjust = 0),  
    axis.text.y = element_text(size = 16),  
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 16, color = "#b35959", face = "bold"), 
    axis.line = element_line(size = 0.5), 
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    panel.background = element_rect(fill = '#f3f6f6'),
    legend.title = element_blank(),
    legend.text = element_text(size = 16), 
    strip.text = element_text(size = 16, face = "bold", color = "black"),  
    legend.position = "none"
  )

p <- ggplot(exprSet2_filtered, aes(x = marker, y = TWAS_Signature1)) +
  geom_violin(aes(fill = cell_type), color = NA,
              alpha = 0.6, width = 0.7,
              trim = TRUE, scale = "width", show.legend = FALSE) + 
  geom_boxplot(aes(fill = cell_type),
               width = 0.7, size = 0.5,
               alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +  
  scale_y_continuous(expand = c(0, 0), limits = c(-3.5, 3.1), 
                     breaks = seq(-3, 3, 1), labels = seq(-3, 3, 1)) +
  scale_fill_manual(values = colors) + 
  scale_color_manual(values = colors) +
  theme +
  facet_wrap(~ AM_group, scales = "free_x") +  
  labs(title = "", 
       x = "",
       y = "TWAS Signature",
       hjust = 0.5)

print(p)


theme <- theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  
    axis.text.x = element_text(angle = 0, size = 16, hjust = 0),  
    axis.text.y = element_text(size = 16), 
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 16, color = "#b35959", face = "bold"), 
    axis.line = element_line(size = 0.5),  
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    panel.background = element_rect(fill = '#f3f6f6'),
    legend.title = element_blank(),
    legend.text = element_text(size = 16), 
    strip.text = element_text(size = 16, face = "bold", color = "black"), 
    legend.position = "none"
  )
colors<-c(brewer.pal(8,"Dark2"),brewer.pal(5,"Set1"))
data$cell_type<-substring(data$cell_type,1,4)
my_comparisons <- list( c("AM-1", "AM-2"), c("AM-2", "AM-3"),c("AM-3", "AM-4"),c("AM-4", "AM-5"))
ggplot(data, aes(x = cell_type, y = log(TWAS_Signature1 + 1))) +
  theme_bw() +
  RotatedAxis() +
  geom_violin(aes(fill = cell_type), color = NA,
              alpha = 0.6, width = 0.7,
              trim = TRUE, scale = "width", show.legend = FALSE) + 
  geom_boxplot(aes(fill = cell_type),
               width = 0.7, size = 0.5,
               alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
  scale_fill_manual(values = colors) + 
  scale_color_manual(values = colors) +
  theme +
  labs(title = "", 
       x = "",
       y = "TWAS Signature",
       hjust = 0.5)+
  stat_compare_means(
    comparisons = my_comparisons, 
    label = "p.signif",  
    size = 6 
  ) +
  labs(
    x = "", y = "TWAS Signature", 
  ) +ylim(0,0.08)+NoLegend()


##################### AddModuleScore_bulk
AddModuleScore_bulk <- function (object, features, nbin = 24, ctrl = 100, 
                                 k = FALSE, assay = NULL, name = "Cluster", seed = 1) 
{
  set.seed(seed = seed)
  
  features <- lapply(X = features, FUN = function(x) {
    missing.features <- setdiff(x = x, y = rownames(x = object))
    if (length(x = missing.features) > 0) {
      warning("The following features are not present in the object: ", 
              paste(missing.features, collapse = ", "), 
              ifelse(test = FALSE, yes = ", attempting to find updated synonyms", 
                     no = ", not searching for symbol synonyms"), 
              call. = FALSE, immediate. = TRUE)
    }
    return(intersect(x = x, y = rownames(x = object)))
  })
  
  cluster.length <- length(x = features)
  
  pool <- rownames(x = object)
  data.avg <- apply(object, 1, mean)
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                         n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(data.cut == 
                                                                              data.cut[features.use[j]])], size = ctrl, replace = FALSE)))
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
                        ncol = ncol(x = object))
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- apply(object[features.use, ], 2, mean)
  }
  features.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, 
                            ncol = ncol(x = object))
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- object[features.use, , drop = FALSE]
    features.scores[i, ] <- apply(data.use, 2, mean)
  }
  features.scores.use <- features.scores - ctrl.scores
  
  colnames(features.scores.use) <- colnames(object)
  rownames(features.scores.use) <- paste0(name, 1:cluster.length)
  return(features.scores.use)
}
# # ########### TWAS signature
result_lung<-get_otters("stage2")
result_lung2<-subset(result_lung,result_lung$OTTERS_p<0.05)
twasgene<-result_lung2$Gene
twasgene<-bitr(twasgene,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
twasgene<-intersect(twasgene$SYMBOL, rownames(bulk_dataset))
twasgene_list <- list(twasgene)

############### POPLAR TWAS_signature
############### Atelizumab/MPDL3280A ~ PFS ; All People ~ OS

data1<-read.csv("anno.csv",header = T)
data2<-read.csv("BYN_n192.csv",header = T)
data3<-read.csv("ss_ss.csv",header = T)
data3<-data3%>%dplyr::select(alias,subjectId,gender,phenotype)

data4<-left_join(data1,data2,by=c("UNI_ID"))
data4<-subset(data4,data4$HIST%in%"NON-SQUAMOUS")
phenotype<-left_join(data4,data3,by=c("trunc_anonymized_sample_ids"="alias"))
phenotype<-phenotype[phenotype$BCOR%in%c("CR","PR","SD"),]
# phenotype<-phenotype[phenotype$ACTARM%in%"MPDL3280A",]  

survival<-phenotype%>%dplyr::select(sampleID=trunc_anonymized_sample_ids,OS_MONTHS,OS_CENSOR)
survival$OS_time<-survival$OS_MONTHS*30
survival<-survival%>%dplyr::select(sampleID,OS_time,OS_status=OS_CENSOR)
bulk_survival<-na.omit(survival)
#################### 


bulk_dataset<-read.csv("POPLAR.csv")
colnames(bulk_dataset) <- gsub("\\.", "-", colnames(bulk_dataset))
bulk_dataset<-bulk_dataset[,bulk_survival$sampleID]

bulk_dataset<-as.matrix(bulk_dataset)

score <- AddModuleScore_bulk(bulk_dataset, twasgene_list, name='TWAS_signature')

TWAS_signature_score<-as.data.frame(t(score))
TWAS_signature_score$trunc_anonymized_sample_ids<-rownames(TWAS_signature_score)
TWAS_signature_score<-left_join(TWAS_signature_score,phenotype,by="trunc_anonymized_sample_ids")
TWAS_signature_score$TWAS_signature<-ifelse(TWAS_signature_score$TWAS_signature1>median(TWAS_signature_score$TWAS_signature1),"High","Low")
TWAS_signature_score_POPLAR<-TWAS_signature_score

my_comparisons <- list( c("High","Low"))
pdf("POPLAR_TWAS_signature_boxplot.pdf",width = 6,height = 4)
ggplot(TWAS_signature_score, aes(x = TWAS_signature, y = OS_MONTHS, fill = TWAS_signature)) +
  geom_boxplot(outlier.shape = NA, alpha = 1, color = "black", width = 0.3, fill = "white") +
  geom_jitter(aes(color = TWAS_signature), position = position_jitter(width = 0.2), size = 4, alpha = 0.8) +
  scale_fill_manual(values = c("High" = "#DA4644", "Low" = "#308ABE")) +
  scale_color_manual(values = c("High" = "#DA4644", "Low" = "#308ABE")) +
  theme_bw() +
  ylab("OS MONTHS") +labs(color="TWAS Signature")+guides(fill="none")+
  theme(
    plot.title = element_blank(),  
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 18, face = "bold"),  
    axis.text.y = element_text(size = 18, face = "bold"), 
    axis.text.x = element_text(size = 18, face = "bold"),  
    panel.grid.major = element_line(color = "gray80"), 
    panel.grid.minor = element_line(color = "gray90"),  
    panel.border = element_blank(), 
    axis.line = element_line(color = "black"), 
    legend.position = "right",  
    legend.text = element_text(size = 16, face = "bold"),  
    legend.title = element_text(size = 16, face = "bold"),  
    plot.margin = margin(t = 50, r = 20, b = 20, l = 20)  
  )
dev.off()

### OAK TWAS_signature
twasgene<-result_lung2$Gene
twasgene<-bitr(twasgene,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
data1<-read.csv("go28915_updated.csv")
data2<-read.csv("BYN_n699.csv")
data3<-read.csv("ss_ss.csv")
biomarker<-read.table("biomarker.txt",header=T,sep="\t")
data4<-left_join(data1,data3,by=c("trunc_anonymized_sample_ids"="alias"))
phenotype<-left_join(data4,biomarker,by=c("trunc_UNI_ID"="trunc_anonymized_patient_id"))
phenotype<-separate(phenotype,col = "full_anonymized_patient_id",into = c("PAT","UNI_ID"))
phenotype<-left_join(phenotype,data2,by="UNI_ID")
phenotype<-phenotype[phenotype$HIST%in%"NON-SQUAMOUS",]
phenotype<-phenotype[phenotype$TISSUE_TYPE.2%in%"lung",]
phenotype<-phenotype[phenotype$BCOR%in%c("CR","PR","SD"),]

bulk_dataset<-readRDS("OAK.rds")

colnames(bulk_dataset) <- gsub("\\.", "-", colnames(bulk_dataset))
bulk_dataset<-bulk_dataset[,phenotype$trunc_anonymized_sample_ids]

score <- AddModuleScore_bulk(bulk_dataset, twasgene_list, name='TWAS_signature')

TWAS_signature_score<-as.data.frame(t(score))
TWAS_signature_score$trunc_anonymized_sample_ids<-rownames(TWAS_signature_score)
TWAS_signature_score<-left_join(TWAS_signature_score,phenotype,by="trunc_anonymized_sample_ids")
TWAS_signature_score$TWAS_signature<-ifelse(TWAS_signature_score$TWAS_signature1>mean(TWAS_signature_score$TWAS_signature1),"High","Low")
TWAS_signature_score_OAK<-TWAS_signature_score
my_comparisons <- list( c("High","Low"))
pdf("OAK_TWAS_signature_boxplot.pdf",width = 6,height = 4)
ggplot(TWAS_signature_score, aes(x = TWAS_signature, y = OS_MONTHS, fill = TWAS_signature)) +
  geom_boxplot(outlier.shape = NA, alpha = 1, color = "black", width = 0.3, fill = "white") +
  geom_jitter(aes(color = TWAS_signature), position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("High" = "#DA4644", "Low" = "#308ABE")) +
  scale_color_manual(values = c("High" = "#DA4644", "Low" = "#308ABE")) +
  theme_bw() +
  ylab("OS MONTHS") +labs(color="TWAS Signature")+guides(fill="none")+
  theme(
    plot.title = element_blank(),  
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 18, face = "bold"),  
    axis.text.y = element_text(size = 18, face = "bold"), 
    axis.text.x = element_text(size = 18, face = "bold"),  
    panel.grid.major = element_line(color = "gray80"), 
    panel.grid.minor = element_line(color = "gray90"),  
    panel.border = element_blank(), 
    axis.line = element_line(color = "black"), 
    legend.position = "right",  
    legend.text = element_text(size = 16, face = "bold"),  
    legend.title = element_text(size = 16, face = "bold"), 
    plot.margin = margin(t = 50, r = 20, b = 20, l = 20)  
  )
dev.off()

library(survival)
library(survminer)
TWAS_signature_score<-TWAS_signature_score_POPLAR
surv_object <- Surv(time = TWAS_signature_score$OS_MONTHS, event = TWAS_signature_score$OS_CENSOR)

cox_model <- coxph(surv_object ~ TWAS_signature1, data = TWAS_signature_score)

summary(cox_model)

TWAS_signature_score$Group <- ifelse(TWAS_signature_score$TWAS_signature1 > mean(TWAS_signature_score$TWAS_signature1), "High", "Low")

p<-ggsurvplot(
  survfit(surv_object ~ Group, data = TWAS_signature_score),
  data=TWAS_signature_score,
  palette =c("#D55E00", "#0072B2"),
  conf.int = TRUE,
  pval = TRUE,
  legend.title = "TWAS Signature",
  legend.labs = c("High", "Low"),
  xlab = "Atelizumab - OS",
  font.x = 20, font.y = 20,font.legend = 18,pval.size = 7, 
  ggtheme = theme_bw(base_size = 25)
)
p

ggexport(
  p, 
  filename = "TWAS_signature_POPLAR_Atelizumab_OS_survivalcurve.pdf",
  width = 6, 
  height = 4.5
)

#####################
load("psetime_res.RData")
data <- readRDS("OAK_ssGSEA.rds")  
OAK  <- readRDS("OAK.rds")   

pt_anno <- psetime_res@pseudotime_anno
center_raw <- pt_anno %>% group_by(time_anno) %>% summarise(mean_time = mean(time), .groups="drop")
m4 <- center_raw %>% filter(grepl("^AM-4", time_anno)) %>% pull(mean_time) %>% mean(na.rm=TRUE)
m1 <- center_raw %>% filter(grepl("^AM-1", time_anno)) %>% pull(mean_time) %>% mean(na.rm=TRUE)
if (!is.na(m4) && !is.na(m1) && m4 > m1) {
  rng <- range(pt_anno$time, na.rm=TRUE)
  pt_anno$time_aligned <- (rng[2] - pt_anno$time) + rng[1]
} else {
  pt_anno$time_aligned <- pt_anno$time
}

pt_long <- psetime_res@pseudotime_matrix_long %>%
  left_join(pt_anno %>% select(time, time_aligned, time_anno),
            by = c("pseudotime" = "time"))

center_tbl <- pt_anno %>%
  group_by(time_anno) %>% summarise(mean_time = mean(time_aligned), .groups="drop") %>%
  mutate(meta = case_when(
    grepl("^AM-2", time_anno) ~ "AM2_like",
    grepl("^AM-5", time_anno) ~ "AM5_like",
    grepl("^AM-4", time_anno) ~ "AM4_like",
    TRUE ~ "Other"
  )) %>%
  group_by(meta) %>%
  summarise(center = mean(mean_time, na.rm=TRUE), .groups="drop") %>%
  filter(meta %in% c("AM2_like","AM5_like","AM4_like"))

c2 <- center_tbl$center[center_tbl$meta=="AM2_like"]
c5 <- center_tbl$center[center_tbl$meta=="AM5_like"]
c4 <- center_tbl$center[center_tbl$meta=="AM4_like"]

rng   <- range(pt_long$time_aligned, na.rm=TRUE)
sigma <- 0.15 * diff(rng)                                  
kernel_w <- function(t, c, s) exp(-0.5 * ((t - c)/s)^2)

pt_long <- pt_long %>%
  mutate(w2 = kernel_w(time_aligned, c2, sigma),
         w5 = kernel_w(time_aligned, c5, sigma),
         w4 = kernel_w(time_aligned, c4, sigma))

trend_245 <- pt_long %>%
  group_by(gene) %>%
  summarise(
    r2 = suppressWarnings(cor(exp, w2, method="spearman", use="pairwise.complete.obs")),
    r5 = suppressWarnings(cor(exp, w5, method="spearman", use="pairwise.complete.obs")),
    r4 = suppressWarnings(cor(exp, w4, method="spearman", use="pairwise.complete.obs")),
    .groups="drop"
  ) %>% filter(!(is.na(r2) & is.na(r5) & is.na(r4)))

thr_contrast <- 0.4
trend_contrast <- trend_245 %>%
  mutate(
    R25   = (r2 + r5)/2,
    R4    = r4,
    delta = R25 - R4,
    group = case_when(
      delta >=  thr_contrast ~ "AM(2+5) Program",
      delta <= -thr_contrast ~ "AM4 Program",
      TRUE ~ "Intermediate"
    )
  )

p_trend_contrast <- ggplot(trend_contrast, aes(x = reorder(gene, delta), y = delta, color = group)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_hline(yintercept = c(thr_contrast, -thr_contrast), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("AM(2+5) Program" = "#E64B35FF",
                                "AM4 Program"     = "#4DBBD5FF",
                                "Intermediate"    = "grey70")) +
  labs(x = "Genes", y = expression(Delta~r=="(AM(2+5) - AM4)"), color = "") +
  theme_classic(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 16),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.82, 0.2),
    legend.text = element_text(size = 16)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3.5)))
p_trend_contrast

genes_AM25_prog <- trend_contrast %>%
  filter(group == "AM(2+5) Program") %>%
  arrange(desc(delta)) %>% pull(gene) %>% unique()

genes_AM4_prog2 <- trend_contrast %>%
  filter(group == "AM4 Program") %>%
  arrange(delta) %>% pull(gene) %>% unique()

expr_z <- scale(t(OAK), center = TRUE, scale = TRUE)  

AM25_use <- intersect(colnames(expr_z), genes_AM25_prog)
AM4_use2 <- intersect(colnames(expr_z), genes_AM4_prog2)


prog_df <- data.frame(
  sample           = rownames(expr_z),
  AM25_prog_score  = rowMeans(expr_z[, AM25_use,  drop=FALSE], na.rm=TRUE),
  AM4_prog_score   = rowMeans(expr_z[, AM4_use2,  drop=FALSE], na.rm=TRUE)
) %>%
  mutate(AM4vsAM25_prog = AM4_prog_score - AM25_prog_score)

meta_weights <- data %>% as.data.frame() %>% tibble::rownames_to_column("sample")
am2_cols <- setdiff(grep("^AM-2", colnames(meta_weights), value = TRUE), "sample")
am5_cols <- setdiff(grep("^AM-5", colnames(meta_weights), value = TRUE), "sample")
am4_cols <- setdiff(grep("^AM-4", colnames(meta_weights), value = TRUE), "sample")

am25_like_tbl <- meta_weights %>%
  mutate(
    AM2_sum = if (length(am2_cols)>0) rowSums(across(all_of(am2_cols)), na.rm=TRUE) else 0,
    AM5_sum = if (length(am5_cols)>0) rowSums(across(all_of(am5_cols)), na.rm=TRUE) else 0,
    AM25_Like_Score = AM2_sum + AM5_sum
  ) %>% select(sample, AM25_Like_Score)

am4_like_tbl <- meta_weights %>%
  mutate(AM4_Like_Score = if (length(am4_cols)>0) rowSums(across(all_of(am4_cols)), na.rm=TRUE) else 0) %>%
  select(sample, AM4_Like_Score)

df_corr2 <- prog_df %>%
  select(sample, AM25_prog_score, AM4_prog_score, AM4vsAM25_prog) %>%
  left_join(am25_like_tbl, by="sample") %>%
  left_join(am4_like_tbl,  by="sample") %>%
  transmute(
    sample,
    AM25_Like_Score   = AM25_Like_Score,
    AM4_Like_Score    = AM4_Like_Score,
    AM25_prog_score   = AM25_prog_score,
    AM4_prog_score    = AM4_prog_score,
    Progression_Score = AM4vsAM25_prog     
  )

BetterVis_Dotplot_DualY <- function(data, x_var, y_left, y_right,
                                    point_size = 5,
                                    smoothline = TRUE,
                                    cor_test = TRUE,
                                    cor_method = "pearson",
                                    x_title = NULL,
                                    y_left_title = NULL, y_right_title = NULL,
                                    legend_position = "bottom",
                                    color_left = "#ACD6EC", color_right = "#F5A889",
                                    label_size = 4.5) {
  library(ggplot2)
  library(ggprism)
  library(scales)

  x_lab <- ifelse(is.null(x_title), x_var, x_title)
  y_left_lab  <- ifelse(is.null(y_left_title), y_left, y_left_title)
  y_right_lab <- ifelse(is.null(y_right_title), y_right, y_right_title)

  range01 <- function(x) (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
  data$y_left_scaled  <- range01(data[[y_left]])
  data$y_right_scaled <- range01(data[[y_right]])

  p <- ggplot(data, aes_string(x = x_var)) +
    geom_point(aes(y = y_left_scaled, color=y_left_lab),
               shape=16, size=point_size, alpha=0.8) +
    {if (smoothline) geom_smooth(aes(y = y_left_scaled, color=y_left_lab),
                                 method="loess", se=FALSE, show.legend=FALSE)} +
    geom_point(aes(y = y_right_scaled, color=y_right_lab),
               shape=16, size=point_size, alpha=0.8) +
    {if (smoothline) geom_smooth(aes(y = y_right_scaled, color=y_right_lab),
                                 method="loess", se=FALSE, show.legend=FALSE)} +

    scale_y_continuous(
      name = y_left_lab,
      sec.axis = sec_axis(~., name=y_right_lab)
    ) +
    scale_color_manual(values = setNames(c(color_left, color_right),
                                         c(y_left_lab, y_right_lab))) +
    
    theme_prism(base_size=18, border=TRUE) +
    labs(x=x_lab) +
    theme(
      legend.position = legend_position,
      legend.title = element_blank(),
      legend.text = element_text(size=14),
      axis.text = element_text(size=14),
      axis.title = element_text(size=16),
      panel.grid.major = element_line(color="grey85", size=0.3, linetype=2),
      panel.grid.minor = element_line(color="grey92", size=0.2, linetype=2),
      plot.title = element_blank()
    )
  if (cor_test) {
    cor_left  <- cor.test(data[[x_var]], data[[y_left]],  method=cor_method)
    cor_right <- cor.test(data[[x_var]], data[[y_right]], method=cor_method)
    
    fmt_p <- function(p) {
      if (p < 0.001) return("p<0.001")
      else return(paste0("p=", formatC(p, format="f", digits=3)))
    }
    label_left  <- paste0("r=", round(cor_left$estimate,2), ", ", fmt_p(cor_left$p.value))
    label_right <- paste0("r=", round(cor_right$estimate,2), ", ", fmt_p(cor_right$p.value))
    
    p <- p +
      annotate("label", x=min(data[[x_var]],na.rm=TRUE), y=1,
               label=label_left, hjust=0, vjust=1,
               size=label_size, color=color_left, fill="white",
               label.size=0.6, linetype=2) +
      annotate("label", x=max(data[[x_var]],na.rm=TRUE), y=1,
               label=label_right, hjust=1, vjust=1,
               size=label_size, color=color_right, fill="white",
               label.size=0.6, linetype=2)
  }
  
  return(p)
}
p_prog_dual <- BetterVis_Dotplot_DualY(
  data = df_corr2,
  x_var = "Progression_Score",
  y_left  = "AM25_Like_Score",
  y_right = "AM4_Like_Score",
  x_title = "Progression",
  y_left_title  = "AM(2+5) ssGSEA Score",
  y_right_title = "AM4 ssGSEA Score",
  point_size = 2
)
p_prog_dual

data1 <- read.csv("go28915_updated.csv")
data2 <- read.csv("BYN_n699.csv")
data3 <- read.csv("ss_ss.csv")
biomarker <- read.table("biomarker.txt", header=TRUE, sep="\t")

data4 <- left_join(data1, data3, by=c("trunc_anonymized_sample_ids"="alias"))
phenotype <- left_join(data4, biomarker, by=c("trunc_UNI_ID"="trunc_anonymized_patient_id"))
phenotype <- tidyr::separate(phenotype, col="full_anonymized_patient_id", into=c("PAT","UNI_ID"))
phenotype <- left_join(phenotype, data2, by = "UNI_ID")

phenotype <- phenotype[phenotype$HIST %in% "NON-SQUAMOUS", ]
phenotype <- phenotype[phenotype$TISSUE_TYPE.2 %in% "lung", ]
phenotype <- phenotype[phenotype$BCOR %in% c("CR", "PR", "SD"), ]


survival <- phenotype %>%
  dplyr::select(sample = trunc_anonymized_sample_ids,
                PFS_MONTHS, PFS_CENSOR) %>%
  transmute(
    sample,
    PFS_time   = PFS_MONTHS,   
    PFS_status = PFS_CENSOR     
  ) %>%
  na.omit()

df_all_treat <- survival %>%
  left_join(phenotype %>% select(sample = trunc_anonymized_sample_ids, ACTARM, gender),
            by = "sample") %>%
  mutate(
    Treatment = case_when(
      ACTARM == "MPDL3280A" ~ "IO",
      ACTARM == "Docetaxel" ~ "Chemo",
      TRUE ~ "Other"
    )
  )

stopifnot("AM4vsAM25_prog" %in% names(prog_df))

cox_chemo_prog_AM25 <- coxph(
  Surv(PFS_time, PFS_status) ~ scale(AM4vsAM25_prog),
  data = subset(df_all_treat %>% left_join(prog_df, by = "sample"),
                ACTARM == "Docetaxel")
)
cox_io_prog_AM25 <- coxph(
  Surv(PFS_time, PFS_status) ~ scale(AM4vsAM25_prog),
  data = subset(df_all_treat %>% left_join(prog_df, by = "sample"),
                ACTARM == "MPDL3280A")
)

cat("\n== Cox: Docetaxel ==\n"); print(summary(cox_chemo_prog_AM25))
cat("\n== Cox: Atezolizumab ==\n"); print(summary(cox_io_prog_AM25))

###########

extract_cox <- function(fit, drug, var_label){
  s <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  tibble::tibble(
    Drug = drug,
    Variable = var_label,
    Estimate = s$estimate,
    lower.limit = s$conf.low,
    upper.limit = s$conf.high,
    pval = s$p.value,
    sig = ifelse(s$p.value < 0.05, "*", "")
  )
}

df_forest_AM25 <- dplyr::bind_rows(
  extract_cox(cox_chemo_prog_AM25, "Docetaxel",    "AM4 − (AM2+AM5)"),
  extract_cox(cox_io_prog_AM25,    "Atezolizumab", "AM4 − (AM2+AM5)")
)

b1_AM25 <- ggplot(df_forest_AM25, aes(x = Estimate, y = Variable)) +
  geom_point(size = 4, aes(color = Variable)) +
  geom_errorbarh(aes(xmin = lower.limit, xmax = upper.limit, color = Variable),
                 size = 1, height = 0.2) +
  geom_vline(xintercept = 1, color = "gray50", linetype = "dashed") +
  geom_text(aes(label = sig, x = pmin(0.12, min(Estimate, na.rm = TRUE))), size = 7) +
  facet_wrap(~ Drug, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("AM4 − (AM2+AM5)" = "#FA7F6F")) +
  scale_x_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5)) +
  xlab("Hazard Ratio (95% CI)") + ylab("") +
  theme_few() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title  = element_text(size = 16, color = "black"),
        strip.text  = element_text(size = 16, face = "bold"),
        legend.position = "none")

b1_AM25
