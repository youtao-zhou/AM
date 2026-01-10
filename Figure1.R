## Figure1

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)
library(forcats)
library(Scissor)
library(tidyverse)
library(ggrepel)
library(cowplot)
library(SeuratObject)
library(scales)
library(grid)
library(ComplexHeatmap)
library(grDevices)
library(colorRamp2)
library(MetBrewer)
library(ReporterScore)
library(stringr)
library(forcats)

files<-as.data.frame(read_excel("HTAN_scRNA_files.xlsx"))

remove_all_na_columns <- function(df) {
  df %>% select_if(~ !all(is.na(.)))
}
files <- remove_all_na_columns(files)
files <- files %>%
  mutate(sample = str_extract(Biospecimen, "^[^_]+_[^_]+"))
files<-subset(files,files$Diagnosis%in%c("Adenocarcinoma NOS","Adenocarcinoma combined with other types of carcinoma"))
files<-subset(files,files$`Data Access`%in%"dbGaP"&files$`File Format`%in%"fastq")
files<-files%>%dplyr::select(Biospecimen,sample)


cases<-as.data.frame(read_excel("Lungcases.xlsx"))
cases<- remove_all_na_columns(cases)
cases<-cases%>%dplyr::select(`HTAN Participant ID`:`Primary Diagnosis`,Gender:`Vital Status`,Stage=`AJCC Pathologic Stage`)
# cases<-merge(cases,files,by.x = "HTAN Participant ID",by.y = "sample")

trnhtan<-as.data.frame(fread("HTAN_TRN_phenotype.csv"))
trnhtan<-subset(trnhtan,trnhtan$Dataset%in%c("HTAPP","HTABU","HTAMSK"))
trnhtan <- trnhtan %>%
  mutate(HTAN_ID = str_extract(sample, "^[^_]+_[^_]+"))
trnhtan<-trnhtan%>%dplyr::select(Group,HTAN_ID,sample,Dataset)


# cases<-cases%>%left_join(trnhtan,by=c("HTAN Participant ID"="HTAN_ID"))
cases<-merge(cases,trnhtan,by.x = "HTAN Participant ID",by.y = "HTAN_ID")
cases<-distinct(.data = cases,.keep_all = T)
cases<-cases%>%dplyr::select(sample=`HTAN Participant ID`,sex=Gender,age=`Age at Diagnosis (years)`,Group,Dataset,Stage)
#######################################

trn<-as.data.frame(fread("TRN_phenotype.csv"))
trncases<-trn%>%dplyr::select(sample,sex,age,Group=disease,Stage=uicc_stage)
trncases<-distinct(.data = trncases,.keep_all = T,sample)
trncases$Group[trncases$Group%in%"lung adenocarcinoma"]<-"tumor"
trncases$Dataset<-"TRN_LUAD"


###################################
TRN_HTAN_phenotype<-rbind(cases,trncases)
TRN_HTAN_phenotype$sex[TRN_HTAN_phenotype$sex%in%NA]<-"unknown"
TRN_HTAN_phenotype$sex[TRN_HTAN_phenotype$age%in%NA]<-"unknown"
TRN_HTAN_phenotype<- TRN_HTAN_phenotype %>%
  distinct()

TRN_HTAN_phenotype <- TRN_HTAN_phenotype %>%
  mutate(
    sex = ifelse(sex == "unknown", "unknown", sex),
    age_group = cut(age, breaks = c(-Inf, 20, 40, 60, Inf), labels = c("20", "40", "60", "80")),
    age_group = factor(age_group, levels = c("20", "40", "60", "80")),
    Dataset = factor(Dataset)  # 确保 Dataset 列是因子类型
  )

TRN_HTAN_phenotype <- TRN_HTAN_phenotype %>%
  group_by(Dataset, age_group) %>%
  mutate(
    x_jitter = as.numeric(Dataset) + runif(n(), min = -0.3, max = 0.3),
    y_jitter = as.numeric(age_group) + runif(n(), min = -0.3, max = 0.3)
  ) %>%
  ungroup()

pdf("phenotype.pdf",width = 12, height = 6)
ggplot(TRN_HTAN_phenotype, aes(x = x_jitter, y = y_jitter)) +
  geom_point(aes(shape = Group, color = sex), alpha = 0.7, size = 5) +
  scale_shape_manual(values = c("tumor" = 16, "normal" = 17)) + 
  scale_color_manual(values = c("male" = "#6BC4A8", "female" = "#AE8CC0", "unknown" = "#00CDD1")) + 
  scale_x_continuous(breaks = 1:length(levels(TRN_HTAN_phenotype$Dataset)), labels = levels(TRN_HTAN_phenotype$Dataset)) +
  scale_y_continuous(breaks = 1:4, labels = c("20", "40", "60", "80")) +
  labs(
    x = "",
    y = "Age",
    color = "Sex",
    shape = "Group",
    size = "Sample Count"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(linetype = "dashed"),  
    axis.text.x = element_text(angle = 0, size = 18, face = "bold"), 
    axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face = "bold"),  
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 20),  
    legend.title = element_text(size = 18),  
    legend.text = element_text(size = 18),  
    legend.key.size = unit(2, "lines"), 
    legend.position = "right",  
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    plot.title = element_blank() 
  )

dev.off()

text_size_adjust <- 14 

TRN_df <- TRN_HTAN_phenotype %>%
  mutate(
    Dataset = factor(Dataset),
    Stage_simple = case_when(
      Stage %in% c("Stage IA","Stage IA1","Stage IA2","I") ~ "Stage I",
      Stage %in% c("Stage IIA","Stage IIB","II")           ~ "Stage II",
      Stage %in% c("III","III or IV")                      ~ "Stage III",
      Stage %in% c("Stage IVC","IV")                       ~ "Stage IV",
      TRUE ~ NA_character_
    ),
    age_group = ifelse(is.na(age_group) | age_group == "", NA, as.character(age_group)),
    age_group = factor(age_group, levels = sort(unique(as.numeric(age_group))), ordered = TRUE),
    sex = ifelse(is.na(sex) | sex == "", NA, sex)
  ) %>%
  arrange(Dataset, Stage_simple, desc(age))

pal_dataset <- c("#F3BDA5" ,"#B7E4C7", "#D8B4E2", "#1f77b4")        
pal_sex     <- c("#F9D68D" ,"#FF847C", "#E0E0E0")                     
pal_stage   <- c("#7FBEEB" ,"#B19CD9", "#F89CAE", "#5FC8D1")          
pal_agegp   <- c("#8DD1C6","#8498AB","#F5BE8F","#5CB0C3") 

p_age <- ggplot(TRN_df, aes(x = Dataset, y = age, fill = Dataset)) +
  geom_boxplot(width = 0.7, color = "grey20", linewidth = 0.5, outlier.alpha = 0) +
  scale_fill_manual(values = pal_dataset, guide = "none") +
  labs(x = NULL, y = "Age") +
  theme_classic(base_size = text_size_adjust) +
  theme(axis.title.x = element_blank())

df_agegp <- TRN_df %>%
  filter(!is.na(age_group)) %>%
  count(Dataset, age_group) %>%
  group_by(Dataset) %>% mutate(prop = n / sum(n)) %>% ungroup()

p_agegrp <- ggplot(df_agegp, aes(x = Dataset, y = prop, fill = age_group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68,
           color = "black", linewidth = 0.5) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = pal_agegp, name = "Age group") +
  labs(x = NULL, y = "Proportion") +
  theme_classic(base_size = text_size_adjust) +
  theme(axis.title.x = element_blank())

df_sex <- TRN_df %>%
  filter(!is.na(sex)) %>%
  count(Dataset, sex) %>%
  group_by(Dataset) %>% mutate(prop = n / sum(n)) %>% ungroup()

p_sex <- ggplot(df_sex, aes(x = Dataset, y = prop, fill = sex)) +
  geom_col(width = 0.85, color = "black", linewidth = 0.5) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = pal_sex, name = "Sex") +
  labs(x = NULL, y = "Proportion") +
  theme_classic(base_size = text_size_adjust) +
  theme(axis.title.x = element_blank())

df_stage <- TRN_df %>%
  filter(!is.na(Stage_simple)) %>%
  count(Dataset, Stage_simple) %>%
  group_by(Dataset) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
  mutate(Stage_simple = factor(Stage_simple, levels = c("Stage I","Stage II","Stage III","Stage IV")))

p_stage <- ggplot(df_stage, aes(x = Dataset, y = prop, fill = Stage_simple)) +
  geom_col(width = 0.85, color = "black", linewidth = 0.5) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = pal_stage, name = "Stage") +
  labs(x = NULL, y = "Proportion") +
  theme_classic(base_size = text_size_adjust) +
  theme(axis.title.x = element_blank())

legend_title_size <- 12   
legend_text_size  <- 12   
legend_key_pt     <- 12   
legend_nrow       <- 2   


p_agegrp <- p_agegrp + guides(
  fill = guide_legend(nrow = legend_nrow,
                      keywidth = unit(legend_key_pt, "pt"),
                      keyheight = unit(legend_key_pt, "pt"))
)
p_sex <- p_sex + guides(
  fill = guide_legend(nrow = legend_nrow,
                      keywidth = unit(legend_key_pt, "pt"),
                      keyheight = unit(legend_key_pt, "pt"))
)
p_stage <- p_stage + guides(
  fill = guide_legend(nrow = legend_nrow,
                      keywidth = unit(legend_key_pt, "pt"),
                      keyheight = unit(legend_key_pt, "pt"))
)


p_left  <- p_age / p_agegrp
p_right <- p_sex / p_stage

p_overview_small <- (p_left | p_right) +
  plot_layout(guides = "collect") &
  theme(
    legend.position   = "bottom",
    legend.title      = element_text(size = legend_title_size),
    legend.text       = element_text(size = legend_text_size),
    legend.key.size   = unit(legend_key_pt, "pt"),
    legend.spacing.x  = unit(3, "pt"),
    legend.spacing.y  = unit(2, "pt"),
    legend.box.margin = margin(t = 2, r = 0, b = 0, l = 0)
  )

p_overview_small

#########################
trn<-as.data.frame(fread("HTAN_TRN_phenotype.csv"))
immune_or_not<-readRDS("immune_or_not.rds")
trn<-trn%>%left_join(immune_or_not,by=c("cell_type"="celltype"))
cell_phenotype<-trn%>%dplyr::select(Dataset,Group,cell_type,Group_Immune)

cell_phenotype$Group_Immune <- factor(cell_phenotype$Group_Immune, levels = c("immune", "not-immune"))

proportions <- cell_phenotype %>%
  group_by(Group_Immune, Group, cell_type) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

immune_data <- proportions %>% filter(Group_Immune == "immune")
not_immune_data <- proportions %>% filter(Group_Immune == "not-immune")

cell_colors <- c(brewer.pal(8, "Set2"), brewer.pal(9, "Set3"), brewer.pal(9, "Set1"), brewer.pal(8, "Accent"))

pdf("celltype_proprotion_notimmune.pdf",width = 8,height = 6)
p1 <- ggplot(not_immune_data, aes(x = Group, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.3, color = "black", size = 0.25) +  
  scale_fill_manual(values = cell_colors) +
  labs(x = "", y = "Percentage", fill = "Not-immune subsets") + 
  ggtitle("") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(linetype = "dashed"),  
    panel.border = element_blank(),  
    plot.title = element_text(size = 14),  
    legend.title = element_text(size = 16),  
    legend.text = element_text(size = 16),  
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.spacing.y = unit(0.1, 'cm'),
  )  

print(p1)
dev.off()


pdf("celltype_proprotion_immune.pdf",width = 8,height = 6)
p2 <- ggplot(immune_data, aes(x = Group, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.3, color = "black", size = 0.25) + 
  scale_fill_manual(values = cell_colors) +
  labs(x = "", y = "Percentage", fill = "Immune subsets") +  
  ggtitle("") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(linetype = "dashed"),  
    panel.border = element_blank(),  
    plot.title = element_text(size = 14),  
    legend.title = element_text(size = 16), 
    legend.text = element_text(size = 16),  
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.spacing.y = unit(0.1, 'cm'), 
  )  

print(p2)
dev.off()

###########
tumor  <- read.csv("Tumor_celltypes_bootstrap_results.csv")  %>%
  select(annotation, bp_value, bias_corrected_estimate, CI_lo, CI_hi) %>%
  mutate(group = "Tumor")
normal <- read.csv("Normal_celltypes_bootstrap_results.csv") %>%
  select(annotation, bp_value, bias_corrected_estimate, CI_lo, CI_hi) %>%
  mutate(group = "Normal")
tumor$annotation[tumor$annotation%in%"CD1c-positive myeloid dendritic cell"]<-"CD1c+ myeloid DC"
normal$annotation[normal$annotation%in%"CD1c-positive myeloid dendritic cell"]<-"CD1c+ myeloid DC"

get_sig_sets <- function(tumor, normal, alpha = 0.05){
  sig_t <- tumor  %>% filter(bp_value < alpha)  %>% pull(annotation)
  sig_n <- normal %>% filter(bp_value < alpha)  %>% pull(annotation)
  list(
    tumor_only  = setdiff(sig_t, sig_n),
    normal_only = setdiff(sig_n, sig_t)
  )
}
sets <- get_sig_sets(tumor, normal, alpha = 0.05)

df_plot <- bind_rows(
  tumor  %>% filter(annotation %in% sets$tumor_only)  %>% mutate(side = "Tumor"),
  normal %>% filter(annotation %in% sets$normal_only) %>% mutate(side = "Normal")
) %>%
  transmute(
    side,
    annotation,
    beta = bias_corrected_estimate,
    lo   = CI_lo,
    hi   = CI_hi,
    p    = pmax(bp_value, 1e-300),    
    neglogp = -log10(p)
  ) %>%
  group_by(side) %>%             
  arrange(desc(abs(beta)), .by_group = TRUE) %>%
  mutate(y = factor(annotation, levels = rev(unique(annotation)))) %>%
  ungroup() %>%
  mutate(side = factor(side, levels = c("Tumor","Normal")))  


p1 <- ggplot(df_plot, aes(x = beta, y = y, color = neglogp)) +
  facet_grid(rows = vars(side), scales = "free_y", space = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_segment(aes(x = lo, xend = hi, yend = y), linewidth = 2.1, lineend = "round") +
  geom_point(size = 3, stroke = 0.3, color = "black", show.legend = FALSE,shape=23,fill = "#FFC75F") +
  scale_color_gradientn(
    colours = c("#FF9671", "#FF6F91", "#845EC2"),
    limits  = c(1.5, 3.5),
    breaks  = c(1.5, 2.5, 3.5),
    oob     = scales::squish,
    name    = expression(-log[10](p))
  ) +
  labs(
    x = "Bias Corrected Estimate (95% CI)",  
    y = NULL
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 16),         
    panel.border = element_rect(color="black", fill=NA, linewidth=0.4), 
    axis.line    = element_line(linewidth=0.1),  
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y.left   = element_text(hjust = 1, margin = margin(r = 5)),
    axis.title.x   = element_text(size = 14),
    panel.spacing.y    = unit(6, "pt"),
    strip.background   = element_rect(fill = "white", colour = NA),
    strip.text.y       = element_text(face = "bold",size=15),
    legend.position    = c(0.85, 0.45),
    legend.key.size = unit(0.4, "cm"),
    legend.text     = element_text(size = 13),
    legend.title    = element_text(size=14)
  )+
  coord_cartesian(clip = "off")

p1



#### 


pbmc<-readRDS("Rename_am.rds")
pbmc<-readRDS("Rename_am_normal.rds")

cell_type_counts <- pbmc@meta.data %>%
  group_by(cell_type) %>%
  summarize(count = n()) %>%
  ungroup()
total_cells <- sum(cell_type_counts$count)
cell_type_counts <- cell_type_counts %>%
  mutate(proportion = count / total_cells * 100)

am_proportion_data <- cell_type_counts

cell_colors <- c(brewer.pal(8, "Set2"), brewer.pal(9, "Set3"), brewer.pal(9, "Set1"), brewer.pal(8, "Accent"))

am_proportion_data <- am_proportion_data %>%
  arrange(desc(proportion))

df_tumor<-am_proportion_data
df_normal<-am_proportion_data

tumor_palette  <- readRDS("am_tumor_palette.rds")
normal_palette <- readRDS("am_normal_palette.rds")

piecewise_perc_trans <- function(cap = 15, tail_end = 30, tail_scale = 0.2) {
  stopifnot(tail_end > cap, tail_scale > 0, tail_scale < 1)
  trans_new(
    name = sprintf("piecewise_cap%d_tail%d_s%.2f", cap, tail_end, tail_scale),
    transform = function(x) {
      out <- rep(NA_real_, length(x))
      ok  <- is.finite(x)
      if (any(ok)) {
        idx  <- which(ok)
        x2   <- pmin(x[idx], tail_end)
        bel  <- x2 <= cap
        out[idx[bel]]  <- x2[bel]
        out[idx[!bel]] <- cap + tail_scale * (x2[!bel] - cap)
      }
      out
    },
    inverse = function(y) {
      out <- rep(NA_real_, length(y))
      ok  <- is.finite(y)
      if (any(ok)) {
        idx  <- which(ok)
        y2   <- y[idx]
        bel  <- y2 <= cap
        out[idx[bel]]  <- y2[bel]
        out[idx[!bel]] <- cap + (y2[!bel] - cap) / tail_scale
        out[idx]       <- pmin(out[idx], tail_end)
      }
      out
    }
  )
}


make_bar <- function(df, title = "", axis_side = c("left","right"),
                     palette = NULL, x_max = 30,
                     breaks_show = c(5, 15, 30),
                     cap = 15, tail_end = 30, tail_scale = 0.2) {
  axis_side <- match.arg(axis_side)
  if (!is.null(palette)) stopifnot(all(levels(df$cell_type) %in% names(palette)))
  
  df <- df %>%
    filter(!is.na(proportion)) %>%
    mutate(cell_type = factor(cell_type, levels = balanced_order(df)))  
  
  ggplot(df, aes(x = proportion, y = cell_type, fill = cell_type)) +
    geom_col(width = 0.70, color = "black", linewidth = 0.2, na.rm = TRUE) +
    {
      if (!is.null(palette)) {
        # 关键：用 df 的 levels 顺序来取色，确保严格一一对应
        scale_fill_manual(values = palette[levels(df$cell_type)], guide = "none")
      } else {
        guides(fill = "none")
      }
    } +
    scale_x_continuous(
      limits = c(0, x_max),
      breaks = c(5, 15),
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0),
      trans  = piecewise_perc_trans(cap = cap, tail_end = tail_end, tail_scale = tail_scale)
    ) +
    { if (axis_side == "right") scale_y_discrete(position = "right") else NULL } +
    labs(x = "", y = NULL, title = title) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(color = "black",size = 12),
      axis.text.x        = element_text(color = "black",,size = 12),
      plot.title         = element_text(hjust = 0.5, face = "bold",size = 16),
      plot.margin        = margin(4, 6, 4, 6)
    )
}

xlim_shared <- 30
p_tumor  <- make_bar(df_tumor,  "Tumor",  "left",
                     palette = tumor_palette,
                     x_max = xlim_shared, breaks_show = c(5, 15, 30),
                     cap = 15, tail_end = 30, tail_scale = 0.2)

p_normal <- make_bar(df_normal, "Normal", "right",
                     palette = normal_palette,
                     x_max = xlim_shared, breaks_show = c(5, 15, 30),
                     cap = 15, tail_end = 30, tail_scale = 0.2)

p_final <- p_tumor + p_normal + patchwork::plot_layout(widths = c(1, 1))
p_final

genes_of_interest <- c("FABP4",  
                       "MT1E", 
                       "RNASEK",
                       "CD4", 
                       "APOC2", 
                       "CCL20",
                       "SPP1", 
                       "SPIDR", 
                       "IGHA1",
                       "RAD51B", 
                       "HLA-DRB5","HLA-DQA2","HLA-DQB2",
                       "VCAN")
pdf("am_marker_dotplot.pdf",width = 18,height = 6)
DotPlot(pbmc,features = genes_of_interest)
dev.off()

#########################
sc_dataset <- readRDS("Rename_am.rds")
sc_dataset$cell_type2<-substring(sc_dataset$cell_type,1,4)
tumor_colors<-readRDS("am_tumor_palette.rds")

emb <- as.data.frame(Embeddings(sc_dataset, "umap"))
colnames(emb)[1:2] <- c("UMAP_1","UMAP_2")
df_plot <- emb %>%
  mutate(
    ident      = as.character(Idents(sc_dataset)),
    cell_type2 = as.character(sc_dataset$cell_type2)
  )

am_levels <- c("AM-1","AM-2","AM-3","AM-4","AM-5")
df_am <- df_plot %>% filter(cell_type2 %in% am_levels)

am_labels <- df_am %>%
  group_by(cell_type2) %>%
  summarise(cx = mean(UMAP_1), cy = mean(UMAP_2), .groups = "drop")

xr <- range(df_plot$UMAP_1, na.rm = TRUE)
yr <- range(df_plot$UMAP_2, na.rm = TRUE)
dx <- diff(xr); dy <- diff(yr)
x0 <- xr[1]; y0 <- yr[1]
ax_len <- 0.12 * min(dx, dy)  # 两轴用同一长度
hx <- ax_len; hy <- ax_len

p <- ggplot(df_plot, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = ident), size = 0.6, alpha = 0.9, show.legend = FALSE) +
  scale_color_manual(values = tumor_colors, guide = "none") +
  stat_ellipse(
    data = df_am,
    aes(group = cell_type2),
    type = "norm", level = 0.9,
    linetype = "dashed", linewidth = 0.35, color = "black",
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = am_labels,
    aes(x = cx, y = cy, label = cell_type2),
    inherit.aes = FALSE,
    size = 6, fontface = "bold",
    box.padding = 0.3, point.padding = 0.2,
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
  annotate("text", x = x0 - dx*0.03, y = y0 + hy*0.55, label = "UMAP2", size = 5, angle = 90, vjust = 0.5) +
  
  coord_fixed(ratio = 1) + 
  theme_void() +
  theme(legend.position = "none")


pdf("am_tumor_UMAP.pdf")
p
dev.off()


df_legend <- data.frame(ident = unique(df_plot$ident))
p_legend <- ggplot(df_legend, aes(x = 1, color = ident)) +
  geom_point(aes(y = 1), size = 6) +
  scale_color_manual(values = tumor_colors, name = "") +
  theme_void() + theme(legend.position = "right",legend.text = element_text(size=13))
legend_only <- cowplot::get_legend(p_legend)
pdf("am_tumor_umap_legend.pdf")
grid::grid.newpage(); grid::grid.draw(legend_only)
dev.off()


library(Seurat)
library(ggplot2)


genes_of_interest <- c("FABP4", "MT1E", "APOC2", "RNASEK", "CD4", "CCL20", "CXCL2","CCL4", "CCL3", "CCL18",
                       "SPP1", "SPIDR", "RAD51B", "IGHA1", "HLA-DRB5","HLA-DQA2", "VCAN")
sc_dataset$cell_type<- factor(
  sub("^(AM-\\d+).*", "\\1", sc_dataset$cell_type),
  levels = paste0("AM-", 1:5)
)

custom_markers <- c(
  "AM-1.FABP4"        = "AM-1",
  "AM-1.MT1E"         = "AM-1",
  "AM-1.APOC2- RNASEK+" = "AM-1",
  "AM-1.APOC2+ RNASEK+" = "AM-1",
  "AM-1.HLA"          = "AM-1",
  "AM-2.CD4"          = "AM-2",
  "AM-2.CK"           = "AM-2",
  "AM-3.SPP1"         = "AM-3",
  "AM-4.SPIDR+ RAD51B-" = "AM-4",
  "AM-4.SPIDR- RAD51B+" = "AM-4",
  "AM-4.SPIDR+ RAD51B+" = "AM-4",
  "AM-5.IGHA1"        = "AM-5",
  "AM-5.VCAN"         = "AM-5"
)


custom_colors <- c(
  "AM-1" = "#7AC7C9",
  "AM-2" = "#FBBE09",
  "AM-3" = "#E3A7C5",
  "AM-4" = "#91C79E",
  "AM-5" = "#DD7572"
)

BetterVis_DotPlot_Heatmap_Seurat(sc_dataset,
                                 features = genes_of_interest,cluster.features=F,cluster.idents = F,rotate.axis = F,idents.kmeans = 1,features.kmeans = 1,
                                 show.data.legend = F,show.idents.legend = F,show.split.legend = F,row.names.side = "right",column.names.side = "top",
                                 gene_markers = custom_markers ,marker_colors = custom_colors,annotation_marker_legend = F,
                                 row.split =c(1,1,1,1,1,2,2,3,4,4,4,5,5),
                                 col.split =c(1,1,1,1,2,2,2,2,2,2,3,4,4,5,1,1,5),
                                 data.colors = c("#1289B2","#FCF6AF","#E53849"),dots.size = 4,dots.type = "square root",inner.border = F,border_gp_control = F,
                                 features.names.size = 12,idents.names.size = 12,keep_marker_order=F
)


#########################
twasgene <- readRDS("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am/twasgene.rds")
gene <- twasgene$SYMBOL

am <- readRDS("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am/Rename_am.rds")
am_normal <- readRDS("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am_normal/Rename_am_normal.rds")


am_cells <- subset(am, idents = c("AM-5.IGHA1"))
am_normal_cells <- subset(am_normal, idents = c("AM-2.IGHA1"))

am_expr <- as.data.frame(GetAssayData(am_cells, assay = "RNA"))
am_normal_expr <- as.data.frame(GetAssayData(am_normal_cells, assay = "RNA"))

am_expr<- am_expr[rownames(am_expr) %in% gene, ]
am_normal_expr<- am_normal_expr[rownames(am_normal_expr) %in% gene, ]

combined_expr <- cbind(am_expr, am_normal_expr)

cell_labels <- c(rep("am_tumor", ncol(am_expr)), rep("am_normal", ncol(am_normal_expr)))
cell_info <- data.frame(Cell = colnames(combined_expr), Group = cell_labels)

cell_info <- cell_info %>%column_to_rownames(var = "Cell")

saveRDS(cell_info,"/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am/GRSA/cell_info.rds")
saveRDS(combined_expr,"/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am/GRSA/combined_expr.rds")



cell_info<-readRDS("cell_info.rds")
combined_expr<-readRDS("combined_expr.rds")
reporter_res_gene <- GRSA(combined_expr, "Group", cell_info,
                          mode = "directed",
                          method = "wilcox.test", perm = 999,
                          type = "hsa", feature = "gene"
)


p1 <- plot_report_bar(reporter_res_gene, rs_threshold = 1)
p1

reporter_res_gene2 <- modify_description(reporter_res_gene, pattern = " - Homo sapiens (human)")
p2 <- plot_report_bar(reporter_res_gene2, rs_threshold = 1)
p2

rs <- reporter_res_gene2$reporter_s
vs <- attr(rs, "vs_group")
g1 <- "AM Normal"  # am_normal
g2 <- "AM Tumor"   # am_tumor



must_show <- c(
  "Phagosome",
  "Endocytosis",
  "Cell adhesion molecules",
  "Th1 and Th2 cell differentiation",
  "Th17 cell differentiation",
  "Antigen processing and presentation",
  "Natural killer cell mediated cytotoxicity",
)


df <- rs %>%
  filter(!is.na(ReporterScore)) %>%  
  filter(Description%in%must_show) %>% 
  transmute(
    Description,
    ReporterScore,
    Group = if_else(ReporterScore >= 0,
                    paste("Enriched in", g2),
                    paste("Enriched in", g1))
  )
df$Description[df$Description%in%"Intestinal immune network for IgA production"]<-"Intestinal IgA production"
df$Description[df$Description%in%"Kaposi sarcoma-associated herpesvirus infection"]<-"Kaposi sarcoma-associated herpesvirus"

topN <- 8
df_use <- bind_rows(
  df %>% group_by(Group) %>%
    slice_max(order_by = abs(ReporterScore), n = topN, with_ties = FALSE) %>%
    ungroup(),
  df %>% filter(Description %in% must_show)
) %>%
  distinct(Description, .keep_all = TRUE) %>%
  mutate(
    is_key = Description %in% must_show,
    label  = if_else(is_key, paste0("★ ", Description), Description),
    label  = str_wrap(label, 38)
  ) %>%
  arrange(ReporterScore) %>%
  mutate(label = fct_inorder(label))

cols <- setNames(
  c("#4DBBD5FF", "#E64B35FF"),  
  c(paste("Enriched in", g1), paste("Enriched in", g2))
)
thr <- 1


p <- ggplot(df_use, aes(x = ReporterScore, y = label, color = Group)) +
  geom_segment(aes(x = 0, xend = ReporterScore, yend = label),
               linewidth = 2, lineend = "round", color = "grey82") +
  geom_point(size = 5) +
  geom_vline(xintercept = c(-thr, thr), linetype = "dashed", color = "grey70") +
  scale_color_manual(values = cols, name = NULL) +  
  labs(x = "ReporterScore", y = NULL) +
  theme_light(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 14, colour = "black", face = "bold"),
    plot.title = element_blank(),  
    legend.position = c(0.95, 0.05),   
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", colour = "black", size = 0.3,linetype = "dashed"), 
    legend.key.size = unit(0.5, "cm"),
    legend.text  = element_text(colour = "black",size = 11)
  )
p

ggsave("AM_GRSA.pdf", plot = p, device = cairo_pdf, width = 7, height =5)

