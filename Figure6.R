### Figure 6

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(colorRamp2)
library(ggrepel)
library(ggforce)
library(grid)
library(future)
library(Seurat)
library(igraph)
library(ggalluvial)
library(patchwork)
library(CellChat)


data<-read.csv("goBP_am_tumor_top10.csv")

data_proc <- data %>%
  mutate(score = -log10(p.adjust)) %>%
  select(Description, cell_type, score)

go_mat <- data_proc %>%
  pivot_wider(names_from = cell_type, values_from = score, values_fill = 0)

mat <- as.data.frame(go_mat)
rownames(mat) <- mat$Description
mat <- mat[,-1]

pca <- prcomp(t(mat), scale.=TRUE)
pca_df <- as.data.frame(pca$x[, 1:2])
pca_df$cell_type <- rownames(pca_df)

col_palette <- c(
  "AM-1.FABP4"       = "#E64B35FF",
  "AM-5.IGHA1"       = "#8491B4FF",
  "AM-4.SPIDR+ RAD51B-" = "#4DBBD5FF",
  "AM-4.SPIDR+ RAD51B+" = "#4DBBD5FF",
  "AM-4.SPIDR- RAD51B+" = "#4DBBD5FF",
  "AM-1.MT1E"        = "#E64B35FF",
  "AM-1.APOC2- RNASEK+" = "#E64B35FF",
  "AM-1.APOC2+ RNASEK+" = "#E64B35FF",
  "AM-1.HLA"         = "#E64B35FF",
  "AM-5.VCAN"        = "#8491B4FF"
)

pca_df$group <- dplyr::case_when(
  grepl("AM-1", pca_df$cell_type) ~ "AM-1",
  grepl("AM-5", pca_df$cell_type) ~ "AM-5",
  grepl("AM-4", pca_df$cell_type) ~ "AM-4",
  TRUE ~ NA_character_
)

p <- ggplot(subset(pca_df, !is.na(group)), aes(PC1, PC2)) +
  geom_mark_ellipse(aes(fill = group, label = group),
                    alpha = 0.08, color = NA,
                    expand = unit(4, "mm"),
                    label.fontsize = 14, label.fontface = "bold",
                    label.fill = NA) +
  geom_segment(data = data.frame(
    x = c(mean(pca_df$PC1[pca_df$group=="AM-1"]),
          mean(pca_df$PC1[pca_df$group=="AM-5"])),
    y = c(mean(pca_df$PC2[pca_df$group=="AM-1"]),
          mean(pca_df$PC2[pca_df$group=="AM-5"])),
    xend = c(mean(pca_df$PC1[pca_df$group=="AM-5"]),
             mean(pca_df$PC1[pca_df$group=="AM-4"])),
    yend = c(mean(pca_df$PC2[pca_df$group=="AM-5"]),
             mean(pca_df$PC2[pca_df$group=="AM-4"]))
  ),
  aes(x = x, y = y, xend = xend, yend = yend),
  color = "#7A7A7A", linewidth = 1, alpha = 0.6) +
  geom_point(aes(color = cell_type), size = 5, alpha = 0.9) +
  geom_text_repel(data = subset(pca_df, group %in% c("AM-1","AM-5","AM-4")),
                  aes(label = cell_type, color = cell_type),
                  size = 4.5, fontface = "bold",
                  box.padding = 0.4, point.padding = 0.4,
                  max.overlaps = Inf, segment.color = "grey80") +
  scale_color_manual(values = col_palette) +
  scale_fill_manual(values = c("AM-1"="#E64B35FF",
                               "AM-5"="#8491B4FF",
                               "AM-4"="#4DBBD5FF")) +
  
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title = element_text(face="bold", size=16),
    axis.text = element_text(size=16, color="black"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size=1, color="black")
  )


p 

################
am_prefix <- function(v) sub("^(AM-[0-9]+).*", "\\1", v)
RA <- regulonActivity_byGroup_Scaled
am_group <- am_prefix(colnames(RA))
regulonActivity_AM <- sapply(unique(am_group), function(am){
  cols <- which(am_group == am)
  rowMeans(RA[, cols, drop = FALSE], na.rm = TRUE)
})
regulonActivity_AM <- as.data.frame(regulonActivity_AM)
zRA <- t(scale(t(as.matrix(regulonActivity_AM)))); zRA[!is.finite(zRA)] <- 0

species <- "Homo sapiens"
msig_h <- msigdbr(species = species, category = "H")

G_imm_clean <- msig_h %>%
  filter(gs_name %in% c("HALLMARK_INTERFERON_GAMMA_RESPONSE",
                        "HALLMARK_TNFA_SIGNALING_VIA_NFKB")) %>%
  group_by(gs_name) %>%
  summarise(genes = list(unique(gene_symbol)), .groups = "drop")

G_met_clean <- msig_h %>%
  filter(gs_name %in% c("HALLMARK_GLYCOLYSIS",
                        "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                        "HALLMARK_FATTY_ACID_METABOLISM")) %>%
  group_by(gs_name) %>%
  summarise(genes = list(unique(gene_symbol)), .groups = "drop")

fisher_one <- function(targets, geneset){
  a <- length(intersect(targets, geneset))
  b <- length(setdiff(geneset, targets))
  c <- length(setdiff(targets, geneset))
  d <- length(setdiff(universe, union(targets, geneset)))
  ft <- suppressWarnings(fisher.test(matrix(c(a,b,c,d),2), alternative="greater"))
  c(p = ft$p.value, OR = unname(ft$estimate))
}

padj <- function(p) p.adjust(p, "BH")

score_regulon <- function(gs_tbl){
  do.call(rbind, lapply(names(reg_targets), function(tf){
    tg <- reg_targets[[tf]]
    res <- lapply(gs_tbl$genes, function(gs) fisher_one(tg, toupper(gs)))
    res <- do.call(rbind, res); colnames(res) <- c("p","OR")
    data.frame(regulon=tf, gs=gs_tbl$gs_name, p=res[,"p"], OR=res[,"OR"])
  }))
}

pick_top <- function(E, K=25){
  E %>% group_by(regulon) %>%
    summarise(FDR = padj(p) |> min(na.rm=TRUE),
              OR = median(OR, na.rm=TRUE), .groups="drop") %>%
    filter(is.finite(FDR), FDR<0.05, OR>1) %>%
    mutate(score = (-log10(FDR)) * log2(OR)) %>%
    arrange(desc(score)) %>% slice_head(n=K) %>% pull(regulon)
}

reg_targets <- lapply(regulons, function(v) unique(toupper(v)))
universe <- unique(unlist(reg_targets)); universe <- universe[!is.na(universe)]

E_imm <- score_regulon(G_imm_clean)
E_met <- score_regulon(G_met_clean)

R_imm <- intersect(pick_top(E_imm, 25), rownames(zRA))
R_met <- intersect(pick_top(E_met, 25), rownames(zRA))

axis_scores_mph <- rbind(
  Immune    = colMeans(zRA[R_imm, , drop=FALSE]),
  Metabolic = colMeans(zRA[R_met, , drop=FALSE])
)
round(axis_scores_mph, 3)

msig_go <- msigdbr(species = species, category = "C5", subcategory = "GO:BP")
pat_rep <- c("MACROPHAGE_ACTIVATION", "MACROPHAGE_MIGRATION", "MACROPHAGE_CHEMOTAXIS",
             "MACROPHAGE_DIFFERENTIATION", "PHAGOCYTOSIS", "EFFEROCYTOSIS",
             "APOPTOTIC_CELL_CLEARANCE", "WOUND_HEALING", "RESPONSE_TO_WOUNDING",
             "TISSUE_REGENERATION", "TISSUE_REMODELING", "REMODELING", "SCAR",
             "EXTRACELLULAR_MATRIX_ORGANIZATION", "COLLAGEN_FIBRIL_ORGANIZATION",
             "ELASTIC_FIBER", "ANGIOGENESIS", "VASCULATURE_DEVELOPMENT",
             "BLOOD_VESSEL_MORPHOGENESIS", "FIBROBLAST_PROLIFERATION",
             "REGULATION_OF_FIBROBLAST_PROLIFERATION", "TGF_BETA")

G_rep_mph <- msig_go %>%
  filter(str_detect(gs_name, str_c(pat_rep, collapse="|"))) %>%
  filter(!str_detect(gs_name, "EPITHELIAL_MESENCHYMAL_TRANSITION|EMT")) %>%
  group_by(gs_name) %>% summarise(genes=list(unique(gene_symbol)), .groups="drop")

E_repM <- score_regulon(G_rep_mph)
R_repM <- intersect(pick_top(E_repM, 25), rownames(zRA))
Repair_new <- colMeans(zRA[R_repM, , drop=FALSE], na.rm=TRUE)

axis_final <- rbind(
  Immune    = axis_scores_mph["Immune", ],
  Metabolic = axis_scores_mph["Metabolic", ],
  Repair    = Repair_new
)
round(axis_final, 3)

AXp <- t(axis_final) %>% as.data.frame()
AXp$AM <- rownames(AXp)
colnames(AXp)[1:3] <- c("Immune","Metabolic","Repair")

min_val <- min(as.matrix(AXp[, 1:3]), na.rm = TRUE)
shift_extra <- 0.2 * abs(min_val)
AXp[, 1:3] <- AXp[, 1:3] - min_val + shift_extra

pal <- MetBrewer::met.brewer("Hiroshige", 5)
names(pal) <- AXp$AM


pt_size <- 6
label_size <- 5
base_font <- 16

p1 <- ggplot(AXp, aes(x = Immune, y = Metabolic, color = AM)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "#C4D7E0", alpha = 0.15) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "#E6EAD0", alpha = 0.15) +
  geom_point(size = pt_size) +
  geom_text_repel(aes(label = AM), size = label_size, fontface = "bold",
                  box.padding = 0.4, max.overlaps = 100, segment.color = "grey60") +
  scale_color_manual(values = pal) +
  labs(
    title = "",
    x = "Immune Regulon Score",
    y = "Metabolic Regulon Score"
  ) +
  theme_classic(base_size = base_font) +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black",size=15)
  )

p2 <- ggplot(AXp, aes(x = Immune, y = Repair, color = AM)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "#D8E2DC", alpha = 0.15) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "#F8EDEB", alpha = 0.15) +
  geom_point(size = pt_size) +
  geom_text_repel(aes(label = AM), size = label_size, fontface = "bold",
                  box.padding = 0.4, max.overlaps = 100, segment.color = "grey60") +
  scale_color_manual(values = pal) +
  labs(
    title = "",
    x = "Immune Regulon Score",
    y = "Repair Regulon Score"
  ) +
  theme_classic(base_size = base_font) +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black",size=15)
  )

p_final <- p1 + p2 +
  plot_annotation(title = NULL) & 
  theme(
    plot.margin = margin(5, 5, 5, 5),       
    plot.title = element_blank()           
  )

p_final
ggsave("AM5_regulon_score.pdf",p_final,width = 8,height = 4.5)

###############################
suppressMessages(options(stringsAsFactors = FALSE))
suppressMessages(options(futrue.globlas.Maxsize=2*1024**3))

groupSize<-readRDS("groupSize.rds")
cellchat <- readRDS("cellchat_am_tumor_computeCommunProb.rds")
cellchat <- filterCommunication(cellchat , min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
df.net <- subsetCommunication(cellchat)

am5_comm <- subsetCommunication(cellchat,
                                sources.use = "AM-5.IGHA1",
                                targets.use = NULL)

am5_comm_rev <- subsetCommunication(cellchat,
                                    sources.use = NULL,
                                    targets.use = "AM-5.IGHA1")
am5_all <- bind_rows(
  am5_comm %>% mutate(direction = "Outgoing"),
  am5_comm_rev %>% mutate(direction = "Incoming")
)


am5_summary <- am5_all %>%
  group_by(direction, target, source) %>%
  summarise(mean_prob = mean(prob), n_pair = n(), .groups="drop") %>%
  mutate(partner = ifelse(direction=="Outgoing", target, source)) %>%
  arrange(desc(mean_prob))


top_pairs <- am5_all %>%
  group_by(interaction_name_2) %>%
  summarise(mean_prob = mean(prob)) %>%
  arrange(desc(mean_prob)) %>%
  slice_head(n = 10)
top_pairs

path_colors <- c(
  "ECM_adhesion"   = "#F1A340",   
  "Antigen_presentation" = "#92C5DE", 
  "Inflammatory"   = "#B2182B",   
  "Repair_Regulation" = "#4DAF4A"
)

top_pairs$module <- dplyr::case_when(
  grepl("FN1|COL1A", top_pairs$interaction_name_2) ~ "ECM_adhesion",
  grepl("HLA", top_pairs$interaction_name_2) ~ "Antigen_presentation",
  grepl("MIF|ANXA1|RETN", top_pairs$interaction_name_2) ~ "Inflammatory",
  grepl("GRN|PTPRC", top_pairs$interaction_name_2) ~ "Repair_Regulation",
  TRUE ~ "Other"
)

ggplot(top_pairs, aes(x = reorder(interaction_name_2, mean_prob),
                      y = mean_prob, fill = module)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = path_colors) +
  labs(x = NULL, y = "Mean communication probability",
       fill = "Pathway category") +
  theme_classic(base_size = 14) +
  theme(axis.text.y = element_text(size = 10, face = "bold"),
        legend.position = "right")

am5_out <- subsetCommunication(cellchat, sources.use = "AM-5.IGHA1", targets.use = NULL) %>%
  mutate(direction = "Outgoing", partner = target)

am5_in  <- subsetCommunication(cellchat, sources.use = NULL, targets.use = "AM-5.IGHA1") %>%
  mutate(direction = "Incoming", partner = source)

am5_all <- bind_rows(am5_out, am5_in)
am5_all <- am5_all %>% filter(!is.na(partner))
am5_all <- am5_all %>%
  semi_join(top_pairs, by = "interaction_name_2")

module_map <- function(x){
  case_when(
    grepl("FN1|COL1A", x, ignore.case = TRUE) ~ "ECM_adhesion",
    grepl("^HLA-| HLA-", x, ignore.case = TRUE) ~ "Antigen_presentation",
    grepl("MIF|ANXA1|RETN", x, ignore.case = TRUE) ~ "Inflammatory",
    grepl("GRN|PTPRC|MRC1|SORT1", x, ignore.case = TRUE) ~ "Repair_Regulation",
    TRUE ~ "Other"
  )
}

top_pairs2 <- top_pairs %>%
  mutate(module = module_map(interaction_name_2))

top_pairs2<-as.data.frame(top_pairs2)

am5_all <- am5_all %>%
  left_join(top_pairs2 %>% dplyr::select(interaction_name_2, module), by = "interaction_name_2")

am5_partner_pair <- am5_all %>%
  group_by(direction, partner, interaction_name_2, module) %>%
  summarise(mean_prob = mean(prob, na.rm = TRUE), .groups = "drop")

am5_partner_module <- am5_partner_pair %>%
  group_by(direction, partner, module) %>%
  summarise(prob_sum = sum(mean_prob), .groups = "drop") %>%
  group_by(direction, partner) %>%
  mutate(frac = prob_sum / sum(prob_sum)) %>%
  ungroup()

partner_order <- c("AM-1.FABP4","AM-1.MT1E","AM-1.APOC2- RNASEK+","AM-1.APOC2+ RNASEK+","AM-1.HLA",
                   "AM-5.IGHA1","AM-5.VCAN",
                   "AM-4.SPIDR+ RAD51B-","AM-4.SPIDR+ RAD51B+","AM-4.SPIDR- RAD51B+")


col_in  <- "#AFCDF2"
col_out <- "#E6B2C6"
partner_colors <- c(
  "#F1A340", "#92C5DE", "#B2182B", "#4DAF4A",
  "#E69F00", "#56B4E9", "#009E73", "#D55E00",
  "#CC79A7", "#0072B2", "#999999", "#A65628", "#984EA3"
)
path_colors <- c(
  "ECM_adhesion" = "#F1A340",
  "Antigen_presentation" = "#92C5DE",
  "Inflammatory" = "#C9C7E3",
  "Repair_Regulation" = "#98C570",
  "Other" = "grey75"
)

df_in  <- subset(am5_partner_pair, direction == "Incoming")
df_out <- subset(am5_partner_pair, direction == "Outgoing")

annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) {
  layer(
    data = data, stat = "identity", position = "identity",
    geom = ggplot2::GeomCustomAnn,
    inherit.aes = FALSE,
    params = list(grob = grob, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  )
}

p <- ggplot(am5_partner_pair,
            aes(x = partner, y = interaction_name_2)) +
  geom_point(aes(size = mean_prob, fill = module),
             shape = 21, color = "black", stroke = 0.2) +
  scale_size_continuous(name = "Mean prob", range = c(2, 10)) +
  scale_fill_manual(values = path_colors, name = "Module") +
  facet_grid(rows = vars(direction), switch = "y") +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1, size = 13, colour = "black", face = "plain"),
    axis.text.y = element_text(size = 13, colour = "black", face = "plain"),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.3, "lines"),
    panel.grid = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 70)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 7, stroke = 0.3))
  )

p <- p +
  annotation_custom2(
    grid::rectGrob(gp = gpar(fill = col_in, col = NA)),
    xmin = -3.8, xmax = -3.3, ymin = -Inf, ymax = Inf, data = df_in
  ) +
  annotation_custom2(
    grid::textGrob("Incoming", rot = 90, y = unit(0.5, "npc"),
                   gp = gpar(fontsize = 15, fontface = "bold", col = "black")),
    xmin = -3.65, xmax = -3.55, ymin = -Inf, ymax = Inf, data = df_in
  ) +
  annotation_custom2(
    grid::rectGrob(gp = gpar(fill = col_out, col = NA)),
    xmin = -3.8, xmax = -3.3, ymin = -Inf, ymax = Inf, data = df_out
  ) +
  annotation_custom2(
    grid::textGrob("Outgoing", rot = 90, y = unit(0.5, "npc"),
                   gp = gpar(fontsize = 15, fontface = "bold", col = "black")),
    xmin = -3.65, xmax = -3.55, ymin = -Inf, ymax = Inf, data = df_out
  )

partners <- levels(factor(am5_partner_pair$partner))
p <- p +
  annotation_custom2(
    grid::grobTree(
      grid::rectGrob(
        x = seq(0.05, 0.95, length.out = length(partners)),
        y = unit(0, "npc"),
        width = unit(1.05 * 0.9 / length(partners), "npc"),
        height = unit(0.03, "npc"),
        gp = grid::gpar(fill = partner_colors[seq_along(partners)], col = NA)
      )
    ),
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    data = df_out
  )

y_levels <- levels(factor(am5_partner_pair$interaction_name_2))
y_band_df <- expand.grid(
  direction = c("Incoming", "Outgoing"),
  interaction_name_2 = y_levels,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
) %>%
  mutate(
    y_idx = as.integer(factor(interaction_name_2, levels = y_levels)),
    ymin  = y_idx - 0.45,
    ymax  = y_idx + 0.45,
    fill_y = partner_colors[(y_idx - 1) %% length(partner_colors) + 1]
  )

p <- p +
  geom_rect(
    data = y_band_df,
    aes(ymin = ymin, ymax = ymax, fill = I(fill_y)),
    xmin = 0.3, xmax = 0.5,
    inherit.aes = FALSE, show.legend = FALSE
  )
p
ggsave("am5_bubble_cellchat.pdf", p, width = 12, height = 8)



############
load("/GPUFS/gyfyy_jxhe_1/User/zyt/scPagwas/scPagwas_lc/am/am_scPagwas.RData")
scPagwas_Gene<-Pagwas@assays$RNA@var.features
twasgene<-readRDS("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am/twasgene.rds")
twasgene$TWAS_Gene<-twasgene$SYMBOL
scPagwas_Gene<-data.frame("SYMBOL"=scPagwas_Gene)
scPagwas_Gene$scPagwas_Gene<-scPagwas_Gene$SYMBOL
venndata<-full_join(scPagwas_Gene,twasgene,by="SYMBOL")
venndata <- venndata %>%column_to_rownames("SYMBOL") %>%dplyr::select(-ENSEMBL)
venndata$TWAS_Gene <- ifelse(is.na(venndata$TWAS_Gene), 0, 1)
venndata$scPagwas_Gene <- ifelse(is.na(venndata$scPagwas_Gene), 0, 1)
intersectgene<-c("IGFBP1","CH25H","FTCD","NOTCH4","PSORS1C1")

library(UpSetR)

pdf("UpsetR_TWAS_scPagwas_intersect.pdf",width = 8)  # 缩小 PDF 大小
p <- upset(venndata, 
           nset = 7, 
           nintersects = 10, 
           order.by = c('degree', 'freq'), 
           decreasing = c(TRUE, TRUE),
           mb.ratio = c(0.7, 0.3),
           point.size = 5,    # 增大点的大小
           line.size = 3,     # 增大连线大小
           set_size.numbers_size = 10, # 增大集合大小的字体
           mainbar.y.label = "Intersection size", 
           sets.x.label = "", 
           main.bar.color = "#F3BCBF", 
           sets.bar.color = "#BCDCF3",
           text.scale = c(2.2,2.2,2.2,1.7,1.7,1.7), # 统一放大所有字体
           queries = list(list(query = intersects, 
                               params = list("scPagwas_Gene", "TWAS_Gene"), 
                               active = TRUE, 
                               color = "#d66a35", 
                               query.name = "scPagwas_Gene vs TWAS_Gene")))

p
grid.text(paste(intersectgene, collapse = "\n"), 
          x = unit(0.8, "npc"),  # 调整 x 位置
          y = unit(0.5, "npc"),  # 调整 y 位置
          just = "left",         # 左对齐
          gp = gpar(fontsize = 12, col = "#F3BCBF"))  # 增大字体至 12

dev.off()

###############
#################### Dotplot scPagwas_am intersect twasgene
seurat_object<-readRDS("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am/Rename_am.rds")
load("/GPUFS/gyfyy_jxhe_1/User/zyt/scPagwas/scPagwas_lc/am/am_scPagwas.RData")
intersectgene<-Pagwas@assays$RNA@var.features
intersectgene<-c("IGFBP1","CH25H","FTCD","NOTCH4","PSORS1C1")
dot_plot <- DotPlot(seurat_object, features = intersectgene)

dot_plot_data <- dot_plot$data
dot_plot_data <- dot_plot_data %>%
  mutate(size_category = cut(pct.exp, breaks = c(0, 1, 2, 10, 40), 
                             labels = c("0.25", "0.50", "0.75", "1.0"), include.lowest = TRUE))
custom_dot_plot <- ggplot(dot_plot_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = size_category, color = avg.exp.scaled)) +
  scale_color_gradient(low = "lightgrey", high = "blue", limits = c(-1.1, 2.5)) +
  scale_size_manual(values = c("0.25" = 2, "0.50" = 4, "0.75" = 6, "1.0" = 8)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  labs(size = "Percent Expressed", color = "Average Expression")

print(custom_dot_plot)




###########
##########
pathway_scores<-as.data.frame(pathway_scores)
DefaultAssay(seurat_object) <- "RNA"
ch25h_expr <- FetchData(seurat_object, vars = "CH25H")
ch25h_expr$cell <- rownames(ch25h_expr)

pathway_scores$cell <- rownames(pathway_scores)
df <- inner_join(ch25h_expr, pathway_scores, by = "cell")

meta <- seurat_object@meta.data
meta$cell <- rownames(meta)
df <- inner_join(df, meta[, c("cell", "cell_type")], by = "cell")

df_am5 <- df %>% filter(cell_type == "AM-5.IGHA1")
dim(df_am5)
df_am5<-subset(df_am5,df_am5$CH25H>0)

pathways <- colnames(pathway_scores)
pathways <- pathways[!pathways %in% "cell"] 

cor_res <- sapply(pathways, function(pwy) {
  cor(df_am5$CH25H, df_am5[[pwy]], method = "pearson")
})

pval_res <- sapply(pathways, function(pwy) {
  cor.test(df_am5$CH25H, df_am5[[pwy]], method = "pearson")$p.value
})

cor_df <- data.frame(
  Pathway = pathways,
  Spearman_rho = round(cor_res, 3),
  P_value = signif(pval_res, 3)
) %>%
  arrange(desc(Spearman_rho))

cor_df[cor_df$P_value<0.05,]


####
df_am5 <- df %>% filter(cell_type == "AM-5.IGHA1")
df_am5 <- df_am5 %>%
  mutate(CH25H_group = ifelse(CH25H > 0, "CH25H+", "CH25H-"))

table(df_am5$CH25H_group)
pw <- c("NFkB", "JAK-STAT", "PI3K", "TNFa", "Trail", "VEGF", "WNT")

test_res <- lapply(pw, function(p) {
  w <- wilcox.test(df_am5[[p]] ~ df_am5$CH25H_group)
  data.frame(
    Pathway = p,
    CH25H_pos_mean = mean(df_am5[df_am5$CH25H_group == "CH25H+", p]),
    CH25H_neg_mean = mean(df_am5[df_am5$CH25H_group == "CH25H-", p]),
    delta_mean = mean(df_am5[df_am5$CH25H_group == "CH25H+", p]) -
      mean(df_am5[df_am5$CH25H_group == "CH25H-", p]),
    p_value = w$p.value
  )
}) %>% bind_rows() %>%
  mutate(FDR = p.adjust(p_value, method="BH"))
test_res[test_res$p_value<0.05,]
ggplot(df_am5, aes(CH25H_group, TNFa, fill = CH25H_group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3) +
  scale_fill_manual(values = c("CH25H+" = "#E64B35FF", "CH25H-" = "#4DBBD5FF")) +
  stat_compare_means(
    comparisons = list(c("CH25H+", "CH25H-")),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE,
    size = 6,    
    vjust = 0.3
  ) +
  labs(title = "TNFa", x = "", y = "AM-5 PROGENy Score") +
  guides(fill = "none") +
  theme_classic(base_size = 18) +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 16,face = "bold", color = "black"          ),
    axis.title.x = element_text(size = 16,color = "black"),
    axis.title.y = element_text(size = 16,color = "black"),
    axis.text.x = element_text(size = 16,color = "black"),
    axis.text.y = element_text(size = 16,color = "black")
  )

#########################

ch25h<-read.csv("CH25H.csv")
ch25h<-na.omit(ch25h)
colnames(ch25h)<-c("OS","Status","CH25H")
ch25h$Group <- ifelse(
  ch25h$CH25H >= median(ch25h$CH25H, na.rm = TRUE),
  "High", "Low"
)
surv_object <- Surv(time = ch25h$OS, event = ch25h$Status)
p <- ggsurvplot(
  survfit(surv_object ~ Group, data = ch25h),
  data = ch25h,
  palette = c("#D55E00", "#0072B2"),
  conf.int = TRUE,
  pval = TRUE,
  legend.title = "CH25H",
  legend.labs = c("High", "Low"),
  xlab = "Overall Survival (days)",
  ylab = "Survival Probability",
  font.x = 20, font.y = 20, font.legend = 18, pval.size = 7,
  ggtheme = theme_bw(base_size = 25)
)
p
ggexport(
  p, 
  filename = "CH25H_survivalcurve.pdf",
  width = 6, 
  height = 4.5
)

######


