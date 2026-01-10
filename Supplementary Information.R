### Supplementary Information
library(RColorBrewer)
library(scPagwas)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(tidyverse)
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(KernSmooth)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(scRNAseq)
library(data.table)
library(pheatmap)
library(future)
library(CellChat)
library(ChIPseeker)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)
library(AnnotationHub)
library(annotate)
library(bumphunter)
library(minfi)
library(nor1mix)
library(illuminaio)
library(GEOquery)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(forcats)
library(regioneR)
library(annotatr)
library(fs)
library(GIFT)
library(Jupiter)
library(ggrepel)
library(scRank)

load("Tumor_scPagwas.RData")

pdf("scPagwas_barplot_tumor.pdf",w=14,h=8)
plot_bar_positie_nagtive(seurat_obj=Pagwas,
                         var_ident="positiveCells",
                         var_group="cell_type", 
                         p_thre = 0.01,
                         vec_group_colors=NULL,
                         f_color=colorRampPalette(brewer.pal(n=10, name="RdYlBu")),
                         do_plot = T,
                         fontsize_title = 15, fontsize_axistitle_x = 22, fontsize_axistitle_y = 22, 
                         fontsize_axistext_x = 22, fontsize_axistext_y = 22, fontsize_legendtitle = 18, 
                         fontsize_legendtext = 18 )
dev.off()


Bootstrap_estimate_Plot(bootstrap_results=Pagwas@misc$bootstrap_results,
                        figurenames = "estimateplot.pdf",
                        width = 9,
                        height = 7,
                        do_plot=T)

load("Normal_scPagwas.RData")

pdf("scPagwas_barplot_normal.pdf",w=14,h=8)
plot_bar_positie_nagtive(seurat_obj=Pagwas,
                         var_ident="positiveCells",
                         var_group="cell_type", 
                         p_thre = 0.01,
                         vec_group_colors=NULL,
                         f_color=colorRampPalette(brewer.pal(n=10, name="RdYlBu")),
                         do_plot = T,
                         fontsize_title = 15, fontsize_axistitle_x = 22, fontsize_axistitle_y = 22, 
                         fontsize_axistext_x = 22, fontsize_axistext_y = 22, fontsize_legendtitle = 18, 
                         fontsize_legendtext = 18 )
dev.off()

Bootstrap_estimate_Plot(bootstrap_results=Pagwas@misc$bootstrap_results,
                        figurenames = "estimateplot.pdf",
                        width = 9,
                        height = 7,
                        do_plot=T)


############
loom <- open_loom('out_SCENIC_50000.loom')

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

rownames(regulonAUC)
names(regulons)
seurat.data<-readRDS("Rename_am.rds")

regulonAUC_colnames <- colnames(regulonAUC)
sub_seurat_data <- seurat.data[, match(regulonAUC_colnames, colnames(seurat.data))]
seurat.data<-sub_seurat_data
sub_regulonAUC <- regulonAUC[,match(colnames(seurat.data),colnames(regulonAUC))]
dim(sub_regulonAUC)

identical(colnames(sub_regulonAUC), colnames(seurat.data))

cellClusters <- data.frame(row.names = colnames(seurat.data), 
                           seurat_clusters = as.character(seurat.data@meta.data$seurat_clusters))
cellTypes <- data.frame(row.names = colnames(seurat.data), 
                        celltype = seurat.data@meta.data$cell_type)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 

save(sub_regulonAUC,cellTypes,cellClusters,seurat.data,
     file = 'for_rss_and_visual.Rdata')

###########

seurat.data<-readRDS("Rename_am.rds")
load('for_rss_and_visual.Rdata')
regulonsToPlot = c('ARID5B(+)','ZNF77(+)','NR2F2(+)','ZNF283(+)','ZNF404(+)')
regulonsToPlot %in% row.names(sub_regulonAUC)
seurat.data@meta.data = cbind(seurat.data@meta.data ,t(assay(sub_regulonAUC[regulonsToPlot,])))


selectedResolution <- "celltype"
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])

sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)

regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)

rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellTypes[colnames(sub_regulonAUC), selectedResolution]) 
rss=na.omit(rss) 

rss=regulonActivity_byGroup_Scaled
head(rss)
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss),
                 cluster =   colnames(rss)[i],
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top10 <- df %>% group_by(cluster) %>% top_n(5, fc)
top10 <- top10 %>% distinct(path, .keep_all = TRUE)
rowcn = data.frame(path = top10$cluster)
n = rss[top10$path,]
pheatmap(n,
         annotation_row = rowcn,
         show_rownames = T)


p<-pheatmap(n,cluster_rows = F,cluster_cols = F,fontsize = 12,angle_col = 315)
ggsave("pheatmap_TF_Tumor.pdf",p,w=6,h=10)

#########
coef <- fread("cytosig_scores.Coef", header = F, sep = "\t")
pvalue <- fread("cytosig_scores.Pvalue", header = F, sep = "\t")
stderr <- fread("cytosig_scores.StdErr", header = F, sep = "\t")
zscore <- fread("cytosig_scores.Zscore", header = F, sep = "\t")

coef <- as.data.frame(coef)
rownames(coef) <- coef[[1]]  
coef <- coef[, -1, drop = FALSE] 

pvalue <- as.data.frame(pvalue)
rownames(pvalue) <- pvalue[[1]]
pvalue <- pvalue[, -1, drop = FALSE]

stderr <- as.data.frame(stderr)
rownames(stderr) <- stderr[[1]]
stderr <- stderr[, -1, drop = FALSE]

zscore <- as.data.frame(zscore)
rownames(zscore) <- zscore[[1]]
zscore <- zscore[, -1, drop = FALSE]

seurat_obj<-readRDS("F:/胸外/twaslc/am/Rename_am.rds")
colnames(coef) <- colnames(seurat_obj)

coef_with_type <- as.data.frame(t(coef))  
coef_with_type$CellType <- seurat_obj@meta.data$cell_type

coef_mean <- coef_with_type %>%
  group_by(CellType) %>%
  summarise_all(mean)

coef_mean <- as.data.frame(t(coef_mean[,-1]))  
colnames(coef_mean) <- coef_with_type$CellType[!duplicated(coef_with_type$CellType)]  

pdf("cytosig_celltype.pdf",width=8,height=10)
pheatmap(coef_mean,
         cluster_rows = TRUE,       
         cluster_cols = FALSE,      
         scale = "row",             
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         fontsize_row = 15,         
         fontsize_col = 16,        
         fontsize_legend = 20,      
         main = "Tumor - CytoSig (Cell Type Averages)",fontsize= 20)

dev.off()

############
seurat_object<-readRDS("Rename_am.rds")
data.input <- seurat_object@assays$RNA@counts 
meta <- seurat_object@meta.data
unique(meta$cell_type)

cell.use <- rownames(meta) 
data.input <- data.input[, cell.use]
meta = meta[cell.use, ]
unique(meta$cell_type)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "cell_type") 
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents)) 
saveRDS(groupSize,"groupSize.rds")
groupSize
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat,features = NULL)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
saveRDS(cellchat,"cellchat_am_tumor_computeCommunProb.rds")
##################################
groupSize<-readRDS("groupSize.rds")
cellchat <- readRDS("cellchat_am_tumor_computeCommunProb.rds")
cellchat <- filterCommunication(cellchat , min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
df.net <- subsetCommunication(cellchat)


pdf("netVisual_circle_tumor.pdf",height = 8)
netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = TRUE, 
                 label.edge = FALSE, 
                 title.name = "Tumor - Number of interactions",
                 vertex.label.cex = 1.2,  
                 edge.width.max = 15,vertex.size = 15)
dev.off()

###############
lcgwas<-fread("hg37_lcgwas_GCST90043864.tsv")
SNP<-get_chr_pos(SNP = lcgwas$variant_id,ref_genome = "GRCh38",dbSNP = 155)
SNP<-SNP%>%left_join(lcgwas,by=c("SNP"="variant_id"))
lcgwas<-SNP%>%dplyr::select(CHR:SNP,effect_allele:CONVERGE)
lcgwas$Zscore<-lcgwas$beta/lcgwas$standard_error
fwrite(lcgwas,"lcgwas_hg38.txt",sep = "\t")

###############################

result <- FindAllMarkers(a, only.pos = TRUE, min.pct = 0.25)
group_of_cluster3<-FindMarkers(object = a,ident.1 = "neutrophil",test.use = "wilcox", only.pos = FALSE,min.pct = 0.1)

aa<-fread("differential_genes_normal.csv")
group_of_cluster4<-subset(aa,aa$group%in%"malignant cell")
group_of_cluster4$logfoldchange<-as.numeric(group_of_cluster4$logfoldchange)
group_of_cluster4<-subset(group_of_cluster4,group_of_cluster4$logfoldchange>1.25)
genetumor<-intersect(group_of_cluster4$gene,twasgene$ENSEMBL)
bitr(genetumor,"ENSEMBL","SYMBOL",OrgDb = "org.Hs.eg.db")
#####################################################   manhattan #########################

Gene<-bitr(result_lung$Gene,"ENSEMBL","SYMBOL",OrgDb = "org.Hs.eg.db")
manhattan_data<-merge(result_lung,Gene,by.x="Gene",by.y="ENSEMBL")
manhattan_data$BP<-as.numeric(manhattan_data$GeneEnd)-as.numeric(manhattan_data$GeneStart)+1
manhattan_data<-manhattan_data%>%dplyr::select(SNP=SYMBOL,CHR=CHROM,BP,P=OTTERS_p,z_minp)%>%distinct(SNP,.keep_all = T)
data<-manhattan_data
max_bp_per_chr <- data %>% 
  group_by(CHR) %>% 
  summarise(max_bp = max(BP))
cumulative_max_bp <- cumsum(max_bp_per_chr$max_bp)

data$all <- 0
for (chr in 1:22) {
  chr_indices <- which(data$CHR == chr)
  if (chr > 1) {
    prev_sum <- cumulative_max_bp[chr - 1]
  } else {
    prev_sum <- 0
  }
  
  data$all[c(chr_indices)] <- data$BP[c(chr_indices)] + prev_sum
}
data$BPcum <- data$all
data<-data%>%mutate(is_highlight=ifelse(SNP%in%twasgene$SYMBOL,"yes","no"))%>%
  mutate(is_annotate=ifelse(SNP%in%twasgene$SYMBOL,"yes","no"))%>%mutate(is_color=ifelse(z_minp<0,"yes","no"))
data$CHR <- factor(data$CHR, levels = as.character(1:22))
my_colors <- c(brewer.pal(n = 9, name = "Set3"),brewer.pal(n=8,"Set2"),brewer.pal(n=5,"Set2"),"#C1007A","#436FB8")

X_axis <-  data %>% dplyr::group_by(CHR) %>% dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
X_axis<-X_axis%>%arrange(CHR)
X_axis$center[7]<-10134437.5
X_axis$center[8]<-11152787.5
X_axis$center[9]<-12046052.5
X_axis$center[10]<-	13471218.0
X_axis$center[11]<-14340894.5
X_axis$center[14]<-18549284.0
X_axis$CHR <- factor(X_axis$CHR, levels = as.character(1:22))
pdf("manhanttan_TWAS.pdf",width = 10)
p1 <- ggplot(data, aes(x = BPcum, y = -log10(P))) +
  geom_point(aes(color = data$CHR), alpha = 0.8, size = 1.3) +
  scale_color_manual(values = my_colors) +
  scale_x_continuous(label = X_axis$CHR, breaks = X_axis$center) +
  scale_y_continuous(expand = c(0, 0)) + ylim(0, 25) +
  geom_point(data = subset(data, is_highlight == "yes"), color = "#1D91C0", size = 2) +
  geom_label_repel(data = subset(data, is_annotate == "yes"), aes(label = SNP, color = is_color), size = 5) +
  theme_bw() + 
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 16, face = "bold"),  
    axis.text.y = element_text(size = 16, face = "bold"), 
    axis.title.x = element_text(size = 18, face = "bold"),  
    axis.title.y = element_text(size = 16, face = "bold"),  
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold")  
  ) +
  geom_hline(yintercept = c(-log10(p_lung)), color = "red", size = 1.2, linetype = "dotted") +
  xlab("Chromosome")  
p1
dev.off()

###################################### TWAS gene chromosomal functional annotation
twasgene<-readRDS("twasgene.rds")
annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data1<-as.data.frame(fread("SDPR.txt"))
head(data1)
data1$CHROM<-paste0("chr",data1$CHROM)
data1<-subset(data1,data1$TargetID%in%twasgene$ENSEMBL)
head(annotation)

gr <- GRanges(
  seqnames = data1$CHROM,
  ranges = IRanges(start = data1$POS, end = data1$POS),
  gene_id = data1$TargetID  
)

mnemonics <- import("E088_25_imputed12marks_hg38lift_mnemonics.bed.gz")

overlaps <- findOverlaps(gr, mnemonics)

query_idx <- queryHits(overlaps)   
subject_idx <- subjectHits(overlaps) 

cpg_overlap <- as.data.frame(gr[query_idx])

mnemonics_overlap <- as.data.frame(mnemonics[subject_idx])[, c("seqnames", "start", "end", "name")]
colnames(mnemonics_overlap) <- c("annot_seqnames", "annot_start", "annot_end", "annot_state")

result <- cbind(cpg_overlap, mnemonics_overlap)
head(result)

result_lung<-readRDS("F:/胸外/twaslc/twaslca/result_lung_twas.rds")
result_lung$chr<-paste0("chr",result_lung$CHROM)
result_lung<-result_lung%>%dplyr::select(Gene,OTTERS_p,z_from=chr,CHROM)%>%filter(Gene%in%unique(result$gene_id))
result<-left_join(result_lung,result,by=c("Gene"="gene_id"))
colnames(result)
sorted_CHR_num <- unique(result$z_from[order(as.numeric(result$CHROM))])
result$z_from <- factor(result$z_from, levels = sorted_CHR_num)
result$annot_number <- as.numeric(sub("^(\\d+)_.*", "\\1", result$annot_state))
sorted_states <- unique(result$annot_state[order(result$annot_number)])
result$annot_state <- factor(result$annot_state, levels = sorted_states)
plot_data <- result %>%
  dplyr::group_by(annot_state, z_from) %>%
  dplyr::summarise(
    count = n(),  
    dplyr::across(c(width, strand, annot_number), dplyr::first) 
  ) %>%
  dplyr::ungroup()

plot_data$count<-log10(plot_data$count)
outcome_colors <- c(brewer.pal(10,"Set3"),brewer.pal(8,"Set2"))
p <- ggplot(plot_data, aes(x = annot_state, y = count, fill = z_from)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(name = "") +
  scale_x_discrete(name = "") +
  scale_fill_manual(
    values = outcome_colors,
    breaks = chr_order,       
    labels = labels_map,      
    name = "chromosome"             
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 16),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = c(0.8, 0.8),
    legend.justification = c(1, 0.5),
    legend.background = element_blank(),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  ) +
  guides(
    fill = guide_legend(
      nrow = 4, byrow = TRUE,
      override.aes = list(alpha = 1, size = 4),
      keyheight = unit(4, "mm"),
      keywidth  = unit(6, "mm")
    )
  )
p

###################
hits_mem_treg <- top_hits %>%
  filter(module %in% c("Memory","Treg")) %>%
  arrange(factor(module, levels = c("Memory","Treg")),
          desc(abs(rho)))

hits_mem_treg

plots_mem_treg <- purrr::map2(
  hits_mem_treg$cytokine,
  hits_mem_treg$module,
  ~ plot_cyto_module_median(
    df_merge = df_merge,
    gene_use = .x,
    t_module = .y
  )
)

n_panels <- length(plots_mem_treg)
ncol <- 3
nrow <- ceiling(n_panels / ncol)

p_mem_treg_grid <- ggarrange(
  plotlist = plots_mem_treg,
  nrow = nrow,
  ncol = ncol,
  align = "hv"
)
p_mem_treg_grid

######################
coef <- fread("cytosig_scores.Coef", header = F, sep = "\t")
pvalue <- fread("cytosig_scores.Pvalue", header = F, sep = "\t")
stderr <- fread("cytosig_scores.StdErr", header = F, sep = "\t")
zscore <- fread("cytosig_scores.Zscore", header = F, sep = "\t")

coef <- as.data.frame(coef)
rownames(coef) <- coef[[1]]  
coef <- coef[, -1, drop = FALSE] 

pvalue <- as.data.frame(pvalue)
rownames(pvalue) <- pvalue[[1]]
pvalue <- pvalue[, -1, drop = FALSE]

stderr <- as.data.frame(stderr)
rownames(stderr) <- stderr[[1]]
stderr <- stderr[, -1, drop = FALSE]

zscore <- as.data.frame(zscore)
rownames(zscore) <- zscore[[1]]
zscore <- zscore[, -1, drop = FALSE]

seurat_obj<-readRDS("Rename_am.rds")
colnames(coef) <- colnames(seurat_obj)

coef_with_type <- as.data.frame(t(coef))  
coef_with_type$CellType <- seurat_obj@meta.data$cell_type

coef_mean <- coef_with_type %>%
  group_by(CellType) %>%
  summarise_all(mean)

coef_mean <- as.data.frame(t(coef_mean[,-1])) 
colnames(coef_mean) <- coef_with_type$CellType[!duplicated(coef_with_type$CellType)]  

seurat_obj<-readRDS("Rename_am.rds")
seurat_obj@meta.data$CellTypeGroup <- sub("^(AM-[1-5]).*", "\\1", seurat_obj@meta.data$cell_type)

coef_with_type <- as.data.frame(t(coef)) 
coef_with_type$CellTypeGroup <- seurat_obj@meta.data$CellTypeGroup

coef_mean <- coef_with_type %>%
  group_by(CellTypeGroup) %>%
  summarise_all(mean)

coef_mean <- as.data.frame(t(coef_mean[,-1]))  
colnames(coef_mean) <- unique(coef_with_type$CellTypeGroup)  

pdf("cytosig_AM1-5.pdf",height=10,width=8)
pheatmap(coef_mean,
         cluster_rows = TRUE,      
         cluster_cols = FALSE,     
         scale = "row",            
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         fontsize_row = 18,        
         fontsize_col = 18,        
         fontsize_legend = 20,     
         main = "Tumor - CytoSig (Cell Type Averages)",fontsize= 20)
dev.off()

###################
TWAS_path <- '/GPUFS/gyfyy_jxhe_1/User/zyt/TWAS/Lung/stage1'
GWAS_file <- "/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/GIFT/lcgwas_hg38_outcome.gz"
GWAS_data <- fread(GWAS_file)
GWAS_data <- as.data.frame(GWAS_data)
GWAS_data$samplesize.outcome <- 500000
GWAS_LD_file <- "/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/GIFT/hg38_1000genomes/hg38"
GWAS_N <- 500000

result_heng1<-as.data.frame(readRDS("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/GIFT/result_lung001.rds"))
twasgene<-readRDS("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/am/twasgene.rds")
result_heng1<-subset(result_heng1,result_heng1$Gene%in%twasgene$ENSEMBL)
genes<- bitr(result_heng1$Gene,"ENSEMBL","SYMBOL","org.Hs.eg.db")
result_heng1<-left_join(genes,result_heng1,c("ENSEMBL"="Gene"))
result_heng1 <- result_heng1[order(as.numeric(as.character(result_heng1$CHROM)), result_heng1$GeneStart), ]

split_genes_by_chromosome <- function(df, chrom_col, gene_col, max_group_size = 7) {
  all_splits <- list()
  unique_chromosomes <- unique(df[[chrom_col]])
  for (chr in unique_chromosomes) {
    genes_in_chr <- df[df[[chrom_col]] == chr, gene_col]
    n <- length(genes_in_chr)
    if (n <= max_group_size) {
      split_list <- list(genes_in_chr)
    } else {
      num_groups <- ceiling(n / max_group_size)
      group_sizes <- rep(n %/% num_groups, num_groups)
      remainder <- n %% num_groups
      
      for (i in 1:remainder) {
        group_sizes[i] <- group_sizes[i] + 1
      }
      
      split_list <- list()
      start <- 1
      for (size in group_sizes) {
        end <- min(start + size - 1, n)
        split_list[[length(split_list) + 1]] <- genes_in_chr[start:end]
        start <- end + 1
      }
    }
    all_splits[[chr]] <- split_list
  }
  
  return(all_splits)
}
GIFT_list<-split_genes_by_chromosome(df =result_heng1,chrom_col = "CHROM",gene_col = "ENSEMBL" )

res_list <- list()
for (chr in names(GIFT_list)) {
  grouped_genes <- GIFT_list[[chr]]
  for (genes in grouped_genes) {
    valid_genes <- genes[genes %in% result_heng1$ENSEMBL]
    
    if (length(valid_genes) > 0) {
      result <- tryCatch({
        JUP_GIFT(
          TWAS_path = TWAS_path,
          GWAS_data = GWAS_data,
          GWAS_LD_file = GWAS_LD_file,
          gene_list = valid_genes,
          GWAS_N = GWAS_N,
          parallel = TRUE
        )
      }, error = function(e) {
        message(sprintf("Error in processing genes %s: %s", paste(valid_genes, collapse = ", "), e$message))
        return(NULL)
      })
  
      if (!is.null(result)) {
        result$CHROM <- chr  
        res_list[[length(res_list) + 1]] <- result
      }
    }
  }
}

res <- do.call(rbind, res_list)
saveRDS(res,"res_20241226.rds")


data <- readRDS("res_20241226.rds")
twasgene <- readRDS("../am/twasgene.rds")
protein <- read.table("../am/Protein_gene.txt")
data <- left_join(twasgene, data, c("ENSEMBL" = "gene"))
data$Protein <- ifelse(data$SYMBOL %in% protein$V1, "YES", "NO")
data <- data[!is.na(data$p_ACAT), ]
data <- data[!is.na(data$CHROM), ]


resultlung <- readRDS("/GPUFS/gyfyy_jxhe_1/User/zyt/HTAN/GIFT/result_lung001.rds")
resultlung <- subset(resultlung, resultlung$Gene %in% twasgene$ENSEMBL)
resultlung <- resultlung %>% dplyr::select(Gene, GeneStart, GeneEnd, z_minp)

data <- left_join(data, resultlung, c("ENSEMBL" = "Gene"))


manhattan_data <- data %>%
  dplyr::select(SNP = SYMBOL, CHR = CHROM, BP = GeneStart, P = p_ACAT, z_minp, Protein) %>%
  distinct(SNP, .keep_all = TRUE)

manhattan_data$BP <- as.numeric(manhattan_data$BP)

manhattan_data$CHR <- factor(manhattan_data$CHR, levels = as.character(sort(unique(as.numeric(manhattan_data$CHR)))))

max_bp_per_chr <- manhattan_data %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP, na.rm = TRUE))

cumulative_max_bp <- cumsum(max_bp_per_chr$max_bp)

max_bp_per_chr <- max_bp_per_chr %>%
  mutate(cumulative_start = lag(cumulative_max_bp, default = 0)) 

manhattan_data <- manhattan_data %>%
  left_join(max_bp_per_chr, by = "CHR") %>%
  mutate(BPcum = BP + cumulative_start)  

manhattan_data <- manhattan_data %>%
  mutate(
    is_highlight = ifelse(Protein == "YES", "yes", "no"),
    is_annotate = ifelse(Protein == "YES", "yes", "no"),
    is_color = ifelse(z_minp < 0, "yes", "no") 
  )

X_axis <- manhattan_data %>%
  group_by(CHR) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

my_colors <- c(brewer.pal(n = 9, name = "Set3"), brewer.pal(n = 8, name = "Set2"), 
               brewer.pal(n = 5, name = "Set2"), "#C1007A", "#436FB8")

manhattan_data<-subset(manhattan_data,manhattan_data$P>0)
pdf("GIFT_manhattan_TWAS.pdf",width=10,height=6)
p1 <- ggplot(manhattan_data, aes(x = BPcum, y = -log10(P))) +
  geom_point(aes(color = CHR), alpha = 0.8, size =6) + 
  scale_color_manual(values = my_colors) + 
  scale_x_continuous(
    labels = X_axis$CHR,  
    breaks = X_axis$center  
  ) +
  scale_y_continuous(expand = c(0, 0)) + ylim(0, 150) + 
  geom_point(data = subset(manhattan_data, is_highlight == "yes"), color = "#1D91C0", size = 2) +  
  geom_label_repel(data = subset(manhattan_data, is_annotate == "yes"), 
                   aes(label = SNP, color = is_color), 
                   size = 8,  
                   box.padding = 0.5, 
                   max.overlaps = Inf,  
                   segment.size = 1,  
                   segment.color = "gray") +  
  theme_bw() + 
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1, face="bold", size=12),
    axis.text.y = element_text(angle = 0, hjust = 1, face="bold", size=12),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15)
  ) +
  xlab("Chromosome") + 
  ylab("-log10(P_ACAT)") 
print(p1)
dev.off()
###############
seuratObj <- readRDS("Rename_am.rds")
all_cells <- Cells(seuratObj)
set.seed(123)  
selected_cells <- sample(all_cells, 2000)
subset_seuratObj <- subset(seuratObj, cells = selected_cells)
seuratObj <- subset_seuratObj 
gene_list <- c("CD279","CD274","CH25H","IGFBP1","PSORS1C1","CTLA4","LAG3","TIGIT","PSMB2", "LMNA", "CPS1", "CD47", "MPIG6B", "CSNK2B", "APOM", 
               "PRRC2A", "LST1", "PPP1R11", "CLIC1","UMAD1", "MSRA", 
               "KRT5", "NSF", "NMT1", "PEPD", "IGHA1")
for (gene in gene_list) {
  obj <- CreateScRank(input = seuratObj,
                      species = 'human',
                      cell_type = 'cell_type',
                      target = gene)
  obj <- scRank::Constr_net(obj)
  obj <- scRank::rank_celltype(obj)
  saveRDS(obj, paste0("obj_", gene, ".rds"))
}


obj_PD1 <- readRDS("obj_PD1.rds")
obj_PDL1 <- readRDS("obj_PDL1.rds")
obj_CH25H <- readRDS("obj_CH25H.rds")
obj_IGFBP1 <- readRDS("obj_IGFBP1.rds")
obj_PSORS1C1 <- readRDS("obj_PSORS1C1.rds")
obj_CTLA4<-readRDS("obj_CTLA4.rds")
obj_LAG3<-readRDS("obj_LAG3.rds")
obj_TIGIT<-readRDS("obj_TIGIT.rds")
obj_APOM<-readRDS("obj_APOM.rds")
obj_CD47<-readRDS("obj_CD47.rds")
obj_CLIC1<-readRDS("obj_CLIC1.rds")
obj_CPS1<-readRDS("obj_CPS1.rds")
obj_CSNK2B<-readRDS("obj_CSNK2B.rds")
obj_IGHA1<-readRDS("obj_IGHA1.rds")
obj_KRT5<-readRDS("obj_KRT5.rds")
obj_LMNA<-readRDS("obj_LMNA.rds")
obj_LST1<-readRDS("obj_LST1.rds")
obj_MPIG6B<-readRDS("obj_MPIG6B.rds")
obj_MSRA<-readRDS("obj_MSRA.rds")
obj_NMT1<-readRDS("obj_NMT1.rds")
obj_NSF<-readRDS("obj_NSF.rds")
obj_PEPD<-readRDS("obj_PEPD.rds")
obj_PPP1R11<-readRDS("obj_PPP1R11.rds")
obj_PRRC2A<-readRDS("obj_PRRC2A.rds")
obj_UMAD1<-readRDS("obj_UMAD1.rds")
obj_PSMB2<-readRDS("obj_PSMB2.rds")

obj_names <- c("obj_PD1", "obj_PDL1", "obj_CH25H", "obj_IGFBP1", "obj_PSORS1C1",
               "obj_CTLA4", "obj_LAG3", "obj_TIGIT", "obj_APOM", "obj_CD47", 
               "obj_CLIC1", "obj_CPS1", "obj_CSNK2B", "obj_IGHA1", "obj_KRT5", 
               "obj_LMNA", "obj_LST1", "obj_MPIG6B", "obj_MSRA", "obj_NMT1", 
               "obj_NSF", "obj_PEPD", "obj_PPP1R11", "obj_PRRC2A", "obj_UMAD1", 
               "obj_PSMB2")
for (i in obj_names){
  obj <- get(i)
  pdf(paste0("DimPlot_",i,".pdf"))
  p<-plot_dim(obj,point_size=2,box_padding=1.2,legend_size = 7, 
              legend_text_size=12,legend_title_size=12,label_text = 5.5, x_title_size = 18, y_title_size = 18,
              title=sub("^obj_", "", i),title_size=18)
  print(p)
  dev.off()
}



get_rank_info <- function(obj, name) {
  rank_info <- as.data.frame(obj@cell_type_rank)
  rank_info$cell_type <- rownames(rank_info)
  rank_info$obj_name <- name
  return(rank_info)
}

rank_PD1 <- get_rank_info(obj_PD1, "PD-1")
rank_PDL1 <- get_rank_info(obj_PDL1, "PD-L1")
rank_CH25H <- get_rank_info(obj_CH25H, "CH25H")
rank_IGFBP1 <- get_rank_info(obj_IGFBP1, "IGFBP1")
rank_PSORS1C1 <- get_rank_info(obj_PSORS1C1, "PSORS1C1")
rank_CTLA4 <- get_rank_info(obj_CTLA4, "CTLA4")
rank_LAG3 <- get_rank_info(obj_LAG3, "LAG3")
rank_TIGIT <- get_rank_info(obj_TIGIT, "TIGIT")
rank_APOM <- get_rank_info(obj_APOM, "APOM")
rank_CD47 <- get_rank_info(obj_CD47, "CD47")
rank_CLIC1 <- get_rank_info(obj_CLIC1, "CLIC1")
rank_CPS1 <- get_rank_info(obj_CPS1, "CPS1")
rank_CSNK2B <- get_rank_info(obj_CSNK2B, "CSNK2B")
rank_IGHA1 <- get_rank_info(obj_IGHA1, "IGHA1")
rank_KRT5 <- get_rank_info(obj_KRT5, "KRT5")
rank_LMNA <- get_rank_info(obj_LMNA, "LMNA")
rank_LST1 <- get_rank_info(obj_LST1, "LST1")
rank_MPIG6B <- get_rank_info(obj_MPIG6B, "MPIG6B")
rank_MSRA <- get_rank_info(obj_MSRA, "MSRA")
rank_NMT1 <- get_rank_info(obj_NMT1, "NMT1")
rank_NSF <- get_rank_info(obj_NSF, "NSF")
rank_PEPD <- get_rank_info(obj_PEPD, "PEPD")
rank_PPP1R11 <- get_rank_info(obj_PPP1R11, "PPP1R11")
rank_PRRC2A <- get_rank_info(obj_PRRC2A, "PRRC2A")
rank_UMAD1 <- get_rank_info(obj_UMAD1, "UMAD1")
rank_PSMB2<-get_rank_info(obj_PSMB2,"PSMB2")

all_rank_info <- bind_rows(rank_PD1, rank_PDL1, rank_CH25H, rank_IGFBP1, rank_PSORS1C1,rank_CTLA4,rank_LAG3,rank_TIGIT)
all_rank_info$cell_type[all_rank_info$cell_type%in%"AM-1.MT1B"]<-"AM-1.MT1E"

colors<-c(brewer.pal(8,"Set2"),brewer.pal(5,"Set3"))


all_rank_info<-readRDS("F:/胸外/twaslc/am/scRank/all_rank_info.rds")
pdf("scRank_IGHA1_barplot.pdf",width=12,height = 8)
ggplot(all_rank_info, aes(x = rank, y = reorder(cell_type, rank), fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  facet_wrap(~ obj_name, nrow = 1) +
  theme_minimal() +
  theme(
    text = element_text(size = 18), 
    axis.title.x = element_blank(),   
    axis.title.y = element_blank(),   
    axis.text.y = element_text(size = 18), 
    axis.text.x = element_text(size = 18),
    legend.position = "none"         
  )
dev.off()

all_rank_info <- bind_rows(all_rank_info,rank_APOM, rank_CD47, rank_CLIC1,
                           rank_CPS1, rank_CSNK2B, rank_IGHA1, rank_KRT5, rank_LMNA, rank_LST1,
                           rank_MPIG6B, rank_MSRA, rank_NMT1, rank_NSF, rank_PEPD, rank_PPP1R11,
                           rank_PRRC2A, rank_UMAD1,rank_PSMB2
)


all_rank_info_sorted <- all_rank_info[order(all_rank_info$rank), ]

all_rank_info_sorted <- all_rank_info_sorted %>%
  mutate(Target = ifelse(obj_name %in% c("TIGIT", "LAG3", "PD-1", "PD-L1", "CTLA4"),
                         "Immune Checkpoint", "TWAS Target"))
heatmap_data <- all_rank_info_sorted %>%
  select(obj_name, cell_type, rank, Target) %>%
  pivot_wider(names_from = cell_type, values_from = rank)

row_annotation <- heatmap_data %>%
  select(obj_name, Target) %>%
  column_to_rownames("obj_name")  

heatmap_data <- heatmap_data %>%
  arrange(factor(Target, levels = c("Immune Checkpoint", "TWAS Target"))) 

heatmap_matrix <- heatmap_data %>%
  select(-Target) %>%
  column_to_rownames("obj_name") %>%
  as.matrix()

color_palette <- colorRampPalette(c("#05A8C5", "#FFFFFF", "#E34A2B"))(100)

annotation_colors <- list(
  Target = c("Immune Checkpoint" = "#F948F3", "TWAS Target" = "#5DE0D0")
)
pdf("scRank_heatmap_celltype_obj.pdf",width=8,height=8)
pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,                
  cluster_cols = FALSE,                
  annotation_row = row_annotation,     
  annotation_colors = annotation_colors, 
  color = color_palette,                
  show_rownames = TRUE,                 
  show_colnames = TRUE,                
  main = "",
  fontsize = 12,
  fontsize_row = 15,                   
  fontsize_col = 15,                    
  legend = TRUE,                       
  angle_col = 270                      
)
dev.off()