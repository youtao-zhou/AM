# AM2 AM4 heatmap
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

options(stringsAsFactors = FALSE)
set.seed(1)

AM_ORDER <- paste0("AM-", 1:5)
AM2_COLOR <- "#2F6FB0"
AM4_COLOR <- "#C64A3B"
GREY_COLOR <- "#9CA3AA"
DARK_COLOR <- "#22272B"

zscore_columns <- function(x) {
  x <- as.matrix(x)
  centered <- sweep(x, 2, colMeans(x, na.rm = TRUE), "-")
  sds <- sqrt(colMeans(centered^2, na.rm = TRUE))
  sweep(centered, 2, sds, "/")
}

SIGNATURES <- list(
  AntigenPres_MHCII = c(
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
    "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "CD74", "CIITA"
  ),
  AntigenPres_MHCI = c(
    "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5"
  ),
  Costimulation = c("CD80", "CD86", "CD40", "ICOSLG", "TNFSF9", "TNFSF4", "CD58"),
  Checkpoint_ligand = c("CD274", "PDCD1LG2", "LGALS9", "VSIR", "CD276", "IDO1"),
  Tcell_recruit = c("CXCL9", "CXCL10", "CXCL11", "CXCL16", "CCL5", "CCL19", "CXCL13"),
  IFN_I = c("ISG15", "IFI6", "IFI27", "MX1", "OAS1", "IFIT1", "IFIT3", "IRF7"),
  IFN_II = c("STAT1", "GBP1", "GBP5", "IDO1", "CXCL9", "CIITA"),
  Inflammatory = c("IL1B", "TNF", "IL6", "CXCL8", "CCL3", "CCL4", "CXCL2", "NFKB1"),
  Complement = c("C1QA", "C1QB", "C1QC", "C3", "CFB", "C2"),
  Reparative_M2 = c("MRC1", "CD163", "STAB1", "F13A1", "FOLR2", "LYVE1", "SELENOP", "GAS6"),
  Efferocytosis = c("MERTK", "AXL", "TYRO3", "GAS6", "PROS1", "CD36", "MFGE8"),
  ProTumor_GF = c("SPP1", "VEGFA", "PDGFB", "INHBA", "MMP9", "MMP12", "MMP14", "GRN"),
  TGFb_signal = c("TGFB1", "TGFB2", "TGFBI", "LTBP1", "SMAD3"),
  ECM_matrix = c("FN1", "SPARC", "COL1A1", "LAMC1", "THBS1", "TIMP1", "CTSB"),
  Angiogenesis = c("VEGFA", "VCAN", "FLT1", "ANGPT2", "CD93"),
  Lipid_metab = c("APOC1", "APOE", "FABP4", "FABP5", "PLIN2", "LPL", "TREM2", "SPP1"),
  OXPHOS = c("NDUFA4", "COX5A", "COX6C", "ATP5F1A", "UQCRB", "NDUFB2"),
  Glycolysis = c("LDHA", "PKM", "ENO1", "GAPDH", "SLC2A1", "PGK1"),
  Proliferation = c("MKI67", "TOP2A", "PCNA", "CCNB1", "CDK1", "UBE2C"),
  Stress_hypoxia = c("HIF1A", "VEGFA", "BNIP3", "DDIT4", "HSPA1A", "FOS", "JUN"),
  Phagocytosis = c("CD68", "MARCO", "MSR1", "CD36", "FCGR1A", "FCGR3A"),
  Interferon_stim_chemok = c("CXCL10", "CXCL11", "CCL8", "CCL2"),
  AM_alveolar_identity = c("MARCO", "MCEMP1", "FABP4", "PPARG", "SIGLEC1"),
  DC_like = c("FLT3", "CLEC9A", "XCR1", "BATF3", "ZBTB46"),
  Stemness_AM4 = c("SPIDR", "RAD51B", "ATG7", "NFAT5", "TPRG1")
)

sc <- readRDS(am_data("HTAN/am/object/Rename_am.rds"))
DefaultAssay(sc) <- "RNA"
sc$fine <- as.character(sc$cell_type)

sig_present <- lapply(SIGNATURES, intersect, y = rownames(sc))
sig_present <- sig_present[lengths(sig_present) >= 3]

sc <- AddModuleScore(sc, features = sig_present, name = "F1F_", seed = 1)
score_cols <- setNames(paste0("F1F_", seq_along(sig_present)), names(sig_present))

md <- sc@meta.data
fine_levels <- unique(as.character(md$fine))
fine_sig <- sapply(names(sig_present), function(s) {
  tapply(md[[score_cols[[s]]]], factor(md$fine, levels = fine_levels), mean, na.rm = TRUE)
})
fine_sig <- as.matrix(fine_sig)
fine_sig <- fine_sig[rowSums(is.na(fine_sig)) == 0, , drop = FALSE]

fine_macro <- sub("^((AM-[0-9]+)).*$", "\\1", rownames(fine_sig))
macro_sig <- t(vapply(
  AM_ORDER,
  function(am) colMeans(fine_sig[fine_macro == am, , drop = FALSE]),
  numeric(ncol(fine_sig))
))
rownames(macro_sig) <- AM_ORDER
signature_z <- zscore_columns(macro_sig)

EXPRESSION_MODULES <- list(
  "Alveolar homeostasis" = c("AM_alveolar_identity", "Phagocytosis", "Efferocytosis"),
  "Lipid metabolism" = "Lipid_metab",
  "MHC-II antigen presentation" = "AntigenPres_MHCII",
  "IFN / inflammatory response" = c("IFN_I", "IFN_II", "Inflammatory", "Interferon_stim_chemok"),
  "T-cell regulation" = c("Costimulation", "Checkpoint_ligand", "Tcell_recruit")
)
expression_module_raw <- sapply(names(EXPRESSION_MODULES), function(m) {
  rowMeans(signature_z[, EXPRESSION_MODULES[[m]], drop = FALSE])
})
rownames(expression_module_raw) <- AM_ORDER
expression_module_z <- zscore_columns(expression_module_raw)

GO_MODULES <- list(
  "Immune-cell recruitment" = c(
    "chemotaxis", "leukocyte migration", "granulocyte migration",
    "neutrophil migration", "myeloid leukocyte migration"
  ),
  "Cell-matrix adhesion" = c("cell-matrix adhesion", "cell-substrate (?:adhesion|junction)"),
  "Tissue morphogenesis / remodeling" = c(
    "GTPase", "Ras protein signal", "cell morphogenesis",
    "plasma membrane tubulation", "actin filament organization", "cell polarity"
  ),
  "Vesicle trafficking / autophagy" = c(
    "Golgi", "endosomal transport", "vesicle organization", "cytosolic transport",
    "macroautophagy", "protein localization to organelle", "intracellular transport"
  ),
  "Stress adaptation / wound repair" = c(
    "wound healing", "reactive oxygen species", "endoplasmic reticulum stress",
    "double-strand break repair", "histone modification", "ubiquitination",
    "protein folding", "cellular catabolic process", "apoptotic signaling"
  )
)

read_sheet <- function(sheet, tissue) {
  x <- read_excel(am_data("HTAN/am/enrichment/am_enrichment.xlsx"), sheet = sheet)
  x$tissue <- tissue
  x
}
enrichment <- bind_rows(read_sheet("am_tumor", "Tumor"), read_sheet("am_normal", "Normal"))
enrichment$AM <- sub("^((AM-[0-9]+)).*$", "\\1", enrichment$cell_type)
enrichment$p.adjust <- as.numeric(enrichment$p.adjust)
enrichment$evidence_score <- pmin(-log10(pmax(enrichment$p.adjust, 1e-300)), 15)

match_first_rule <- function(description, patterns) {
  hit <- rep("", length(description))
  for (p in patterns) {
    cur <- grepl(p, description, ignore.case = TRUE, perl = TRUE)
    hit[cur & hit == ""] <- p
  }
  hit
}

go_evidence <- bind_rows(lapply(names(GO_MODULES), function(m) {
  keep <- match_first_rule(enrichment$Description, GO_MODULES[[m]]) != ""
  out <- enrichment[keep, , drop = FALSE]
  out$module <- m
  out
}))

go_macro_df <- crossing(
  distinct(enrichment, tissue, cell_type, AM),
  module = names(GO_MODULES)
) %>%
  left_join(
    go_evidence %>%
      group_by(tissue, cell_type, AM, module) %>%
      summarise(evidence_score = max(evidence_score), .groups = "drop"),
    by = c("tissue", "cell_type", "AM", "module")
  ) %>%
  mutate(evidence_score = replace_na(evidence_score, 0)) %>%
  group_by(tissue, AM, module) %>%
  summarise(evidence_score = mean(evidence_score), .groups = "drop") %>%
  group_by(AM, module) %>%
  summarise(evidence_score = mean(evidence_score), .groups = "drop") %>%
  pivot_wider(names_from = module, values_from = evidence_score) %>%
  arrange(factor(AM, levels = AM_ORDER))

go_macro <- as.matrix(go_macro_df[, names(GO_MODULES), drop = FALSE])
rownames(go_macro) <- go_macro_df$AM
go_module_z <- zscore_columns(go_macro)

MODULE_ORDER <- c(
  "MHC-II antigen presentation",
  "T-cell regulation",
  "Immune-cell recruitment",
  "Cell-matrix adhesion",
  "Tissue morphogenesis / remodeling",
  "Vesicle trafficking / autophagy",
  "Stress adaptation / wound repair"
)
module_z <- cbind(expression_module_z, go_module_z)[AM_ORDER, MODULE_ORDER, drop = FALSE]

pca <- prcomp(module_z, center = TRUE, scale. = FALSE, rank. = 2)
coords <- pca$x[, 1:2, drop = FALSE]
variance <- pca$sdev^2 / sum(pca$sdev^2)

if (coords["AM-4", "PC1"] < coords["AM-2", "PC1"]) coords[, "PC1"] <- -coords[, "PC1"]
if (coords["AM-2", "PC2"] < coords["AM-3", "PC2"]) coords[, "PC2"] <- -coords[, "PC2"]

MODULE_LABELS <- c(
  "MHC-II antigen presentation"       = "MHC-II antigen presentation",
  "T-cell regulation"                 = "T-cell regulation",
  "Immune-cell recruitment"           = "Immune-cell recruitment",
  "Cell-matrix adhesion"              = "Cell-matrix adhesion",
  "Tissue morphogenesis / remodeling" = "Tissue remodeling",
  "Vesicle trafficking / autophagy"   = "Intracellular trafficking",
  "Stress adaptation / wound repair"  = "Tissue repair"
)

heatmap_long <- as.data.frame(module_z, check.names = FALSE) %>%
  mutate(AM = rownames(.)) %>%
  pivot_longer(all_of(MODULE_ORDER), names_to = "module", values_to = "z") %>%
  mutate(
    AM = factor(AM, levels = AM_ORDER),
    module = factor(module,
      levels = rev(MODULE_ORDER),
      labels = unname(MODULE_LABELS[rev(MODULE_ORDER)])
    )
  )

HEAT_COLORS <- c(
  "#4575B4", "#74ADD1", "#ABD9E9", "#FFFFBF",
  "#FFFFBF", "#FEE090", "#F46D43", "#D73027"
)

HEAT_VALUES <- c(0.000, 0.300, 0.375, 0.450, 0.550, 0.625, 0.700, 1.000)

p_heat <- ggplot(heatmap_long, aes(AM, module, fill = z)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_gradientn(
    colors = HEAT_COLORS, values = HEAT_VALUES,
    limits = c(-2, 2), oob = squish, breaks = c(-2, 0, 2),
    name = "Zscore"
  ) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), labels = scales::label_wrap(16)) +
  guides(fill = guide_colorbar(
    title.position = "top", title.hjust = 0.5,
    barwidth = grid::unit(4, "mm"), barheight = grid::unit(32, "mm"),
    ticks = FALSE, frame.colour = NA
  )) +
  labs(x = NULL, y = NULL) +
  theme(
    text = element_text(color = "black"),
    axis.text.x = element_text(
      size = 12, color = "black",
      margin = margin(b = 14)
    ),
    axis.text.y = element_text(
      size = 12, color = "black",
      lineheight = 0.9, margin = margin(r = 6)
    ),
    axis.ticks = element_blank(), axis.line = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    plot.margin = margin(t = 10, r = 8, b = 0, l = 2)
  )
p_heat
