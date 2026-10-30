# AM2 AM4 functional analysis
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

ROOT <- am_data("HTAN")
RDS_FILE <- file.path(ROOT, "am/object/Rename_am.rds")
ENRICH_FILE <- file.path(ROOT, "am/enrichment/am_enrichment.xlsx")
OUT_DIR <- Sys.getenv(
  "AM_FUNCTIONAL_OUT_DIR",
  unset = file.path(ROOT, "am/AM2_AM4_PCA")
)
REBUILD_SIGNATURE_MATRIX <- tolower(
  Sys.getenv("AM_FUNCTIONAL_REBUILD_SIGNATURES", unset = "true")
) %in% c("1", "true", "yes", "y")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SIGNATURE_MATRIX_FILE <- file.path(OUT_DIR, "AM_fine_subtype_signature_scores.csv")
AM_ORDER <- paste0("AM-", 1:5)
AM2_COLOR <- "#2F6FB0"
AM4_COLOR <- "#C64A3B"
GREY_COLOR <- "#9CA3AA"
DARK_COLOR <- "#22272B"

stopifnot(file.exists(RDS_FILE), file.exists(ENRICH_FILE))

message("AM-2/AM-4 functional-analysis output: ", OUT_DIR)
message("Rebuild signature matrix: ", REBUILD_SIGNATURE_MATRIX)

zscore_columns <- function(x) {
  x <- as.matrix(x)
  means <- colMeans(x, na.rm = TRUE)
  centered <- sweep(x, 2, means, FUN = "-")
  sds <- sqrt(colMeans(centered^2, na.rm = TRUE))
  if (any(!is.finite(sds) | sds == 0)) {
    bad <- colnames(x)[!is.finite(sds) | sds == 0]
    stop("Cannot z-score invariant column(s): ", paste(bad, collapse = ", "))
  }
  sweep(centered, 2, sds, FUN = "/")
}

write_matrix_csv <- function(x, file, row_name = "row") {
  out <- data.frame(x, check.names = FALSE)
  out <- cbind(setNames(data.frame(rownames(out)), row_name), out)
  write.csv(out, file, row.names = FALSE, quote = TRUE)
}

SIGNATURES <- list(
  AntigenPres_MHCII = c(
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
    "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "CD74", "CIITA"
  ),
  AntigenPres_MHCI = c(
    "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9",
    "NLRC5"
  ),
  Costimulation = c(
    "CD80", "CD86", "CD40", "ICOSLG", "TNFSF9", "TNFSF4", "CD58"
  ),
  Checkpoint_ligand = c(
    "CD274", "PDCD1LG2", "LGALS9", "VSIR", "CD276", "IDO1"
  ),
  Tcell_recruit = c(
    "CXCL9", "CXCL10", "CXCL11", "CXCL16", "CCL5", "CCL19", "CXCL13"
  ),
  IFN_I = c("ISG15", "IFI6", "IFI27", "MX1", "OAS1", "IFIT1", "IFIT3", "IRF7"),
  IFN_II = c("STAT1", "GBP1", "GBP5", "IDO1", "CXCL9", "CIITA"),
  Inflammatory = c("IL1B", "TNF", "IL6", "CXCL8", "CCL3", "CCL4", "CXCL2", "NFKB1"),
  Complement = c("C1QA", "C1QB", "C1QC", "C3", "CFB", "C2"),
  Reparative_M2 = c(
    "MRC1", "CD163", "STAB1", "F13A1", "FOLR2", "LYVE1", "SELENOP", "GAS6"
  ),
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

build_signature_matrix <- function() {
  message("1/5 Reading Seurat object ...")
  sc <- readRDS(RDS_FILE)
  if (!"RNA" %in% Assays(sc)) {
    stop("RNA assay is absent from ", RDS_FILE)
  }
  if (!"cell_type" %in% colnames(sc@meta.data)) {
    stop("cell_type is absent from Seurat metadata")
  }
  DefaultAssay(sc) <- "RNA"
  sc$fine <- as.character(sc$cell_type)

  signatures_present <- lapply(SIGNATURES, intersect, y = rownames(sc))
  retained <- lengths(signatures_present) >= 3
  if (any(!retained)) {
    message(
      "Signatures omitted (<3 detected genes): ",
      paste(names(signatures_present)[!retained], collapse = ", ")
    )
  }
  signatures_present <- signatures_present[retained]

  gene_audit <- bind_rows(lapply(names(SIGNATURES), function(signature) {
    requested <- SIGNATURES[[signature]]
    detected <- intersect(requested, rownames(sc))
    data.frame(
      signature = signature,
      requested_n = length(requested),
      detected_n = length(detected),
      retained = length(detected) >= 3,
      detected_genes = paste(detected, collapse = ";"),
      missing_genes = paste(setdiff(requested, detected), collapse = ";")
    )
  }))
  write.csv(
    gene_audit,
    file.path(OUT_DIR, "AM_functional_signature_gene_audit.csv"),
    row.names = FALSE,
    quote = TRUE
  )

  required_for_final <- c(
    "AntigenPres_MHCII", "Costimulation", "Checkpoint_ligand", "Tcell_recruit"
  )
  missing_required <- setdiff(required_for_final, names(signatures_present))
  if (length(missing_required) > 0) {
    stop(
      "Required final-panel signatures have <3 detected genes: ",
      paste(missing_required, collapse = ", ")
    )
  }

  message("2/5 Calculating AddModuleScore signatures ...")

  sc <- AddModuleScore(
    object = sc,
    features = signatures_present,
    name = "AMFUNC_",
    seed = 1
  )
  score_columns <- paste0("AMFUNC_", seq_along(signatures_present))
  names(score_columns) <- names(signatures_present)
  if (!all(score_columns %in% colnames(sc@meta.data))) {
    stop("AddModuleScore did not create the expected score columns")
  }

  md <- sc@meta.data
  fine_levels <- unique(as.character(md$fine))
  signature_matrix <- sapply(names(signatures_present), function(signature) {
    tapply(
      md[[score_columns[[signature]]]],
      factor(md$fine, levels = fine_levels),
      mean,
      na.rm = TRUE
    )
  })
  signature_matrix <- as.matrix(signature_matrix)
  signature_matrix <- signature_matrix[rowSums(is.na(signature_matrix)) == 0, , drop = FALSE]
  colnames(signature_matrix) <- names(signatures_present)
  rownames(signature_matrix) <- fine_levels[match(rownames(signature_matrix), fine_levels)]

  macro <- sub("^((AM-[0-9]+)).*$", "\\1", rownames(signature_matrix))
  if (any(!macro %in% AM_ORDER)) {
    stop(
      "Could not map fine subtype(s) to AM-1...AM-5: ",
      paste(rownames(signature_matrix)[!macro %in% AM_ORDER], collapse = ", ")
    )
  }

  write_matrix_csv(
    signature_matrix,
    SIGNATURE_MATRIX_FILE,
    row_name = "fine_subtype"
  )
  signature_matrix
}

read_signature_matrix <- function() {
  if (!file.exists(SIGNATURE_MATRIX_FILE)) {
    stop(
      "AM_FUNCTIONAL_REBUILD_SIGNATURES=false, but cached matrix is absent: ",
      SIGNATURE_MATRIX_FILE
    )
  }
  x <- read.csv(SIGNATURE_MATRIX_FILE, check.names = FALSE)
  rownames(x) <- x[[1]]
  as.matrix(x[, -1, drop = FALSE])
}

fine_signature_matrix <- if (REBUILD_SIGNATURE_MATRIX) {
  build_signature_matrix()
} else {
  message("1-2/5 Reusing cached signature matrix ...")
  read_signature_matrix()
}

fine_macro <- sub("^((AM-[0-9]+)).*$", "\\1", rownames(fine_signature_matrix))
if (any(!fine_macro %in% AM_ORDER)) {
  stop("The fine-subtype matrix contains labels outside AM-1...AM-5")
}

macro_signature_matrix <- vapply(
  AM_ORDER,
  function(am) colMeans(fine_signature_matrix[fine_macro == am, , drop = FALSE]),
  numeric(ncol(fine_signature_matrix))
)
macro_signature_matrix <- t(macro_signature_matrix)
rownames(macro_signature_matrix) <- AM_ORDER
colnames(macro_signature_matrix) <- colnames(fine_signature_matrix)
signature_z <- zscore_columns(macro_signature_matrix)

write_matrix_csv(
  macro_signature_matrix,
  file.path(OUT_DIR, "AM_major_class_signature_scores.csv"),
  row_name = "AM"
)
write_matrix_csv(
  signature_z,
  file.path(OUT_DIR, "AM_major_class_signature_scores_z.csv"),
  row_name = "AM"
)

EXPRESSION_MODULES <- list(
  "Alveolar homeostasis" = c(
    "AM_alveolar_identity", "Phagocytosis", "Efferocytosis"
  ),
  "Lipid metabolism" = "Lipid_metab",
  "MHC-II antigen presentation" = "AntigenPres_MHCII",
  "IFN / inflammatory response" = c(
    "IFN_I", "IFN_II", "Inflammatory", "Interferon_stim_chemok"
  ),
  "T-cell regulation" = c(
    "Costimulation", "Checkpoint_ligand", "Tcell_recruit"
  )
)

expression_module_raw <- sapply(names(EXPRESSION_MODULES), function(module) {
  components <- EXPRESSION_MODULES[[module]]
  missing <- setdiff(components, colnames(signature_z))
  if (length(missing) > 0) {
    stop(
      "Missing expression signature(s) for ", module, ": ",
      paste(missing, collapse = ", ")
    )
  }
  rowMeans(signature_z[, components, drop = FALSE])
})
rownames(expression_module_raw) <- AM_ORDER
expression_module_z <- zscore_columns(expression_module_raw)

GO_MODULES <- list(
  "Immune-cell recruitment" = c(
    "chemotaxis",
    "leukocyte migration",
    "granulocyte migration",
    "neutrophil migration",
    "myeloid leukocyte migration"
  ),
  "Cell-matrix adhesion" = c(
    "cell-matrix adhesion",
    "cell-substrate (?:adhesion|junction)"
  ),
  "Tissue morphogenesis / remodeling" = c(
    "GTPase",
    "Ras protein signal",
    "cell morphogenesis",
    "plasma membrane tubulation",
    "actin filament organization",
    "cell polarity"
  ),
  "Vesicle trafficking / autophagy" = c(
    "Golgi",
    "endosomal transport",
    "vesicle organization",
    "cytosolic transport",
    "macroautophagy",
    "protein localization to organelle",
    "intracellular transport"
  ),
  "Stress adaptation / wound repair" = c(
    "wound healing",
    "reactive oxygen species",
    "endoplasmic reticulum stress",
    "double-strand break repair",
    "histone modification",
    "ubiquitination",
    "protein folding",
    "cellular catabolic process",
    "apoptotic signaling"
  )
)

message("3/5 Reading Normal and Tumor GO enrichment ...")
read_enrichment_sheet <- function(sheet, tissue) {
  x <- read_excel(ENRICH_FILE, sheet = sheet)
  required <- c("cell_type", "Description", "ID", "p.adjust")
  missing <- setdiff(required, colnames(x))
  if (length(missing) > 0) {
    stop(
      "Missing column(s) in enrichment sheet ", sheet, ": ",
      paste(missing, collapse = ", ")
    )
  }
  x$tissue <- tissue
  x
}

enrichment <- bind_rows(
  read_enrichment_sheet("am_tumor", "Tumor"),
  read_enrichment_sheet("am_normal", "Normal")
)
enrichment$AM <- sub("^((AM-[0-9]+)).*$", "\\1", enrichment$cell_type)
enrichment$p.adjust <- as.numeric(enrichment$p.adjust)
if (any(!enrichment$AM %in% AM_ORDER) || any(!is.finite(enrichment$p.adjust))) {
  stop("Enrichment data contain invalid AM labels or adjusted P values")
}
enrichment$evidence_score <- pmin(
  -log10(pmax(enrichment$p.adjust, 1e-300)),
  15
)

match_first_rule <- function(description, patterns) {
  hit <- rep("", length(description))
  for (pattern in patterns) {
    current <- grepl(pattern, description, ignore.case = TRUE, perl = TRUE)
    hit[current & hit == ""] <- pattern
  }
  hit
}

go_evidence <- bind_rows(lapply(names(GO_MODULES), function(module) {
  patterns <- GO_MODULES[[module]]
  matched_rule <- match_first_rule(enrichment$Description, patterns)
  keep <- matched_rule != ""
  out <- enrichment[keep, , drop = FALSE]
  out$module <- module
  out$matched_rule <- matched_rule[keep]
  out
}))

fine_groups <- distinct(enrichment, tissue, cell_type, AM)
go_fine <- crossing(
  fine_groups,
  module = names(GO_MODULES)
) %>%
  left_join(
    go_evidence %>%
      group_by(tissue, cell_type, AM, module) %>%
      summarise(evidence_score = max(evidence_score), .groups = "drop"),
    by = c("tissue", "cell_type", "AM", "module")
  ) %>%
  mutate(evidence_score = replace_na(evidence_score, 0))

go_macro_long <- go_fine %>%
  group_by(tissue, AM, module) %>%
  summarise(evidence_score = mean(evidence_score), .groups = "drop") %>%
  group_by(AM, module) %>%
  summarise(evidence_score = mean(evidence_score), .groups = "drop")

go_macro_df <- go_macro_long %>%
  pivot_wider(names_from = module, values_from = evidence_score) %>%
  arrange(factor(AM, levels = AM_ORDER))
go_macro <- as.matrix(go_macro_df[, names(GO_MODULES), drop = FALSE])
rownames(go_macro) <- go_macro_df$AM
go_module_z <- zscore_columns(go_macro)

write_matrix_csv(
  expression_module_z,
  file.path(OUT_DIR, "AM_expression_modules_z_all.csv"),
  row_name = "AM"
)
write_matrix_csv(
  go_macro,
  file.path(OUT_DIR, "AM_GO_module_evidence_scores.csv"),
  row_name = "AM"
)
write_matrix_csv(
  go_module_z,
  file.path(OUT_DIR, "AM_GO_modules_z.csv"),
  row_name = "AM"
)
write.csv(
  go_evidence,
  file.path(OUT_DIR, "AM_GO_matched_term_evidence.csv"),
  row.names = FALSE,
  quote = TRUE
)

MODULE_ORDER <- c(
  "MHC-II antigen presentation",
  "T-cell regulation",
  "Immune-cell recruitment",
  "Cell-matrix adhesion",
  "Tissue morphogenesis / remodeling",
  "Vesicle trafficking / autophagy",
  "Stress adaptation / wound repair"
)
IMMUNE_ENDPOINT <- MODULE_ORDER[1:3]
REPAIR_ENDPOINT <- MODULE_ORDER[4:7]

all_module_z <- cbind(expression_module_z, go_module_z)
module_z <- all_module_z[AM_ORDER, MODULE_ORDER, drop = FALSE]
if (any(!is.finite(module_z))) {
  stop("The final 5 x 7 module matrix contains non-finite values")
}

if (any(apply(module_z[, IMMUNE_ENDPOINT, drop = FALSE], 2, which.max) != 2)) {
  stop("At least one immune module is no longer maximal in AM-2")
}
if (any(apply(module_z[, REPAIR_ENDPOINT, drop = FALSE], 2, which.max) != 4)) {
  stop("At least one repair/remodeling module is no longer maximal in AM-4")
}

write_matrix_csv(
  t(module_z),
  file.path(OUT_DIR, "AM2_AM4_functional_module_scores_z.csv"),
  row_name = "Functional module"
)

pca <- prcomp(module_z, center = TRUE, scale. = FALSE, rank. = 2)
coords <- pca$x[, 1:2, drop = FALSE]
loadings <- pca$rotation[, 1:2, drop = FALSE]
variance <- pca$sdev^2 / sum(pca$sdev^2)

if (coords["AM-4", "PC1"] < coords["AM-2", "PC1"]) {
  coords[, "PC1"] <- -coords[, "PC1"]
  loadings[, "PC1"] <- -loadings[, "PC1"]
}
if (coords["AM-2", "PC2"] < coords["AM-3", "PC2"]) {
  coords[, "PC2"] <- -coords[, "PC2"]
  loadings[, "PC2"] <- -loadings[, "PC2"]
}

endpoint_composites <- data.frame(
  AM = AM_ORDER,
  `Immune regulation / T-cell coordination` =
    rowMeans(module_z[, IMMUNE_ENDPOINT, drop = FALSE]),
  `Repair / remodeling` =
    rowMeans(module_z[, REPAIR_ENDPOINT, drop = FALSE]),
  check.names = FALSE
)

pca_coordinates <- data.frame(
  AM = AM_ORDER,
  PC1 = coords[AM_ORDER, "PC1"],
  PC2 = coords[AM_ORDER, "PC2"],
  variance_explained_PC1 = variance[1],
  variance_explained_PC2 = variance[2],
  check.names = FALSE
)
pca_loadings <- data.frame(
  `Functional module` = rownames(loadings),
  PC1 = loadings[, "PC1"],
  PC2 = loadings[, "PC2"],
  check.names = FALSE
)

write.csv(
  endpoint_composites,
  file.path(OUT_DIR, "AM2_AM4_functional_endpoint_scores.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  pca_coordinates,
  file.path(OUT_DIR, "AM2_AM4_functional_PCA_coordinates.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  pca_loadings,
  file.path(OUT_DIR, "AM2_AM4_functional_PCA_loadings.csv"),
  row.names = FALSE,
  quote = TRUE
)

definitions <- bind_rows(
  lapply(MODULE_ORDER, function(module) {
    is_expression <- module %in% names(EXPRESSION_MODULES)
    components <- if (is_expression) {
      EXPRESSION_MODULES[[module]]
    } else {
      GO_MODULES[[module]]
    }
    data.frame(
      module = module,
      domain = if (module %in% IMMUNE_ENDPOINT) {
        "AM-2 immune function"
      } else {
        "AM-4 repair / remodeling"
      },
      source = if (is_expression) {
        "scRNA gene-signature score"
      } else {
        "Normal + Tumor top-GO enrichment"
      },
      components_or_matching_rules = paste(components, collapse = "; ")
    )
  })
)
write.csv(
  definitions,
  file.path(OUT_DIR, "AM2_AM4_functional_module_definitions.csv"),
  row.names = FALSE,
  quote = TRUE
)

excluded_reasons <- c(
  "Alveolar homeostasis" =
    "Excluded: AM identity/clearance rather than the AM-2 immune endpoint",
  "Lipid metabolism" =
    "Excluded: maximum occurs outside AM-2/AM-4",
  "IFN / inflammatory response" =
    "Excluded: maximum occurs in AM-3 rather than AM-2/AM-4"
)
selection_audit <- bind_rows(
  lapply(MODULE_ORDER, function(module) {
    values <- module_z[, module]
    data.frame(
      module = module,
      decision = "Retained",
      reason = "Contributes directly to the AM-2 immune or AM-4 repair endpoint",
      highest_AM = names(which.max(values)),
      highest_z = max(values)
    )
  }),
  lapply(names(excluded_reasons), function(module) {
    values <- all_module_z[, module]
    data.frame(
      module = module,
      decision = "Excluded from selection-focused panel",
      reason = unname(excluded_reasons[[module]]),
      highest_AM = names(which.max(values)),
      highest_z = max(values)
    )
  })
)
write.csv(
  selection_audit,
  file.path(OUT_DIR, "AM2_AM4_functional_module_selection_audit.csv"),
  row.names = FALSE,
  quote = TRUE
)

message("4/5 Drawing heatmap and PCA ...")
theme_set(theme_classic(base_family = "Arial", base_size = 10))

heatmap_long <- as.data.frame(module_z, check.names = FALSE) %>%
  mutate(AM = rownames(.)) %>%
  pivot_longer(
    cols = all_of(MODULE_ORDER),
    names_to = "module",
    values_to = "z"
  ) %>%
  mutate(
    AM = factor(AM, levels = AM_ORDER),
    module = factor(module, levels = rev(MODULE_ORDER))
  )

HEAT_COLORS <- c(
  "#3B6FA1", "#8FAFC6", "#E7EDF0", "#F8F7F3",
  "#F8F7F3", "#F2E3DE", "#D98E7A", "#C45A45"
)
HEAT_VALUES <- c(0.000, 0.300, 0.375, 0.450, 0.550, 0.625, 0.700, 1.000)

p_heat <- ggplot(heatmap_long, aes(x = AM, y = module, fill = z)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_gradientn(
    colors = HEAT_COLORS,
    values = HEAT_VALUES,
    limits = c(-2, 2),
    oob = squish,
    breaks = c(-2, 0, 2),
    name = "Relative module activity (z-score)"
  ) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  guides(fill = guide_colorbar(
    title.position = "bottom",
    title.hjust = 0.5,
    barwidth = grid::unit(38, "mm"),
    barheight = grid::unit(3.4, "mm"),
    ticks.colour = "#555B61",
    frame.colour = NA
  )) +
  labs(title = "Selected functional programs", x = NULL, y = NULL) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    axis.text.x = element_text(size = 10, color = "#5E646B", margin = margin(b = 7)),
    axis.text.y = element_text(size = 9.2, color = "#30343A", margin = margin(r = 5)),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 8),
    plot.margin = margin(t = 31, r = 8, b = 0, l = 2)
  )

domain_strip <- data.frame(
  x = 1,
  module = factor(rev(MODULE_ORDER), levels = rev(MODULE_ORDER)),
  domain = c(rep("AM-4 repair / remodeling", 4), rep("AM-2 immune function", 3))
)
p_strip <- ggplot(domain_strip, aes(x = x, y = module, fill = domain)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_manual(
    values = c(
      "AM-2 immune function" = AM2_COLOR,
      "AM-4 repair / remodeling" = AM4_COLOR
    ),
    guide = "none"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(t = 74, r = 0, b = 47, l = 0))

heatmap_panel <- p_strip + p_heat +
  plot_layout(widths = c(0.035, 1))

coords_df <- data.frame(
  AM = factor(AM_ORDER, levels = AM_ORDER),
  PC1 = coords[AM_ORDER, "PC1"],
  PC2 = coords[AM_ORDER, "PC2"]
)
coords_df$focus <- ifelse(
  coords_df$AM == "AM-2",
  "AM-2",
  ifelse(coords_df$AM == "AM-4", "AM-4", "Other")
)

x_range <- diff(range(coords_df$PC1))
y_range <- diff(range(coords_df$PC2))
xpad <- x_range * 0.20
ypad <- y_range * 0.28

label_offsets <- data.frame(
  AM = factor(AM_ORDER, levels = AM_ORDER),
  dx = c(-0.08, -0.12, -0.06, 0.12, 0.10),
  dy = c(-0.24, -0.30, -0.26, -0.25, 0.13),
  hjust = c(1, 1, 1, 0, 0),
  vjust = c(1, 1, 1, 1, 0)
)
labels_df <- left_join(coords_df, label_offsets, by = "AM") %>%
  mutate(
    label_x = PC1 + dx,
    label_y = PC2 + dy,
    label_color = ifelse(
      AM == "AM-2",
      AM2_COLOR,
      ifelse(AM == "AM-4", AM4_COLOR, "#626970")
    ),
    label_face = ifelse(AM %in% c("AM-2", "AM-4"), "bold", "plain")
  )

p_pca <- ggplot(coords_df, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, color = "#D8DBDE", linewidth = 0.45) +
  geom_vline(xintercept = 0, color = "#D8DBDE", linewidth = 0.45) +
  geom_point(
    data = filter(coords_df, AM %in% c("AM-2", "AM-4")),
    aes(color = focus),
    size = 9,
    alpha = 0.12,
    stroke = 0
  ) +
  geom_point(
    aes(fill = focus, size = focus),
    shape = 21,
    color = "white",
    stroke = 0.8
  ) +
  scale_color_manual(
    values = c("AM-2" = AM2_COLOR, "AM-4" = AM4_COLOR, "Other" = GREY_COLOR),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c("AM-2" = AM2_COLOR, "AM-4" = AM4_COLOR, "Other" = GREY_COLOR),
    guide = "none"
  ) +
  scale_size_manual(
    values = c("AM-2" = 4.7, "AM-4" = 4.7, "Other" = 3.6),
    guide = "none"
  ) +
  geom_text(
    data = labels_df,
    aes(
      x = label_x, y = label_y, label = AM,
      hjust = hjust, vjust = vjust,
      color = focus, fontface = label_face
    ),
    size = 3.4,
    show.legend = FALSE
  ) +
  annotate(
    "segment",
    x = coords["AM-2", "PC1"] - 0.18,
    y = coords["AM-2", "PC2"] + 0.76,
    xend = coords["AM-2", "PC1"] - 0.03,
    yend = coords["AM-2", "PC2"] + 0.10,
    color = AM2_COLOR,
    linewidth = 0.55
  ) +
  annotate(
    "text",
    x = coords["AM-2", "PC1"] - 0.18,
    y = coords["AM-2", "PC2"] + 0.82,
    label = "IMMUNE REGULATION /\nT-CELL COORDINATION",
    color = AM2_COLOR,
    fontface = "bold",
    size = 2.9,
    hjust = 1,
    vjust = 0
  ) +
  annotate(
    "segment",
    x = coords["AM-4", "PC1"] + 0.12,
    y = coords["AM-4", "PC2"] + 0.72,
    xend = coords["AM-4", "PC1"] + 0.03,
    yend = coords["AM-4", "PC2"] + 0.10,
    color = AM4_COLOR,
    linewidth = 0.55
  ) +
  annotate(
    "text",
    x = coords["AM-4", "PC1"] + 0.12,
    y = coords["AM-4", "PC2"] + 0.78,
    label = "REPAIR / REMODELING",
    color = AM4_COLOR,
    fontface = "bold",
    size = 2.9,
    hjust = 0,
    vjust = 0
  ) +
  coord_cartesian(
    xlim = c(min(coords_df$PC1) - xpad, max(coords_df$PC1) + xpad),
    ylim = c(min(coords_df$PC2) - ypad * 0.65, max(coords_df$PC2) + ypad),
    clip = "off"
  ) +
  labs(
    title = "Functional-state space",
    subtitle = "PCA of the same 7-module matrix",
    x = sprintf("PC1 (%.1f%%)", variance[1] * 100),
    y = sprintf("PC2 (%.1f%%)", variance[2] * 100)
  ) +
  theme_classic(base_family = "Arial", base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 8.7, color = "#6D7379", hjust = 0),
    axis.line = element_line(color = "#72777C", linewidth = 0.45),
    axis.ticks = element_line(color = "#72777C", linewidth = 0.4),
    axis.text = element_text(color = "#555B61", size = 8.5),
    axis.title = element_text(size = 10),
    plot.margin = margin(t = 31, r = 12, b = 5, l = 8)
  )

combined <- heatmap_panel | p_pca
combined <- combined +
  plot_layout(widths = c(1.48, 1)) +
  plot_annotation(
    title = "AM-2 and AM-4 define distinct functional endpoints",
    theme = theme(
      plot.title = element_text(
        family = "Arial", face = "bold", size = 15,
        color = DARK_COLOR, hjust = 0.5,
        margin = margin(b = 7)
      )
    )
  )

heatmap_only <- heatmap_panel +
  plot_annotation(
    title = "Selected functional programs",
    theme = theme(
      plot.title = element_text(
        family = "Arial", face = "bold", size = 14, hjust = 0.5
      )
    )
  )

pca_only <- p_pca +
  plot_annotation(
    title = "AM functional-state space",
    theme = theme(
      plot.title = element_text(
        family = "Arial", face = "bold", size = 14, hjust = 0.5
      )
    )
  )

ggsave(
  file.path(OUT_DIR, "AM2_AM4_functional_heatmap_PCA.pdf"),
  combined,
  width = 13.2,
  height = 5.65,
  units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(OUT_DIR, "AM2_AM4_functional_heatmap_PCA.svg"),
  combined,
  width = 13.2,
  height = 5.65,
  units = "in",
  device = grDevices::svg
)
ggsave(
  file.path(OUT_DIR, "AM2_AM4_functional_heatmap_PCA.png"),
  combined,
  width = 13.2,
  height = 5.65,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(OUT_DIR, "AM2_AM4_functional_module_heatmap.pdf"),
  heatmap_only,
  width = 7.8,
  height = 5.4,
  units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(OUT_DIR, "AM2_AM4_major_class_functional_PCA.pdf"),
  pca_only,
  width = 5.6,
  height = 5.4,
  units = "in",
  device = cairo_pdf
)

expected_outputs <- c(
  "AM_fine_subtype_signature_scores.csv",
  "AM_major_class_signature_scores.csv",
  "AM_major_class_signature_scores_z.csv",
  "AM_expression_modules_z_all.csv",
  "AM_GO_module_evidence_scores.csv",
  "AM_GO_modules_z.csv",
  "AM_GO_matched_term_evidence.csv",
  "AM2_AM4_functional_module_scores_z.csv",
  "AM2_AM4_functional_endpoint_scores.csv",
  "AM2_AM4_functional_PCA_coordinates.csv",
  "AM2_AM4_functional_PCA_loadings.csv",
  "AM2_AM4_functional_module_definitions.csv",
  "AM2_AM4_functional_module_selection_audit.csv",
  "AM2_AM4_functional_heatmap_PCA.pdf",
  "AM2_AM4_functional_heatmap_PCA.svg",
  "AM2_AM4_functional_heatmap_PCA.png",
  "AM2_AM4_functional_module_heatmap.pdf",
  "AM2_AM4_major_class_functional_PCA.pdf"
)
if (REBUILD_SIGNATURE_MATRIX) {
  expected_outputs <- c(
    expected_outputs,
    "AM_functional_signature_gene_audit.csv"
  )
}
output_paths <- file.path(OUT_DIR, expected_outputs)
if (any(!file.exists(output_paths)) || any(file.info(output_paths)$size == 0)) {
  stop(
    "Missing or empty output(s): ",
    paste(expected_outputs[!file.exists(output_paths) |
      is.na(file.info(output_paths)$size) |
      file.info(output_paths)$size == 0], collapse = ", ")
  )
}

summary_lines <- c(
  "# AM-2/AM-4 functional-analysis reproducibility summary",
  "",
  paste0("- Source Seurat object: `", RDS_FILE, "`"),
  paste0("- Source GO enrichment: `", ENRICH_FILE, "`"),
  paste0(
    "- Signature matrix: ", nrow(fine_signature_matrix), " fine subtypes x ",
    ncol(fine_signature_matrix), " retained signatures"
  ),
  "- Major-class aggregation: equal mean of fine subtypes within AM-1...AM-5",
  "- GO aggregation: strongest matching top term per fine subtype; equal fine-subtype mean within tissue; equal Normal/Tumor mean",
  "- Standardization: population z-score across the five major AM classes (ddof = 0)",
  "- Displayed matrix: 5 AM classes x 7 modules",
  sprintf("- PCA variance explained: PC1 %.1f%%; PC2 %.1f%%", variance[1] * 100, variance[2] * 100),
  "- PCA input is exactly the same standardized 5 x 7 matrix shown in the heatmap",
  "- PCA signs: PC1 oriented toward AM-4; PC2 oriented toward AM-2",
  "- Heatmap scale: continuous -2 to 2, with a broad near-white interval from approximately -0.5 to +0.5",
  "",
  "## Endpoint composites",
  "",
  paste(capture.output(print(endpoint_composites, row.names = FALSE)), collapse = "\n")
)
writeLines(
  summary_lines,
  file.path(OUT_DIR, "AM2_AM4_functional_analysis_summary.md"),
  useBytes = TRUE
)

message("5/5 Complete.")
message(sprintf(
  "PCA variance: PC1 %.1f%%; PC2 %.1f%%",
  variance[1] * 100,
  variance[2] * 100
))
message("Figures and intermediate tables written to: ", OUT_DIR)
