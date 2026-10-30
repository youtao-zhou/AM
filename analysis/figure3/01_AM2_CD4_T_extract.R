# 01 AM2 CD4 T extract
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
library(data.table)
library(CellChat)

setwd(am_data("HTAN/am/cellchat/data"))

cellchat <- readRDS("CellChat_full_AM_full100k_T.rds")
info <- fread("AM_T_edges_with_coverage.csv.gz")
selection <- fread("AM2_CD4_representative_selection.csv")
support <- fread("AM2_CD4_patient_support.csv.gz")

T_states <- c(
  "CD4 Tn/Tcm", "CD4 Tem", "CD4 Tfh", "CD8 Tm", "CD8 CTL",
  "CD8 Tex", "Treg", "T-IFN", "Cycling T"
)
keys <- c("source", "target", "interaction_name")

edges <- as.data.table(as.table(cellchat@net$prob))
setnames(edges, c(keys, "prob"))
edges[, (keys) := lapply(.SD, as.character), .SDcols = keys]
edges[, p_cellchat := as.vector(cellchat@net$pval)]
edges[, p_plus1 := pmin(1, (1 + round(p_cellchat * 2000)) / 2001)]
edges[, q_model := p.adjust(p_plus1, "BH")]

cd4 <- edges[(source == "AM-2.CD4" & target %in% T_states) |
  (target == "AM-2.CD4" & source %in% T_states)]
ann <- setdiff(names(info), c(keys, "prob", "p_cellchat", "p_plus1", "q_model"))
cd4 <- merge(cd4, info[AM == "AM-2.CD4", c(keys, ann), with = FALSE], by = keys)
cd4 <- merge(cd4, support[, .(
  source, target, interaction_name, patient_equal_support,
  patient_equal_coexpression, patient_equal_supported_coexpression
)],
by = keys, all.x = TRUE
)
stopifnot(nrow(cd4) == 31284)

cd4[, paired_support := !is.na(patient_equal_support) &
  n_supported_patients >= 3 & patient_equal_support >= 0.2]
cd4[, support_status := fcase(is.na(patient_equal_support), "No matched-sample estimate",
  paired_support, "Supported",
  default = "Limited coverage"
)]
cd4[, main_supported := prob > 0 & q_model <= 0.05 & paired_support]
fwrite(cd4, "AM2_CD4_T_all_pairs.csv.gz")

plot_data <- merge(cd4, selection[, .(direction, interaction_name)],
  by = c("direction", "interaction_name")
)
fwrite(plot_data, "AM2_CD4_T_plot_data.csv")

old_pairs <- data.table(
  direction = rep(c("T to AM", "AM to T"), each = 6),
  ligand = c(
    "IL16", "HLA-DRB5", "HLA-DRB1", "HLA-DRA", "HLA-DPB1", "HLA-DPA1",
    "HLA-B", "HLA-A", "FN1", "LGALS9", "THBS1", "CLEC2C"
  ),
  receptor = c(rep("CD4", 6), "CD8A", "CD8A", "CD44", "CD44", "CD47", "KLRB1")
)
fwrite(
  merge(cd4, old_pairs, by = c("direction", "ligand", "receptor")),
  "Original_Figure4_pairs_updated.csv"
)
