# Fig4g 01 extract cells
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ.get("AM_REPO_DIR", ".")) / "python"))
from paths import data_path, code_path
import pathlib, time
import numpy as np, pandas as pd, fast_matrix_market as fmm

root = data_path("HTAN/am/GSE243013")
dd = root / "data"
T_GENES = [
    "CTLA4",
    "CD28",
    "CCR8",
    "TNFRSF4",
    "TNFRSF18",
    "TNFRSF9",
    "ICOS",
    "IL10",
    "LAYN",
    "PDCD1",
    "HAVCR2",
    "LAG3",
    "TIGIT",
    "TOX",
    "ENTPD1",
    "CXCL13",
    "GZMB",
    "PRF1",
    "IFNG",
    "GZMK",
    "IL7R",
    "TCF7",
]
t0 = time.time()


def say(s):
    print(f"[{time.time() - t0:5.0f}s] {s}", flush=True)


meta = pd.read_csv(
    dd / "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz", low_memory=False
)
bc = pd.read_csv(dd / "GSE243013_barcodes.csv.gz")["barcode"].values
genes = pd.read_csv(dd / "GSE243013_genes.csv.gz")["geneSymbol"].values
assert (
    len(bc) == len(meta) and (meta["cellID"].values == bc).all()
), "metadata 与 barcode 顺序不一致"
sub = meta["sub_cell_type"].astype(str)
is_mac = sub.str.startswith("Mφ").values
is_T = sub.str.match(r"^(CD4T_|CD8T_)").values
gidx = {g: i for i, g in enumerate(genes)}
want = ["CD80"] + [g for g in T_GENES if g in gidx]
say(f"metadata: macrophages={is_mac.sum()}, T cells={is_T.sum()}, genes={len(want)}")

clin_cols = [
    "sampleID",
    "gender",
    "age",
    "smoking_history",
    "cancer_type",
    "pre_treatment_staging",
    "anti-PD1_therapy",
    "chemotherapy",
    "targeted_therapy",
    "cycles",
    "pathological_response",
    "pathological_response_rate",
    "radiological_response",
]
meta.drop_duplicates("sampleID")[clin_cols].to_csv(
    dd / "Fig4g_sample_clinical.csv", index=False
)

fmm.PARALLELISM = 32
(data, (r, c)), shape = fmm.read_coo(
    str(dd / "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz")
)
say(f"matrix read: shape={shape}, nnz={len(data)}")
gmap = np.full(shape[1], -1, np.int32)
gmap[[gidx[g] for g in want]] = np.arange(len(want))
keep_cell = is_mac | is_T
m = keep_cell[r] & (gmap[c] >= 0)
r, c, data = r[m], gmap[c[m]], data[m]
del m
Y = np.zeros((shape[0], len(want)), np.float32)
np.add.at(Y, (r, c), data)
say("subset done")


def table(mask, cols):
    df = pd.DataFrame(
        Y[mask][:, [want.index(g) for g in cols]].astype(np.int32), columns=cols
    )
    for k, v in reversed(
        list(
            {
                "cellID": meta["cellID"],
                "sampleID": meta["sampleID"],
                "sub_cell_type": meta["sub_cell_type"],
                "total_counts": meta["total_counts"],
            }.items()
        )
    ):
        df.insert(0, k, v.values[mask])
    return df


table(is_mac, ["CD80"]).to_csv(dd / "Fig4g_macrophage_cells.csv.gz", index=False)
table(is_T, [g for g in want if g != "CD80"]).to_csv(
    dd / "Fig4g_T_cells.csv.gz", index=False
)
say(
    "written: Fig4g_macrophage_cells.csv.gz, Fig4g_T_cells.csv.gz, Fig4g_sample_clinical.csv"
)
