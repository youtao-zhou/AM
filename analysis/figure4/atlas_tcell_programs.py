# atlas tcell programs
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ.get("AM_REPO_DIR", ".")) / "python"))
from paths import data_path, code_path
import pathlib, datetime, numpy as np, pandas as pd
from scipy import sparse, stats
from scipy.spatial import cKDTree

ROOT = pathlib.Path(data_path(""))
out = ROOT / "vibe/219_fig4f_atlas_reference"
dd = out / "data"
dd.mkdir(parents=True, exist_ok=True)
t0 = datetime.datetime.now()


def log_phase(p, note=""):
    with open(out / "run.log", "a") as fh:
        fh.write(
            f"{datetime.datetime.now().astimezone().isoformat(timespec='seconds')}\t{p}\t{(datetime.datetime.now()-t0).seconds}\t{note}\n"
        )


def zc(x):
    x = np.asarray(x, float)
    return (x - x.mean(0)) / np.maximum(x.std(0), 1e-10)


def ring(xy, r=1):
    ns = cKDTree(xy).query_ball_point(xy, r + 0.01)
    rr, cc = [], []
    for i, js in enumerate(ns):
        for j in js:
            if j != i:
                rr.append(i)
                cc.append(j)
    a = sparse.csr_matrix((np.ones(len(rr)), (rr, cc)), shape=(len(xy), len(xy)))
    n = np.asarray(a.sum(1)).ravel()
    return sparse.diags(1 / np.maximum(n, 1)) @ a, n


dom = pd.read_csv(ROOT / "sp/BANKSY/data_domain_am_cd8_26.csv.gz")
gs = pd.read_csv(ROOT / "vibe/82_am2_visium_t_states/data/program_gene_membership.csv")
gs = gs[(gs.version == "predefined") & gs.include & gs.program_eligible]
progs = {p: g.gene.tolist() for p, g in gs.groupby("program")}
progs["panT"] = ["CD3D", "CD3E", "TRAC", "CD2"]
rows = []
for sample, d in dom.groupby("sample", sort=True):
    ef = ROOT / f"vibe/82_am2_visium_t_states/spatial/{sample}_expression.csv.gz"
    wf = ROOT / f"vibe/15_tumor_states_spatial_screen/weights/{sample}.csv"
    if not (ef.exists() and wf.exists()):
        continue
    ee = pd.read_csv(ef, index_col="barcode")
    d = d.set_index("barcode")
    w = d[["AM2", "AM4"]].rename(columns={"AM2": "AM-2", "AM4": "AM-4"})
    ix = d.index.intersection(ee.index).intersection(w.index)
    d, ee, w = d.loc[ix], ee.loc[ix], w.loc[ix]
    xy = np.c_[d.array_col / 2, d.array_row * np.sqrt(3) / 2]
    A, n = ring(xy)
    lib = ee.library_size.to_numpy()
    X = ee.drop(columns="library_size")
    gx = pd.DataFrame(
        zc(np.log1p(X.to_numpy() / np.maximum(lib[:, None], 1) * 1e4)),
        index=ix,
        columns=X.columns,
    )
    sc = pd.DataFrame(
        {
            p: gx[[g for g in gg if g in gx]].mean(1)
            for p, gg in progs.items()
            if sum(g in gx for g in gg) >= 3
        },
        index=ix,
    )
    n4 = A @ w["AM-4"].to_numpy()
    n2 = A @ w["AM-2"].to_numpy()
    ok = (n >= 3) & np.isfinite(n4) & np.isfinite(n2)
    tumor = d.compartment.eq("Tumor").to_numpy()
    for top in (0.2, 0.1):
        for region in ("all", "tumor"):
            m = ok & (tumor if region == "tumor" else True)
            if m.sum() < 200:
                continue
            g4 = m & (n4 >= np.quantile(n4[m], 1 - top))
            g2 = m & (n2 >= np.quantile(n2[m], 1 - top))
            both = g4 & g2
            g4 &= ~both
            g2 &= ~both
            for grp, sel in (("AM-4-rich", g4), ("AM-2-rich", g2)):
                r = dict(
                    sample=sample,
                    patient=sample.split("_")[0],
                    top=top,
                    region=region,
                    group=grp,
                    n_spots=int(sel.sum()),
                    n_overlap_removed=int(both.sum()),
                )
                r.update({p: float(sc.loc[sel, p].mean()) for p in sc.columns})
                rows.append(r)
    log_phase("section", sample)
sec = pd.DataFrame(rows)
sec.to_csv(dd / "section_group_scores.csv", index=False)
empty = sec[sec.n_spots == 0]
log_phase(
    "empty_groups",
    f"{len(empty)} section-groups with 0 spots after overlap removal: "
    + ";".join(
        f"{a}/{b}/{c}/{e}"
        for a, b, c, e in empty[["sample", "top", "region", "group"]].itertuples(
            index=False
        )
    ),
)
sec = sec[sec.n_spots > 0]
pcols = [c for c in sec.columns if c in progs]


def wmean(g):
    return pd.Series(
        {
            **{c: np.average(g[c], weights=g.n_spots) for c in pcols},
            "n_spots": g.n_spots.sum(),
            "n_sections": len(g),
        }
    )


pat = sec.groupby(["patient", "top", "region", "group"]).apply(wmean).reset_index()
pat.to_csv(dd / "patient_group_scores.csv", index=False)
tests = []
for (top, region), g in pat.groupby(["top", "region"]):
    wide = g.pivot(index="patient", columns="group", values=pcols)
    for p in pcols:
        for ver in ("score", "minus_panT"):
            if ver == "minus_panT" and p == "panT":
                continue
            a = wide[(p, "AM-4-rich")]
            b = wide[(p, "AM-2-rich")]
            if ver == "minus_panT":
                a = a - wide[("panT", "AM-4-rich")]
                b = b - wide[("panT", "AM-2-rich")]
            k = a.notna() & b.notna()
            a, b = a[k], b[k]
            tests.append(
                dict(
                    top=top,
                    region=region,
                    version=ver,
                    program=p,
                    n=len(a),
                    median_AM2_rich=b.median(),
                    median_AM4_rich=a.median(),
                    n_AM4_lower=int((a < b).sum()),
                    P=stats.wilcoxon(a, b).pvalue if len(a) >= 5 else np.nan,
                )
            )
pd.DataFrame(tests).to_csv(dd / "tests.csv", index=False)
log_phase("done", f"patients={pat.patient.nunique()}")
print(pd.DataFrame(tests).query("version=='score'").round(4).to_string())
