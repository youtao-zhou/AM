# 02 run cNMF
from __future__ import annotations
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ.get("AM_REPO_DIR", ".")) / "python"))
from paths import data_path, code_path

import json
import os
import time
import warnings
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "48")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "48")
os.environ.setdefault("MKL_NUM_THREADS", "48")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from joblib import Parallel, delayed
from scipy.io import mmread
from scipy.stats import fisher_exact, mannwhitneyu
from sklearn.decomposition import NMF
from sklearn.metrics import roc_auc_score

ROOT = Path(data_path("HTAN"))
OUT = ROOT / "am/cNMF/data"
SETS = ROOT / "vibe/41_am_unsupervised_method_sweep/evaluation_gene_sets.csv"

N_HVG = 2000
K_RANGE = list(range(5, 21))
SEEDS = [0, 1, 2, 3, 4]
TARGET, CONTRAST = "AM-4.CD80", "AM-5.IGHA1"
T0 = time.time()


def log_phase(phase: str, msg: str) -> None:
    el = f"{time.time() - T0:.1f}s"
    with open(OUT / "run.log", "a") as fh:
        fh.write(f"{time.strftime('%Y-%m-%dT%H:%M:%S%z')}\t{phase}\t{el}\t{msg}\n")
    print(f"[{el}] {phase} :: {msg}", flush=True)


def load():
    genes = [g.strip() for g in open(OUT / "tumorobject_genes.txt")]
    bcs = [b.strip() for b in open(OUT / "tumorobject_barcodes.txt")]
    x = (
        sp.csr_matrix(mmread(OUT / "tumorobject_counts.mtx"))
        .T.tocsr()
        .astype(np.float64)
    )
    obs = pd.read_csv(OUT / "tumorobject_cell_metadata.csv")
    if x.shape != (len(bcs), len(genes)):
        raise RuntimeError(f"shape {x.shape} != {(len(bcs), len(genes))}")
    if list(obs["barcode"].astype(str)) != bcs:
        raise RuntimeError("metadata order mismatch")
    expm1 = np.expm1(x.data)
    scale = np.minimum.reduceat(expm1, x.indptr[:-1])
    rec = np.rint(expm1 / np.repeat(scale, np.diff(x.indptr)))
    counts = sp.csr_matrix((rec, x.indices.copy(), x.indptr.copy()), shape=x.shape)
    counts.eliminate_zeros()
    return counts, obs, genes


def label_blind_hvg(counts, genes):
    lib = np.asarray(counts.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0
    norm = counts.multiply((1e4 / lib)[:, None]).tocsr()
    np.log1p(norm.data, out=norm.data)
    mean = np.asarray(norm.mean(axis=0)).ravel()
    sq = np.asarray(norm.multiply(norm).mean(axis=0)).ravel()
    var = np.maximum(sq - mean**2, 0.0)
    disp = np.full_like(mean, -np.inf)
    ok = mean > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        disp[ok] = np.log1p(var[ok] / mean[ok])
    return np.argsort(-disp)[:N_HVG]


def fit(mat, k, seed):
    m = NMF(
        n_components=k,
        init="random",
        solver="cd",
        max_iter=400,
        random_state=seed,
        tol=1e-4,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        w = m.fit_transform(mat)
    return (
        w.astype(np.float32),
        m.components_.astype(np.float32),
        float(m.reconstruction_err_),
    )


def unit(a):
    n = np.linalg.norm(a, axis=1, keepdims=True)
    n[n == 0] = 1.0
    return a / n


def stability(specs):
    v = []
    for i in range(len(specs)):
        for j in range(i + 1, len(specs)):
            c = unit(specs[i]) @ unit(specs[j]).T
            v.append(0.5 * (c.max(axis=1).mean() + c.max(axis=0).mean()))
    return float(np.mean(v)) if v else float("nan")


def auc_safe(y, s):
    return float(roc_auc_score(y, s)) if len(np.unique(y)) > 1 else float("nan")


def evaluate(usage, obs):
    frac = usage / np.maximum(usage.sum(1, keepdims=True), 1e-12)
    state = obs["cell_type"].to_numpy()
    is_htan = (obs["Dataset2"].to_numpy() == "HTAN").astype(int)
    y_state = (state == TARGET).astype(int)
    sub = (is_htan == 1) & np.isin(state, [TARGET, CONTRAST])
    y_within = (state[sub] == TARGET).astype(int)

    best = None
    for p in range(frac.shape[1]):
        v = frac[:, p]
        a_state = auc_safe(y_state, v)
        if best is None or a_state > best["AUC_target_state"]:
            means = pd.Series(v).groupby(state).mean()
            if len(np.unique(y_within)) == 2:
                _, pv = mannwhitneyu(v[sub][y_within == 1], v[sub][y_within == 0])
            else:
                pv = float("nan")
            best = {
                "program": p,
                "AUC_target_state": a_state,
                "AUC_is_HTAN": auc_safe(is_htan, v),
                "AUC_within_HTAN": auc_safe(y_within, v[sub]),
                "within_HTAN_p": float(pv),
                "target_is_top_state": bool(means.idxmax() == TARGET),
                "target_rank_among_states": int(
                    (means.sort_values(ascending=False).index == TARGET).argmax() + 1
                ),
                "n_within_test": int(sub.sum()),
            }
    return best


def enrich(top, universe, sets_df):
    uni, ts = set(universe), set(top)
    rows = []
    for (dom, name), g in sets_df.groupby(["domain", "gene_set"], sort=False):
        gs = set(g["gene"].astype(str)) & uni
        if len(gs) < 3:
            continue
        a = len(ts & gs)
        b = len(ts) - a
        c = len(gs) - a
        d = len(uni) - a - b - c
        orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append(
            {
                "domain": dom,
                "gene_set": name,
                "set_size": len(gs),
                "overlap": a,
                "odds_ratio": float(orr),
                "p_value": float(p),
            }
        )
    o = pd.DataFrame(rows).sort_values("p_value")
    m = len(o)
    o["FDR"] = np.minimum.accumulate(
        (o["p_value"].to_numpy() * m / np.arange(1, m + 1))[::-1]
    )[::-1]
    return o


def main():
    log_phase("run_started", "full tumor-object label-blind NMF")
    counts, obs, genes = load()
    log_phase("phase1_loaded", f"{counts.shape[0]} cells x {counts.shape[1]} genes")

    hidx = label_blind_hvg(counts, genes)
    hvg = [genes[i] for i in hidx]
    pd.DataFrame({"gene": hvg}).to_csv(OUT / "label_blind_HVGs.csv", index=False)

    sub = counts[:, hidx]
    lib = np.asarray(sub.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0
    norm = sub.multiply((1e4 / lib)[:, None]).tocsr()
    np.log1p(norm.data, out=norm.data)
    mat = np.asarray(norm.todense(), dtype=np.float32)
    log_phase("phase2_hvg_ready", f"{mat.shape[0]} x {mat.shape[1]}")

    rank_rows, sweep, keep = [], [], {}
    for k in K_RANGE:
        res = Parallel(n_jobs=len(SEEDS), backend="loky")(
            delayed(fit)(mat, k, s) for s in SEEDS
        )
        errs = [r[2] for r in res]
        rank_rows.append(
            {
                "K": k,
                "mean_error": float(np.mean(errs)),
                "sd_error": float(np.std(errs)),
                "stability": stability([r[1] for r in res]),
            }
        )
        ev = evaluate(res[0][0], obs)
        ev["K"] = k
        sweep.append(ev)
        keep[k] = (res[0][0], res[0][1])
        log_phase(
            "phase_rank",
            f"K={k} stab={rank_rows[-1]['stability']:.3f} "
            f"AUCstate={ev['AUC_target_state']:.3f} "
            f"AUChtan={ev['AUC_is_HTAN']:.3f} "
            f"AUCwithin={ev['AUC_within_HTAN']:.3f} "
            f"rank={ev['target_rank_among_states']}",
        )

    rank = pd.DataFrame(rank_rows)
    e = rank["mean_error"].to_numpy()
    ks = rank["K"].to_numpy()
    es = (e - e.min()) / (e.max() - e.min())
    kss = (ks - ks.min()) / (ks.max() - ks.min())
    rank["knee_score"] = kss * (es[-1] - es[0]) - (es - es[0]) * (kss[-1] - kss[0])
    chosen = int(
        rank.nlargest(3, "knee_score")
        .sort_values("stability", ascending=False)["K"]
        .iloc[0]
    )
    rank["selected"] = rank["K"] == chosen
    rank.to_csv(OUT / "nmf_rank_selection_random_init.csv", index=False)
    log_phase("phase3_rank_selected", f"label-blind K = {chosen}")

    sw = pd.DataFrame(sweep)
    sw["AUC_gap_state_minus_cohort"] = sw["AUC_target_state"] - sw["AUC_is_HTAN"]
    sw.to_csv(OUT / "confounding_sweep_all_K.csv", index=False)
    sw[sw.K == chosen].to_csv(OUT / "within_HTAN_test.csv", index=False)

    usage, spectra = keep[chosen]
    frac = usage / np.maximum(usage.sum(1, keepdims=True), 1e-12)
    state = obs["cell_type"].to_numpy()
    rows = []
    for p in range(frac.shape[1]):
        g = pd.Series(frac[:, p]).groupby(state).agg(["mean", "median", "size"])
        z = (g["mean"] - g["mean"].mean()) / (g["mean"].std(ddof=0) or 1.0)
        for s in g.index:
            rows.append(
                {
                    "program": p,
                    "cell_type": s,
                    "cells": int(g.loc[s, "size"]),
                    "mean_usage": float(g.loc[s, "mean"]),
                    "usage_zscore": float(z.loc[s]),
                }
            )
    ubs = pd.DataFrame(rows)
    ubs.to_csv(OUT / "program_usage_by_state.csv", index=False)

    cand_p = int(sw.loc[sw.K == chosen, "program"].iloc[0])
    tg = [hvg[i] for i in np.argsort(-spectra[cand_p])[:100]]
    pd.DataFrame({"program": cand_p, "rank": range(1, 101), "gene": tg}).to_csv(
        OUT / "am4_program_top100_genes.csv", index=False
    )
    enr = enrich(tg, hvg, pd.read_csv(SETS))
    enr.insert(0, "program", cand_p)
    enr.to_csv(OUT / "program_gene_set_enrichment.csv", index=False)

    dep = obs.groupby("cell_type")[["total_counts", "n_genes_by_counts"]].median()
    dep.to_csv(OUT / "depth_summary_by_state.csv")

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    a = ax[0]
    a.plot(rank["K"], rank["stability"], "s-", color="tab:green")
    a.axvline(chosen, ls="--", c="k")
    a.set_xlabel("NMF rank")
    a.set_ylabel("Cross-seed stability")
    a.set_title(f"A  Label-blind rank selection (K={chosen})")

    a = ax[1]
    a.plot(sw["K"], sw["AUC_target_state"], "o-", label=f"predict {TARGET}")
    a.plot(sw["K"], sw["AUC_is_HTAN"], "s-", label="predict is_HTAN")
    a.axhline(0.5, ls="--", c="k", lw=1)
    a.set_ylim(0.4, 1.02)
    a.set_xlabel("NMF rank")
    a.set_ylabel("AUC")
    a.set_title("B  State vs cohort")
    a.legend(fontsize=8)

    a = ax[2]
    a.plot(sw["K"], sw["AUC_within_HTAN"], "o-", color="tab:red")
    a.axhline(0.5, ls="--", c="k", lw=1)
    a.set_ylim(0.4, 1.02)
    a.set_xlabel("NMF rank")
    a.set_ylabel(f"AUC {TARGET} vs {CONTRAST}")
    a.set_title("C  Within HTAN tumor (same cohort)")

    fig.suptitle(
        f"Full tumor AM object ({mat.shape[0]} cells, current grouping)", fontsize=12
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(OUT / "Figure_tumorobject_NMF.png", dpi=200)
    fig.savefig(OUT / "Figure_tumorobject_NMF.pdf")
    plt.close(fig)

    concl = {
        "cells": int(mat.shape[0]),
        "selected_K": chosen,
        "K_with_target_as_top_state": int(sw["target_is_top_state"].sum()),
        "K_total": len(K_RANGE),
        "median_AUC_target_state": float(sw["AUC_target_state"].median()),
        "median_AUC_is_HTAN": float(sw["AUC_is_HTAN"].median()),
        "median_AUC_within_HTAN": float(sw["AUC_within_HTAN"].median()),
        "K_state_AUC_above_cohort_AUC": int(
            (sw["AUC_gap_state_minus_cohort"] > 0).sum()
        ),
        "K_within_HTAN_above_0.7": int((sw["AUC_within_HTAN"] > 0.7).sum()),
        "repair_FDR_lt_005": int(((enr.domain == "Repair") & (enr.FDR < 0.05)).sum()),
        "immune_FDR_lt_005": int(((enr.domain == "Immune") & (enr.FDR < 0.05)).sum()),
        "project_AM4_Repair_overlap": (
            int(enr.loc[enr.gene_set == "Project_AM4_Repair", "overlap"].iloc[0])
            if (enr.gene_set == "Project_AM4_Repair").any()
            else None
        ),
    }
    with open(OUT / "conclusion.json", "w") as fh:
        json.dump(concl, fh, indent=2)
    log_phase("phase4_done", json.dumps(concl))
    print("PY_ALL_DONE")


if __name__ == "__main__":
    main()
