# response trend permutation
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ.get("AM_REPO_DIR", ".")) / "python"))
from paths import data_path, code_path
import csv, pathlib, itertools, collections, math
import numpy as np

out_dir = pathlib.Path(data_path("vibe/238_fig4d_trend_redraw"))
rows = list(
    csv.DictReader(open(data_path("tsa/data/AM4_post21_NR_MPR_pCR_patient_data.csv")))
)
y = np.array([float(r["am4_pct_mac"]) for r in rows])
code = {"NR": 0, "MPR": 1, "PCR": 2}
o = np.array([code[r["response"]] for r in rows])
s = [r["delivery_source"] for r in rows]


def rank(v):
    v = np.asarray(v, float)
    order = v.argsort()
    r = np.empty(len(v))
    r[order] = np.arange(1, len(v) + 1)
    for u in np.unique(v):
        m = v == u
        r[m] = r[m].mean()
    return r


ry = rank(y)
ro_of = {k: rank(o)[o == k][0] for k in (0, 1, 2)}
n = len(y)
my = ry.mean()
mo = rank(o).mean()
sy = np.sqrt(((ry - my) ** 2).sum())
so = np.sqrt(((rank(o) - mo) ** 2).sum())


def rho_from_S(S):
    return (S - n * my * mo) / (sy * so)


obs = rho_from_S(float((ry * np.array([ro_of[k] for k in o])).sum()))
dists = []
for b in sorted(set(s)):
    idx = [i for i in range(n) if s[i] == b]
    labs = [o[i] for i in idx]
    cnt = collections.Counter()
    for perm in set(itertools.permutations(labs)) if len(idx) <= 8 else None or []:
        cnt[round(sum(ry[i] * ro_of[k] for i, k in zip(idx, perm)), 6)] += 1
    if len(idx) > 8:
        from sympy.utilities.iterables import multiset_permutations

        for perm in multiset_permutations(labs):
            cnt[round(sum(ry[i] * ro_of[k] for i, k in zip(idx, perm)), 6)] += 1
    dists.append(cnt)
    print(b, len(idx), sum(cnt.values()), "arrangements")
tot = collections.Counter({0.0: 1})
for d in dists:
    new = collections.Counter()
    for a, ca in tot.items():
        for c, cc in d.items():
            new[round(a + c, 6)] += ca * cc
    tot = new
N = sum(tot.values())
ext = sum(c for S, c in tot.items() if abs(rho_from_S(S)) >= abs(obs) - 1e-9)
print(
    f"rho = {obs:.4f}; exact two-sided within-batch permutation P = {ext}/{N} = {ext / N:.4f}"
)
open(out_dir / "exact_trend_P.txt", "w").write(
    f"rho\t{obs:.6f}\nP_exact\t{ext / N:.6f}\nextreme\t{ext}\ntotal\t{N}\n"
)
