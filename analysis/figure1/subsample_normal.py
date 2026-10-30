# subsample normal
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ.get("AM_REPO_DIR", ".")) / "python"))
from paths import data_path, code_path
import pathlib, numpy as np, scipy.io, scipy.sparse as sp, anndata as ad

out = pathlib.Path(data_path("vibe/191_scpagwas_hg19/normal_allcelltype/export"))
out.mkdir(parents=True, exist_ok=True)
a = ad.read_h5ad(data_path("HTAN/cellhint/TRN_HTAN_normal.h5ad"), backed="r")
rng = np.random.default_rng(191)
idx = np.sort(rng.choice(a.n_obs, size=min(100000, a.n_obs), replace=False))
X = a.X[idx]
X = X if sp.issparse(X) else sp.csr_matrix(X)
scipy.io.mmwrite(str(out / "matrix.mtx"), X.T.tocoo())
a.var_names.to_series().to_csv(out / "genes.tsv", index=False, header=False)
a.obs.iloc[idx][["cell_type", "sample", "Dataset"]].to_csv(out / "cells.csv")
print("exported", X.shape)
