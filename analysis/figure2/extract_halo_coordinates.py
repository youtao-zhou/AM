# extract halo coordinates
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ.get("AM_REPO_DIR", ".")) / "python"))
from paths import data_path, code_path
from pathlib import Path
import zipfile, time, json
import pandas as pd, numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

ROOT = Path(data_path(""))
O = ROOT / "vibe/275_source_data_missing_panels"
T = time.monotonic()
cols = [
    "Analysis Region",
    "Object Id",
    "XMin",
    "XMax",
    "YMin",
    "YMax",
    "all macrophage",
    "tumor cells",
    "AM2 marophagee",
    "AM4 marophagee",
    "stromal cells",
    "CD8a positive cells",
]
records = []
for panel, zname in [
    ("AM2", "Halo_archive_20260930_oldall2.zip"),
    ("AM4", "Halo_archive_20260930_oldall.zip"),
]:
    with zipfile.ZipFile(ROOT / "tsa" / zname) as z:
        names = [
            n
            for n in z.namelist()
            if "D2406578-" + panel in n and n.endswith("_object_results.csv")
        ]
        for n in names:
            job = n.split("_job")[1].split("Indica")[0]
            key = f"F2h_{panel}_job{job}"
            with z.open(n) as f:
                d = pd.read_csv(f, usecols=lambda x: x in cols)
            d["x_px"] = (d.XMin + d.XMax) / 2
            d["y_px"] = (d.YMin + d.YMax) / 2
            d.to_csv(O / (key + "_INTERNAL.csv.gz"), index=False, compression="gzip")
            counts = {
                c: int(d[c].sum())
                for c in [
                    "all macrophage",
                    "tumor cells",
                    f"{panel} marophagee",
                    "stromal cells",
                ]
                if c in d
            }
            records.append(
                dict(
                    key=key,
                    panel=panel,
                    job=job,
                    source_zip=zname,
                    member=n,
                    rows=len(d),
                    counts=counts,
                )
            )
            h, _, _ = np.histogram2d(d.x_px, d.y_px, bins=220)
            fig, ax = plt.subplots(figsize=(4, 4))
            ax.imshow(gaussian_filter(h.T, 1) > 0.3, cmap="Greys")
            ax.set_title(key)
            ax.axis("off")
            fig.savefig(O / (key + "_geometry.png"), dpi=120)
            plt.close(fig)
            print(
                key,
                len(d),
                counts,
                "seconds",
                round(time.monotonic() - T, 1),
                flush=True,
            )
            (O / "F2h_source_versions_INTERNAL.json").write_text(
                json.dumps(records, ensure_ascii=False, indent=2)
            )
