# atlas patient tables
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ.get("AM_REPO_DIR", ".")) / "python"))
from paths import data_path, code_path
import pathlib, datetime
import pandas as pd
from scipy.stats import wilcoxon

ROOT = pathlib.Path(data_path(""))
SP = ROOT / "sp"
out_dir = ROOT / "vibe/220_fig2e_4f_atlas_redraw"


def log_phase(p):
    s = f"{datetime.datetime.now().astimezone().isoformat(timespec='seconds')}\t{p}"
    print(s)
    open(out_dir / "run.log", "a").write(s + "\n")


pm = pd.read_csv(
    ROOT / "vibe/16_am24_relative_tumor_position/tables/patient_metrics.csv"
)
sm = pd.read_csv(
    ROOT / "vibe/16_am24_relative_tumor_position/tables/section_metrics.csv"
)
cur = pd.read_csv(SP / "AM24_tumor_fraction_patient_data.csv")
reb = pm[(pm.route == "rebuilt") & (pm.setting == "all")].set_index("patient")
assert (cur.set_index("patient").delta_pct - reb.delta_pct).abs().max() < 1e-9
leg = pm[(pm.route == "legacy") & (pm.setting == "all")].copy()
nsec = cur.set_index("patient").n_sections
leg["n_sections"] = leg.patient.map(nsec)
assert leg.n_sections.notna().all() and len(leg) == 20
leg = leg.sort_values("delta_pct", ascending=False).reset_index(drop=True)
leg.insert(0, "display_order", range(1, len(leg) + 1))
leg[
    [
        "display_order",
        "patient",
        "n_sections",
        "AM2_tumor_pct",
        "AM4_tumor_pct",
        "delta_pct",
    ]
].to_csv(SP / "AM24_tumor_fraction_patient_data_atlas.csv", index=False)
long = pd.concat(
    [
        leg[["patient", "n_sections"]].assign(
            AM="AM-2", tumor_fraction_percent=leg.AM2_tumor_pct
        ),
        leg[["patient", "n_sections"]].assign(
            AM="AM-4", tumor_fraction_percent=leg.AM4_tumor_pct
        ),
    ]
)
long.to_csv(SP / "AM24_tumor_fraction_plot_data_atlas.csv", index=False)
sm[(sm.route == "legacy") & (sm.setting == "all")][
    ["sample", "patient", "AM2_tumor_pct", "AM4_tumor_pct", "delta_pct"]
].to_csv(SP / "AM24_tumor_fraction_section_data_atlas.csv", index=False)
p = wilcoxon(leg.AM4_tumor_pct, leg.AM2_tumor_pct).pvalue
assert abs(p - 3.948e-4) < 1e-6 and (leg.delta_pct > 0).sum() == 17
pd.DataFrame(
    [
        dict(
            route="atlas",
            setting="all",
            metric="delta_pct",
            n=20,
            effect=leg.delta_pct.median(),
            positive=int((leg.delta_pct > 0).sum()),
            p=p,
        )
    ]
).to_csv(SP / "AM24_tumor_fraction_statistics_atlas.csv", index=False)
log_phase(
    f"2e atlas data: median delta {leg.delta_pct.median():.2f}, 17/20, P = {p:.2g}"
)

pg = pd.read_csv(ROOT / "vibe/219_fig4f_atlas_reference/data/patient_group_scores.csv")
pg = pg[(pg.top == 0.2) & (pg.region == "all")]
progs = ["CD8_CTL", "CD4_Tem", "CD8_Tm", "CD4_TnTcm", "Treg", "panT"]
lf = pg.melt(
    id_vars=["patient", "group", "n_spots", "n_sections"],
    value_vars=progs,
    var_name="program",
    value_name="score",
)
lf["n_spots"] = lf.n_spots.astype(int)
lf["n_sections"] = lf.n_sections.astype(int)
assert lf.patient.nunique() == 25 and len(lf) == 300
lf[["patient", "group", "n_spots", "n_sections", "program", "score"]].to_csv(
    SP / "AM4_AM2_visium_T_programme_scores_patient_data_atlas.csv", index=False
)
log_phase("4f atlas data written (25 patients)")
