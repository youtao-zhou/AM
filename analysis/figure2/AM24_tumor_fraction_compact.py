# AM24 tumor fraction compact
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ.get("AM_REPO_DIR", ".")) / "python"))
from paths import data_path, code_path
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=data_path("sp"), help="sp directory"
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output directory; default sp/Figure/AM24_tumor_fraction_compact",
    )
    parser.add_argument("--lang", choices=["both", "cn", "en"], default="both")
    parser.add_argument(
        "--cjk-font", type=Path, default=None, help="Optional Chinese font file"
    )
    args = parser.parse_args()
    root = args.root.resolve()
    out = (
        args.out.resolve() if args.out else root / "Figure/AM24_tumor_fraction_compact"
    )
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()

    def log_phase(phase):
        s = f'{time.strftime("%Y-%m-%dT%H:%M:%S%z")}\t{phase}\t{time.monotonic()-started:.2f}s'
        print(s, flush=True)
        with (out / "run.log").open("a") as f:
            f.write(s + "\n")

    log_phase("read_packaged_data")
    paths = {
        k: root / f"AM24_tumor_fraction_{k}.csv"
        for k in ["patient_data", "section_data", "plot_data", "statistics"]
    }
    hashes = {k: hashlib.sha256(p.read_bytes()).hexdigest() for k, p in paths.items()}
    dat = {k: pd.read_csv(p, float_precision="round_trip") for k, p in paths.items()}
    d = dat["patient_data"].sort_values("patient").reset_index(drop=True)
    sections = dat["section_data"]
    frozen = dat["statistics"].iloc[0]
    cols = ["AM2_tumor_pct", "AM4_tumor_pct", "delta_pct"]
    reconstructed = sections.groupby("patient")[cols].mean().sort_index()
    wide = d.set_index("patient")[cols].sort_index()
    p = float(wilcoxon(d.AM4_tumor_pct - d.AM2_tumor_pct).pvalue)
    checks = {
        "20_paired_patients": len(d) == 20 and d.patient.nunique() == 20,
        "21_unique_sections": len(sections) == 21
        and sections["sample"].nunique() == 21,
        "same_patient_IDs": reconstructed.index.tolist() == wide.index.tolist(),
        "section_to_patient_values": np.allclose(
            reconstructed.to_numpy(), wide.to_numpy(), rtol=0, atol=1e-12
        ),
        "counts_match_sections": np.array_equal(
            d.n_sections.to_numpy(),
            sections.groupby("patient").size().reindex(d.patient).to_numpy(),
        ),
        "delta_identity": np.allclose(
            d.delta_pct, d.AM4_tumor_pct - d.AM2_tumor_pct, rtol=0, atol=1e-12
        ),
        "frozen_nominal_P": np.isclose(p, float(frozen.p), rtol=0, atol=1e-12),
        "frozen_median_delta": np.isclose(
            d.delta_pct.median(), float(frozen.effect), rtol=0, atol=1e-12
        ),
        "all_finite": np.isfinite(wide.to_numpy()).all(),
        "bounded_percentages": d[cols[:2]].ge(0).all().all()
        and d[cols[:2]].le(100).all().all(),
    }
    long = d.melt(
        id_vars=["patient", "n_sections"],
        value_vars=cols[:2],
        var_name="AM",
        value_name="tumor_fraction_percent",
    )
    long["AM"] = long.AM.map({"AM2_tumor_pct": "AM-2", "AM4_tumor_pct": "AM-4"})
    keys = ["patient", "AM"]
    old = dat["plot_data"].sort_values(keys).reset_index(drop=True)
    fresh = long.sort_values(keys).reset_index(drop=True)
    checks["40_unique_plot_rows"] = len(old) == 40 and not old.duplicated(keys).any()
    checks["plot_patient_AM_IDs"] = old[keys].equals(fresh[keys])
    checks["plot_values_match"] = np.allclose(
        old.tumor_fraction_percent, fresh.tumor_fraction_percent, rtol=0, atol=1e-12
    )
    assert all(checks.values()), checks
    langs = ["cn", "en"] if args.lang == "both" else [args.lang]
    cn = None
    if "cn" in langs:
        font = args.cjk_font or Path(
            "/usr/share/fonts/google-noto-cjk/NotoSansCJK-Regular.ttc"
        )
        if font.exists():
            cn = fm.FontProperties(fname=str(font))
        else:
            for family in [
                "Noto Sans CJK SC",
                "Microsoft YaHei",
                "SimHei",
                "PingFang SC",
            ]:
                try:
                    cn = fm.FontProperties(
                        fname=fm.findfont(family, fallback_to_default=False)
                    )
                    break
                except ValueError:
                    pass
        if cn is None:
            raise RuntimeError(
                "No Chinese font found. Supply --cjk-font FONT or use --lang en."
            )
    en = fm.FontProperties(family="DejaVu Sans")
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.linewidth": 1.1,
        }
    )
    colors = ["#237F96", "#C74940"]
    values = [d.AM2_tumor_pct.to_numpy(), d.AM4_tumor_pct.to_numpy()]
    jitter = np.random.default_rng(2201).uniform(-0.105, 0.105, size=(2, len(d)))
    log_phase("validated_draw_figures")
    for lang in langs:
        fp = cn if lang == "cn" else en
        fig, ax = plt.subplots(figsize=(3.7, 3.8))
        fig.subplots_adjust(left=0.25, right=0.96, top=0.97, bottom=0.21)
        bp = ax.boxplot(
            values,
            positions=[0, 1],
            widths=0.45,
            patch_artist=True,
            showfliers=False,
            whis=1.5,
            medianprops={"color": "#222222", "linewidth": 1.8},
            whiskerprops={"color": "#444444", "linewidth": 1.1},
            capprops={"color": "#444444", "linewidth": 1.1},
            boxprops={"linewidth": 1.1, "edgecolor": "#444444"},
        )
        for box, color in zip(bp["boxes"], colors):
            box.set_facecolor(color)
            box.set_alpha(0.22)
        for i, (v, color) in enumerate(zip(values, colors)):
            ax.scatter(
                i + jitter[i],
                v,
                s=24,
                c=color,
                alpha=0.8,
                edgecolors="white",
                linewidths=0.35,
                zorder=4,
            )
        ax.plot([0, 0, 1, 1], [35, 36, 36, 35], color="#333333", lw=1)
        ax.text(0.5, 36.5, f"P = {p:.4f}", ha="center", va="bottom", fontsize=11)
        ax.set_xticks([0, 1], ["AM-2", "AM-4"], fontsize=12)
        ax.set_xlim(-0.55, 1.55)
        ax.set_ylim(0, 40)
        ax.set_yticks([0, 10, 20, 30, 40])
        ax.set_ylabel(
            (
                "该类 AM 位于肿瘤区的比例（%）"
                if lang == "cn"
                else "Estimated AM fraction\nin tumor domain (%)"
            ),
            fontproperties=fp,
            fontsize=11,
            labelpad=7,
        )
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="x", length=0, pad=8)
        ax.tick_params(axis="y", direction="out", length=3)
        fig.text(
            0.60,
            0.068,
            "20 位患者 · 配对检验" if lang == "cn" else "20 patients · paired test",
            ha="center",
            fontproperties=fp,
            fontsize=9,
            color="#444444",
        )
        for ext in ["pdf", "png", "svg"]:
            file = out / f"AM24_tumor_fraction_compact_{lang}.{ext}"
            fig.savefig(file, dpi=220, bbox_inches="tight", pad_inches=0.10)
            checks["created_" + file.name] = (
                file.exists() and file.stat().st_size > 1000
            )
        plt.close(fig)
    for key, path in paths.items():
        checks["unchanged_" + key] = (
            hashlib.sha256(path.read_bytes()).hexdigest() == hashes[key]
        )
    pd.DataFrame([{"check": k, "passed": bool(v)} for k, v in checks.items()]).to_csv(
        out / "validation_checks.csv", index=False
    )
    pd.DataFrame(
        [{"file": str(paths[k]), "sha256": v} for k, v in hashes.items()]
    ).to_csv(out / "input_manifest.csv", index=False)
    assert all(checks.values()), checks
    result = {
        "n_patients": len(d),
        "n_sections": len(sections),
        "nominal_p": p,
        "median_paired_difference_pp": float(d.delta_pct.median()),
        "checks_passed": len(checks),
        "status": "complete",
        "note": "Exploratory paired test; reproduction from frozen intermediate data, not a refit of spatial models.",
    }
    (out / "reproduction_status.json").write_text(json.dumps(result, indent=2) + "\n")
    log_phase("complete")


if __name__ == "__main__":
    main()
