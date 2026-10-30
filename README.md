# AM

Code for the multiomics analysis of alveolar macrophages in Zhou et al.
This version follows Figures 1–5 of the manuscript dated 28 October 2026.

## Run a figure

Use R 4.4.1 and run commands from the repository root. Package versions used
for the plotting scripts are listed in [environment/r-packages.csv](environment/r-packages.csv).
Inputs are supplied separately; see [data/README.md](data/README.md).

```bash
export AM_DATA_DIR=/path/to/AM-data
Rscript Figure1.R       # all available Figure 1 panels
Rscript Figure1.R b     # one panel
Rscript Figure5.R d e   # shared scripts run once
```

Results are written to `output/Figure1/` through `output/Figure5/`.
Set `AM_OUTPUT_DIR` to use another output directory. With no `AM_DATA_DIR`,
scripts look in `data/`.

## Files

| Path | Contents |
| --- | --- |
| `Figure1.R`–`Figure5.R` | Figure entry points |
| `figures/` | Panel-level plotting and statistics |
| `analysis/` | Upstream processing and analysis scripts |
| `scripts.csv` | Panel-to-script mapping |
| `data/inputs.csv` | Required plotting inputs and checksums |
| `docs/figures.md` | Panel coverage and figure-specific notes |

The 37 plotting scripts cover 42 of 51 main-figure panels. They use saved
analysis inputs and recalculate statistics where specified in the code.
Manual illustrations and image assembly are not reproduced by these scripts.
The upstream scripts require additional matrices and objects and have not
all been rerun for this release; see [analysis/README.md](analysis/README.md).

Supplementary figures and the old repository's Figure 6 are outside this release.
