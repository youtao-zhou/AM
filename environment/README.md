# Environment

The plotting scripts were run with R 4.4.1. `r-packages.csv` records the installed
versions of directly referenced packages. Rows marked `plot` are used by the
figure entry points; `analysis` rows are additional upstream dependencies.
A blank version means the package was not installed in the validation environment.
`showtext` and `sysfonts` are optional font helpers; the plotting scripts have
fallbacks when they are unavailable.

Install packages needed for the scripts you intend to run. Some upstream
dependencies, including Seurat, CellChat, scPagwas, SCENIC and Bioconductor
annotation packages, require their respective installation procedures.

Python is used only by upstream analyses. `requirements-analysis.txt` lists
direct dependencies; those workflows were not executed as part of the plotting
validation and their versions are not pinned here.

PDF export requires an R build with Cairo support. Plots use Arial when present,
with Liberation Sans as a fallback; several spatial plots use Noto Sans.
Font availability can change text placement.
