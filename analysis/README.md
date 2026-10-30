# Upstream analyses

These scripts prepare the inputs used by `figures/`. They are grouped by the
current figure number; some filenames retain earlier panel names.

| Directory | Analyses |
| --- | --- |
| `figure1/` | scPagwas, AM annotation, GRSA, cNMF, PCoA and functional modules |
| `figure2/` | Atlas-reference tables, HALO coordinates and tumor proximity |
| `figure3/` | SCENIC activity, PROGENy, T-cell processing, CellChat and spatial statistics |
| `figure4/` | Spatial statistics, response permutations, TWAS and CD80/T-cell tables |
| `figure5/` | Cohort scoring, clinical models and tissue signatures |

Run scripts from the repository root with `AM_DATA_DIR` set. They retain the
analysis dataset layout and write intermediate tables there. Use a working
copy of the data when rebuilding inputs. Python dependencies are listed in
`environment/requirements-analysis.txt`; R dependencies are listed in
`environment/r-packages.csv`.

These workflows were syntax-checked but were not all executed for this release.
They require the original matrices, metadata, reference objects and compatible
package versions. In particular:

- scPagwas uses hg19 GWAS coordinates and saved tumor/normal inputs.
- The normal subsample contains up to 100,000 cells with seed 191.
- Historical AM plots may require the 13-subcluster labels and missing palettes.
- `response_trend_permutation.py` enumerates within-batch response-label
  permutations; the full permutation space is 37,837,800 arrangements.
- The GSE243013 cell extraction script was written for roughly 60 GB of memory.
- `figure2/atlas_patient_tables.py` also prepares inputs for Figure 4f, after
  `figure4/atlas_tcell_programs.py`.

Use the figure entry points to redraw existing results without rerunning these analyses.
