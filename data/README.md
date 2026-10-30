# Data

Keep data outside the repository and set `AM_DATA_DIR`, or place them here.
The relative paths in [inputs.csv](inputs.csv) are resolved under that directory.
The manifest lists the files read by the plotting scripts, with SHA256 checksums.
It does not list every input required by the upstream analyses.

The expected layout is:

```text
AM-data/
  Figure1/
    plot_inputs/
    source_data/
  HTAN/
  sp/
  tsa/
  twaslca/
  reference/
```

`Figure1/` contains the saved plotting tables from the Figure 1 code bundle.
The remaining directories retain the analysis dataset layout. `reference/` contains
reference statistics used by the Figure 4c and Figure 5 consistency checks.
Data files, expression matrices, cell coordinates and patient tables are not
included in this code repository.

Run `Rscript R/check_inputs.R` to list missing plotting inputs before drawing figures.
