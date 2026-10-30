# Validation

Validated against the saved manuscript inputs on 30 October 2026 with R 4.4.1.

- All 37 plotting scripts completed successfully.
- 122 output comparisons passed: CSV values at tolerance 1e-12, PDF pages
  rasterized at 72 dpi, PNG pixel comparisons, and input-checksum comparisons.
- Filesystem tracing found no writes to the input data directories.
- Figure 1b and Figure 5g also ran from a relocated repository with spaces in
  its path, using copied inputs and no original project-path configuration.
- Requests for unavailable panels or missing input files fail before execution.

Upstream workflows were syntax-checked, not fully rerun. The manual and
unresolved panels are listed in [figures.md](figures.md).
