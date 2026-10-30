# Figure notes

| Figure | Available panels | Missing or manually assembled |
| --- | --- | --- |
| 1 | b, d, e, g, h, i, j | a, c, f |
| 2 | a, b, c, d, e, i, j, k | f, g, h |
| 3 | a, b, c, d, e, f, g, i, j, k, l | h |
| 4 | b, c, d, e, f, g, h, i, j | a |
| 5 | b, c, d, e, f, g, h | a |

## Figure 1

- **a:** illustration; the editable source was not recovered.
- **c:** tumor UMAP code is in `analysis/figure1/tumor_umap.R`. The normal
  UMAP script and two palette objects are missing.
- **f:** Ro/e calculations are available in `analysis/figure1/STARTRAC_am.R`;
  the final color and row-order settings were not recovered.
- **h, i:** plots use the historical 13-subcluster results. The current object
  has 11 labels and cannot be substituted for that analysis.
- **i:** the saved script labels two modules “Intracellular trafficking” and
  “Tissue repair”; the assembled figure uses “Vesicle trafficking” and “Wound repair”.

## Figure 2

- **a, b:** frozen scores for 319 OAK patients. The script defaults to plot mode;
  `--mode=full` also requires the expression matrix.
- **c:** saved results for 26 sections, based on the historical 31-cell-type reference.
- **d:** the spatial gradient is recalculated. The historical mixed-model fit
  is available only as a coefficient table; its fitting code was not recovered.
- **e:** atlas-reference results for 20 paired patients.
- **f, g:** final image crops, arrows and assembly require manual work.
- **h:** coordinates are available for UT_REP01. The original gating version
  and screenshot crop remain unresolved. This image is separate from the
  10-patient quantitative cohort used in panel i.
- **j, k:** single-cell tumor-contact analysis. Panel k uses nominal P values.

## Figure 3

- **a–g:** saved SCENIC, PROGENy, CytoSig, UMAP and CellChat results. PROGENy
  values are standardized across pathways within each AM state. CellChat
  probabilities are normalized within each T-cell group.
- **h:** fluorescence composition, crops, arrows and assembly require manual work.
- **i:** saved coordinates and cell labels; scale bar is 200 pixels.
- **j, k:** 21 patients. The P = 0.032 and n = 21 text on panel j was added
  during assembly and is absent from the plotting script.
- **l:** retains the existing NR/MPR/pCR labels and stratified comparisons.

## Figure 4

- **a:** manually assembled image panel.
- **b:** saved coordinates and cell labels; scale bar is 200 pixels.
- **c, d:** 21 patients. The exact trend P in d is 1940153/37837800 ≈ 0.0513;
  its two-decimal display is 0.05, not a P < 0.05 result.
- **f:** atlas-reference analysis in 25 patients.
- **g:** the main figure uses the full Manhattan output, including noncoding genes.
- **h:** smoothing affects display only. Statistics use unsmoothed TWAS scores.
  The saved script includes a P15 title omitted from the assembled figure.
- **i, j:** 52 LUAD patients. CD80 positivity is based on RNA expression and
  macrophage annotations; it is not CD68 protein gating.

## Figure 5

- **a:** illustration; the editable source was not recovered.
- **b:** patient scores are supplied as inputs; correlations and pooled estimates
  are recalculated.
- **c–g:** 421 randomized OAK/POPLAR patients and OAK-reference scores. Panel e
  uses the 390 response-evaluable patients. Pooled Cox models are stratified by trial.
- **c:** the saved script includes continuous interaction P values omitted
  from the assembled figure.
- **f:** the final OAK/POPLAR legend differs from the saved plotting output.
- **g:** both treatment-by-score interactions enter the same model.
- **h:** 52 LUAD patients, 50 with regression measurements. The tissue score
  is `z(AM2 signature) - z(log2(CD80 O/E + 0.01))`, distinct from the bulk score.

The code retains the current analysis labels. The previously identified N12/P292
and NR/non-MPR classification questions are not resolved by this release.
The missing historical single-cell TWAS scoring for Supplementary Figure 5f
also remains outside its scope.
