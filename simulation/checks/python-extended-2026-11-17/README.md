# Extended synthetic-data pilot — Study 1 (17 November 2026)

## Scope and result
600 independently generated datasets, 150 per truth (linear, diminishing, plateau, genuine downturn); each N=474. Each dataset is fitted twice with Gaussian kernel ridge regression, widths 1 and 2 (1,200 fits), ridge=0.05. Derivatives are smoothed using LOWESS span=0.8, it=0. Figure 1 of Eldor et al. (2023) is digitised from the previously supplied curve trace. Original PDF has been consulted directly for Study 1 sample size, summary moments and plotted quantity.

The mean realised x-y correlation was 0.2594–0.2596 by scenario, matching the published 0.2595. Synthetic predictor mean ~3.33, SD ~0.74 vs published 3.335, 0.751; outcome mean ~3.683, SD ~0.731 vs published 3.695, 0.752. Tenure is simulated as a correlated normal covariate but its published mean and SD are **not** matched. The original sample's precise distribution, covariance structure and latent function are unknown.

For width 2, 30.7% of linear datasets produced negative upper endpoint slopes and 19.3% had an integrated apparent fall >=0.0184; median full-curve RMSE 0.293. Genuine downturn produced 65.3%, 54.7%, and 0.471, respectively. Linear had lower median full-curve RMSE than downturn at both widths.

## Critical qualifications
- **Not the official KRLS implementation**; no equivalence to the original Stata analysis or R KRLS package has been established. Gaussian kernel-ridge is a related but not identical estimator.
- A simulated negative slope is not evidence of a true negative population relationship. Conversely, matching a published derivative curve does not recover its underlying data or uniquely identify the true function.
- Outcome signal strengths are calibrated to the same realised correlation, but this does not ensure equivalent local signal-to-noise ratios, effect curvature or identical empirical distributions.
- The integrated fall is measured over each synthetic sample's observed predictor support, not exactly the full digitised figure domain.
- Whole-curve RMSE is in derivative units and is **not** normalized by the published curve's amplitude; matching threshold choices are not preregistered. Shape-RMSE uses a post-hoc multiplicative scale adjustment and should be viewed as exploratory.
- Results depend on ridge regularization, kernel width, LOWESS settings, discrete scale rounding, and tail sample sizes. This experiment tests two widths, one LOWESS span, and one ridge parameter, not a comprehensive sensitivity grid.
- The comparison includes simulated points within observed x support; the digitised published line has uncertainty and includes a LOWESS line whose full reproduction requires original plotting conventions.
- **Do not treat these rates as posterior probabilities or as an exact replication of the original paper.**

## Reproduction
`OPENBLAS_NUM_THREADS=1 REPS=150 python run.py` (requires numpy, pandas, scipy, statsmodels). The random seed is fixed at 314159. The script saves `per_run.csv`, `summary.csv`, `example_curves.csv`. The digitised target `digitised_all.npy` is included.

## Next step
Run an official KRLS simulation with matched original covariate distributions, outcome moments, exact kernel and regularisation, and a prespecified whole-curve matching distance. Await requested original participant-level data for direct model validation.
