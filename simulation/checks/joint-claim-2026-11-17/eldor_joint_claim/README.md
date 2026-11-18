# Exploratory joint-pattern test — Study 1

17 November 2026. Re-analysis of **previously saved** Python pilot data (600 independently generated datasets, 150 per scenario, each analysed at two Gaussian kernel widths). No additional model fits were performed for this follow-up. The original paper's Figure 1 plots marginal derivatives, not the outcome curve.

## Test

The *same* smoothed derivative curve must meet all three criteria: (i) positive derivative at psychological safety x=4, (ii) negative derivative at the maximum observed x, and (iii) implied upper-end fall of at least 0.0184 outcome-scale points (the previously digitised published fall). Rates at 0.0092 and 0.0368 are provided as threshold sensitivity checks. Binomial Wilson 95% Monte Carlo confidence intervals are descriptive, not uncertainty intervals for the original paper.

## Major limitations

* The rule was specified after the initial simulation had been run and inspected. It is **exploratory and not preregistered**. The 0.0184 cutoff is drawn from the digitised target and is uncertain.
* This rule **does not test statistical significance** or require a specific zero-crossing location, monotonic derivative between x=4 and the endpoint, or a close match to the entire published derivative curve. It is a composite *shape-and-magnitude indicator*, not a full formal test of the paper's evidential claim.
* The fitted estimator is Python Gaussian kernel ridge with fixed ridge 0.05, NOT official KRLS, and the kernel widths use raw score units. The other reviewer's official-package analysis is more relevant for package-specific false-positive rates.
* The simulations were calibrated to the same realised x-y correlation, not to equivalent local signal strength or all original covariate distributions. The true downturn scenario was a particular quadratic shape, not the full space of possible downturns.
* These are conditional Monte Carlo rates, not probabilities that any hypothesis is true. No causal interpretation is possible.
* Plateau and downturn rates are close, with overlapping Monte Carlo confidence intervals. Do not assert a meaningful difference between them.

## Reproduce

Place this directory alongside `eldor_python_extended/`, then run `python eldor_joint_claim/analyse.py`. Input: `eldor_python_extended/per_run.csv`. Output: `joint_claim_rates.csv`.

## Further analysis before publication

Specify an inferential decision rule before simulations, incorporate the zero-crossing and peak locations, evaluate whole-curve shape after sensible amplitude calibration, and run with the official R KRLS package under published settings. Obtain the raw participant-level data for direct uncertainty and influence analysis.
