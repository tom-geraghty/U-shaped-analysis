# Outside check: Study 1 in plain kernel ridge regression (Python), 17 Nov 2026

Supplied by Tom Geraghty on 17 Nov 2026 as a zip (`eldor_python_extended_2026-11-17.zip`). Who ran it is not recorded here. The files in this folder are as supplied, unchanged; this note is ours. The supplied README describes the run and its own limits.

## Re-run

`OPENBLAS_NUM_THREADS=1 REPS=150 python3 run.py` (numpy, pandas, scipy, statsmodels; about 10 seconds) reproduced `per_run.csv` on 17 Nov 2026 to within 7e-13 on every number, with every negative-endpoint flag and every "fall at least 0.0184" flag identical. `digitised_all.npy` is byte-identical to `../../../digitisation/digitised_all.npy`.

## How it relates to our simulation

It is not the official KRLS package. Three differences matter:

1. **Kernel width on raw scores.** The script uses exp(-d²/(2w²)) on unstandardised safety and tenure. KRLS standardises each column and uses exp(-d²/σ), with default σ = number of columns = 2. With safety's SD of 0.751, width w = 1 equals σ ≈ 3.5 for safety (about 1.8 times the default) and σ = 2 for tenure; w = 2 equals σ ≈ 14 for safety (about 7 times the default) and σ = 8 for tenure (4 times). Approximate, since the two columns are scaled differently.
2. **Fixed penalty** (ridge 0.05 on the unstandardised outcome) in place of KRLS's leave-one-out choice. A lighter penalty means less pull towards the mean.
3. Derivative of the fitted function on raw scores, not the KRLS pointwise derivative.

Results (straight line / real fall, line negative at the top): 40.7% / 67.3% at w = 1; 30.7% / 65.3% at w = 2. Fall at the top at least 0.0184: 26.0% / 62.0% and 19.3% / 54.7%.

Against ours (`../../fall-size/`, `../../results/`, official package, default settings): 64.5% / 78.5% negative at the top; 34% / 67% with a fall at least the published one. Our kernel-width grid (14 Nov) gave 28% for the straight line at four times the default width, close to this run's 30.7% at a width of roughly four to seven times the default. So the two agree once width and penalty are allowed for: false falls are common, and rarer as the kernel widens or the penalty lightens. This is the width limit stated in draft section 3.5, from a second implementation.

## Not to lean on

The supplied README reports that the straight line has the lowest whole-curve RMSE against Figure 1. That RMSE is in raw slope units and is not scaled to the figure. The script's own scale-adjusted measure shows that none of the four truths matches Figure 1's shape: the best-fitting multiplier is near zero for three of them (median −0.07 to 0.03) and about 0.3 for the straight line (whose fitted slopes therefore run about three times the published ones), and the scale-adjusted RMSE is about the same for all four (median 0.24 to 0.25), close to what a flat line at zero would give. So the RMSE ranking says little about shape. Our 14 Nov reading (`eldor-figure-digitisation`, section 4) found that a straight line does not reproduce Figure 1's peak or low crossing. Do not cite the RMSE result as support for a straight line.

It does not replace the second reviewer's re-run with the official package.
