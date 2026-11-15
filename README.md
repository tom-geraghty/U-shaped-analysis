# Does performance fall at high psychological safety? Code and outputs

Code, simulation results and digitised figure readings behind a critique of:

Eldor, L., Hodor, M., & Cappelli, P. (2023). The limits of psychological safety: Nonlinear relationships with performance. *Organizational Behavior and Human Decision Processes*, 177, 104255. https://doi.org/10.1016/j.obhdp.2023.104255

The critique (Geraghty, in preparation) argues that the published figures cannot show whether performance falls at high psychological safety. This repository holds everything needed to check the two parts of that argument that rest on our own computation: the simulation and the digitised figures. None of it uses the original study data, which are not available at the link the paper gives.

## The simulation (`simulation/`)

**Question.** Eldor et al. read a fall in performance from a smoothed line through the pointwise slopes that kernel regularised least squares (KRLS; Hainmueller & Hazlett, 2014) estimates. If performance in fact rose in a straight line, with no fall at all, how often would that line still end below zero at the highest scores?

**Method.** Made-up datasets copy Study 1 (N = 474; Table 1 of the paper): safety with mean 3.335 and SD 0.751, rounded to sevenths and clipped to 1 to 5; tenure correlated −.17 with safety; performance generated from a known curve, rescaled to mean 3.695 and SD 0.752, rounded to quarters and clipped to 1 to 5. KRLS is fitted on safety and tenure with the method's original defaults (kernel width equal to the number of covariates, penalty by leave-one-out, exact solution). The statistic is a lowess smooth of the pointwise slopes (span 0.8, Stata's default; the paper does not give its span), read at the highest safety score.

**Results** (200 datasets per curve; Wilson 95% intervals):

| Curve that generated the data | Smoothed line negative at the top |
|---|---|
| Straight line (no fall) | 129/200 = 64.5% (58% to 71%) |
| Rises, then flat above the mean | 177/200 = 88.5% (83% to 92%) |
| Rises to a peak at 4.75, then falls | 157/200 = 78.5% (72% to 84%) |

The second and third rows come from `study1_matched_signal.R`, in which the noise is set for each curve so that all three give the same correlation between safety and performance (about .26). The first run, without that matching, is in `perrun_S1_seed101.csv` and `perrun_S1_seed102.csv` (straight line 64.5%, flattening 84.5%, fall 79.0%).

**To re-run** (about 10 minutes on two cores):

```
apt-get update && apt-get install -y --no-install-recommends r-base-core
git clone https://github.com/cran/KRLS && R CMD INSTALL KRLS
Rscript simulation/study1_official_perrun.R S1 100 101
Rscript simulation/study1_official_perrun.R S1 100 102
Rscript simulation/study1_matched_signal.R S1 100 101
Rscript simulation/study1_matched_signal.R S1 100 102
python3 simulation/perrun_intervals.py perrun_S1_seed101.csv perrun_S1_seed102.csv
```

The same seeds reproduce the saved files exactly. The scripts set the kernel width explicitly, because KRLS 1.7 (2026) changed its default; results are identical under KRLS 1.0-0 (2017), the release current when the paper was written.

**Checks** (`simulation/checks/`): the same run under KRLS 1.0-0 (identical results), and on a safety grid of 0.05, which is where the dots in the paper's Figure 1 sit (straight line 66%, flattening 82%, fall 79%).

**Limits.** Study 1 only. The figures depend on the kernel width and the smoothing span, which the paper does not report: at four times the default width the straight-line figure falls to about 28% (`archive-2026-11-14/logs/grid_sigma.log`). The made-up data follow the paper's published descriptive statistics, not its data. A negative end to the smoothed line is a reading of a plot, not a test.

## The digitised figures (`digitisation/`)

Readings of the lines and dots in the paper's Figures 1 to 8, taken from the 300 dpi images embedded in the published PDF. The scripts that produced them are in `archive-2026-11-14/code/` (`digitise.py`, `digitise_all.py`, `digitise_bands.py`, `digitise_fig1.py`). They need the paper's PDF, from which the figures are extracted with `pdfimages -png`. The PDF and its images are not included here.

## Archive (`archive-2026-11-14/`)

The earlier code and logs, unchanged: a first simulation in Python, its repetition in the KRLS package, the settings grids for kernel width and penalty, the runs matching Figure 1, and per-study runs for Studies 2 to 5 that the critique does not cite. Its own README lists the files.

## How this was produced

The simulation, digitisation and analysis were carried out with an AI model (Claude, Anthropic) and checked by re-running them in separate sessions. They have not yet been independently re-run by a person. Corrections are welcome as issues.

## Licence

Code: MIT. Results and digitised readings: CC BY 4.0.
