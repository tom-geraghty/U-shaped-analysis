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
git clone https://github.com/cran/KRLS && git -C KRLS checkout f71767c && R CMD INSTALL KRLS   # KRLS 1.7-1
Rscript simulation/study1_official_perrun.R S1 100 101
Rscript simulation/study1_official_perrun.R S1 100 102
Rscript simulation/study1_matched_signal.R S1 100 101
Rscript simulation/study1_matched_signal.R S1 100 102
python3 simulation/perrun_intervals.py perrun_S1_seed101.csv perrun_S1_seed102.csv
```

The same seeds reproduce the saved files exactly. The scripts set the kernel width explicitly, because KRLS 1.7 (2026) changed its default; results are identical under KRLS 1.0-0 (2017), the release current when the paper was written.

**Checks** (`simulation/checks/`): the same run under KRLS 1.0-0 (identical results), and on a safety grid of 0.05, which is where the dots in the paper's Figure 1 sit (straight line 66%, flattening 82%, fall 79%).

**Limits.** Study 1 only (the fall-size runs below cover all five studies). The figures depend on the kernel width and the smoothing span, which the paper does not report: at four times the default width the straight-line figure falls to about 28% (`archive-2026-11-14/logs/grid_sigma.log`). The made-up data follow the paper's published descriptive statistics, not its data. A negative end to the smoothed line is a reading of a plot, not a test.

## How large a fall? (`simulation/fall-size/`, added 17 Nov 2026)

**Question.** The published figures imply falls at the top of about 0.01 to 0.08 scale points, each with a fall at the bottom of similar or greater size on the lines the paper reads (`digitisation/published_falls.py`). For Study 4 that line (Fig. 5, solid) omits the two units with safety scores below 3; the paper's dashed line over the full sample stays at about +0.2 below a score of 4, so on the full sample Study 4 has no fall at the bottom. That dashed line was read by eye from the figure and is not digitised here, because it runs among the dots and under the solid line above about 6.6. How often do data with and without a real fall produce falls like these?

**Method.** Stand-ins for all five studies, built from each study's published descriptives (Study 3's are not published and are assumed). Same KRLS settings and smoothing as above. In each dataset the smoothed line is added up to rebuild the outcome curve, as for the published figures, and the fall at the top and bottom is measured in scale points. 200 datasets per curve per study.

**Results** (share of datasets with a fall at the top at least as large as the published one):

| Study (published fall at the top) | Straight line | Real fall | Flattening |
|---|---|---|---|
| 1 (0.018) | 34% | 67% | 78% |
| 2 in-role (0.080) | 15% | 32% | 45% |
| 3 (0.025) | 24% | 65% | 66% |
| 4 (0.011) | 21% | 49% | 74% |
| 5 service (0.071) | 32% | 48% | 62% |

Where the line fell at the top, a fall at the bottom at least as large appeared in 35% to 69% of straight-line datasets, 10% to 15% under a real fall and 3% to 8% under a flattening. The published falls at the bottom are larger than any of these curves usually produce (3% to 18% of straight-line datasets reach them). Full output: `simulation/fall-size/fallsize_summary.log`. The comparisons for Studies 2 to 5 are exploratory: their stand-ins leave out the grouping and repeated measures of the real studies, Study 2's outcome scale is our construction, and Study 3's descriptives are assumed. A skewed-tenure version of Study 1 changes nothing beyond chance (`simulation/checks/tenure-skew-2026-11-18/`).

**To re-run** (about an hour run one after another; Study 5 takes most of it, and the runs can go in parallel):

```
cd simulation/fall-size
Rscript fallsize_study1.R 100 201 && Rscript fallsize_study1.R 100 202
mv fallsize_seed201.csv fallsize_S1_seed201.csv && mv fallsize_seed202.csv fallsize_S1_seed202.csv
for s in S2 S3 S4 S5; do Rscript fallsize_all.R $s 100 301; Rscript fallsize_all.R $s 100 302; done
Rscript fallsum.R S1 S2 S3 S4 S5
```

A fresh clone reproduced `fallsize_S4_seed301.csv` byte for byte on 17 Nov 2026. Exact equality was checked on the platform in `simulation/sessionInfo.txt` (R 4.3.3, reference BLAS); another linear-algebra library may change the last digits.

## Re-running everything: what you need

- R (tested on 4.3.3) and KRLS 1.7-1, pinned as above. Python 3 with numpy for the digitised falls and intervals.
- Nothing else for the simulations: they use only the paper's published statistics, which are written into the scripts. Please check those against Tables 1 to 10 of the paper; they are the assumptions most worth a second look.
- The paper's PDF, only if you want to redo the digitisation from the images. We cannot include it.
- Figure 5 (Study 4) was read on 17 Nov 2026 with `digitisation/digitise_fig5.py` (solid line, the paper's restricted sample). The published falls in `fallsum.R` are typed in from `digitisation/published_falls_output.txt`.

## The digitised figures (`digitisation/`)

Readings of the lines and dots in the paper's Figures 1 to 8, taken from the 300 dpi images embedded in the published PDF. The scripts that produced them are in `archive-2026-11-14/code/` (`digitise.py`, `digitise_all.py`, `digitise_bands.py`, `digitise_fig1.py`). They need the paper's PDF, from which the figures are extracted with `pdfimages -png`. The PDF and its images are not included here. `published_falls.py` recomputes the implied falls from the saved readings without the PDF (`published_falls_output.txt`).

## Archive (`archive-2026-11-14/`)

The earlier code and logs, unchanged: a first simulation in Python, its repetition in the KRLS package, the settings grids for kernel width and penalty, the runs matching Figure 1, and per-study runs for Studies 2 to 5 that the critique does not cite. Its own README lists the files.

## How this was produced

The simulation, digitisation and analysis were carried out with an AI model (Claude, Anthropic) and checked by re-running them in separate sessions. They have not yet been independently re-run by a person. Corrections are welcome as issues.

## Licence

Code: MIT. Results and digitised readings: CC BY 4.0.
