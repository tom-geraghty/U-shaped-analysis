# Sensitivity check: skewed tenure in the Study 1 stand-in (18 Nov 2026)

Raised by the second reviewer's assumptions audit (18 Nov): our stand-ins draw tenure as a standard normal, while Study 1 reports tenure with mean 5.504 and SD 7.306 years (positively skewed, bounded at zero). KRLS standardises each column, so the mean and SD do not matter, but the skew and the bound could change the kernel's geometry.

`fallsize_skewt.R` is `../../fall-size/fallsize_all.R` with one change: when `SKEWT=1`, tenure is passed through a Gaussian copula to a gamma distribution with Study 1's mean and SD. Its realised correlation with safety is about −.15 (target −.17; the copula weakens it slightly).

    SKEWT=1 Rscript fallsize_skewt.R S1 100 201; SKEWT=1 Rscript fallsize_skewt.R S1 100 202
    SKEWT=0 Rscript fallsize_skewt.R S1 100 201   # normal tenure, same seed, for comparison

| Truth | Normal tenure (100 per truth): line negative at top / fall ≥ 0.018 | Skewed tenure (200 per truth) |
|---|---|---|
| Straight line | 65% / 33% | 59% / 31% |
| Rises then flat above the mean | 82% / 72% | 84% / 74% |
| Paper's reading (peak then fall) | 90% / 69% | 83% / 63% |

No conclusion changes. Within Monte Carlo error at these run counts (about ±7 points). Not yet re-run by a second person.
