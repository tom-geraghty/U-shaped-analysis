# Test 3 — synthetic upper-tail influence sensitivity (18 November 2026)

## Design
Study 1 synthetic data: N=474; safety approximately mean 3.335, SD 0.751, discretised to seven-item increments and bounded 1–5; tenure correlated with safety at approximately −0.1724; performance bounded 1–5 and rounded to quarter points. Two data-generating processes: a plateau after psychological safety 4.75, and a genuine upper-tail decline totalling 0.0184 performance points between 4.75 and 5.0. Noise SD 0.72. There are 1,200 attempted replicates per truth; datasets with fewer than five observations above 4.75 or unidentified baseline model are excluded (four plateau, eight downturn).

For each eligible dataset, estimate an OLS segmented regression controlling for tenure, with the knot fixed at its *true* value 4.75. The estimated upper-tail slope is the sum of the pre-knot slope and hinge coefficient. Then perturb 1, 3 or 5 randomly selected observations above 4.75 by ±0.25 performance points (clipped to 1–5), remove randomly selected upper-tail observations, or remove those with the highest safety scores. The same baseline dataset is reused across interventions. A new random selection is drawn for each k. Evaluate changes in the sign of the fitted upper-tail slope and changes in one-sided 5% statistical significance.

## Findings
For the true 0.0184-point decline, randomly deleting 1, 3 and 5 upper-tail observations reverses the estimated slope sign in 9.1%, 13.9% and 22.2% of eligible datasets. Deleting the five highest-safety observations reverses sign in 27.2% of *estimable* modified datasets; 0.76% of these modified datasets lose full model rank. Raising five random upper-tail performance ratings by 0.25 points reverses sign in 16.1%; lowering them reverses sign in 17.5%. The median number of upper-tail observations is 13. The plateau condition shows a similar pattern (21.3% sign reversals for deleting five random upper-tail observations).

**Important:** These are changes in the sign of the estimated coefficient, **not** evidence that a statistically significant effect is reversed. Significance changes are much rarer (see summary.csv). Because the fitted tail coefficient is often close to zero, even small changes can flip its sign. The experiment does not reproduce the published KRLS fit or demonstrate that the original observations were influential.

## Limitations
- The fixed-knot segmented model is deliberately simpler than KRLS; the original turning point is not known to be 4.75 in individual-level outcome data.
- The true outcome distribution and residual structure are only approximated, not recovered from the paper.
- One-sided p values are exploratory and do not correct for the selection of a turning point in actual analyses.
- The percentage of sign flips is conditional on the chosen perturbations, simulated distributions and eligible samples; it is not a probability that the published result is wrong.
- The published 0.0184 value is reconstructed by integrating a digitised LOWESS curve of marginal effects; it is not a directly observed performance decline.

## Reproduction
Run `python run.py` with NumPy, SciPy and pandas. Outputs `per_run.csv` and `summary.csv`. The script has a fixed RNG seed. Results are from the Python execution in this ChatGPT session.
