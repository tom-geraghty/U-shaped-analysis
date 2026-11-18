# Study 1: oracle-turning-point power experiment (18 November 2026)

## Question
If the true performance response increases to psychological safety 4.75, then falls by the approximately 0.0184 outcome points reconstructed from Eldor et al.'s published Study 1 Figure 1, how often would a statistical test detect a negative slope?

## Design
30,000 independently simulated datasets: 5,000 at each of six true upper-tail declines (0, 0.0184, 0.08, 0.20, 0.50, 1.00). Each dataset contains 474 observations. Predictor generated as a clipped/rounded Gaussian with mean target 3.335 and SD target 0.751; tenure is a correlated normal covariate (correlation approximately -0.1724), not a faithful tenure distribution. Outcome is a clipped/rounded 1–5 scale in quarter-point steps, with residual SD 0.72. The baseline conditional mean increases at slope 0.27 until x=4.75; from there it decreases linearly by the specified amount by x=5. This is a *kinked* curve, not the smooth KRLS function or the original data-generating process.

Fit OLS with intercept, x, hinge max(x-4.75,0), and tenure. Test the linear combination slope_x + slope_hinge < 0 with a one-sided t-test at alpha=0.05, treating the knot as *known in advance*. This oracle-knot design avoids turning-point selection and is favourable to detection relative to a situation where the analyst must discover the knot. However, a parametric test correctly specifying a kink is not generally a strict power upper bound on KRLS; the comparison is illustrative.

## Results
See oracle_power_results.csv. Under a true plateau, 4.82% were significantly negative; for a fall of 0.0184, 5.78% were significantly negative. The latter is effectively indistinguishable from the nominal false-positive rate with this design. At falls of 0.08, 0.20, 0.50, 1.00, detection rose to 9.20%, 18.36%, 58.66%, 96.76%. Median count x>4.75 was 12–13.

## Important limitations
- This is **not official R KRLS**, and is not a replication of the authors' inferential method.
- The upper-end 0.0184 is obtained by integrating a digitised *marginal-effect* line, not a directly observed performance difference or a confidence interval.
- Outcome and covariate distributions, conditional mean, residual variance, knot, and kink shape are assumed; only a subset of Study 1's published moments is approximately matched. Realised x-y correlations are reported in the CSV and change with fall size.
- The test is one-sided with a **known** knot. It does not account for data-driven search for a turning point or simultaneous uncertainty over a derivative curve.
- The 0.0184 value is extremely small relative to the 0.752 published outcome SD; lack of detection does not show absence of a genuine decline.
- These results estimate power **conditional on the simulated model**. They do not estimate power of the actual published KRLS analysis or establish the true original relationship.
- Simulations were run in Python by ChatGPT (GPT-6); no human independent replication.

## Reproduce
`python run.py` (NumPy, SciPy, pandas); seed 20261118. Results are deterministic for the specified environment.
