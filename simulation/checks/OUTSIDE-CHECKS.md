# Outside checks: what each one adds

Four packages supplied by Tom Geraghty on 17 and 18 Nov 2026. Their READMEs say they were run in Python by ChatGPT; they are AI checks, not the independent human re-run the paper still needs. Each folder holds the files as supplied, unchanged. This note and `python-extended-2026-11-17/NOTE.md` are ours.

All four were re-run here on 18 Nov 2026 (Python 3, numpy, scipy, pandas, statsmodels) and reproduced: power test byte for byte; influence test and the Python pilot to within 1e-12; joint-claim rates byte for byte. Each takes under ten seconds.

None of them uses the official KRLS package. They use Study 1 stand-ins built from the paper's published moments, much as ours do.

## Power test (`power-test-2026-11-18/`) — the useful one

**Question.** If performance really did fall by the 0.0184 points read from Figure 1, beyond a score of 4.75, could Study 1's data detect it?

**Design.** 30,000 datasets. Performance rises to 4.75 and then falls by a set amount. An ordinary regression with a hinge at 4.75, the turning point given in advance (which favours detection), tests whether the slope above 4.75 is below zero, one-sided, at 5%.

**Result.** With no fall, 4.8% of datasets are significant (the false-positive rate). With a fall of 0.0184, 5.8%: no better than chance. A fall of 0.08 is detected 9% of the time, 0.20 18%, 0.50 59%, and 1.00 97%. About 13 people score above 4.75.

**What it adds.** Even a real fall of the size the figure implies could not be told apart from no fall with these data, even with the turning point known and the right shape assumed. That is section 3.2's point (no test of a downturn is reported) put as a number: no test with these data could have found it. A fall would have to be around half a point to be found more often than not.

**Checks on its assumptions.** Slope 0.27 below the knot and residual SD 0.72 agree with Table 1 (r = .26, outcome SD .752: slope about .26, residual SD about .73). The scores are on sevenths and clipped at 5, as in ours. Limits are in its README: a kinked line, not KRLS; assumed distributions; a known knot.

## Influence test (`influence-test-2026-11-18/`) — modest

Perturbs or deletes one to five of the roughly 13 people above 4.75 and records how often the fitted slope there changes sign. Deleting five at random flips the sign in about 22% of datasets, under a real 0.0184 fall and under a flat line alike; significance rarely changes, because the slope is near zero to begin with. It shows that the top-end estimate rests on a handful of people, which the paper's own counts already show (`eldor-paper-scaffold`, T2). Not needed for the paper.

## Joint-claim test (`joint-claim-2026-11-17/`) — adds nothing new

Re-analyses the Python pilot's saved runs (no new fits; `eldor_python_extended/per_run.csv` is an identical copy of `../python-extended-2026-11-17/per_run.csv`, kept so the script runs as supplied). It counts datasets where the line is positive at 4, negative at the top, and implies a fall of at least 0.0184: straight line 19% to 25%, real fall 47% to 53%, flattening 37% to 49%. The rule was set after seeing the results, and the estimator is not KRLS (see `python-extended-2026-11-17/NOTE.md` on its kernel width). It agrees with our official-package fall-size result (`../fall-size/`): a fall that size is about twice as common under a real fall, and common without one.

## Python pilot (`python-extended-2026-11-17/`)

See its `NOTE.md`.
