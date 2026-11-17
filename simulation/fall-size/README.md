# Size of the fall at the top and bottom (all five studies)

Added 17 Nov 2026. Measures how large a fall the footnote-4 smoothed line implies, in outcome-scale points,
in simulated datasets with a known truth, for comparison with the falls read from the published figures.

- `fallsize_study1.R`: Study 1 stand-in (seeds 201, 202; 100 datasets per truth each).
- `fallsize_all.R`: same method for any study stand-in (`Rscript fallsize_all.R S2 100 301`). Studies 2 to 5 used seeds 301 and 302.
- `fallsum.R`: summary against the published falls (`Rscript fallsum.R S1 S2 S3 S4 S5`); output in `fallsize_summary.log`.

Method: authors' KRLS package defaults (sigma = number of columns, exact, leave-one-out lambda); lowess of the pointwise
derivatives with span 0.8; fall = area below zero of that line beyond its last zero crossing (top) or before its first
(bottom). Stand-ins as in `../study1_official_perrun.R`. Study 3's descriptives are not published and are assumed.
KRLS 1.7-1 warns that the organisation dummies in the Study 3 and 4 stand-ins look categorical; they are passed as numeric
columns, as in the earlier runs. Not yet re-run by a second reviewer.
