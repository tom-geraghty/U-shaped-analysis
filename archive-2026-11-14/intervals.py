# Wilson 95% intervals for the saved simulation percentages (no new runs).
# Counts are rebuilt as round(percent x runs) because the logs keep only rounded percentages,
# so each count may be off by one run. Usage: python3 intervals.py
import math
def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z*z/n
    c = (p + z*z/(2*n)) / d; h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
    return c - h, c + h
rows = [  # label, percent, runs, log file
 ("Study 1, straight line, plot negative in top 2.5%",            87, 200, "ext_S1.log"),
 ("Study 1, straight line, smoothed line (f=.8) negative at top",  64, 200, "ext_S1.log"),
 ("Study 1, paper's reading, smoothed line (f=.8) negative at top",82, 200, "ext_S1.log"),
 ("Grid, width = default (2), straight line, line(.8)",            62, 150, "grid_sigma.log"),
 ("Grid, width = 1, straight line, line(.8)",                      73, 150, "grid_sigma.log"),
 ("Grid, width = 4 (2x default), straight line, line(.8)",         52, 150, "grid_sigma.log"),
 ("Grid, width = 8 (4x default), straight line, line(.8)",         28, 150, "grid_sigma.log"),
 ("Grid, width = default, paper's reading, line(.8)",              81, 150, "grid_sigma.log"),
 ("Grid, width = 8, paper's reading, line(.8)",                    77, 150, "grid_sigma.log"),
]
for lab, pct, n, f in rows:
    lo, hi = wilson(round(pct*n/100), n)
    print(f"{lab}: {pct}% of {n} -> {lo:.0%} to {hi:.0%}   [{f}]")
