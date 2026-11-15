# Exact counts and Wilson 95% intervals from the per-run CSVs (official KRLS, Study 1, default settings).
# Usage: python3 perrun_intervals.py perrun_S1_seed101.csv perrun_S1_seed102.csv
import csv, math, sys
def wilson(k, n, z=1.96):
    p = k/n; d = 1 + z*z/n; c = (p + z*z/(2*n))/d; h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return c-h, c+h
rows = [r for f in sys.argv[1:] for r in csv.DictReader(open(f))]
stats = [("top", "plot negative in top 2.5% of safety scores"),
         ("e8", "footnote-4 smoothed line (f=.8) negative at top"),
         ("e5", "smoothed line (f=.5) negative at top"),
         ("bot", "plot negative in bottom 2.5%")]
for truth in dict.fromkeys(r["truth"] for r in rows):
    R = [r for r in rows if r["truth"] == truth]; n = len(R); print(f"\n{truth} ({n} datasets)")
    for key, lab in stats:
        k = sum(float(r[key]) < 0 for r in R); lo, hi = wilson(k, n)
        print(f"  {lab}: {k}/{n} = {k/n:.0%}  (95% interval {lo:.0%} to {hi:.0%})")
