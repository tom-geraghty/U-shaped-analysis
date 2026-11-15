import numpy as np, sys
from krls_sim import krls_fit, deriv_at
MU, SD, LO, HI, STEP = 3.335, 0.751, 1, 5, 1/7
N, RHO = 474, -0.17
def lowess(x, y, grid, frac):
    out = []; k = int(np.ceil(frac * len(x)))
    for g in grid:
        d = np.abs(x - g); h = max(np.sort(d)[k - 1], 1e-9)
        w = np.clip(1 - (d / h) ** 3, 0, None) ** 3
        A = np.column_stack([np.ones_like(x), x - g]) * np.sqrt(w)[:, None]
        out.append(np.linalg.lstsq(A, y * np.sqrt(w), rcond=None)[0][0])
    return np.array(out)
def study(name, truth, r, reps=200):
    rng = np.random.default_rng(20261114); res = []
    for _ in range(reps):
        z = rng.standard_normal(N)
        x = np.clip(np.round((MU + SD * z) / STEP) * STEP, LO, HI)
        t = RHO * z + np.sqrt(1 - RHO ** 2) * rng.standard_normal(N)
        m = truth(x); m = (m - m.mean()) / m.std()
        y = r * m + np.sqrt(1 - r ** 2) * rng.standard_normal(N)
        fit = krls_fit(np.column_stack([x, t]), y); z1 = fit["Z"][:, 0]
        pw = deriv_at(fit, fit["Z"], 0)
        g = np.linspace(x.min(), x.max(), 120)
        zg = (g - x.mean()) / x.std(ddof=1)
        cur = deriv_at(fit, np.column_stack([zg, np.zeros_like(zg)]), 0)
        row = dict(top=pw[x >= np.quantile(x, .975)].mean(), bot=pw[x <= np.quantile(x, .025)].mean(),
                   mid=pw[(x > np.quantile(x, .25)) & (x < np.quantile(x, .75))].mean(),
                   share_neg_top=(pw[x >= np.quantile(x, .975)] < 0).mean(),
                   cur_end=cur[-1], cur_start=cur[0], corr=np.corrcoef(x, y)[0, 1])
        for f in (0.8, 0.5, 0.3):
            sm = lowess(x, pw, g, f); row[f"lo{f}_end"] = sm[-1]; row[f"lo{f}_start"] = sm[0]; row[f"lo{f}_pk"] = (x < g[sm.argmax()]).mean()
        res.append(row)
    A = lambda k: np.array([q[k] for q in res])
    print(f"\n{name}  (mean corr of safety and performance {A('corr').mean():.2f}; {reps} runs)")
    print(f"  average marginal effect: middle half {A('mid').mean():+.3f}; top 2.5% {A('top').mean():+.3f}; bottom 2.5% {A('bot').mean():+.3f}")
    print(f"  NEGATIVE at the top of the observed range, share of runs:")
    print(f"    mean of the dots in the top 2.5%: {100*(A('top')<0).mean():.0f}%   | curve at mean tenure, at the maximum score: {100*(A('cur_end')<0).mean():.0f}%")
    print(f"    smoothed line through the dots at the maximum score: bandwidth 0.8: {100*(A('lo0.8_end')<0).mean():.0f}%, 0.5: {100*(A('lo0.5_end')<0).mean():.0f}%, 0.3: {100*(A('lo0.3_end')<0).mean():.0f}%")
    print(f"  NEGATIVE at the bottom: dots in bottom 2.5%: {100*(A('bot')<0).mean():.0f}% | smoothed line at the minimum (bw 0.8): {100*(A('lo0.8_start')<0).mean():.0f}%")
    print(f"  smoothed line (bw 0.8) peaks at percentile: median {100*np.median(A('lo0.8_pk')):.0f}")
    return res
study("STRAIGHT LINE truth", lambda x: x, .26)
study("PAPER'S OWN QUADRATIC as truth (min at 2.48)", lambda x: (x - 2.48) ** 2, .28)
study("RISE THEN FLAT truth (flat above the mean)", lambda x: np.minimum(x, MU), .26)
study("GENUINE DOWNTURN truth (peak at 4.0, falls after)", lambda x: -(x - 4.0) ** 2, .26)
