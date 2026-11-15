import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from krls_sim import krls_fit, deriv_at

MU, SD, LO, HI, STEP = 3.335, 0.751, 1, 5, 1/7      # Study 1 safety scores (Table 1)
PM, PS_ = 3.695, 0.752                              # Study 1 performance mean, SD
N, R, RHO = 474, .26, -0.17                         # Study 1's N and correlation

def gen(rng):
    z = rng.standard_normal(N)
    x = np.clip(np.round((MU + SD * z) / STEP) * STEP, LO, HI)
    t = RHO * z + np.sqrt(1 - RHO ** 2) * rng.standard_normal(N)
    m = x.copy(); a, b = m.mean(), m.std()
    y = PM + PS_ * (R * (m - a) / b + np.sqrt(1 - R ** 2) * rng.standard_normal(N))
    return x, t, y, (a, b)

def lowess(x, y, grid, frac=0.8):
    out = []; k = int(np.ceil(frac * len(x)))
    for g in grid:
        d = np.abs(x - g); h = np.sort(d)[k - 1]
        w = np.clip(1 - (d / h) ** 3, 0, None) ** 3
        A = np.column_stack([np.ones_like(x), x - g]) * np.sqrt(w)[:, None]
        out.append(np.linalg.lstsq(A, y * np.sqrt(w), rcond=None)[0][0])
    return np.array(out)

rng = np.random.default_rng(20261114); runs = []
for i in range(200):
    x, t, y, ab = gen(rng)
    fit = krls_fit(np.column_stack([x, t]), y)
    pw = deriv_at(fit, fit["Z"], 0)
    g = np.linspace(x.min(), x.max(), 200); sm = lowess(x, pw, g)
    q = np.quantile(x, .975)
    runs.append(dict(x=x, t=t, y=y, ab=ab, pw=pw, g=g, sm=sm, top=pw[x >= q].mean(), q=q, ntop=(x >= q).sum(), peak=g[sm.argmax()]))
tops = np.array([r["top"] for r in runs])
r = runs[int(np.argmin(np.abs(tops - np.median(tops))))]
print(f"top-2.5% mean effect: median {np.median(tops):+.3f}; negative in {100*(tops<0).mean():.0f}% of runs")
print(f"chosen run: top mean {r['top']:+.3f}, cut {r['q']:.2f}, n {r['ntop']}, smoothed peak at {r['peak']:.2f}, line end {r['sm'][-1]:+.3f}, dots min {r['pw'].min():+.2f} max {r['pw'].max():+.2f}")
BLUE, INK, MUTE, DOT, SURF = "#2a78d6", "#0b0b0b", "#52514e", "#9a9993", "#fcfcfb"
plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans", "axes.edgecolor": "#c9c8c2", "axes.labelcolor": MUTE,
                     "xtick.color": MUTE, "ytick.color": MUTE})
fig, ax = plt.subplots(2, 1, figsize=(6.4, 9.2), facecolor=SURF)
a, b = r["ab"]; g = np.linspace(1, 5, 200)
truth = PM + PS_ * R * (g - a) / b
slope = np.full_like(g, PS_ * R / b)
for A in ax:
    A.set_facecolor(SURF); A.spines[["top", "right"]].set_visible(False); A.grid(axis="y", color="#e6e5df", lw=0.8); A.set_axisbelow(True)
    A.set_xlim(0.9, 5.1); A.set_xlabel("Psychological safety score (1 to 5)")
ax[0].scatter(r["x"], r["y"], s=9, color=DOT, alpha=.55, lw=0)
ax[0].plot(g, truth, color=BLUE, lw=2)
ax[0].set_ylabel("Performance rating"); ax[0].set_ylim(1.4, 6.2)
ax[0].set_title("1. The made-up data: a straight line, no downturn", loc="left", color=INK, fontsize=12.5, fontweight="bold", pad=10)
ax[0].annotate("True average performance:\nthe same gain at every level of safety", xy=(4.3, truth[np.searchsorted(g, 4.3)]), xytext=(1.15, 5.45),
               color=INK, fontsize=10, arrowprops=dict(arrowstyle="-", color=MUTE, lw=.8))
ax[1].axhline(0, color=MUTE, lw=1)
ax[1].scatter(r["x"], r["pw"], s=9, color=DOT, alpha=.7, lw=0)
ax[1].plot(r["g"], r["sm"], color=INK, lw=2)
ax[1].plot(g, slope, color=BLUE, lw=2, ls=(0, (5, 3)))
ax[1].set_ylabel("Marginal effect of safety on performance")
ax[1].set_ylim(-0.5, 0.85)
ax[1].set_title("2. What the paper's method draws from the same data", loc="left", color=INK, fontsize=12.5, fontweight="bold", pad=10)
ax[1].annotate("True marginal effect (dashed):\nthe same everywhere", xy=(1.6, slope[0]), xytext=(1.05, 0.62), color=INK, fontsize=10,
               arrowprops=dict(arrowstyle="-", color=MUTE, lw=.8))
ax[1].annotate(f"Estimated (dots, smoothed line):\na hump. Dots go negative at both ends.\nThe {r['ntop']} highest scorers average {r['top']:+.2f}", xy=(4.9, r["top"]), xytext=(2.35, -0.38),
               color=INK, fontsize=10, arrowprops=dict(arrowstyle="-", color=MUTE, lw=.8))
fig.text(0.02, 0.012, "Simulation, 14 Nov 2026. 474 made-up people with Study 1's scale, sample size and correlation (r = .26).\n"
         "Estimator: KRLS with the published defaults, our own code. One run of 200: the median run for the highest scorers.", fontsize=8.2, color=MUTE, va="bottom")
fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig("krls-simulation-fig1.png", dpi=150, facecolor=SURF)
