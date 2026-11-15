import numpy as np, json
from scipy import ndimage as ndi
from digitise import load
a = load("figs/fig1-000.png")
X0, PXU = 108.5, 199.125          # column of x = 1.00; pixels per unit (ticks at 108.5, 507, 706, 905)
Y06, PYU = 49.0, 600.0            # row of y = .6; pixels per unit (ticks at 49, 169, 289, 409 for .6, .4, .2, 0)
tox = lambda c: 1 + (c - X0) / PXU
toy = lambda r: 0.6 - (r - Y06) / PYU
H, W, _ = a.shape
inside = np.zeros((H, W), bool); inside[20:588, 98:930] = True
dark = (a.max(2) < 70) & inside
grey = (np.abs(a[:, :, 0] - 112) < 40) & (np.ptp(a, axis=2) < 25) & inside & ~dark
# ---- the line: mean row of dark pixels per column
cols = np.where(dark.any(0))[0]
line_x = tox(cols); line_y = np.array([toy(np.where(dark[:, c])[0].mean()) for c in cols])
thick = np.array([dark[:, c].sum() for c in cols])
print(f"line: {len(cols)} columns, x from {line_x.min():.3f} to {line_x.max():.3f}; median thickness {np.median(thick):.0f}px")
# ---- the dots: disc template matching
lab, n = ndi.label(grey); areas = ndi.sum(grey, lab, range(1, n + 1))
single = np.median(areas[(areas > 30) & (areas < 140)]); r = np.sqrt(single / np.pi)
print(f"grey components {n}; typical single-dot area {single:.0f}px, radius {r:.1f}px")
R = int(np.ceil(r)); yy, xx = np.mgrid[-R:R + 1, -R:R + 1]; disc = (xx ** 2 + yy ** 2 <= r ** 2).astype(float)
score = ndi.correlate((grey | dark).astype(float), disc, mode="constant") / disc.sum()
gscore = ndi.correlate(grey.astype(float), disc, mode="constant") / disc.sum()
peak = (score == ndi.maximum_filter(score, size=2 * R - 1)) & (score > 0.9) & (gscore > 0.55)
pr, pc = np.where(peak)
dots_x, dots_y = tox(pc), toy(pr)
print(f"dots found: {len(pr)} (of 474 observations; the rest are hidden under others)")
est_hidden = grey.sum() / single
print(f"grey area is worth about {est_hidden:.0f} separate dots")
np.savez("figs/fig1_digitised.npz", line_x=line_x, line_y=line_y, dots_x=dots_x, dots_y=dots_y)
# ---- features of the line
def cross(xs, ys, lo, hi):
    m = (xs >= lo) & (xs <= hi); xs, ys = xs[m], ys[m]; s = np.where(np.diff(np.sign(ys)) != 0)[0]
    return [round(float(xs[i] - ys[i] * (xs[i + 1] - xs[i]) / (ys[i + 1] - ys[i])), 3) for i in s]
k = line_y.argmax()
print(f"line: starts {line_y[0]:+.3f} at x={line_x[0]:.2f}; peak {line_y[k]:+.3f} at x={line_x[k]:.2f}; ends {line_y[-1]:+.3f} at x={line_x[-1]:.2f}")
print("line crosses zero at x =", cross(line_x, line_y, 1, 5))
print(f"dots: max {dots_y.max():+.3f} at x={dots_x[dots_y.argmax()]:.2f}; min {dots_y.min():+.3f} at x={dots_x[dots_y.argmin()]:.2f}")
top = dots_x >= 4.74
print(f"dots at 4.75 or above: {top.sum()} visible; values {np.round(np.sort(dots_y[top]), 2).tolist()}")
print(f"distinct x positions of dots (rounded to 0.05): {len(np.unique(np.round(dots_x / 0.05) * 0.05))}; x of dots min {dots_x.min():.2f} max {dots_x.max():.2f}")
# ---- integrate the line: the fitted performance curve it implies (up to a constant)
xs = np.linspace(line_x.min(), line_x.max(), 801); ys = np.interp(xs, line_x, line_y)
F = np.concatenate([[0], np.cumsum((ys[1:] + ys[:-1]) / 2 * np.diff(xs))])
i_min, i_max = F.argmin(), F.argmax()
print(f"implied performance (relative): lowest at x={xs[i_min]:.2f}, highest at x={xs[i_max]:.2f}; rise between them {F[i_max]-F[i_min]:.3f} points")
print(f"  change from x=1.00 to the low point: {F[i_min]-F[0]:+.3f}; from the high point to x=5.00: {F[-1]-F[i_max]:+.3f}")
for x0 in (1, 2, 2.5, 3, 3.335, 4, 4.5, 4.75, 5): print(f"   x={x0}: {np.interp(x0, xs, F) - np.interp(3.335, xs, F):+.3f}", end=";")
print()
