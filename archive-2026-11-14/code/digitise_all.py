import numpy as np
from scipy import ndimage as ndi
from digitise import load
# panel: file, column range, (col_a, x_a, col_b, x_b), (row_a, y_a, row_b, y_b), axis row, axis col
P = {
 "Fig 1 (Study 1, in-role)":        ("figs/fig1-000.png", (0, 945),    (108.5, 1, 905, 5),      (49, .6, 409, 0),      593, 93),
 "Fig 2a (Study 2, service)":       ("figs/fig2-000.png", (0, 708),    (81, 1, 678, 7),         (47, .03, 417, -.02),  445, 70),
 "Fig 2b (Study 2, in-role)":       ("figs/fig2-000.png", (708, 1417), (889, 2, 1387, 7),       (29, .6, 409, -.2),    445, 778),
 "Fig 3 (Study 3, in-role)":        ("figs/fig3-000.png", (0, 945),    (170, 2.5, 905, 4),      (39, .3, 577, -.2),    593, 93),
 "Fig 7a (Study 5, service)":       ("figs/fig7-000.png", (0, 767),    (118, 2, 735, 4.5),      (61, .6, 449, -.4),    482, 75),
 "Fig 7b (Study 5, productivity)":  ("figs/fig7-000.png", (767, 1535), (886, 2, 1502, 4.5),     (58, 100, 461, -100),  482, 843),
}
def digitise(name):
    path, (c0, c1), (ca, xa, cb, xb), (ra, ya, rb, yb), arow, acol = P[name]
    a = load(path); H, W, _ = a.shape
    tox = lambda c: xa + (c - ca) * (xb - xa) / (cb - ca); toy = lambda r: ya + (r - ra) * (yb - ya) / (rb - ra)
    inside = np.zeros((H, W), bool); inside[8:arow - 5, acol + 5:c1 - 6] = True
    dark = (a.max(2) < 70) & inside
    grey = (np.abs(a[:, :, 0] - 112) < 45) & (np.ptp(a, axis=2) < 25) & inside & ~dark
    # line: thin dark structure; drop columns where the dark run is much thicker than the line (text, overlaps)
    cols = np.where(dark.any(0))[0]; th = np.array([dark[:, c].sum() for c in cols]); keep = th <= np.median(th) + 3
    lx = tox(cols[keep]); ly = np.array([toy(np.where(dark[:, c])[0].mean()) for c in cols[keep]])
    # dots
    lab, n = ndi.label(grey); areas = ndi.sum(grey, lab, range(1, n + 1))
    single = np.median(areas[(areas > 25) & (areas < 160)]); r = np.sqrt(single / np.pi); R = int(np.ceil(r))
    yy, xx = np.mgrid[-R:R + 1, -R:R + 1]; disc = (xx ** 2 + yy ** 2 <= r ** 2).astype(float)
    score = ndi.correlate((grey | dark).astype(float), disc, mode="constant") / disc.sum()
    gs = ndi.correlate(grey.astype(float), disc, mode="constant") / disc.sum()
    pk = (score == ndi.maximum_filter(score, size=2 * R + 1)) & (score > 0.9) & (gs > 0.55)
    pr, pc = np.where(pk)
    # merge detections closer than one radius
    pts = []; 
    for rr, cc in sorted(zip(pr, pc)):
        if all((rr - q[0]) ** 2 + (cc - q[1]) ** 2 > r ** 2 for q in pts): pts.append((rr, cc))
    pts = np.array(pts); dx, dy = tox(pts[:, 1]), toy(pts[:, 0])
    return dict(lx=lx, ly=ly, dx=dx, dy=dy, area_dots=grey.sum() / single, px=abs((xb - xa) / (cb - ca)))
def cross(xs, ys):
    s = np.where(np.diff(np.sign(ys)) != 0)[0]
    return [float(xs[i] - ys[i] * (xs[i + 1] - xs[i]) / (ys[i + 1] - ys[i])) for i in s if ys[i] != ys[i + 1]]
if __name__ == "__main__":
    out = {}
    for name in P:
        d = digitise(name); out[name] = d; lx, ly = d["lx"], d["ly"]; k = ly.argmax()
        xs = np.linspace(lx.min(), lx.max(), 1201); ys = np.interp(xs, lx, ly)
        F = np.concatenate([[0], np.cumsum((ys[1:] + ys[:-1]) / 2 * np.diff(xs))])
        cr = cross(xs, ys)
        print(f"\n{name}: line from x={lx.min():.2f} ({ly[0]:+.3g}) to x={lx.max():.2f} ({ly[-1]:+.3g}); peak {ly[k]:+.3g} at x={lx[k]:.2f}")
        print("   zero crossings at x =", [round(c, 2) for c in cr])
        hi = F.argmax(); lo = F[:hi + 1].argmin()
        print(f"   implied outcome: low point x={xs[lo]:.2f}, high point x={xs[hi]:.2f}; rise {F[hi]-F[lo]:.3g}; fall from high point to the top of the range {F[-1]-F[hi]:+.3g}; change from the bottom of the range to the low point {F[lo]-F[0]:+.3g}")
        dx, dy = d["dx"], d["dy"]
        print(f"   dots: {len(dx)} separate detections (grey area worth about {d['area_dots']:.0f}); effects from {dy.min():+.3g} to {dy.max():+.3g}; x from {dx.min():.2f} to {dx.max():.2f}")
        tol = 1.2 * d["px"]; base = np.floor(dx.min())
        g = []
        for step, nm in ((1/7, "sevenths"), (1/3, "thirds"), (1/6, "sixths"), (0.05, "0.05"), (0.25, "quarters"), (0.2, "fifths"), (0.1, "tenths")):
            rr = np.abs((dx - base) / step - np.round((dx - base) / step)) * step
            g.append(f"{nm} {100*(rr < tol).mean():.0f}% (chance {100*min(1, 2*tol/step):.0f}%)")
        print("   share of dots on a grid of: " + "; ".join(g))
    np.save("figs/digitised_all.npy", out, allow_pickle=True)
