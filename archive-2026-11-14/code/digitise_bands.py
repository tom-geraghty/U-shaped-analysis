"""Read the two predicted-performance lines and their 95% bands off Figs 4, 6 and 8 at chosen safety scores."""
import numpy as np
from digitise import load
def gridrows(a, c0, c1, r0, r1):
    g = a[r0:r1, c0:c1]; m = (g.mean(2) > 228) & (g.mean(2) < 249) & (g[:, :, 2] - g[:, :, 0] >= 4)
    cnt = m.sum(1); rows = np.where(cnt > 0.25 * (c1 - c0))[0]
    out = []; 
    for r in rows:
        if not out or r - out[-1][-1] > 2: out.append([r])
        else: out[-1].append(r)
    return [float(np.mean(o)) + r0 for o in out]
P = {
 "Fig 4 (Study 3)":            dict(path="figs/fig4-000.png", cols=(120, 1205), rows=(60, 690), x=(184, 2.5, 1183.5, 5.0), y=(146, 6, 647, 2)),
 "Fig 6 (Study 4)":            dict(path="figs/fig6-000.png", cols=(95, 930),   rows=(40, 528), x=(103, 2.5, 905, 7.0),    y=(76.5, 6, 517, -2)),
 "Fig 8a (Study 5, service)":  dict(path="figs/fig8-000.png", cols=(45, 760),   rows=(30, 426), x=(108, 2.0, 628, 4.5),    y=None),
 "Fig 8b (Study 5, product.)": dict(path="figs/fig8-000.png", cols=(848, 1525), rows=(30, 426), x=(854, 2.0, 1503, 5.0),   y=(62, 600, 416, 400)),
}
def read(name, xs):
    p = P[name]; a = load(p["path"]); (c0, c1), (r0, r1) = p["cols"], p["rows"]
    ca, xa, cb, xb = p["x"]; tocol = lambda x: ca + (x - xa) * (cb - ca) / (xb - xa); tox = lambda c: xa + (c - ca) * (xb - xa) / (cb - ca)
    if p["y"] is None:
        gr = gridrows(a, c0, c1, r0, r1); print("   gridline rows found:", [round(v, 1) for v in gr]); ra, ya, rb, yb = gr[0], 10, gr[-1], 8
    else: ra, ya, rb, yb = p["y"]
    toy = lambda r: ya + (r - ra) * (yb - ya) / (rb - ra)
    g = a.mean(2); flat = np.ptp(a, axis=2) < 14
    light = flat & (g > 213) & (g < 229); darkb = flat & (g > 183) & (g < 199); over = flat & (g >= 199) & (g <= 213)
    black = a.max(2) < 75; grey = flat & (g > 100) & (g < 150)
    band = light | darkb | over | grey | black
    # extent of the drawn curves
    cols_drawn = [c for c in range(c0, c1) if (darkb[r0:r1, c] | light[r0:r1, c]).sum() > 3]
    print(f"   curves drawn from x={tox(min(cols_drawn)):.2f} to x={tox(max(cols_drawn)):.2f}")
    res = {}
    for x in xs:
        c = int(round(tocol(x))); c = min(max(c, min(cols_drawn) + 2), max(cols_drawn) - 2)
        def span(mask):
            rr = np.where(mask[r0:r1, c - 1 : c + 2].any(1))[0]
            return (toy(rr.max() + r0), toy(rr.min() + r0)) if len(rr) else (np.nan, np.nan)
        dlo, dhi = span(darkb | over)          # +1SD band (darker)
        llo, lhi = span(light | over)          # -1SD band (lighter)
        bl = [np.where(black[r0:r1, cc])[0] for cc in range(c - 9, c + 10)]; bl = [b.mean() + r0 for b in bl if len(b)]
        hi_line = toy(np.median(bl)) if bl else np.nan
        # the solid grey line of the -1SD curve: grey pixels well inside the light band
        rr = np.where(grey[r0:r1, c - 1 : c + 2].any(1))[0] + r0
        inner = [r for r in rr if toy(r) < lhi - 0.04 * abs(lhi - llo) - 1e-9 and toy(r) > llo + 0.04 * abs(lhi - llo)]
        # remove rows that belong to the darker band's outline
        cand = [r for r in inner if not (abs(toy(r) - dlo) < 0.02 * abs(ya - yb) or abs(toy(r) - dhi) < 0.02 * abs(ya - yb))] or inner
        lo_line = toy(np.median(cand)) if cand else np.nan
        res[x] = (lo_line, llo, lhi, hi_line, dlo, dhi)
        ov = max(0.0, min(lhi, dhi) - max(llo, dlo))
        print(f"   x={x:.2f}: -1SD line {lo_line:7.2f} [{llo:7.2f}, {lhi:7.2f}] | +1SD line {hi_line:7.2f} [{dlo:7.2f}, {dhi:7.2f}] | bands overlap: {'yes' if ov > 0 else 'no'}")
    return res
if __name__ == "__main__":
    for name, xs in (("Fig 4 (Study 3)", [2.5, 3.0, 3.38, 3.71, 3.97, 4.5, 4.95]),
                     ("Fig 6 (Study 4)", [2.5, 3.4, 4.0, 4.5, 5.0, 5.29, 5.5, 6.0, 6.5, 7.0]),
                     ("Fig 8a (Study 5, service)", [2.0, 2.97, 3.39, 3.82, 4.09, 4.40, 4.95]),
                     ("Fig 8b (Study 5, product.)", [2.05, 2.97, 3.39, 3.82, 4.09, 4.40, 4.95])):
        print(name); read(name, xs)
