import numpy as np
from statistics import NormalDist as ND
from krls_sim import krls_fit, deriv_at

# (a) implementation checks
rng = np.random.default_rng(1)
n = 300; X = rng.standard_normal((n, 2)); y = 0.5 * X[:, 0] + rng.standard_normal(n)
fit = krls_fit(X, y)
def fhat(Zq):
    D2 = ((Zq[:, None, :] - fit["Z"][None, :, :]) ** 2).sum(-1)
    return (np.exp(-D2 / fit["s2"]) @ fit["c"]) * fit["ys"] + y.mean()
Zq = np.column_stack([np.linspace(-2.5, 2.5, 11), np.zeros(11)]); h = 1e-5
num = (fhat(Zq + [h, 0]) - fhat(Zq - [h, 0])) / (2 * h) / fit["xs"][0]
print("derivative formula vs finite difference, max abs gap:", float(np.abs(num - deriv_at(fit, Zq, 0)).max()))
print("fitted value far from the data (z=8) minus mean of y:", float(fhat(np.array([[8.0, 0.0]]))[0] - y.mean()))

def run(name, n, gen, reps=200, seed=20261114):
    rng = np.random.default_rng(seed)
    ch, cl, pk, top, bot, pct_hi = [], [], [], [], [], []
    for _ in range(reps):
        x1, x2, yy = gen(n, rng)
        fit = krls_fit(np.column_stack([x1, x2]), yy)
        z1 = fit["Z"][:, 0]
        pw = deriv_at(fit, fit["Z"], 0)
        top.append(pw[z1 >= np.quantile(z1, .975)].mean()); bot.append(pw[z1 <= np.quantile(z1, .025)].mean())
        g = np.linspace(z1.min(), z1.max(), 401)                 # only over the observed range
        dg = deriv_at(fit, np.column_stack([g, np.zeros_like(g)]), 0)
        k = dg.argmax(); pk.append((z1 < g[k]).mean())
        neg = np.where(dg[k:] < 0)[0]
        if len(neg): ch.append(g[k + neg[0]]); pct_hi.append((z1 < g[k + neg[0]]).mean())
        else: ch.append(np.nan); pct_hi.append(np.nan)
    top, bot, pk, pct_hi = map(np.array, (top, bot, pk, pct_hi))
    print(f"\n{name} ({reps} runs)")
    print(f"  marginal-effect curve peaks at sample percentile (median): {100*np.median(pk):.0f}")
    print(f"  turns negative above the peak, inside the observed range, in {100*np.isfinite(pct_hi).mean():.0f}% of runs; median percentile of crossing {100*np.nanmedian(pct_hi):.1f} (IQR {100*np.nanpercentile(pct_hi,25):.1f} to {100*np.nanpercentile(pct_hi,75):.1f})")
    print(f"  mean pointwise effect negative in top 2.5%: {100*(top<0).mean():.0f}% of runs; in bottom 2.5%: {100*(bot<0).mean():.0f}%")

def mk(truth, r, rho=-0.17, scale=None):
    def gen(n, rng):
        z = rng.standard_normal(n)
        if scale:                                   # bounded, coarse survey score: mean, sd, lo, hi, step
            mu, sd, lo, hi, step = scale
            x = np.clip(np.round((mu + sd * z) / step) * step, lo, hi)
        else:
            x = z
        x2 = rho * z + np.sqrt(1 - rho ** 2) * rng.standard_normal(n)
        m = truth(x); m = (m - m.mean()) / m.std()
        return x, x2, r * m + np.sqrt(1 - r ** 2) * rng.standard_normal(n)
    return gen

S1 = (3.335, 0.751, 1, 5, 1/7)      # seven-item mean on 1 to 5
run("Straight line, Study 1's N, r and scale (bounded 1 to 5, steps of 1/7)", 474, mk(lambda x: x, .26, scale=S1))
run("Study 1's OWN quadratic as the truth (convex: minimum at 2.48, ever steeper above it)", 474, mk(lambda x: (x - 2.48) ** 2, .28, scale=S1))
run("Same convex truth, unbounded normal scores", 474, mk(lambda z: (0.855 + 0.751 * z) ** 2, .28))
S5 = (3.390, 0.425, 1, 5, 0.01)
run("Straight line, Study 5 service: N=800, r=.18", 800, mk(lambda x: x, .18, rho=-0.08, scale=S5))
S2 = (4.020, 1.306, 1, 7, 1/7)
run("Straight line, Study 2 in-role: N=301, r=.19, scale to 7", 301, mk(lambda x: x, .19, rho=0.0, scale=S2))
