"""KRLS edge behaviour: what do pointwise marginal effects look like when the truth has NO downturn?
Implements Hainmueller & Hazlett (2014) defaults: standardised X and y, Gaussian kernel
k = exp(-||xi-xj||^2 / sigma2), sigma2 = number of covariates, lambda by leave-one-out error.
"""
import numpy as np, sys
from scipy.optimize import minimize_scalar

def krls_fit(X, y):
    n, d = X.shape
    xm, xs = X.mean(0), X.std(0, ddof=1)
    ym, ys = y.mean(), y.std(ddof=1)
    Z = (X - xm) / xs; t = (y - ym) / ys
    s2 = d
    D2 = ((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1)
    K = np.exp(-D2 / s2)
    w, Q = np.linalg.eigh(K)
    Qt = Q.T @ t
    def looe(loglam):
        lam = np.exp(loglam)
        inv = 1.0 / (w + lam)
        c = Q @ (inv * Qt)
        gdiag = (Q ** 2) @ inv
        return float(((c / gdiag) ** 2).sum())
    res = minimize_scalar(looe, bounds=(np.log(1e-4), np.log(n)), method="bounded")
    lam = np.exp(res.x)
    c = Q @ (Qt / (w + lam))
    return dict(Z=Z, c=c, s2=s2, lam=lam, xs=xs, ys=ys, xm=xm, K=K)

def deriv_at(fit, Zq, dim=0):
    """d fhat / d x_dim at query points Zq (standardised), returned in original units."""
    Z, c, s2 = fit["Z"], fit["c"], fit["s2"]
    D2 = ((Zq[:, None, :] - Z[None, :, :]) ** 2).sum(-1)
    Kq = np.exp(-D2 / s2)
    diff = Zq[:, None, dim] - Z[None, :, dim]
    d = (-2.0 / s2) * (Kq * diff * c[None, :]).sum(1)
    return d * fit["ys"] / fit["xs"][dim]

def one(n, r, truth, rng, rho_x=0.0, discrete=None):
    x1 = rng.standard_normal(n)
    x2 = rho_x * x1 + np.sqrt(1 - rho_x ** 2) * rng.standard_normal(n)   # the control (tenure)
    if truth == "linear":
        m = x1
    elif truth == "plateau":            # rises, then flat above the mean: no downturn anywhere
        m = np.minimum(x1, 0.0)
    elif truth == "null":
        m = np.zeros(n)
    m = (m - m.mean()) / (m.std() if m.std() > 0 else 1)
    y = r * m + np.sqrt(1 - r ** 2) * rng.standard_normal(n) if truth != "null" else rng.standard_normal(n)
    X = np.column_stack([x1, x2])
    fit = krls_fit(X, y)
    pw = deriv_at(fit, fit["Z"], 0)                       # pointwise, at each observation (what the paper plots)
    z1 = fit["Z"][:, 0]
    top = z1 >= np.quantile(z1, 0.975); bot = z1 <= np.quantile(z1, 0.025)
    mid = (z1 > np.quantile(z1, 0.25)) & (z1 < np.quantile(z1, 0.75))
    # where does the derivative curve (control held at its mean) cross zero going down, above its peak?
    g = np.linspace(-3, 3, 601); Zq = np.column_stack([g, np.zeros_like(g)])
    dg = deriv_at(fit, Zq, 0)
    pk = dg.argmax(); cross_hi = np.nan; cross_lo = np.nan
    neg = np.where(dg[pk:] < 0)[0]
    if len(neg): cross_hi = g[pk + neg[0]]
    negl = np.where(dg[:pk] < 0)[0]
    if len(negl): cross_lo = g[negl[-1]]
    return dict(lam=fit["lam"], top=pw[top].mean(), bot=pw[bot].mean(), mid=pw[mid].mean(),
                peak=g[pk], cross_hi=cross_hi, cross_lo=cross_lo)

if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    from statistics import NormalDist as ND
    scen = [("S1-like: linear truth, r=.26, N=474", 474, .26, "linear", -0.17),
            ("S2-like: linear truth, r=.19, N=301", 301, .19, "linear", 0.0),
            ("S5-like: linear truth, r=.18, N=800", 800, .18, "linear", -0.08),
            ("Plateau truth (rise then flat), r=.26, N=474", 474, .26, "plateau", -0.17),
            ("Strong linear truth, r=.60, N=474", 474, .60, "linear", -0.17),
            ("No relationship at all, N=474", 474, 0.0, "null", -0.17)]
    for name, n, r, truth, rho in scen:
        rng = np.random.default_rng(20261114)
        R = [one(n, r, truth, rng, rho) for _ in range(reps)]
        A = lambda k: np.array([q[k] for q in R])
        ch, cl = A("cross_hi"), A("cross_lo")
        print(f"\n{name}  ({reps} runs)")
        print(f"  mean marginal effect: middle half {A('mid').mean():+.3f} | top 2.5% {A('top').mean():+.3f} | bottom 2.5% {A('bot').mean():+.3f}")
        print(f"  runs with NEGATIVE mean marginal effect in the top 2.5%: {100*(A('top')<0).mean():.0f}% | in the bottom 2.5%: {100*(A('bot')<0).mean():.0f}%")
        print(f"  curve turns negative above its peak within +3 SD in {100*np.isfinite(ch).mean():.0f}% of runs; median crossing z={np.nanmedian(ch):.2f} (percentile {100*ND().cdf(np.nanmedian(ch)):.1f}); IQR z {np.nanpercentile(ch,25):.2f} to {np.nanpercentile(ch,75):.2f}")
        print(f"  curve negative below its peak within -3 SD in {100*np.isfinite(cl).mean():.0f}% of runs; median crossing z={np.nanmedian(cl):.2f}")
        print(f"  median peak of the marginal-effect curve z={np.median(A('peak')):.2f}; median lambda {np.median(A('lam')):.2f}")
