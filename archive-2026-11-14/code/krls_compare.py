import numpy as np, csv
from krls_sim import krls_fit, deriv_at
def fit_fixed(X, y, s2, lam):
    xs = X.std(0, ddof=1); ys = y.std(ddof=1); Z = (X - X.mean(0)) / xs; t = (y - y.mean()) / ys
    K = np.exp(-((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1) / s2)
    return dict(Z=Z, c=np.linalg.solve(K + lam * np.eye(len(y)), t), s2=s2, xs=xs, ys=ys)
for nm in ("line", "quad"):
    rows = list(csv.DictReader(open(f"rcheck/{nm}.csv")))
    X = np.array([[float(r["safety"]), float(r["tenure"])] for r in rows]); y = np.array([float(r["y"]) for r in rows])
    dR = np.array([float(r["d_safety"]) for r in csv.DictReader(open(f"rcheck/{nm}_old.csv"))])
    lamR = float(open(f"rcheck/{nm}_old.lambda").read())
    a = fit_fixed(X, y, 2.0, lamR); dA = deriv_at(a, a["Z"], 0)
    b = krls_fit(X, y); dB = deriv_at(b, b["Z"], 0)
    print(f"{nm}: same sigma and lambda -> max abs difference in pointwise derivatives {np.abs(dA - dR).max():.2e} (range of derivatives {dR.min():+.3f} to {dR.max():+.3f})")
    print(f"      own lambda search -> ours {b['lam']:.3f} vs official {lamR:.3f}; derivatives correlate {np.corrcoef(dB, dR)[0,1]:.5f}; max abs difference {np.abs(dB - dR).max():.4f}")
