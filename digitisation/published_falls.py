"""Recompute the falls the published figures imply, from the digitised lines in digitised_all.npy.

Adds the plotted slope up along the line (trapezium rule) to rebuild the outcome curve, then reports the rise,
the fall from its highest point to the top of the range, and the fall from the bottom of the range to its lowest
point. Needs only numpy; the paper's PDF is not required. Run: python3 digitisation/published_falls.py
Study 4 (Fig. 5) is not in digitised_all.npy; its values (top 0.011, bottom 0.07) are not reproducible from this repo yet.
"""
import os
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(here, "digitised_all.npy"), allow_pickle=True).item()
for name, v in d.items():
    xs, ys = np.asarray(v["lx"], float), np.asarray(v["ly"], float)
    o = np.argsort(xs); xs, ys = xs[o], ys[o]
    F = np.concatenate([[0], np.cumsum((ys[1:] + ys[:-1]) / 2 * np.diff(xs))])
    hi = int(np.argmax(F)); lo = int(np.argmin(F[: hi + 1]))
    print(f"{name}: range {xs[0]:.2f} to {xs[-1]:.2f}; rise {F[hi] - F[lo]:.3f}; "
          f"fall at the top {F[hi] - F[-1]:.3f}; fall at the bottom {F[0] - F[lo]:.3f}")
