"""Digitise the solid line of Fig. 5 (Study 4) from the paper's PDF.

Needs the PDF: pdfimages -png -f 10 -l 10 <paper.pdf> p10   (gives p10-000.png, 945 x 687, 300 dpi)
Then: python3 digitisation/digitise_fig5.py p10-000.png
Axis calibration from the tick marks: x 2 at column 108, 7 at 905; y 1 at row 38.5, -0.5 at row 578; axis row 593, axis column 93.
The solid black line is the paper's restricted sample (two units below a safety score of 3 omitted); the dashed grey line,
over the full sample, is not read here because it runs under the solid line above about 6.6.
Writes fig5_digitised.npz (lx, ly) next to this script.
"""
import os, sys
import numpy as np
from PIL import Image

a = np.asarray(Image.open(sys.argv[1]).convert("RGB")).astype(int); H, W, _ = a.shape
tox = lambda c: 2 + (c - 108) * 5 / (905 - 108)
toy = lambda r: 1 + (r - 38.5) * (-1.5) / (578 - 38.5)
inside = np.zeros((H, W), bool); inside[8:593 - 5, 93 + 5:W - 6] = True
dark = (a.max(2) < 70) & inside
cols = np.where(dark.any(0))[0]; th = np.array([dark[:, c].sum() for c in cols]); keep = th <= np.median(th) + 3
lx = tox(cols[keep]); ly = np.array([toy(np.where(dark[:, c])[0].mean()) for c in cols[keep]])
o = np.argsort(lx)
np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig5_digitised.npz"), lx=lx[o], ly=ly[o])
