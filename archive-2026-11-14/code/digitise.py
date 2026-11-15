"""Digitise a Stata marginal-effect figure: calibrate axes from tick marks, then read the line and the dots."""
import numpy as np, sys, json
from PIL import Image
from scipy import ndimage as ndi

def load(path): return np.asarray(Image.open(path).convert("RGB")).astype(int)

def runs(mask):
    idx = np.where(mask)[0]
    if len(idx) == 0: return []
    cuts = np.where(np.diff(idx) > 1)[0]
    starts = np.r_[idx[0], idx[cuts + 1]]; ends = np.r_[idx[cuts], idx[-1]]
    return [(int(s), int(e)) for s, e in zip(starts, ends)]

def axes(a):
    dark = a.max(2) < 60
    H, W = dark.shape
    # x axis: the row with the longest dark run; y axis: the column with the longest dark run
    row = int(np.argmax(dark.sum(1))); col = int(np.argmax(dark.sum(0)))
    # thickness
    r1 = row
    while r1 + 1 < H and dark[r1 + 1].sum() > 0.5 * dark[row].sum(): r1 += 1
    c0 = col
    while c0 - 1 >= 0 and dark[:, c0 - 1].sum() > 0.5 * dark[:, col].sum(): c0 -= 1
    # x ticks: dark pixels in the rows just below the axis line
    xt = [ (s + e) / 2 for s, e in runs(dark[r1 + 3 : r1 + 7].all(0)) ]
    # y ticks: dark pixels in the columns just left of the axis line
    yt = [ (s + e) / 2 for s, e in runs(dark[:, c0 - 7 : c0 - 3].all(1)) ]
    return row, col, xt, yt

if __name__ == "__main__":
    a = load(sys.argv[1]); row, col, xt, yt = axes(a)
    print("image", a.shape, "| x-axis row", row, "| y-axis col", col)
    print("x tick columns:", [round(v, 1) for v in xt])
    print("y tick rows:", [round(v, 1) for v in yt])
