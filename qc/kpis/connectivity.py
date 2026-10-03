from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from skimage.measure import euler_number


def connectivity(ctx) -> dict:
    lab, n = ndi.label(ctx.pore)  # 4-connectivity: conservative for 2-D percolation
    tot = ctx.pore.sum()
    out = {}
    if n == 0 or tot == 0:
        out["pore_percolating_frac_v"] = out["pore_percolating_frac_h"] = 0.0
    else:
        sizes = np.bincount(lab.ravel())
        v = np.intersect1d(lab[0], lab[-1])
        hh = np.intersect1d(lab[:, 0], lab[:, -1])
        v, hh = v[v > 0], hh[hh > 0]
        out["pore_percolating_frac_v"] = float(sizes[v].sum() / tot)
        out["pore_percolating_frac_h"] = float(sizes[hh].sum() / tot)
    out["pore_euler_density_mm2"] = float(euler_number(ctx.pore, connectivity=1) / ctx.area_mm2)
    return out
