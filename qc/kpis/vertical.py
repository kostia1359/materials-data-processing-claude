from __future__ import annotations

import numpy as np


def vertical(ctx, n_bands: int = 5) -> dict:
    """Least-squares slope of phase fraction vs depth (pp per 10 um, positive = increases downward).
    Diagnostic only (weight 0): confounded by shading and unknown frame position."""
    h = ctx.shape[0]
    centres, fp, fs = [], [], []
    for k in range(n_bands):
        sl = slice(k * h // n_bands, (k + 1) * h // n_bands)
        centres.append((sl.start + sl.stop) / 2 * ctx.px_um)
        fp.append(ctx.pore[sl].mean() * 100)
        fs.append(ctx.si[sl].mean() * 100)
    c = np.array(centres)
    return {
        "vertical_pore_slope": float(np.polyfit(c, fp, 1)[0] * 10),
        "vertical_si_slope": float(np.polyfit(c, fs, 1)[0] * 10),
    }
