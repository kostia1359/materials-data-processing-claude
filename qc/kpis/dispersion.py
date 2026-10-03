from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage.morphology import disk


def dispersion(ctx, centroids_px: np.ndarray | None) -> dict:
    q = int(round(ctx.px_len(ctx.cfg["quadrat_px"])))
    h, w = ctx.shape
    fr = []
    for r in range(0, h - q + 1, q):
        for c in range(0, w - q + 1, q):
            fr.append(ctx.si[r : r + q, c : c + q].mean())
    fr = np.array(fr)
    cv = float(fr.std() / fr.mean()) if len(fr) >= 2 and fr.mean() > 0 else np.nan
    # Clark-Evans on sized-particle centroids, with Donnelly's (1978) edge correction of the expected
    # nearest-neighbour distance: the 50-um-tall frame is small enough that the uncorrected ratio is
    # biased upward by ~10-20 % (synthetic Poisson truth gave R = 1.26 uncorrected; DECISIONS.md).
    R = np.nan
    if centroids_px is not None and len(centroids_px) >= 5:
        pts = centroids_px * ctx.px_um
        d, _ = cKDTree(pts).query(pts, k=2)
        n = len(pts)
        perim = 2 * (ctx.shape[0] + ctx.shape[1]) * ctx.px_um
        expected = 0.5 * np.sqrt(ctx.area_um2 / n) + (0.0514 + 0.041 / np.sqrt(n)) * perim / n
        R = float(d[:, 1].mean() / expected)
    # agglomerates: closing merges near-touching particles
    closed = ndi.binary_closing(ctx.si, structure=disk(2).astype(bool))
    lab, n = ndi.label(closed)
    agg = np.nan
    if n > 0:
        sizes = np.bincount(lab.ravel())[1:]
        sizes = sizes[sizes >= ctx.px_area(ctx.cfg["si_min_px_count"])]
        if len(sizes):
            agg = float(sizes[sizes > 3 * np.median(sizes)].sum() / sizes.sum())
    return {"si_quadrat_cv": cv, "si_clark_evans_R": R, "si_agglomerate_frac": agg}
