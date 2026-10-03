"""Pore chords and local thickness."""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi


def chords(mask: np.ndarray, axis: int = 1) -> np.ndarray:
    """Run lengths of True along rows (axis=1) or columns (axis=0); runs touching the frame are dropped."""
    m = mask if axis == 1 else mask.T
    h, w = m.shape
    p = np.zeros((h, w + 2), np.int8)
    p[:, 1:-1] = m
    d = np.diff(p, axis=1)
    sr, sc = np.nonzero(d == 1)
    er, ec = np.nonzero(d == -1)
    lengths = ec - sc
    ok = (sc > 0) & (ec < w)
    return lengths[ok].astype(np.float32)


def local_thickness_map(pore: np.ndarray, px_um: float, factor: int = 2) -> np.ndarray:
    """Local thickness (diameter of the largest inscribed disk covering each pixel), via openings
    by disks computed with two EDTs per radius on a x`factor` downsampled mask."""
    small = pore[::factor, ::factor]
    edt = ndi.distance_transform_edt(small)
    lt = np.zeros(small.shape, np.float32)
    rmax = int(np.ceil(edt.max())) if edt.size else 0
    radii = [r for r in [*range(1, 21), *range(22, 41, 2), 45, 50, 60, 70, 80, 100] if r <= rmax]
    for r in radii:
        seeds = edt >= r
        if not seeds.any():
            break
        opened = ndi.distance_transform_edt(~seeds) <= r
        lt[opened & small] = 2 * r * factor * px_um
    lt[small & (lt == 0)] = factor * px_um  # thinner than the smallest disk
    full = np.repeat(np.repeat(lt, factor, 0), factor, 1)[: pore.shape[0], : pore.shape[1]]
    return full


def binned_median(v: np.ndarray) -> float:
    """Median of a quantised variable, interpolated inside the bin that holds the 50 % point."""
    levels, counts = np.unique(v, return_counts=True)
    if len(levels) == 1:
        return float(levels[0])
    edges = np.concatenate([[levels[0] - (levels[1] - levels[0]) / 2], (levels[1:] + levels[:-1]) / 2,
                            [levels[-1] + (levels[-1] - levels[-2]) / 2]])
    cum = np.concatenate([[0], np.cumsum(counts)]) / counts.sum()
    return float(np.interp(0.5, cum, edges))


def pore_chords(ctx) -> tuple[dict, dict]:
    ch = chords(ctx.pore, 1) * ctx.px_um
    cv = chords(ctx.pore, 0) * ctx.px_um
    out = {
        "pore_chord_h_mean_um": float(ch.mean()) if len(ch) else np.nan,
        "pore_chord_v_mean_um": float(cv.mean()) if len(cv) else np.nan,
        "pore_chord_h_p90_um": float(np.percentile(ch, 90)) if len(ch) else np.nan,
    }
    if ctx.lt_map is not None and ctx.pore.any():
        out["pore_lt_d50_um"] = binned_median(ctx.lt_map[ctx.pore])
    else:
        out["pore_lt_d50_um"] = np.nan
    return out, {"pore_chord_h_um": ch, "pore_chord_v_um": cv}


def pore_threshold_slope(bse_s, etd_s, t1, t2, cfg, delta: float = 3.0) -> float:
    """Change of deep-pore fraction per grey level when t1 moves by ±delta (W-class diagnostic)."""
    from ..segment import pore_mask

    fr = []
    for s in (-delta, delta):
        t = t1 + s
        carbon = (bse_s >= t) & (bse_s < t2)
        p, _, _ = pore_mask(bse_s, etd_s, t, carbon, cfg)
        fr.append(p.mean())
    return float((fr[1] - fr[0]) / (2 * delta))
