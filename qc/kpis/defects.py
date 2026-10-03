"""Cracks (Sato ridges, curtaining-direction excluded), largest void, high-Z inclusions."""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from skimage.measure import regionprops
from skimage.morphology import skeletonize


def sato_response(bse_s: np.ndarray, sigmas) -> np.ndarray:
    from skimage.filters import sato

    return sato(255.0 - bse_s, sigmas=sigmas, black_ridges=False).astype(np.float16)


def curtaining_angle(img: np.ndarray, lo: float, hi: float, factor: int = 4) -> tuple[float, float]:
    """Orientation (deg from horizontal, CCW positive, y up) of straight stripes in [lo, hi] (either
    sign), from the angular distribution of Fourier power. Returns (angle, peak/median ratio)."""
    a = img[::factor, ::factor].astype(np.float32)
    n = min(a.shape)
    a = a[:n, :n] - a[:n, :n].mean()
    win = np.outer(np.hanning(n), np.hanning(n))
    P = np.abs(np.fft.fftshift(np.fft.fft2(a * win))) ** 2
    y, x = np.indices(P.shape) - n // 2
    r = np.hypot(x, y)
    keep = (r > n * 0.05) & (r < n * 0.45)
    # spectral direction is perpendicular to stripe direction; image y axis points down
    spec_ang = np.degrees(np.arctan2(-y, x))
    stripe = (spec_ang[keep] + 90 + 90) % 180 - 90
    hist, edges = np.histogram(stripe, bins=180, range=(-90, 90), weights=P[keep])
    centres = 0.5 * (edges[1:] + edges[:-1])
    sel = (np.abs(centres) >= lo) & (np.abs(centres) <= hi)
    k = np.argmax(np.where(sel, hist, -1))
    return float(centres[k]), float(hist[k] / max(np.median(hist), 1e-12))


def _angle_diff(a: float, b: float) -> float:
    d = abs(a - b) % 180
    return min(d, 180 - d)


def crack_mask(resp: np.ndarray, thr: float, cfg: dict, curtain_deg: float | None, px_scale: float = 1.0):
    """Returns (kept skeleton mask, n_cracks, total skeleton length px). px_scale = 25 nm / pixel size."""
    m = resp > thr
    lab, n = ndi.label(m, structure=np.ones((3, 3)))
    if n == 0:
        return np.zeros_like(m), 0, 0
    skel = skeletonize(m)
    sk_len = np.bincount(lab[skel], minlength=n + 1)
    area = np.bincount(lab.ravel(), minlength=n + 1)
    keep = np.zeros(n + 1, bool)
    cand = np.where((sk_len >= cfg["crack_min_length_px"] * px_scale))[0]
    cand = cand[cand > 0]
    if len(cand):
        props = {p.label: p for p in regionprops(np.where(np.isin(lab, cand), lab, 0))}
        for i in cand:
            width = area[i] / max(sk_len[i], 1)
            if width > cfg["crack_max_width_px"] * px_scale:
                continue
            if curtain_deg is not None:
                ang = np.degrees(props[i].orientation)  # skimage: angle between row axis and major axis
                ang = ((90 - ang) + 90) % 180 - 90  # -> from horizontal, CCW positive (y up)
                if _angle_diff(ang, curtain_deg) <= cfg["curtaining_exclusion_deg"]:
                    continue
            keep[i] = True
    kept = skel & keep[lab]
    return kept, int(keep.sum()), int(kept.sum())


def defects(ctx) -> dict:
    out = {}
    if ctx.crack_resp is not None and ctx.crack_threshold is not None:
        kept, n, length = crack_mask(ctx.crack_resp, ctx.crack_threshold, ctx.cfg, ctx.curtain_deg, ctx.px_len(1.0))
        out["crack_density_um_per_mm2"] = length * ctx.px_um / ctx.area_mm2
        out["crack_count"] = n
        ctx.extras["crack_skeleton"] = kept
    else:
        out["crack_density_um_per_mm2"] = np.nan
        out["crack_count"] = np.nan
    lab, n = ndi.label(ctx.pore)
    if n:
        big = np.bincount(lab.ravel())[1:].max()
        out["largest_void_ecd_um"] = float(2 * np.sqrt(big / np.pi) * ctx.px_um)
    else:
        out["largest_void_ecd_um"] = 0.0
    hz = ctx.bse_s > ctx.cfg["hiz_threshold"]
    lab, n = ndi.label(hz)
    cnt = int((np.bincount(lab.ravel())[1:] >= ctx.px_area(ctx.cfg["hiz_min_px"])).sum()) if n else 0
    out["hiZ_inclusion_count"] = cnt
    out["hiZ_inclusion_count_mm2"] = cnt / ctx.area_mm2
    out["area_mm2"] = ctx.area_mm2
    return out
