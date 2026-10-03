"""Two-point statistics: autocovariance lengths, integral range, ImageRep-style SE of fractions."""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi


def autocovariance(ind: np.ndarray) -> np.ndarray:
    f = ind.astype(np.float32) - ind.mean()
    var = float(ind.mean() * (1 - ind.mean()))
    if var <= 0:
        return np.zeros_like(f)
    F = np.fft.rfft2(f)
    c = np.fft.irfft2(np.abs(F) ** 2, s=f.shape) / (f.size * var)
    return c.astype(np.float32)


def corr_length(profile: np.ndarray) -> float:
    below = np.where(profile < 1 / np.e)[0]
    if len(below) == 0:
        return float(len(profile))
    r = below[0]
    if r == 0:
        return 0.0
    # linear interpolation between r-1 and r
    a, b = profile[r - 1], profile[r]
    return float(r - 1 + (a - 1 / np.e) / (a - b)) if a != b else float(r)


def integral_range(c: np.ndarray) -> float:
    """Sum of the normalised autocovariance over the central disc out to the first zero crossing of
    its radial average. (Summing every connected positive pixel lets small positive noise percolate
    across the whole frame and inflates A2 by ~7x on synthetic truth; see VALIDATION.md.)"""
    cs = np.fft.fftshift(c)
    h, w = cs.shape
    rmax = min(h, w) // 4
    cy, cx = h // 2, w // 2
    win = cs[cy - rmax : cy + rmax + 1, cx - rmax : cx + rmax + 1]
    y, x = np.indices(win.shape) - rmax
    r = np.hypot(x, y).astype(int)
    prof = np.bincount(r.ravel(), win.ravel()) / np.maximum(np.bincount(r.ravel()), 1)
    below = np.where(prof[: rmax + 1] <= 0)[0]
    r0 = below[0] if len(below) else rmax
    return float(max(win[r < r0].sum(), 1.0))


def two_point(ctx) -> dict:
    """Lengths are reported in 25-nm reference pixels ("_px") so they survive a pixel-size change."""
    out = {}
    h, w = ctx.shape
    f = ctx.px_um / 0.025
    cp = autocovariance(ctx.pore)
    out["s2_len_pore_h_px"] = corr_length(cp[0, : w // 2]) * f
    out["s2_len_pore_v_px"] = corr_length(cp[: h // 2, 0]) * f
    a2p = integral_range(cp)
    out["s2_integral_range_pore_px2"] = a2p * f * f
    phi = ctx.pore.mean()
    out["pore_frac_se_imagerep"] = float(np.sqrt(phi * (1 - phi) * a2p / ctx.area_px))
    cs = autocovariance(ctx.si)
    out["s2_len_si_px"] = 0.5 * (corr_length(cs[0, : w // 2]) + corr_length(cs[: h // 2, 0])) * f
    a2s = integral_range(cs)
    phi = ctx.si.mean()
    out["si_frac_se_imagerep"] = float(np.sqrt(phi * (1 - phi) * a2s / ctx.area_px))
    return out
