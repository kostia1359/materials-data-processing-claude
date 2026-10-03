from __future__ import annotations

import numpy as np

from .pores import chords


def structure_tensor_maps(bse_s: np.ndarray, sigma: float):
    from skimage.feature import structure_tensor

    return tuple(a.astype(np.float32) for a in structure_tensor(bse_s, sigma=sigma, order="rc"))


def anisotropy(ctx, pore_lists: dict) -> dict:
    ch, cv = pore_lists["pore_chord_h_um"], pore_lists["pore_chord_v_um"]
    out = {"aniso_pore_chord_ratio": float(ch.mean() / cv.mean()) if len(ch) and len(cv) else np.nan}
    kh, kv = chords(ctx.carbon, 1), chords(ctx.carbon, 0)
    out["aniso_carbon_chord_ratio"] = float(kh.mean() / kv.mean()) if len(kh) and len(kv) else np.nan
    if ctx.st is not None:
        arr, arc, acc = (float(a.mean()) for a in ctx.st)
        jxx, jyy, jxy = acc, arr, -arc  # x = column, y = up
        tr = jxx + jyy
        coh = np.sqrt((jxx - jyy) ** 2 + 4 * jxy**2) / tr if tr > 0 else np.nan
        grad_angle = 0.5 * np.degrees(np.arctan2(2 * jxy, jxx - jyy))
        orient = (grad_angle + 90 + 90) % 180 - 90  # structure direction, in [-90, 90)
        out["st_coherence"] = float(coh)
        out["st_orientation_deg"] = float(orient)
    else:
        out["st_coherence"] = out["st_orientation_deg"] = np.nan
    return out
