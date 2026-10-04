from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi


def fractions(ctx) -> dict:
    a = ctx.area_px
    n_si, n_c, n_p = int(ctx.si.sum()), int(ctx.carbon.sum()), int(ctx.pore.sum())
    out = {
        "si_frac_solid": n_si / max(n_si + n_c, 1),
        "si_frac_total": n_si / a,
        "pore_frac_deep": n_p / a,
        "ambiguous_frac": float(ctx.ambiguous.sum()) / a,
    }
    out["si_bse_contrast"] = bse_contrast(ctx)
    return out


def bse_contrast(ctx) -> float:
    """(Si-candidate - carbon) / (carbon - pore floor) on smoothed BSE, each phase's median inside its
    eroded mask. A Z-contrast ratio: invariant to linear brightness/contrast changes, so a shift means
    the bright phase has a different mean Z (e.g. SiOx or Si-C composite vs Si). Post-hoc KPI (Batch_1)."""
    def med(mask, it):
        core = ndi.binary_erosion(mask, iterations=it)
        if core.sum() < 50:
            core = mask
        return float(np.median(ctx.bse_s[core])) if core.any() else np.nan

    s, c, p = med(ctx.si, 2), med(ctx.carbon, 3), med(ctx.pore, 2)
    return (s - c) / (c - p) if np.isfinite(s) and np.isfinite(c) and np.isfinite(p) and c > p else np.nan
