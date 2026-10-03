from __future__ import annotations

import numpy as np


def fractions(ctx) -> dict:
    a = ctx.area_px
    n_si, n_c, n_p = int(ctx.si.sum()), int(ctx.carbon.sum()), int(ctx.pore.sum())
    return {
        "si_frac_solid": n_si / max(n_si + n_c, 1),
        "si_frac_total": n_si / a,
        "pore_frac_deep": n_p / a,
        "ambiguous_frac": float(ctx.ambiguous.sum()) / a,
    }
