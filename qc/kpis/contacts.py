from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from skimage.measure import perimeter_crofton


def contacts(ctx) -> dict:
    per_px = perimeter_crofton(ctx.pore, directions=4)
    out = {"interface_pore_solid_um_per_um2": float(per_px * ctx.px_um / ctx.area_um2)}
    # The smoothed BSE always passes through the carbon grey band between a bright Si-candidate and a
    # dark pore, so the immediate 1-px neighbour is nearly always "carbon". The contact ring is
    # therefore taken 2-4 px out from the particle (beyond the transition rim); see DECISIONS.md.
    inner = ndi.binary_dilation(ctx.si, iterations=1)
    outer = ndi.binary_dilation(ctx.si, iterations=4)
    ring = outer & ~inner & ~ctx.si
    n = ring.sum()
    out["si_contact_pore_frac"] = float((ring & ctx.pore).sum() / n) if n else np.nan
    out["si_contact_carbon_frac"] = float((ring & ctx.carbon).sum() / n) if n else np.nan
    return out
