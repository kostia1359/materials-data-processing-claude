"""Fused phase map with an explicit 'uncertain' class and exact bound maps (brief 10.2).

Evidence per pixel (Dempster-Shafer over the frame {P(ore), C(arbon), S(i)}):
  * BSE: class plausibilities from three Gaussians fitted to the smoothed-histogram modes. Fully trusted
    for Si-vs-carbon and for dark => pore; 'mid-grey => carbon' is discounted by beta (open shallow pores
    read as carbon), the discounted mass going to {P, C}.
  * ETD (two-point normalised on its pore-floor and graphite modes): dark supports {P}, bright supports
    the set {C, S}; 10 % of its mass stays on the whole frame.
  * InLens carries no class mass; its gradient magnitude only gates a 3x3 majority filter.
Decision: pignistic probability >= betp_min and conflict <= conflict_max (or top-two margin >= margin_min), else "uncertain".
Uncertainty is restricted to the pore/carbon decision (Si-vs-carbon is trusted from BSE), which keeps
the conducting sets nested, so D_eff(L_solid) <= D_eff(L_mid) <= D_eff(L_pore) holds exactly.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from .laplace import CARBON, PORE, SI, UNCERTAIN

_FOUR = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], bool)


def fuse_phase_map(bse_s, etd_s, inlens, t1, t2, cfg: dict, si_mask: np.ndarray | None = None, soft_levels: float = 3.0) -> dict:
    """si_mask: the core Si-candidate mask (BSE is fully trusted for Si vs carbon, brief 10.2). Pixels
    outside it are decided pore vs carbon by evidence fusion. BSE pore plausibility is a logistic around
    t1 with width `soft_levels` grey levels (instead of free Gaussians: an EM fit on the comb-LUT
    histogram put 23-29 % of pixels in the 'Si' component, see DECISIONS.md)."""
    fc = cfg["sim"]["fuse"]
    x = bse_s.astype(np.float32)
    if si_mask is None:
        si_mask = x >= t2
    pP = 1.0 / (1.0 + np.exp((x - t1) / soft_levels))
    pC = 1.0 - pP
    beta = fc["beta_bse_carbon"]
    bP, bC, bPC = pP, beta * pC, (1 - beta) * pC
    if etd_s is not None:
        e = etd_s.astype(np.float32)
        floor = float(np.percentile(e, 2))
        graph = float(np.median(e[(x >= t1) & (x < t2)]))
        en = (e - floor) / max(graph - floor, 1e-6)
        dark = 1.0 / (1.0 + np.exp((en - 0.5) / 0.1))
        eP, eCS, eT = 0.9 * dark, 0.9 * (1 - dark), np.full_like(dark, 0.1)
    else:
        eP, eCS, eT = np.zeros_like(x), np.zeros_like(x), np.ones_like(x)
    K = bP * eCS + bC * eP
    norm = np.maximum(1 - K, 1e-6)
    mP = (bP * (eP + eT) + bPC * eP) / norm
    mC = (bC * (eCS + eT) + bPC * eCS) / norm
    mPC = (bPC * eT) / norm
    betP_P, betP_C = mP + mPC / 2, mC + mPC / 2
    arg = np.where(betP_P >= betP_C, PORE, CARBON).astype(np.uint8)
    arg[si_mask] = SI
    conf = np.where(si_mask, 1.0, np.maximum(betP_P, betP_C)).astype(np.float32)
    margin = np.abs(betP_P - betP_C)
    # InLens-gated 3x3 majority smoothing (low gradient = smooth region = allowed to change); Si untouched
    if inlens is not None and fc.get("majority_filter", 0):
        gi = ndi.gaussian_gradient_magnitude(inlens.astype(np.float32), 1.0)
        gn = gi / max(np.percentile(gi, 95), 1e-6)
        counts = np.stack([ndi.uniform_filter((arg == k).astype(np.float32), fc["majority_filter"]) for k in range(3)])
        maj = np.argmax(counts, 0).astype(np.uint8)
        smooth = (gn * fc.get("inlens_gradient_weight", 1.0) < 0.5) & ~si_mask & (maj != SI)
        arg = np.where(smooth, maj, arg)
    arg = _remove_islands(arg, fc["min_island_px"], protect=si_mask)
    uncertain = (((conf < fc["betp_min"]) | (K > fc["conflict_max"])) & (margin < fc["margin_min"]) | (K > fc["conflict_max"])) & (arg != SI)
    L_mid = arg
    L_solid = np.where(uncertain, CARBON, L_mid).astype(np.uint8)
    L_pore = np.where(uncertain, PORE, L_mid).astype(np.uint8)
    L_vis = np.where(uncertain, UNCERTAIN, L_mid).astype(np.uint8)
    return dict(L_mid=L_mid, L_solid=L_solid, L_pore=L_pore, L_vis=L_vis, confidence=conf, conflict=K.astype(np.float32),
                uncertain_frac=float(uncertain.mean()), conflict_mean=float(K[~si_mask].mean()))


def _remove_islands(lab: np.ndarray, min_px: int, protect: np.ndarray | None = None) -> np.ndarray:
    out = lab.copy()
    for k in (PORE, CARBON):
        comp, n = ndi.label(lab == k, structure=_FOUR)
        if n == 0:
            continue
        sizes = np.bincount(comp.ravel())
        small = (sizes < min_px)
        small[0] = False
        m = small[comp]
        if protect is not None:
            m &= ~protect
        if not m.any():
            continue
        counts = np.stack([ndi.uniform_filter(((lab == j) & ~m).astype(np.float32), 5) for j in range(3)])
        rep = np.argmax(counts, 0)
        out[m] = np.where(rep[m] == SI, CARBON, rep[m])
    return out


def downsample_labels(labels: np.ndarray, factor: int = 4) -> tuple[np.ndarray, dict]:
    """Block-majority on labels; ties -> UNCERTAIN; no morphology. Audit Euler number, component count and
    percolating fraction of the pore phase before/after."""
    from skimage.measure import euler_number

    from .laplace import percolation_check

    h, w = (labels.shape[0] // factor) * factor, (labels.shape[1] // factor) * factor
    L = labels[:h, :w]
    blocks = L.reshape(h // factor, factor, w // factor, factor)
    counts = np.stack([(blocks == k).sum(axis=(1, 3)) for k in range(4)])
    best = counts.max(0)
    ties = (counts == best).sum(0) > 1
    ds = np.argmax(counts, 0).astype(np.uint8)
    ds[ties] = UNCERTAIN

    def stats(lab, f):
        p = lab == PORE
        a_mm2 = p.size * (f * 0.025e-3) ** 2
        return dict(euler_density_mm2=float(euler_number(p, connectivity=1) / a_mm2), components=int(ndi.label(p, structure=_FOUR)[1]),
                    percolating_frac_TP=percolation_check(p, "TP")["spanning_frac"],
                    percolating_frac_IP=percolation_check(p, "IP")["spanning_frac"], pore_frac=float(p.mean()))

    before, after = stats(L, 1), stats(ds, factor)
    audit = dict(before=before, after=after, tie_frac=float(ties.mean()),
                 d_percolating=max(abs(after["percolating_frac_TP"] - before["percolating_frac_TP"]),
                                   abs(after["percolating_frac_IP"] - before["percolating_frac_IP"])))
    return ds, audit
