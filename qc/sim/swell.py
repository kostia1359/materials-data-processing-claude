"""Geometric Si-swelling scenario (brief 10.3, swell.py).

Each Si-candidate component grows by successive 1-px rings until its area reaches f_A x A0. In
'pore_first' mode a ring consumes pore pixels before carbon pixels; in 'isotropic' mode all ring pixels
are taken in raster order. Carbon pixels taken are 'displaced' matrix. Particles are processed largest
first on a shared occupancy map, so growth that reaches another (original or grown) Si region counts as
a touch / merge. Particles whose growth reaches the frame are excluded (edge_hit).

Outputs are relative quantities: pore closure, constraint index CI = displaced carbon / growth area,
buffer sufficiency = share of particles with CI < 0.5, touching/merging statistics and a 2-D bound on
the electrode-thickness change dH/H (displaced area / frame area, all pushed out of plane).
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from .laplace import CARBON, PORE, SI

_FOUR = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], bool)


def swell(labels: np.ndarray, f_A: float, mode: str = "pore_first", px_um: float = 0.1, hotspot_um: float = 5.0,
          return_maps: bool = False) -> dict:
    lab = labels.copy()
    si0 = lab == SI
    comp, n = ndi.label(si0, structure=_FOUR)
    if n == 0:
        return dict(n_particles=0)
    sizes = np.bincount(comp.ravel())
    order = np.argsort(-sizes[1:]) + 1
    owner = comp.astype(np.int32)  # which particle owns each Si pixel (original or grown)
    slices = ndi.find_objects(comp)
    h, w = lab.shape
    pore_total = int((lab == PORE).sum())
    pore_consumed = 0
    rows = []
    for pid in order:
        a0 = int(sizes[pid])
        if a0 < 4:
            continue
        target = int(round(f_A * a0)) - a0
        sl = slices[pid - 1]
        grow = int(np.ceil(np.sqrt(f_A) * np.sqrt(a0 / np.pi))) + 3
        r0, r1 = max(sl[0].start - grow, 0), min(sl[0].stop + grow, h)
        c0, c1 = max(sl[1].start - grow, 0), min(sl[1].stop + grow, w)
        sub_owner = owner[r0:r1, c0:c1]
        sub_lab = lab[r0:r1, c0:c1]
        region = sub_owner == pid
        taken_pore = taken_carbon = 0
        touched = set()
        edge_hit = False
        while taken_pore + taken_carbon < target:
            ring = ndi.binary_dilation(region, structure=_FOUR) & ~region
            if not ring.any():
                break
            others = ring & (sub_lab == SI) & (sub_owner != pid)
            if others.any():
                touched |= set(np.unique(sub_owner[others]).tolist())
            free = ring & (sub_lab != SI)
            if ((r0 == 0) and free[0].any()) or ((r1 == h) and free[-1].any()) or ((c0 == 0) and free[:, 0].any()) or ((c1 == w) and free[:, -1].any()):
                edge_hit = True
            need = target - taken_pore - taken_carbon
            if mode == "pore_first":
                cand = [np.argwhere(free & (sub_lab == PORE)), np.argwhere(free & (sub_lab == CARBON)),
                        np.argwhere(free & (sub_lab != PORE) & (sub_lab != CARBON))]
                pts = np.concatenate([c for c in cand if len(c)]) if any(len(c) for c in cand) else np.zeros((0, 2), int)
            else:
                pts = np.argwhere(free)
            if len(pts) == 0:
                region |= ring  # grows through other Si: count as merged contact, no new area
                continue
            pts = pts[:need]
            vals = sub_lab[pts[:, 0], pts[:, 1]]
            taken_pore += int((vals == PORE).sum())
            taken_carbon += int((vals != PORE).sum())
            sub_lab[pts[:, 0], pts[:, 1]] = SI
            sub_owner[pts[:, 0], pts[:, 1]] = pid
            region[pts[:, 0], pts[:, 1]] = True
        pore_consumed += taken_pore
        g = max(taken_pore + taken_carbon, 1)
        rows.append(dict(pid=int(pid), a0=a0, ecd_um=2 * np.sqrt(a0 / np.pi) * px_um, ci=taken_carbon / g, displaced=taken_carbon,
                         pore_taken=taken_pore, touched=len(touched), edge_hit=edge_hit))
    keep = [r for r in rows if not r["edge_hit"]]
    if not keep:
        return dict(n_particles=0, n_edge_hit=len(rows))
    ci = np.array([r["ci"] for r in keep])
    a = np.array([r["a0"] for r in keep], float)
    ecd = np.array([r["ecd_um"] for r in keep])
    disp = sum(r["displaced"] for r in keep)
    # merged clusters after swelling
    final_si = lab == SI
    clab, cn = ndi.label(final_si, structure=_FOUR)
    csz = np.bincount(clab.ravel())[1:] if cn else np.array([])
    cecd = 2 * np.sqrt(csz / np.pi) * px_um if len(csz) else np.array([0.0])
    by_size = {}
    for name, lo, hi in (("small_lt2um", 0, 2), ("mid_2_5um", 2, 5), ("large_gt5um", 5, np.inf)):
        m = (ecd >= lo) & (ecd < hi)
        by_size[name] = float(np.average(ci[m], weights=a[m])) if m.any() else np.nan
    maps = {}
    if return_maps:
        ci_map = np.full(labels.shape, np.nan, np.float32)
        ci_of = {r["pid"]: r["ci"] for r in keep}
        for pid, v in ci_of.items():
            ci_map[comp == pid] = v
        maps = dict(swollen=lab, constraint_map=ci_map)
    return dict(
        maps=maps, f_A=f_A, mode=mode, n_particles=len(keep), n_edge_hit=len(rows) - len(keep),
        pore_closure_frac=pore_consumed / max(pore_total, 1),
        constraint_index=float(ci.mean()), constraint_index_aw=float(np.average(ci, weights=a)), constraint_index_by_size=by_size,
        buffer_sufficiency_frac=float((ci < 0.5).mean()),
        si_touch_frac=float(np.mean([r["touched"] > 0 for r in keep])),
        si_merges=int(sum(r["touched"] for r in keep)),
        largest_merged_cluster_ecd_um=float(cecd.max()),
        si_area_in_clusters_gt_5um=float(csz[cecd > hotspot_um].sum() / max(csz.sum(), 1)) if len(csz) else 0.0,
        dH_H_bound=float(disp / labels.size),
    )
