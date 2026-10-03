"""Si-candidate particle sizes and shapes on watershed-split components."""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from skimage.measure import regionprops_table


def area_weighted_quantiles(ecd: np.ndarray, qs=(0.1, 0.5, 0.9)) -> list[float]:
    if len(ecd) == 0:
        return [np.nan] * len(qs)
    order = np.argsort(ecd)
    e = ecd[order]
    w = e**2  # area ∝ ECD²
    cw = np.cumsum(w) / w.sum()
    return [float(np.interp(q, cw, e)) for q in qs]


def si_particles(ctx) -> tuple[dict, dict]:
    cfg = ctx.cfg
    lab = ctx.si_labels
    counts = np.bincount(lab.ravel())
    counts[0] = 0
    n_count = int((counts >= ctx.px_area(cfg["si_min_px_count"])).sum())
    out = {"si_num_density_mm2": n_count / ctx.area_mm2}
    big = np.where(counts >= ctx.px_area(cfg["si_min_px_size"]))[0]
    lists = {"si_ecd_um": np.array([])}
    nan_keys = ["si_d10_aw_um", "si_d50_aw_um", "si_d90_aw_um", "si_span", "si_solidity_med", "si_aspect_med",
                "si_circularity_med", "si_interior_texture", "si_n_sized"]
    if len(big) == 0:
        out.update({k: np.nan for k in nan_keys})
        out["si_n_sized"] = 0
        return out, lists
    keep = np.zeros(len(counts), bool)
    keep[big] = True
    lab_big = np.where(keep[lab], lab, 0)
    props = regionprops_table(lab_big, properties=("label", "area", "solidity", "axis_major_length",
                                                   "axis_minor_length", "perimeter_crofton", "centroid"))
    area = props["area"].astype(float)
    ecd = 2 * np.sqrt(area / np.pi) * ctx.px_um
    d10, d50, d90 = area_weighted_quantiles(ecd)
    minor = np.maximum(props["axis_minor_length"], 1e-6)
    per = np.maximum(props["perimeter_crofton"], 1e-6)
    # interior texture: SD of smoothed BSE in the 1-px-eroded particle, relative to image noise
    interior = ndi.binary_erosion(lab_big > 0)
    ilab = np.where(interior, lab_big, 0)
    sd = np.asarray(ndi.standard_deviation(ctx.bse_s, labels=ilab, index=props["label"]), float)
    tex = np.nanmedian(sd) / max(ctx.noise_sd, 1e-6) if np.isfinite(sd).any() else np.nan
    out.update({
        "si_d10_aw_um": d10, "si_d50_aw_um": d50, "si_d90_aw_um": d90,
        "si_span": (d90 - d10) / d50 if d50 > 0 else np.nan,
        "si_solidity_med": float(np.median(props["solidity"])),
        "si_aspect_med": float(np.median(props["axis_major_length"] / minor)),
        "si_circularity_med": float(np.median(np.clip(4 * np.pi * area / per**2, 0, 1.5))),
        "si_interior_texture": float(tex),
        "si_n_sized": int(len(big)),
    })
    lists["si_ecd_um"] = ecd
    lists["si_centroids_px"] = np.stack([props["centroid-0"], props["centroid-1"]], 1)
    return out, lists
