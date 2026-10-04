"""Acquisition covariates ("gates"). They never enter a material KPI or the classifier; they get
their own baseline z-scores and drive the 'acquisition drift suspected' banner."""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

GATE_KEYS = [
    "n_grey_levels_bse", "comb_step_bse", "noise_sd_bse", "blur_bse", "sat0_bse", "sat255_bse",
    "sat0_etd", "sat255_etd", "sat0_inlens", "sat255_inlens", "bse_c_med", "etd_c_med", "inlens_c_med",
    "bse_si_med", "etd_si_med", "inlens_si_med", "bse_row_gradient_pct", "t1", "t2",
    "reg_shift_etd_px", "reg_shift_inlens_px", "height_px", "px_nm", "masked_rows", "etd_is_se",
]
# gates used by the shadow (acquisition-only) classifier, brief 6.3
SHADOW_GATES = ["inlens_c_med", "etd_c_med", "comb_step_bse", "n_grey_levels_bse", "noise_sd_bse", "height_px", "blur_bse"]


def noise_sd(bse: np.ndarray, carbon: np.ndarray) -> float:
    core = ndi.binary_erosion(carbon, iterations=5)
    if core.sum() < 100:
        core = carbon
    b = bse.astype(np.float32)
    resid = b - ndi.uniform_filter(b, 5)
    return float(resid[core].std())


def comb_step(img: np.ndarray) -> float:
    u = np.unique(img)
    u = u[(u >= 20) & (u <= 200)]
    return float(np.median(np.diff(u))) if len(u) > 2 else np.nan


def registration_shift(a: np.ndarray, b: np.ndarray, size: int = 1024) -> float:
    from skimage.registration import phase_cross_correlation

    h, w = a.shape
    s = min(size, h, w)
    r0, c0 = (h - s) // 2, (w - s) // 2
    pa = a[r0 : r0 + s, c0 : c0 + s].astype(np.float32)
    pb = b[r0 : r0 + s, c0 : c0 + s].astype(np.float32)
    # gradient magnitude makes the comparison robust to the contrast inversion between detectors
    ga = np.hypot(*np.gradient(ndi.gaussian_filter(pa, 1.5)))
    gb = np.hypot(*np.gradient(ndi.gaussian_filter(pb, 1.5)))
    shift, _, _ = phase_cross_correlation(ga, gb, upsample_factor=10)
    return float(np.hypot(*shift))


def compute_gates(sample, seg, cfg) -> dict:
    from skimage.measure import blur_effect

    bse, etd, inl = sample.bse, sample.etd, sample.inlens
    core_c = ndi.binary_erosion(seg.carbon, iterations=5)
    core_s = ndi.binary_erosion(seg.si, iterations=2)
    h = bse.shape[0]
    e = max(h // 8, 1)
    g = {
        "n_grey_levels_bse": sample.n_levels.get("BSE", np.nan),
        "comb_step_bse": comb_step(bse),
        "noise_sd_bse": noise_sd(bse, seg.carbon),
        "blur_bse": float(blur_effect(bse[:, : min(bse.shape[1], 2000)])),
        "bse_row_gradient_pct": float((bse[:e].mean() - bse[-e:].mean()) / max(bse.mean(), 1e-6) * 100),
        "t1": seg.t1, "t2": seg.t2,
        "height_px": int(sample.raw_shape[0]),
        "px_nm": sample.px_nm,
        "masked_rows": int(sample.crop["band_top"] + sample.crop["band_bottom"]),
        "shading_flag": bool(seg.shading_flag),
    }
    for name, img in (("bse", bse), ("etd", etd), ("inlens", inl)):
        if img is None:
            for k in ("sat0", "sat255"):
                g[f"{k}_{name}"] = np.nan
            g[f"{name}_c_med"] = g[f"{name}_si_med"] = np.nan
            continue
        g[f"sat0_{name}"] = float((img == 0).mean())
        g[f"sat255_{name}"] = float((img == 255).mean())
        g[f"{name}_c_med"] = float(np.median(img[core_c])) if core_c.any() else np.nan
        g[f"{name}_si_med"] = float(np.median(img[core_s])) if core_s.any() else np.nan
    g["reg_shift_etd_px"] = registration_shift(bse, etd) if etd is not None else np.nan
    g["reg_shift_inlens_px"] = registration_shift(bse, inl) if inl is not None else np.nan
    flags = []
    for k in ("reg_shift_etd_px", "reg_shift_inlens_px"):
        if np.isfinite(g[k]) and g[k] > cfg["registration_flag_px"]:
            flags.append(f"{k}={g[k]:.1f}")
    if seg.shading_flag:
        flags.append("shading_flag")
    g["gate_flags"] = ";".join(flags)
    return g
