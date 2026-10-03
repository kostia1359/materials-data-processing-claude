"""Exactly three classes from BSE: deep pore / carbon matrix / Si-candidate.

ETD is used only to *confirm* pore floors (dark in ETD); it never moves a threshold.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage as ndi
from skimage.filters import apply_hysteresis_threshold, gaussian, threshold_multiotsu
from skimage.morphology import disk
from skimage.segmentation import watershed


def remove_small(mask: np.ndarray, min_px: int) -> np.ndarray:
    lab, n = ndi.label(mask)
    if n == 0:
        return mask
    sizes = np.bincount(lab.ravel())
    keep = sizes >= min_px
    keep[0] = False
    return keep[lab]


@dataclass
class Segmentation:
    bse_s: np.ndarray  # float32 smoothed BSE
    etd_s: np.ndarray | None
    t1: float
    t2: float
    pore: np.ndarray  # bool, final deep pore
    ambiguous: np.ndarray  # BSE-dark but ETD-bright
    si: np.ndarray  # bool, final Si-candidate
    carbon: np.ndarray  # bool
    si_labels: np.ndarray  # int32, watershed-split Si particles (all sizes >= 1 px of si)
    band_thresholds: list
    shading_flag: bool
    etd_pore_cut: float | None


def smooth(img: np.ndarray, sigma: float) -> np.ndarray:
    return gaussian(img.astype(np.float32), sigma=sigma, preserve_range=True).astype(np.float32)


def pore_mask(bse_s, etd_s, t1, carbon_for_etd, cfg, px_nm: float = 25.0):
    lv = cfg["pore_hysteresis_levels"]
    grown = apply_hysteresis_threshold(-bse_s, low=-(t1 + lv), high=-t1)
    etd_cut = None
    if etd_s is not None and carbon_for_etd.any():
        etd_cut = float(np.percentile(etd_s[carbon_for_etd], cfg["pore_etd_percentile"]))
        dark = etd_s < etd_cut
        pore = grown & dark
        ambiguous = grown & ~dark
    else:
        pore = grown
        ambiguous = np.zeros_like(grown)
    pore = remove_small(pore, max(1, int(round(cfg["pore_min_px"] * (25.0 / px_nm) ** 2))))
    return pore, ambiguous, etd_cut


def si_mask(bse_s, t2, cfg):
    si = bse_s >= t2
    r = cfg["si_opening_radius_px"]
    if r > 0:
        si = ndi.binary_opening(si, structure=disk(r).astype(bool))
    lab, n = ndi.label(si)
    if n == 0:
        return si, lab
    # interior check: median of the 1-px-eroded interior must clear t2 + margin
    interior = ndi.binary_erosion(si)
    ilab = np.where(interior, lab, 0)
    idx = np.arange(1, n + 1)
    med = ndi.median(bse_s, labels=ilab, index=idx)
    med = np.nan_to_num(np.asarray(med, dtype=float), nan=-1.0)
    keep = np.zeros(n + 1, bool)
    keep[1:] = med >= t2 + cfg["si_interior_margin_levels"]
    si = keep[lab]
    return si, lab * si


def split_particles(si: np.ndarray, h: float) -> np.ndarray:
    """Watershed on the EDT, markers from h-maxima (brief KPI 2). Done per component bbox."""
    from skimage.morphology import h_maxima

    lab, n = ndi.label(si)
    out = np.zeros(si.shape, np.int32)
    nxt = 1
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        if sl is None:
            continue
        sub = lab[sl] == i
        area = int(sub.sum())
        if area < 200:
            out[sl][sub] = nxt
            nxt += 1
            continue
        pad = np.pad(sub, 1)
        edt = ndi.distance_transform_edt(pad)
        mk, nm = ndi.label(h_maxima(edt, h))
        if nm <= 1:
            out[sl][sub] = nxt
            nxt += 1
            continue
        ws = watershed(-edt, mk, mask=pad)[1:-1, 1:-1]
        region = out[sl]
        region[sub] = ws[sub] + nxt - 1
        nxt += nm
    return out


def segment(bse: np.ndarray, etd: np.ndarray | None, cfg: dict, t_shift=(0.0, 0.0), split=True, px_nm: float = 25.0) -> Segmentation:
    sig = cfg["smoothing_sigma_px"]
    bse_s = smooth(bse, sig)
    etd_s = smooth(etd, sig) if etd is not None else None
    t1, t2 = (float(t) for t in threshold_multiotsu(bse_s, classes=cfg["otsu_classes"]))
    t1 += t_shift[0]
    t2 += t_shift[1]
    carbon0 = (bse_s >= t1) & (bse_s < t2)
    pore, amb, etd_cut = pore_mask(bse_s, etd_s, t1, carbon0, cfg, px_nm)
    si, _ = si_mask(bse_s, t2, cfg)
    si &= ~pore
    carbon = ~pore & ~si
    si_labels = split_particles(si, cfg["si_watershed_h"]) if split else ndi.label(si)[0].astype(np.int32)
    # shading diagnostic: Otsu per horizontal band
    bands = []
    h = bse_s.shape[0]
    for k in range(5):
        sub = bse_s[k * h // 5 : (k + 1) * h // 5]
        try:
            bands.append([float(t) for t in threshold_multiotsu(sub[::2, ::2], classes=3)])
        except ValueError:
            bands.append([np.nan, np.nan])
    b = np.array(bands)
    drift = np.nanmax(b, axis=0) - np.nanmin(b, axis=0)
    shading = bool(np.nanmax(drift) > cfg["shading_flag_levels"])
    return Segmentation(bse_s, etd_s, t1, t2, pore, amb, si, carbon, si_labels, bands, shading, etd_cut)


def overlay_rgb(bse: np.ndarray, seg: Segmentation, factor: int = 4) -> np.ndarray:
    g = bse[::factor, ::factor].astype(np.float32) / 255.0
    rgb = np.stack([g, g, g], -1)
    p = seg.pore[::factor, ::factor]
    s = seg.si[::factor, ::factor]
    rgb[p] = rgb[p] * 0.35 + np.array([0.15, 0.35, 0.95]) * 0.65
    rgb[s] = rgb[s] * 0.35 + np.array([1.0, 0.55, 0.05]) * 0.65
    return (np.clip(rgb, 0, 1) * 255).astype(np.uint8)
