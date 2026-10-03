"""Synthetic BSE/ETD/InLens triples with known ground truth (brief 6.1).

Truth values are measured on the *drawn* masks (after overlaps), using the same definitions as
the KPI code, so recovery tests compare like with like.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi
from skimage.draw import disk as draw_disk
from skimage.draw import ellipse as draw_ellipse

PX_UM = 0.025


@dataclass
class SynthParams:
    width: int = 7000
    height: int = 2000
    pore_frac: float = 0.15
    pore_aspect: float = 2.0  # horizontal : vertical semi-axis
    pore_minor_um_median: float = 0.35
    si_frac: float = 0.07
    si_d50_um: float = 4.0  # number-median section ECD (log-normal)
    si_sigma_log: float = 0.35
    clustered: bool = False
    poisson_count: bool = False  # draw the particle count from a Poisson law instead of filling to si_frac
    carbon_grey: float = 58
    pore_grey: float = 12
    si_grey: float = 112
    noise_sd: float = 8
    comb_keep_every: int = 3
    gamma: float = 0.9
    seed: int = 0


def _place_pores(rng, mask, p: SynthParams):
    h, w = mask.shape
    target = p.pore_frac * h * w
    while mask.sum() < target:
        b = rng.lognormal(np.log(p.pore_minor_um_median / PX_UM), 0.5)
        b = max(b, 2.0)
        a = b * p.pore_aspect
        r, c = rng.uniform(0, h), rng.uniform(0, w)
        rr, cc = draw_ellipse(r, c, b, a, shape=mask.shape)
        mask[rr, cc] = True


def _place_si(rng, si, pore, p: SynthParams):
    h, w = si.shape
    target = p.si_frac * h * w
    occupied = np.zeros_like(si)
    centres = None
    if p.clustered:
        n_c = max(2, int(h * w * PX_UM**2 / 600))
        centres = np.column_stack([rng.uniform(0, h, n_c), rng.uniform(0, w, n_c)])
    tries = 0
    radii = []
    n_target = None
    if p.poisson_count:
        mean_area = np.pi * (p.si_d50_um / PX_UM / 2) ** 2 * np.exp(2 * p.si_sigma_log**2)
        n_target = rng.poisson(target / mean_area)

    def more():
        return len(radii) < n_target if n_target is not None else si.sum() < target

    while more() and tries < 200000:
        tries += 1
        ecd_um = rng.lognormal(np.log(p.si_d50_um), p.si_sigma_log)
        rad = ecd_um / PX_UM / 2
        if centres is not None:
            k = rng.integers(len(centres))
            r, c = rng.normal(centres[k], 2.0 / PX_UM)
        else:
            r, c = rng.uniform(rad, h - rad), rng.uniform(rad, w - rad)
        if not (rad <= r < h - rad and rad <= c < w - rad):
            continue
        rr, cc = draw_disk((r, c), rad + 3, shape=si.shape)
        if occupied[rr, cc].any():
            continue
        occupied[rr, cc] = True
        rr, cc = draw_disk((r, c), rad, shape=si.shape)
        si[rr, cc] = True
        radii.append(rad)
    pore &= ~occupied  # Si particles sit in solid; keep a carbon rim between Si and pores
    return radii


def _comb_gamma(img: np.ndarray, keep: int, gamma: float) -> np.ndarray:
    x = np.clip(img, 0, 255) / 255.0
    x = x**gamma * 255
    x = np.round(x / keep) * keep
    return np.clip(x, 0, 255).astype(np.uint8)


def make_sample(p: SynthParams):
    rng = np.random.default_rng(p.seed)
    h, w = p.height, p.width
    pore = np.zeros((h, w), bool)
    _place_pores(rng, pore, p)
    si = np.zeros((h, w), bool)
    _place_si(rng, si, pore, p)
    pore &= ~si
    base = np.full((h, w), p.carbon_grey, np.float32)
    base[pore] = p.pore_grey
    base[si] = p.si_grey
    base = ndi.gaussian_filter(base, 0.8)
    bse = base + rng.normal(0, p.noise_sd, (h, w))
    etd_base = np.where(pore, 5.0, np.where(si, 120.0, 70.0)).astype(np.float32)
    etd = ndi.gaussian_filter(etd_base, 0.8) + rng.normal(0, 9, (h, w))
    inl_base = np.where(pore, 40.0, np.where(si, 160.0, 90.0)).astype(np.float32)
    inl = ndi.gaussian_filter(inl_base, 0.8) + rng.normal(0, 12, (h, w))
    imgs = {
        "BSE": _comb_gamma(bse, p.comb_keep_every, p.gamma),
        "ETD": _comb_gamma(etd, 2, 1.0),
        "Inlens": _comb_gamma(inl, 2, 1.0),
    }
    return imgs, {"pore": pore, "si": si}


def truth_values(masks: dict) -> dict:
    from .kpis.pores import chords
    from .kpis.sizes import area_weighted_quantiles

    pore, si = masks["pore"], masks["si"]
    lab, n = ndi.label(si)
    areas = np.bincount(lab.ravel())[1:]
    areas = areas[areas >= 500]
    ecd = 2 * np.sqrt(areas / np.pi) * PX_UM
    ch, cv = chords(pore, 1), chords(pore, 0)
    return {
        "si_frac_total": float(si.mean()),
        "pore_frac_deep": float(pore.mean()),
        "si_d50_aw_um": area_weighted_quantiles(ecd)[1],
        "aniso_pore_chord_ratio": float(ch.mean() / cv.mean()),
    }


def perturb(img: np.ndarray, brightness=1.0, contrast=1.0, gamma=1.0, noise=0.0, seed=0) -> np.ndarray:
    x = img.astype(np.float32)
    m = x.mean()
    x = (x - m) * contrast + m
    x = x * brightness
    x = 255 * (np.clip(x, 0, 255) / 255) ** gamma
    if noise:
        x = x + np.random.default_rng(seed).normal(0, noise, x.shape)
    return np.clip(np.round(x), 0, 255).astype(np.uint8)


def write_triple(folder: Path, sample_id: str, imgs: dict, px_nm: float = 25.0):
    folder.mkdir(parents=True, exist_ok=True)
    res = 0.0254 / (px_nm * 1e-9)  # px per inch
    for det, a in imgs.items():
        rgb = np.repeat(a[..., None], 3, axis=2)
        tifffile.imwrite(folder / f"img_{sample_id}_{det}.tif", rgb, photometric="rgb", compression="lzw",
                         resolution=(res, res), resolutionunit="INCH")


# Variations used for the synthetic multi-batch dataset (pipeline tests / demos only, never real data).
BATCH_VARIANTS = {
    "3": dict(),  # baseline
    "1": dict(si_frac=0.10),  # more Si loading
    "2": dict(pore_aspect=3.0, si_d50_um=5.5, pore_frac=0.12),  # coarser Si, flatter/fewer pores
}


def make_dataset(out: Path, per_batch: int = 4, width: int = 3500, height: int = 1000, seed: int = 100):
    out = Path(out)
    k = 0
    for b, var in BATCH_VARIANTS.items():
        for i in range(per_batch):
            p = SynthParams(width=width, height=height, seed=seed + k, **var)
            imgs, _ = make_sample(p)
            write_triple(out / f"Batch_{b}", f"syn{b}x{i:02d}", imgs)
            k += 1
