"""Whole-image and per-strip KPI evaluation, strip SE and block bootstrap."""
from __future__ import annotations

import numpy as np

from .kpis import anisotropy, connectivity, contacts, correlation, defects, dispersion, fractions, pores, sizes, vertical


def compute_kpis(ctx, with_lists: bool = True) -> tuple[dict, dict]:
    out: dict = {}
    lists: dict = {}
    out.update(fractions.fractions(ctx))
    k, l_si = sizes.si_particles(ctx)
    out.update(k)
    out.update(dispersion.dispersion(ctx, l_si.get("si_centroids_px")))
    k, l_p = pores.pore_chords(ctx)
    out.update(k)
    out.update(anisotropy.anisotropy(ctx, l_p))
    out.update(correlation.two_point(ctx))
    out.update(connectivity.connectivity(ctx))
    out.update(contacts.contacts(ctx))
    out.update(defects.defects(ctx))
    out.update(vertical.vertical(ctx))
    if with_lists:
        lists.update({"si_ecd_um": l_si["si_ecd_um"], **l_p})
    return out, lists


def strip_slices(width: int, n: int) -> list[slice]:
    edges = np.linspace(0, width, n + 1).astype(int)
    return [slice(int(a), int(b)) for a, b in zip(edges[:-1], edges[1:])]


def strip_kpis(ctx, n: int) -> list[dict]:
    rows = []
    for i, sl in enumerate(strip_slices(ctx.shape[1], n)):
        k, _ = compute_kpis(ctx.sub(sl), with_lists=False)
        k["strip"] = i
        rows.append(k)
    return rows


def strip_se(strips: list[dict], keys) -> dict:
    out = {}
    n = len(strips)
    for k in keys:
        v = np.array([s.get(k, np.nan) for s in strips], float)
        v = v[np.isfinite(v)]
        out[k] = float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) >= 2 else np.nan
    return out


def bootstrap_strip_means(strips: list[dict], keys, n_boot: int, seed: int) -> np.ndarray:
    """(n_boot, len(keys)) array of strip-resampled means (block bootstrap over strips)."""
    rng = np.random.default_rng(seed)
    m = np.array([[s.get(k, np.nan) for k in keys] for s in strips], float)
    idx = rng.integers(0, len(strips), size=(n_boot, len(strips)))
    return np.nanmean(m[idx], axis=1)


def bootstrap_ci(strips: list[dict], keys, n_boot: int, seed: int) -> dict:
    b = bootstrap_strip_means(strips, keys, n_boot, seed)
    with np.errstate(all="ignore"):
        lo = np.nanpercentile(b, 2.5, axis=0)
        hi = np.nanpercentile(b, 97.5, axis=0)
    return {k: [float(lo[i]), float(hi[i])] for i, k in enumerate(keys)}
