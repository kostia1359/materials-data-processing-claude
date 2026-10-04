"""Exact reference cases for the Laplace solver (brief 10.3)."""
from __future__ import annotations

import numpy as np

from .laplace import fv_laplace


def straight_channels(h=120, w=120, width=6, pitch=20):
    m = np.zeros((h, w))
    for c0 in range(0, w, pitch):
        m[:, c0 : c0 + width] = 1
    return m


def tilted_channel(theta_deg=30.0, h=200, w=600, width=16):
    """A periodic-free tilted band of conducting pixels; ends clipped by the frame."""
    t = np.tan(np.radians(theta_deg))
    r, c = np.indices((h, w))
    centre = w / 2 + (r - h / 2) * t
    halfw = width / 2 / np.cos(np.radians(theta_deg))
    return (np.abs(c - centre) <= halfw).astype(float)


def serpentine(h=300, w=200, width=8, legs=5, x0=20, x1=180):
    """Single channel of constant width: vertical legs alternating between x0 and x1, joined by
    horizontal runs. Returns (map, centre-line path length / height)."""
    m = np.zeros((h, w))
    seg = h // legs
    length = 0.0
    xs = [x0 if k % 2 == 0 else x1 for k in range(legs)]
    for k in range(legs):
        ya, yb = k * seg, (k + 1) * seg if k < legs - 1 else h
        m[ya:yb, xs[k] : xs[k] + width] = 1
        length += yb - ya
        if k < legs - 1:
            lo, hi = sorted((xs[k], xs[k + 1]))
            m[yb - width : yb, lo : hi + width] = 1
            length += hi - lo
    return m, length / h


def checkerboard(n=1024, cell=64, s1=1.0, s2=0.1):
    r, c = np.indices((n, n))
    return np.where(((r // cell) + (c // cell)) % 2 == 0, s1, s2)


def run_all() -> list[dict]:
    out = []
    m = straight_channels()
    r = fv_laplace(m)
    eps = m.mean()
    out.append(dict(case="straight channels τ = 1", expected=1.0, got=eps / r["D_eff_rel"], tol=1e-6))
    for th in (20.0, 35.0):
        m = tilted_channel(th)
        r = fv_laplace(m)
        # effective flux of a band of width W tilted by θ through height H: D_eff = W_perp/(w) * cos θ ... compare τ
        eps = m.mean()
        out.append(dict(case=f"tilted channel θ={th:.0f}° τ = 1/cos²θ", expected=1 / np.cos(np.radians(th)) ** 2,
                        got=eps / r["D_eff_rel"], tol=0.05))
    m, ratio = serpentine()
    r = fv_laplace(m)
    out.append(dict(case="serpentine τ = (L_path/L)²", expected=ratio**2, got=m.mean() / r["D_eff_rel"], tol=0.15))
    m = checkerboard()
    r = fv_laplace(m)
    out.append(dict(case="Keller checkerboard σ_eff = √(σ1σ2)", expected=np.sqrt(1.0 * 0.1), got=r["D_eff_rel"], tol=0.03))
    for o in out:
        o["rel_err"] = abs(o["got"] - o["expected"]) / o["expected"]
        o["ok"] = bool(o["rel_err"] <= o["tol"])
    return out
