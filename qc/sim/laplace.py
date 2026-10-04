"""Steady-state diffusion / conduction on a pixel map (TauFactor conventions, scipy solver).

5-point finite-volume Laplacian on 4-connected pixels, harmonic-mean face conductances, Dirichlet
planes half a pixel outside the first/last row (ΔC = 1), no-flux sides. D_eff/D = Q·L/(ΔC·A) with unit
reference diffusivity, so for a single conducting phase τ = ε·D/D_eff. Only pixels connected to *both*
electrodes enter the system (isolated clusters carry no steady flux and would make it singular);
when nothing spans, the solve is refused with a reason instead of returning a meaningless number.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy import ndimage as ndi

PORE, CARBON, SI, UNCERTAIN = 0, 1, 2, 3
_FOUR = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], bool)


def _orient(a: np.ndarray, axis: str) -> np.ndarray:
    """Flow is always solved top -> bottom; 'IP' (in-plane) transposes the map."""
    return a if axis == "TP" else a.T


def percolation_check(mask: np.ndarray, axis: str = "TP") -> dict:
    m = _orient(mask, axis)
    lab, n = ndi.label(m, structure=_FOUR)
    if n == 0:
        return dict(spans=False, n_clusters=0, spanning_frac=0.0, spanning_labels=[])
    top = np.unique(lab[0][lab[0] > 0])
    bot = np.unique(lab[-1][lab[-1] > 0])
    span = np.intersect1d(top, bot)
    tot = m.sum()
    sizes = np.bincount(lab.ravel())
    return dict(spans=bool(len(span)), n_clusters=int(n), spanning_frac=float(sizes[span].sum() / tot) if tot else 0.0,
                spanning_labels=span.tolist())


def fv_laplace(cond: np.ndarray, axis: str = "TP", solver: str = "spsolve", tol: float = 1e-8) -> dict:
    """Solve div(c grad C) = 0 with C = 1 above the top row and 0 below the bottom row."""
    c = np.asarray(_orient(cond, axis), float)
    h, w = c.shape
    conducting = c > 0
    perc = percolation_check(conducting, "TP")
    if not perc["spans"]:
        return dict(defined=False, reason=f"no conducting cluster spans the frame ({perc['n_clusters']} clusters)",
                    D_eff_rel=np.nan, n_clusters=perc["n_clusters"], spanning_frac=0.0)
    lab, _ = ndi.label(conducting, structure=_FOUR)
    active = np.isin(lab, perc["spanning_labels"])
    idx = -np.ones((h, w), np.int64)
    n = int(active.sum())
    idx[active] = np.arange(n)
    rows, cols, vals = [], [], []
    diag = np.zeros(n)
    b = np.zeros(n)
    # interior faces (vertical neighbours then horizontal neighbours)
    for (sa, sb) in (((slice(0, h - 1), slice(None)), (slice(1, h), slice(None))),
                     ((slice(None), slice(0, w - 1)), (slice(None), slice(1, w)))):
        a_act = active[sa] & active[sb]
        ia, ib = idx[sa][a_act], idx[sb][a_act]
        ca, cb = c[sa][a_act], c[sb][a_act]
        g = 2 * ca * cb / (ca + cb)
        rows += [ia, ib]
        cols += [ib, ia]
        vals += [-g, -g]
        np.add.at(diag, ia, g)
        np.add.at(diag, ib, g)
    # Dirichlet planes half a pixel outside: face conductance = c / 0.5
    top = active[0]
    it = idx[0][top]
    g_top = 2 * c[0][top]
    np.add.at(diag, it, g_top)
    np.add.at(b, it, g_top * 1.0)
    bot = active[-1]
    ib_ = idx[-1][bot]
    g_bot = 2 * c[-1][bot]
    np.add.at(diag, ib_, g_bot)
    A = sp.coo_matrix((np.concatenate(vals + [diag]), (np.concatenate(rows + [np.arange(n)]), np.concatenate(cols + [np.arange(n)]))),
                      shape=(n, n)).tocsr()
    if solver == "amg_cg":
        import pyamg
        from scipy.sparse.linalg import cg

        ml = pyamg.smoothed_aggregation_solver(A, symmetry="symmetric")
        x, info = cg(A, b, rtol=tol, maxiter=2000, M=ml.aspreconditioner())
        if info != 0:
            return dict(defined=False, reason=f"AMG-CG did not converge (info={info})", D_eff_rel=np.nan)
    else:
        from scipy.sparse.linalg import spsolve

        x = spsolve(A.tocsc(), b)
    q_in = float((g_top * (1.0 - x[it])).sum())
    q_out = float((g_bot * x[ib_]).sum())
    q = 0.5 * (q_in + q_out)
    d_eff = q * h / w  # L = h rows, A = w columns (unit depth), ΔC = 1
    field = np.full((h, w), np.nan)
    field[active] = x
    return dict(defined=True, D_eff_rel=float(d_eff), Q_in=q_in, Q_out=q_out,
                flux_balance=float(abs(q_in - q_out) / max(abs(q_in), 1e-300)), n_unknowns=n,
                n_clusters=perc["n_clusters"], spanning_frac=perc["spanning_frac"], field=_orient(field, axis))


def cond_map(labels: np.ndarray, values: dict) -> np.ndarray:
    out = np.zeros(labels.shape, float)
    for k, v in values.items():
        out[labels == k] = v
    return out


def ionic_index(labels: np.ndarray, D_c: float = 0.05, solver: str = "spsolve") -> dict:
    """Two-conductivity electrolyte transport: pore 1, carbon matrix D_c (stands in for unresolved porosity), Si 0."""
    c = cond_map(labels, {PORE: 1.0, CARBON: D_c, SI: 0.0})
    eps = float((labels == PORE).mean())
    tp = fv_laplace(c, "TP", solver)
    ip = fv_laplace(c, "IP", solver)
    out = dict(D_c=D_c, eps_pore=eps,
               D_eff_rel_TP=tp["D_eff_rel"], D_eff_rel_IP=ip["D_eff_rel"],
               flux_balance_TP=tp.get("flux_balance", np.nan), flux_balance_IP=ip.get("flux_balance", np.nan))
    out["aniso_ratio"] = out["D_eff_rel_IP"] / out["D_eff_rel_TP"] if tp["defined"] and ip["defined"] else np.nan
    out["tau_p"] = eps / out["D_eff_rel_TP"] if tp["defined"] and out["D_eff_rel_TP"] > 0 else np.nan
    out["N_M"] = 1.0 / out["D_eff_rel_TP"] if tp["defined"] and out["D_eff_rel_TP"] > 0 else np.nan
    pore_only = percolation_check(labels == PORE, "TP")
    out["pore_spans_TP"] = pore_only["spans"]
    out["pore_clusters"] = pore_only["n_clusters"]
    out["undefined"] = [f"{ax}: {r['reason']}" for ax, r in (("TP", tp), ("IP", ip)) if not r["defined"]]
    return out


def electronic_index(labels: np.ndarray, sigma_si: float = 0.0, solver: str = "spsolve") -> dict:
    """Electronic network: carbon matrix 1, Si sigma_si (0 = worst case), pore 0."""
    c = cond_map(labels, {PORE: 0.0, CARBON: 1.0, SI: sigma_si})
    tp = fv_laplace(c, "TP", solver)
    carbon_perc = percolation_check(labels == CARBON, "TP")
    si = labels == SI
    ring = ndi.binary_dilation(si, structure=_FOUR) & ~si
    exposed = (ring & (labels == PORE)).sum() / max(ring.sum(), 1)
    return dict(sigma_eff_rel_TP=tp["D_eff_rel"], flux_balance=tp.get("flux_balance", np.nan),
                carbon_spanning_frac=carbon_perc["spanning_frac"], exposed_si_frac=float(exposed),
                undefined=[] if tp["defined"] else [f"electronic TP: {tp['reason']}"])
