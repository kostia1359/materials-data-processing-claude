"""PyBaMM composite graphite-Si DFN: relative rate-capability and plating-indicator shifts (brief 10.3).

Frozen LG-M50-type parameter set (Chen2020_composite); only the negative-electrode microstructure
parameters are mapped from the KPIs. Outputs are compared with the baseline median under identical
assumptions - never read as predictions of the real cell.
"""
from __future__ import annotations

import time

import numpy as np

_CACHE: dict = {}


def _model(cfg):
    import pybamm

    key = "model"
    if key not in _CACHE:
        opts = {"particle phases": ("2", "1"), "open-circuit potential": (("single", "current sigmoid"), "single"),
                "transport efficiency": "tortuosity factor"}
        _CACHE[key] = pybamm.lithium_ion.DFN(opts)
    return _CACHE[key]


def map_parameters(eps: float, f_si: float, r_si_um: float, tau_factor: float, cfg: dict) -> dict:
    """KPI -> parameter mapping. eps = electrode porosity (convention applied by the caller);
    f_si = Si-candidate share of solids; CBD = cbd_frac_of_solids of the solids (inactive)."""
    ec = cfg["sim"]["electrode"]
    solids = 1.0 - eps
    cbd = ec["cbd_frac_of_solids"] * solids
    v_si = f_si * solids
    v_gr = max(solids - v_si - cbd, 0.01)
    return {
        "Negative electrode porosity": eps,
        "Primary: Negative electrode active material volume fraction": v_gr,
        "Secondary: Negative electrode active material volume fraction": v_si,
        "Secondary: Negative particle radius [m]": r_si_um * 1e-6,
        "Negative electrode thickness [m]": ec["L_um"] * 1e-6,
        "Negative electrode tortuosity factor (electrolyte)": tau_factor,
    }


def run_point(params: dict, rates=(1, 3), cfg: dict | None = None) -> dict:
    import pybamm

    model = _model(cfg)
    pv = pybamm.ParameterValues("Chen2020_composite")
    # Bruggeman-style tortuosity factors required by the "tortuosity factor" transport option
    for dom, eps_key in (("Negative electrode", "Negative electrode porosity"), ("Separator", "Separator porosity"),
                         ("Positive electrode", "Positive electrode porosity")):
        b = pv[f"{dom} Bruggeman coefficient (electrolyte)"]
        pv.update({f"{dom} tortuosity factor (electrolyte)": pv[eps_key] ** (1 - b),
                   f"{dom} tortuosity factor (electrode)": 1.0}, check_already_exists=False)
    pv.update(params, check_already_exists=False)
    out = {}
    try:
        solver = pybamm.IDAKLUSolver()
    except Exception:  # noqa: BLE001
        solver = pybamm.CasadiSolver(mode="safe")
    q_ref = None
    for rate in rates:
        # one tuple = one cycle (a plain list would make every step its own cycle)
        exp = pybamm.Experiment([("Discharge at C/10 until 2.5 V", "Rest for 15 minutes", f"Charge at {rate}C until 4.2 V")])
        t0 = time.time()
        try:
            sim = pybamm.Simulation(model, parameter_values=pv, experiment=exp, solver=solver)
            sol = sim.solve()
        except Exception as e:  # noqa: BLE001 - reported as undefined, never silently zero
            out[f"error_{rate}C"] = str(e)[:200]
            continue
        dis, chg = sol.cycles[0].steps[0], sol.cycles[0].steps[-1]
        q_dis = float(abs(dis["Discharge capacity [A.h]"].entries[-1] - dis["Discharge capacity [A.h]"].entries[0]))
        q_cc = float(abs(chg["Discharge capacity [A.h]"].entries[-1] - chg["Discharge capacity [A.h]"].entries[0]))
        q_ref = q_dis
        out[f"Q_CC_{rate}C_over_Q_C10"] = q_cc / q_dis if q_dis > 0 else np.nan
        try:
            phi = chg["Negative electrode surface potential difference at separator interface [V]"].entries
            out[f"min_neg_surface_dphi_sep_{rate}C_V"] = float(np.min(phi))
            cap = chg["Discharge capacity [A.h]"].entries
            cross = np.where(phi < 0)[0]
            out[f"capacity_at_0V_crossing_{rate}C"] = float(abs(cap[cross[0]] - cap[0]) / q_dis) if len(cross) else np.nan
        except Exception:  # noqa: BLE001
            out[f"min_neg_surface_dphi_sep_{rate}C_V"] = np.nan
        try:
            ce = chg["Electrolyte concentration [mol.m-3]"].entries
            out[f"min_electrolyte_conc_{rate}C"] = float(np.min(ce))
        except Exception:  # noqa: BLE001
            pass
        out[f"seconds_{rate}C"] = round(time.time() - t0, 1)
    out["Q_C10_Ah"] = q_ref
    out["N_P"] = np_ratio(pv)
    return out


def np_ratio(pv) -> float:
    """Areal capacity ratio negative / positive from the parameter set (both phases on the negative side)."""
    F = 96485.0
    neg = 0.0
    for ph in ("Primary: ", "Secondary: "):
        neg += pv[f"{ph}Negative electrode active material volume fraction"] * pv[f"{ph}Maximum concentration in negative electrode [mol.m-3]"]
    neg *= pv["Negative electrode thickness [m]"] * F
    pos = pv["Positive electrode active material volume fraction"] * pv["Maximum concentration in positive electrode [mol.m-3]"]
    pos *= pv["Positive electrode thickness [m]"] * F
    return float(neg / pos)


def pybamm_runs(eps: float, f_si: float, r_si_um: float, tau_factor: float, cfg: dict, se: dict | None = None,
                rates=None) -> dict:
    """Point run plus one-at-a-time ±SE perturbations (eps, f_si, r_si, tau) when `se` is given."""
    rates = rates or tuple(cfg["sim"]["cell"]["rates"])
    point = run_point(map_parameters(eps, f_si, r_si_um, tau_factor, cfg), rates, cfg)
    res = {"point": point}
    if se:
        base = dict(eps=eps, f_si=f_si, r_si_um=r_si_um, tau_factor=tau_factor)
        for k, d in se.items():
            if not np.isfinite(d) or d <= 0:
                continue
            for sgn in (-1, 1):
                p = dict(base)
                p[k] = max(p[k] + sgn * d, 1e-3)
                res[f"{k}{'+' if sgn > 0 else '-'}"] = run_point(map_parameters(p["eps"], p["f_si"], p["r_si_um"], p["tau_factor"], cfg), rates[:1], cfg)
    return res
