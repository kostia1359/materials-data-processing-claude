"""Section 10 exact cases and invariants."""
import numpy as np
import pytest

from qc.sim import laplace, swell
from qc.sim.phantoms import checkerboard, run_all, serpentine, straight_channels, tilted_channel


def test_phantoms_exact_cases():
    for r in run_all():
        assert r["ok"], r


def test_refuses_non_spanning():
    m = np.zeros((50, 50))
    m[10:20, 10:20] = 1
    r = laplace.fv_laplace(m)
    assert not r["defined"] and "spans" in r["reason"]


def test_d_eff_below_porosity_and_flux_conserved():
    rng = np.random.default_rng(0)
    lab = rng.choice([0, 1, 2], size=(120, 120), p=[0.4, 0.5, 0.1]).astype(np.uint8)
    c = laplace.cond_map(lab, {0: 1.0, 1: 0.05, 2: 0.0})
    r = laplace.fv_laplace(c)
    eps_eff = (c > 0).mean() * c[c > 0].mean()
    assert r["D_eff_rel"] <= eps_eff + 1e-12
    assert r["flux_balance"] < 1e-8


def test_bounds_nested():
    rng = np.random.default_rng(1)
    mid = rng.choice([0, 1, 2], size=(100, 100), p=[0.3, 0.6, 0.1]).astype(np.uint8)
    unc = (rng.random(mid.shape) < 0.1) & (mid != 2)
    solid = np.where(unc, 1, mid)
    pore = np.where(unc, 0, mid)
    v = [laplace.ionic_index(m, 0.05)["D_eff_rel_TP"] for m in (solid, mid, pore)]
    assert v[0] <= v[1] <= v[2]


def test_flip_symmetry_and_reciprocity():
    rng = np.random.default_rng(2)
    c = rng.choice([1.0, 0.05], size=(80, 90), p=[0.4, 0.6])
    a = laplace.fv_laplace(c, "TP")["D_eff_rel"]
    b = laplace.fv_laplace(c[::-1], "TP")["D_eff_rel"]  # reversing the gradient direction
    d = laplace.fv_laplace(c[:, ::-1], "TP")["D_eff_rel"]  # mirror image
    assert abs(a - b) / a < 1e-8 and abs(a - d) / a < 1e-8
    assert abs(laplace.fv_laplace(c, "IP")["D_eff_rel"] - laplace.fv_laplace(c.T, "TP")["D_eff_rel"]) < 1e-12


def test_swelling_area_conservation():
    lab = np.zeros((80, 80), np.uint8)
    lab[30:40, 30:40] = 2
    for mode in ("pore_first", "isotropic"):
        r = swell.swell(lab, 2.0, mode, return_maps=True)
        assert (r["maps"]["swollen"] == 2).sum() == 200
    lab2 = np.ones((80, 80), np.uint8)  # all carbon: everything displaced, CI = 1
    lab2[30:40, 30:40] = 2
    r = swell.swell(lab2, 2.0)
    assert r["constraint_index"] == 1.0 and r["pore_closure_frac"] == 0


def test_pybamm_tortuosity_factor_equals_bruggeman():
    pybamm = pytest.importorskip("pybamm")
    pv = pybamm.ParameterValues("Chen2020_composite")
    base_opts = {"particle phases": ("2", "1"), "open-circuit potential": (("single", "current sigmoid"), "single")}
    exp = pybamm.Experiment([("Discharge at C/10 until 2.5 V", "Rest for 15 minutes", "Charge at 1C until 4.2 V")])

    def q_cc(model, values):
        sol = pybamm.Simulation(model, parameter_values=values, experiment=exp, solver=pybamm.IDAKLUSolver()).solve()
        st = sol.cycles[0].steps[-1]
        return abs(st["Discharge capacity [A.h]"].entries[-1] - st["Discharge capacity [A.h]"].entries[0])

    q_brugg = q_cc(pybamm.lithium_ion.DFN(base_opts), pv.copy())
    pv2 = pv.copy()
    for dom, e in (("Negative electrode", "Negative electrode porosity"), ("Separator", "Separator porosity"),
                   ("Positive electrode", "Positive electrode porosity")):
        b = pv2[f"{dom} Bruggeman coefficient (electrolyte)"]
        key_e = f"{dom} Bruggeman coefficient (electrode)"
        be = pv2[key_e] if key_e in pv2.keys() else 0.0
        pv2.update({f"{dom} tortuosity factor (electrolyte)": pv2[e] ** (1 - b),
                    f"{dom} tortuosity factor (electrode)": 1.0 if be == 0 else float("nan")}, check_already_exists=False)
    q_tau = q_cc(pybamm.lithium_ion.DFN({**base_opts, "transport efficiency": "tortuosity factor"}), pv2)
    assert abs(q_tau - q_brugg) / q_brugg < 0.005
