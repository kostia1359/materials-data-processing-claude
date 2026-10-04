"""Simulation-layer orchestration.

Per image (cached stage, never touches the KPI tables): fused map -> x4 block-majority maps with audit ->
two-conductivity ionic solves (L_mid TP/IP, L_solid/L_pore TP bounds, D_c sweep), electronic solve,
swelling scenarios, ICL proxy, per-strip spread and checks (bound ordering, flux conservation).
Baseline level: medians, ratios to baseline, PyBaMM runs under both porosity conventions, trust grades.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np

from . import cell, fuse, indices, laplace, swell
from .laplace import CARBON, PORE, UNCERTAIN

SIM_VERSION = 1
QC_DIR = Path(__file__).resolve().parent


def sim_stage_version(cfg: dict) -> str:
    h = hashlib.sha1(str(SIM_VERSION).encode())
    for name in ("fuse.py", "laplace.py", "swell.py", "indices.py", "run.py"):  # per-image code only
        h.update((QC_DIR / name).read_bytes())
    h.update(json.dumps(cfg["sim"], sort_keys=True).encode())
    return h.hexdigest()[:10]


def assumptions_hash(cfg: dict) -> str:
    return hashlib.sha1(json.dumps(cfg["sim"], sort_keys=True).encode()).hexdigest()[:10]


def _resolve(ds: np.ndarray, to: int) -> np.ndarray:
    return np.where(ds == UNCERTAIN, to, ds).astype(np.uint8)


def per_image(sample, seg, cfg: dict, n_strips: int = 5) -> tuple[dict, dict]:
    """Returns (row, maps) for one image. Row values are absolute under the frozen assumptions;
    ratios to baseline are formed later."""
    t0 = time.time()
    sc = cfg["sim"]
    f = sc["downsample"]["factor"]
    px = sample.px_um * f
    F = fuse.fuse_phase_map(seg.bse_s, seg.etd_s, sample.inlens, seg.t1, seg.t2, cfg, si_mask=seg.si)
    ds_mid, audit = fuse.downsample_labels(F["L_mid"], f)
    # block ties: carbon for the mid and solid maps, pore for the pore-bound map (keeps the sets nested)
    maps = {"mid": _resolve(ds_mid, CARBON), "solid": _resolve(fuse.downsample_labels(F["L_solid"], f)[0], CARBON),
            "pore": _resolve(fuse.downsample_labels(F["L_pore"], f)[0], PORE)}
    row = dict(uncertain_frac=F["uncertain_frac"], conflict_mean=F["conflict_mean"], ds_tie_frac=audit["tie_frac"],
               ds_d_percolating=audit["d_percolating"], ds_audit_ok=bool(audit["d_percolating"] < sc["downsample"]["max_percolating_fraction_change"]),
               ds_euler_before=audit["before"]["euler_density_mm2"], ds_euler_after=audit["after"]["euler_density_mm2"],
               ds_components_before=audit["before"]["components"], ds_components_after=audit["after"]["components"])
    tr = sc["transport"]
    undefined = []
    # All sparse solves are independent: run them in a thread pool (SuperLU releases the GIL), then
    # collect. Results are identical to the serial order; only wall time changes.
    from concurrent.futures import ThreadPoolExecutor

    def _tp(name, dc, cols=None):
        m = maps[name] if cols is None else maps[name][:, cols]
        return laplace.fv_laplace(laplace.cond_map(m, {PORE: 1.0, CARBON: dc, 2: 0.0}), "TP", tr["solver"])

    w = maps["mid"].shape[1]
    edges = np.linspace(0, w, n_strips + 1).astype(int)
    pool = ThreadPoolExecutor(max_workers=sc.get("threads", 4))
    f_ion = pool.submit(laplace.ionic_index, maps["mid"], tr["D_c"], tr["solver"])
    f_el = pool.submit(laplace.electronic_index, maps["mid"], tr["sigma_si"], tr["solver"])
    f_sweep = {(name, dc): pool.submit(_tp, name, dc) for dc in tr["D_c_sweep"] for name in ("solid", "mid", "pore")
               if not (dc == tr["D_c"] and name == "mid")}
    f_strip = {(i, name): pool.submit(_tp, name, tr["D_c"], slice(a, b))
               for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])) for name in ("solid", "mid", "pore")}
    ion = f_ion.result()
    undefined += ion["undefined"]
    row.update({k: ion[k] for k in ("eps_pore", "D_eff_rel_TP", "D_eff_rel_IP", "aniso_ratio", "tau_p", "N_M",
                                    "flux_balance_TP", "flux_balance_IP", "pore_spans_TP", "pore_clusters")})
    sweep = {}
    for dc in tr["D_c_sweep"]:
        for name in ("solid", "mid", "pore"):
            if dc == tr["D_c"] and name == "mid":
                v = ion["D_eff_rel_TP"]
            else:
                v = f_sweep[(name, dc)].result().get("D_eff_rel", np.nan)
            sweep[f"{name}@{dc}"] = v
    row["D_eff_rel_TP_solid"] = sweep[f"solid@{tr['D_c']}"]
    row["D_eff_rel_TP_pore"] = sweep[f"pore@{tr['D_c']}"]
    row["bounds_ordered"] = bool(all(sweep[f"solid@{dc}"] <= sweep[f"mid@{dc}"] * (1 + 1e-9) <= sweep[f"pore@{dc}"] * (1 + 1e-9) for dc in tr["D_c_sweep"]))
    row["D_c_sweep"] = sweep
    el = f_el.result()
    undefined += el["undefined"]
    row.update({k: el[k] for k in ("sigma_eff_rel_TP", "carbon_spanning_frac", "exposed_si_frac")})
    row["icl_raw"] = indices.icl_index(maps["mid"], px)
    sw = sc["swelling"]
    for fa in sw["f_A_scenarios"]:
        r = swell.swell(maps["mid"], fa, "pore_first", px, sw["cluster_hotspot_um"], return_maps=(fa == max(sw["f_A_scenarios"])))
        if r.get("maps"):
            maps["swollen"] = r["maps"]["swollen"]
            maps["constraint_map"] = r["maps"]["constraint_map"]
        for k, v in r.items():
            if k in ("maps", "f_A", "mode"):
                continue
            if isinstance(v, dict):
                for kk, vv in v.items():
                    row[f"swell{fa}_{k}_{kk}"] = vv
            else:
                row[f"swell{fa}_{k}"] = v
    fa = max(sw["f_A_scenarios"])
    for name in ("solid", "pore"):
        r = swell.swell(maps[name], fa, "pore_first", px, sw["cluster_hotspot_um"])
        row[f"swell{fa}_buffer_sufficiency_frac_{name}"] = r.get("buffer_sufficiency_frac", np.nan)
        row[f"swell{fa}_constraint_index_{name}"] = r.get("constraint_index", np.nan)
    # per-strip spread (empirical REV) and per-strip checks
    strip_vals, strip_ok, strip_flux = [], [], []
    for i in range(n_strips):
        vals = []
        for name in ("solid", "mid", "pore"):
            rr = f_strip[(i, name)].result()
            vals.append(rr.get("D_eff_rel", np.nan))
            if name == "mid":
                strip_flux.append(rr.get("flux_balance", np.nan))
        strip_vals.append(vals[1])
        if all(np.isfinite(v) for v in vals):  # strips where a solve was refused (nothing spans) are not compared
            strip_ok.append(bool(vals[0] <= vals[1] * (1 + 1e-9) <= vals[2] * (1 + 1e-9)))
    sv = np.array([v if v is not None else np.nan for v in strip_vals], float)
    row["strip_D_eff_rel_TP"] = strip_vals
    row["strip_n_undefined"] = int((~np.isfinite(sv)).sum())
    fin = sv[np.isfinite(sv)]
    row["strip_D_eff_rel_TP_cv"] = float(np.std(fin, ddof=1) / np.mean(fin)) if len(fin) > 1 else np.nan
    pool.shutdown()
    row["strip_bounds_ordered"] = bool(all(strip_ok))
    row["strip_max_flux_imbalance"] = float(np.nanmax(strip_flux)) if strip_flux else np.nan
    row["undefined"] = undefined
    row["seconds"] = round(time.time() - t0, 1)
    maps["vis"] = fuse.downsample_labels(F["L_vis"], f)[0]
    return row, maps


# ---------------------------------------------------------------- baseline-relative layer

REL_KEYS = ["D_eff_rel_TP", "D_eff_rel_IP", "aniso_ratio", "tau_p", "sigma_eff_rel_TP", "icl_raw",
            "swell2.43_pore_closure_frac", "swell2.43_constraint_index_aw", "swell2.43_buffer_sufficiency_frac",
            "swell2.43_si_touch_frac", "swell2.43_largest_merged_cluster_ecd_um", "swell2.43_dH_H_bound"]


def baseline_summary(rows: dict, batches: dict, cfg: dict) -> dict:
    """rows: sample_id -> per-image row. Medians and 10-90 % bands over the baseline (batch 3)."""
    base = [r for s, r in rows.items() if batches.get(s) == "3"]
    out = {"n": len(base), "median": {}, "p10": {}, "p90": {}, "strip_cv_median": {}}
    keys = set(REL_KEYS) | {k for r in base for k in r if k.startswith("swell") and isinstance(r[k], (int, float))} | {"eps_pore"}
    for k in keys:
        v = np.array([r.get(k, np.nan) for r in base], float)
        v = v[np.isfinite(v)]
        if len(v):
            out["median"][k], out["p10"][k], out["p90"][k] = float(np.median(v)), float(np.percentile(v, 10)), float(np.percentile(v, 90))
    sweep = {}
    for key in (base[0]["D_c_sweep"] if base else {}):
        v = np.array([r["D_c_sweep"].get(key, np.nan) for r in base], float)
        v = v[np.isfinite(v)]
        if len(v):
            sweep[key] = dict(median=float(np.median(v)), p10=float(np.percentile(v, 10)), p90=float(np.percentile(v, 90)))
    out["D_c_sweep"] = sweep
    out["strip_cv_median"]["D_eff_rel_TP"] = float(np.nanmedian([r.get("strip_D_eff_rel_TP_cv", np.nan) for r in base])) if base else np.nan
    return out


def relative(row: dict, kpis: dict, base_sim: dict, base_kpi_med: dict, cfg: dict, run_cell: bool = True, se: dict | None = None) -> dict:
    """Ratios to baseline + PyBaMM under both porosity conventions; trust grades."""
    sc = cfg["sim"]
    med = base_sim["median"]
    rel = {}
    for k in REL_KEYS:
        v, b = row.get(k, np.nan), med.get(k, np.nan)
        rel[k] = float(v / b) if np.isfinite(v) and np.isfinite(b) and b != 0 else np.nan
    tp = transport = indices.transport_proxies(row.get("D_eff_rel_TP", np.nan), cfg)
    tb = indices.transport_proxies(med.get("D_eff_rel_TP", np.nan), cfg)
    for k in tp:
        rel[k] = tp[k] / tb[k] if np.isfinite(tp[k]) and np.isfinite(tb[k]) and tb[k] else np.nan
    rel["plating_risk_index"] = indices.plating_risk_index(row.get("D_eff_rel_TP", np.nan), med.get("D_eff_rel_TP", np.nan))
    rel["cli_index"] = indices.cli_index(
        dict(kpis, buffer_sufficiency_frac=row.get("swell2.43_buffer_sufficiency_frac")),
        dict(base_kpi_med, buffer_sufficiency_frac=med.get("swell2.43_buffer_sufficiency_frac")), cfg)
    bounds = {"D_eff_rel_TP": [row.get("D_eff_rel_TP_solid"), row.get("D_eff_rel_TP_pore")],
              "buffer_sufficiency_frac@2.43": [row.get("swell2.43_buffer_sufficiency_frac_solid"), row.get("swell2.43_buffer_sufficiency_frac_pore")]}
    for k in ("D_eff_rel_TP",):
        b = med.get(k)
        if b:
            bounds[f"{k}_ratio"] = [v / b if v is not None and np.isfinite(v) else None for v in bounds[k]]
    conv = {}
    f_si = kpis.get("si_frac_solid", np.nan)
    r_si = kpis.get("si_d50_aw_um", np.nan) / 2 * sc["cell"]["si_radius_factor"]
    tau_rel = rel.get("tau_p", np.nan)
    for c in sc["electrode"]["porosity_convention"]:
        eps = indices.porosity(kpis.get("pore_frac_deep", np.nan), base_kpi_med.get("pore_frac_deep"), c)
        eps_b = indices.porosity(base_kpi_med.get("pore_frac_deep", np.nan), base_kpi_med.get("pore_frac_deep"), c)
        e = indices.energy_density(eps, f_si, cfg)
        eb = indices.energy_density(eps_b, base_kpi_med.get("si_frac_solid", np.nan), cfg)
        d = dict(eps=eps, energy_density_mAh_cm3=e, energy_density_ratio=e / eb if eb else np.nan)
        if run_cell and np.isfinite(tau_rel) and np.isfinite(r_si) and np.isfinite(f_si):
            tau_factor = eps_b ** (1 - 1.5) * tau_rel  # baseline median maps to Bruggeman 1.5
            try:
                res = cell.pybamm_runs(eps, f_si, r_si, tau_factor, cfg, se=se if c == sc["electrode"]["porosity_convention"][0] else None)
                d["pybamm"] = res
            except Exception as ex:  # noqa: BLE001
                d["pybamm_error"] = str(ex)[:200]
        conv[c] = d
    grades = {k: indices.GRADES.get(k, "B") for k in list(rel) + ["energy_density_mAh_cm3", "Q_CC_3C_over_Q_C10"]}
    grades["cli_index"] = "C"
    return dict(relative_to_baseline=rel, bounds=bounds, conventions=conv, grade=grades,
                undefined=row.get("undefined", []), assumptions_hash=assumptions_hash(cfg),
                checks=dict(bounds_ordered=row.get("bounds_ordered"), strip_bounds_ordered=row.get("strip_bounds_ordered"),
                            flux_balance_TP=row.get("flux_balance_TP"), downsample_audit_ok=row.get("ds_audit_ok")))
