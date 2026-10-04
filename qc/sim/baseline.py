"""build-baseline side of the simulation layer: per-image stage for all samples, baseline summary,
relative indices + PyBaMM (with ±SE runs), rank robustness across conventions, output tables."""
from __future__ import annotations

import hashlib
import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from ..analyze import _clean, _unclean, analyze_sim
from . import phantoms
from .run import assumptions_hash, baseline_summary, relative

KPI_FOR_SIM = ["si_frac_solid", "si_d50_aw_um", "pore_frac_deep", "si_agglomerate_frac", "si_contact_carbon_frac"]


def _sim_worker(args):
    files, cfg, cache = args
    return files.sample_id, analyze_sim(files, cfg, cache)


def _rel_worker(args):
    sid, row, kp, base_sim, base_med, cfg, se, cache = args
    code = (Path(__file__).parent / "cell.py").read_bytes() + (Path(__file__).parent / "run.py").read_bytes() + (Path(__file__).parent / "indices.py").read_bytes()
    key = hashlib.sha1(json.dumps(_clean(dict(row=row, kp=kp, b=base_sim, m=base_med, se=se, a=cfg["sim"])), sort_keys=True, default=str).encode()
                       + code).hexdigest()[:12]
    p = Path(cache) / sid / "sim_rel.json"
    if p.exists():
        d = json.loads(p.read_text())
        if d.get("key") == key:
            return sid, _unclean(d["rel"])
    rel = relative(row, kp, base_sim, base_med, cfg, run_cell=True, se=se)
    p.write_text(json.dumps(dict(key=key, rel=_clean(rel))))
    return sid, rel


def _se_for_cell(r: pd.Series, row: dict) -> dict:
    tau_se = row.get("strip_D_eff_rel_TP_cv", np.nan) / np.sqrt(5)  # relative SE of the transport index
    return dict(eps=r.get("pore_frac_deep_se", np.nan), f_si=r.get("si_frac_solid_se", np.nan),
                r_si_um=r.get("si_d50_aw_um_se", np.nan) / 2 * 1.273, tau_factor=tau_se)


def run(samples, df: pd.DataFrame, cfg: dict, out: Path, workers: int | None = None, log=print) -> dict:
    cache = out / "cache"
    workers = workers or int(os.environ.get("QC_WORKERS", min(3, os.cpu_count() or 1)))
    rows = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for sid, row in ex.map(_sim_worker, [(s, cfg, cache) for s in samples]):
            rows[sid] = row
            log(f"  sim {sid} ({row.get('seconds', '?')} s)")
    batches = dict(zip(df["sample_id"], df["batch"].astype(str)))
    base_sim = baseline_summary(rows, batches, cfg)
    base_med = df[df["batch"] == "3"][KPI_FOR_SIM].median().to_dict()
    kp = df.set_index("sample_id")
    jobs = []
    for sid, row in rows.items():
        r = kp.loc[sid]
        jobs.append((sid, row, {k: r.get(k) for k in KPI_FOR_SIM}, base_sim, base_med, cfg, _se_for_cell(r, row), str(cache)))
    rels = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for sid, rel in ex.map(_rel_worker, jobs):
            rels[sid] = rel
    sim_rows, bound_rows = [], []
    for sid, row in rows.items():
        rel = rels[sid]
        flat = dict(sample_id=sid, batch=batches.get(sid))
        flat.update({k: v for k, v in row.items() if not isinstance(v, (list, dict))})
        flat.update({f"rel_{k}": v for k, v in rel["relative_to_baseline"].items()})
        for c, d in rel["conventions"].items():
            tag = "off" if c.startswith("offset") else "asis"
            flat[f"eps_{tag}"] = d["eps"]
            flat[f"energy_ratio_{tag}"] = d.get("energy_density_ratio")
            pt = (d.get("pybamm") or {}).get("point", {})
            for k in ("Q_CC_1C_over_Q_C10", "Q_CC_3C_over_Q_C10", "min_neg_surface_dphi_sep_3C_V", "min_electrolyte_conc_3C", "N_P"):
                flat[f"{k}_{tag}"] = pt.get(k)
        sim_rows.append(flat)
        b = dict(sample_id=sid, batch=batches.get(sid), D_eff_rel_TP_solid=row.get("D_eff_rel_TP_solid"),
                 D_eff_rel_TP_mid=row.get("D_eff_rel_TP"), D_eff_rel_TP_pore=row.get("D_eff_rel_TP_pore"),
                 buffer_suff_solid=row.get("swell2.43_buffer_sufficiency_frac_solid"), buffer_suff_mid=row.get("swell2.43_buffer_sufficiency_frac"),
                 buffer_suff_pore=row.get("swell2.43_buffer_sufficiency_frac_pore"))
        for k, v in (row.get("D_c_sweep") or {}).items():
            b[f"sweep_{k}"] = v
        conv0 = next(iter(rel["conventions"].values()))
        for name, res in (conv0.get("pybamm") or {}).items():
            if name != "point":
                b[f"pybamm_{name}_Q1C"] = res.get("Q_CC_1C_over_Q_C10")
        bound_rows.append(b)
    sim_df = pd.DataFrame(sim_rows)
    sim_df.to_csv(out / "sim.csv", index=False)
    pd.DataFrame(bound_rows).to_csv(out / "sim_bounds.csv", index=False)
    robust = rank_robustness(sim_df, pd.DataFrame(bound_rows))
    ph = phantoms.run_all()
    summary = dict(base_sim, assumptions_hash=assumptions_hash(cfg), rank_robustness=robust, phantoms=ph,
                   base_kpi_median=base_med, checks=dict(
                       bounds_ordered_all=bool(sim_df["bounds_ordered"].all()), strip_bounds_ordered_all=bool(sim_df["strip_bounds_ordered"].all()),
                       max_flux_imbalance=float(np.nanmax(sim_df[["flux_balance_TP", "flux_balance_IP", "strip_max_flux_imbalance"]].to_numpy(float))),
                       downsample_audit_ok_all=bool(sim_df["ds_audit_ok"].all()),
                       strips_undefined_total=int(sim_df["strip_n_undefined"].sum()) if "strip_n_undefined" in sim_df else None))
    (out / "sim_baseline.json").write_text(json.dumps(_clean(summary), indent=1))
    log(f"simulation layer: {len(rows)} samples; bounds ordered: {summary['checks']['bounds_ordered_all']}; "
        f"max flux imbalance {summary['checks']['max_flux_imbalance']:.1e}; phantoms ok: {all(p['ok'] for p in ph)}")
    return summary


def rank_robustness(sim_df: pd.DataFrame, bounds: pd.DataFrame, min_rho: float = 0.8) -> dict:
    from scipy.stats import spearmanr

    out = {}
    cols = [c for c in bounds.columns if c.startswith("sweep_mid@")]
    if len(cols) >= 2:
        ref = "sweep_mid@0.05" if "sweep_mid@0.05" in cols else cols[0]
        for c in cols:
            if c != ref:
                rho = spearmanr(bounds[ref], bounds[c]).correlation
                out[f"D_eff_TP rank {ref.split('@')[1]} vs {c.split('@')[1]}"] = dict(rho=float(rho), stable=bool(rho >= min_rho))
    for k in ("energy_ratio", "Q_CC_3C_over_Q_C10", "Q_CC_1C_over_Q_C10"):
        a, b = f"{k}_off" if k == "energy_ratio" else f"{k}_off", f"{k}_asis"
        if a in sim_df and b in sim_df and sim_df[[a, b]].notna().sum().min() > 3:
            rho = spearmanr(sim_df[a], sim_df[b], nan_policy="omit").correlation
            out[f"{k} rank: porosity offset vs as-is"] = dict(rho=float(rho), stable=bool(rho >= min_rho))
    return out
