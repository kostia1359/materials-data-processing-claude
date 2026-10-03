"""Baseline statistics (batch 3 only), robust z, verdict, conformal rank, defect rule."""
from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd

from .gates import GATE_KEYS
from .kpis.registry import DEFECT_COUNTS, KPI_NAMES, SHORTLIST, meaning, rclass, unit

BASELINE_BATCH = "3"


def baseline_version(ids) -> str:
    return hashlib.sha1(",".join(sorted(ids)).encode()).hexdigest()[:8]


def robust_stats(values: np.ndarray, within_se: np.ndarray | None, rel_floor: float) -> dict:
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return dict(n=0, med=np.nan, mad=np.nan, within_se=np.nan, scale=np.nan, min=np.nan, max=np.nan)
    med = float(np.median(v))
    mad = float(1.4826 * np.median(np.abs(v - med)))
    wse = float(np.nanmedian(within_se)) if within_se is not None and np.isfinite(within_se).any() else 0.0
    scale = max(mad, wse, rel_floor * abs(med))
    if scale <= 0:
        scale = 1e-9
    return dict(n=int(len(v)), med=med, mad=mad, within_se=wse, scale=float(scale), min=float(v.min()), max=float(v.max()))


def baseline_stats(kpis: pd.DataFrame, cfg: dict, exclude: str | None = None) -> dict:
    """kpis: one row per sample with KPI columns, <kpi>_se columns and gate columns."""
    base = kpis[(kpis["batch"].astype(str) == BASELINE_BATCH)]
    if exclude is not None:
        base = base[base["sample_id"] != exclude]
    out = {"kpi": {}, "gate": {}, "ids": sorted(base["sample_id"].tolist())}
    for k in KPI_NAMES + list(DEFECT_COUNTS.values()):
        if k not in base:
            continue
        se = base[f"{k}_se"].to_numpy(float) if f"{k}_se" in base else None
        out["kpi"][k] = robust_stats(base[k].to_numpy(float), se, cfg["scale_rel_floor"])
    for k in GATE_KEYS:
        if k in base:
            out["gate"][k] = robust_stats(base[k].to_numpy(float), None, cfg["scale_rel_floor"])
    if "area_mm2" in base:
        out["area_mm2_median"] = float(base["area_mm2"].median()) if len(base) else np.nan
    out["n"] = int(len(base))
    out["version"] = baseline_version(out["ids"])
    return out


def zscores(x: dict, stats: dict, section: str = "kpi") -> dict:
    z = {}
    for k, s in stats[section].items():
        v = x.get(k, np.nan)
        if v is None or not np.isfinite(v) or not np.isfinite(s["med"]):
            z[k] = np.nan
        else:
            z[k] = float((v - s["med"]) / s["scale"])
    return z


def aggregate_score(z: dict, keys=SHORTLIST) -> float:
    a = np.sort(np.abs([z.get(k, np.nan) for k in keys if np.isfinite(z.get(k, np.nan))]))[::-1]
    return float(a[:3].mean()) if len(a) else np.nan


def defect_alarm(x: dict, stats: dict) -> list[str]:
    """Defect KPI exceeds the baseline maximum by > 3 Poisson SD (counts scaled to this sample's area)."""
    hits = []
    area = x.get("area_mm2", np.nan)
    a_ref = stats.get("area_mm2_median", np.nan)
    for dens, cnt in DEFECT_COUNTS.items():
        if cnt not in stats["kpi"] or not np.isfinite(x.get(cnt, np.nan)):
            continue
        mx = stats["kpi"][cnt]["max"]  # largest event count in any baseline image
        if not np.isfinite(mx):
            continue
        expected = mx * (area / a_ref if np.isfinite(a_ref) and a_ref > 0 else 1.0)
        observed = x[cnt]
        if observed > expected + 3 * np.sqrt(max(expected, 1.0)):
            hits.append(dens)
    return hits


def verdict(x: dict, stats: dict, cfg: dict, s_baseline: list[float] | None = None) -> dict:
    z = zscores(x, stats)
    zs = {k: z.get(k, np.nan) for k in SHORTLIST}
    absz = {k: abs(v) for k, v in zs.items() if np.isfinite(v)}
    zc = cfg["zones"]
    n_rej = sum(v > zc["reject"] for v in absz.values())
    n_25 = sum(v > zc["investigate"] for v in absz.values())
    n_2 = sum(v > zc["investigate_count_z"] for v in absz.values())
    defects = defect_alarm(x, stats)
    reasons = []
    if n_rej:
        reasons.append(f"{n_rej} KPI(s) beyond {zc['reject']} robust sigma")
    if n_25 >= zc["multi_count_n"]:
        reasons.append(f"{n_25} KPIs beyond {zc['multi_count_z']} robust sigma")
    if defects:
        reasons.append("defect count above baseline maximum: " + ", ".join(defects))
    if reasons:
        decision = "REJECT"
    elif n_25 or n_2 >= zc["investigate_count_n"]:
        decision = "INVESTIGATE"
        reasons.append(f"{n_25} KPI(s) beyond {zc['investigate']} sigma, {n_2} beyond {zc['investigate_count_z']} sigma")
    else:
        decision = "ACCEPT"
    S = aggregate_score(z)
    n = stats.get("n", 0)
    p_floor = 1.0 / (n + 1)
    if s_baseline:
        p_rank = (1 + sum(si >= S for si in s_baseline)) / (len(s_baseline) + 1)
        p_floor = 1.0 / (len(s_baseline) + 1)
    else:
        p_rank = np.nan
    drivers = []
    for k in sorted(absz, key=lambda k: -absz[k]):
        s = stats["kpi"][k]
        direction = "higher" if zs[k] > 0 else "lower"
        drivers.append(dict(kpi=k, value=x.get(k), baseline_median=s["med"], baseline_scale=s["scale"], z=zs[k],
                            direction=direction, robustness=rclass(k), unit=unit(k), meaning=meaning(k, direction)))
    return dict(decision=decision, reasons=reasons, aggregate_score_S=S, conformal_rank_p=p_rank, p_floor=p_floor,
                n_kpis_beyond_2_5=int(n_25), n_kpis_beyond_2=int(n_2), defect_alarms=defects, drivers=drivers, z_all=z)


def gate_status(gates: dict, stats: dict, cfg: dict) -> dict:
    z = zscores(gates, stats, "gate")
    flags = [f"{k} z={v:+.1f}" for k, v in z.items() if np.isfinite(v) and abs(v) > cfg["gate_caution_z"]]
    if gates.get("gate_flags"):
        flags += [f for f in str(gates["gate_flags"]).split(";") if f]
    status = "caution" if any("z=" in f for f in flags) else "ok"
    return dict(status=status, flags=flags, z=z)


def loo_baseline_scores(kpis: pd.DataFrame, cfg: dict) -> dict:
    """S_i for each baseline image scored against the baseline without it (for conformal rank)."""
    base = kpis[kpis["batch"].astype(str) == BASELINE_BATCH]
    out = {}
    for _, row in base.iterrows():
        st = baseline_stats(kpis, cfg, exclude=row["sample_id"])
        if st["n"] == 0:
            continue
        out[row["sample_id"]] = aggregate_score(zscores(row.to_dict(), st))
    return out
