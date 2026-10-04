"""Section 0.3 development loop: plausibility checks and KPI snapshot regression.

A change is insignificant when |Δ| < 0.1 × the current baseline scale_k and, for fraction KPIs, also
< 1 % relative. Anything else is a behaviour change and fails the run unless accepted with --accept.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .kpis.registry import KPI_NAMES, unit

# brief 0.3.2(d): anything outside is a bug until proven material
PLAUSIBLE = {
    "si_frac_solid": (0.03, 0.15),
    "pore_frac_deep": (0.08, 0.25),
    "aniso_pore_chord_ratio": (1.1, 2.5),
    "si_d50_aw_um": (2.0, 8.0),
}


def plausibility(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in df.iterrows():
        for k, (lo, hi) in PLAUSIBLE.items():
            v = r.get(k, np.nan)
            ok = bool(np.isfinite(v) and lo <= v <= hi)
            rows.append(dict(sample_id=r["sample_id"], batch=r["batch"], kpi=k, value=v, lo=lo, hi=hi, ok=ok))
    return pd.DataFrame(rows)


def snapshots(out: Path) -> list[Path]:
    d = Path(out) / "snapshots"
    files = list(d.glob("kpis_*.csv")) if d.exists() else []
    return sorted(files, key=lambda p: int(re.search(r"kpis_(\d+)", p.name).group(1)))


def compare(old: pd.DataFrame, new: pd.DataFrame, stats: dict, abs_tol_scale: float = 0.1, rel_tol_frac: float = 0.01):
    old = old.set_index("sample_id")
    new = new.set_index("sample_id")
    common = [s for s in new.index if s in old.index]
    changes = []
    baseline_driven = []
    for sid in common:
        thr_changed = ("crack_threshold" in old.columns and "crack_threshold" in new.columns
                       and not np.isclose(float(old.at[sid, "crack_threshold"]), float(new.at[sid, "crack_threshold"])))
        for k in KPI_NAMES:
            if k not in new.columns or k not in old.columns:
                continue
            if thr_changed and k.startswith("crack_"):
                baseline_driven.append(dict(sample_id=sid, kpi=k, old=float(old.at[sid, k]), new=float(new.at[sid, k])))
                continue
            a, b = float(old.at[sid, k]), float(new.at[sid, k])
            if not np.isfinite(a) and not np.isfinite(b):
                continue
            if np.isfinite(a) != np.isfinite(b):
                changes.append(dict(sample_id=sid, kpi=k, old=a, new=b, delta_scale=np.inf, rel=np.inf))
                continue
            sc = (stats.get("kpi", {}).get(k) or {}).get("scale", np.nan)
            d = abs(b - a)
            ds = d / sc if np.isfinite(sc) and sc > 0 else (0.0 if d == 0 else np.inf)
            rel = d / abs(a) if a != 0 else (0.0 if d == 0 else np.inf)
            bad = ds >= abs_tol_scale or (unit(k) == "fraction" and rel >= rel_tol_frac)
            if bad:
                changes.append(dict(sample_id=sid, kpi=k, old=a, new=b, delta_scale=ds, rel=rel))
    return dict(common=common, new_samples=[s for s in new.index if s not in old.index],
                missing=[s for s in old.index if s not in new.index], changes=pd.DataFrame(changes),
                baseline_driven=pd.DataFrame(baseline_driven))


def run(out: Path, accept: bool = False, log=print) -> int:
    out = Path(out)
    new = pd.read_csv(out / "kpis.csv", dtype={"batch": str})
    stats = json.loads((out / "baseline_stats.json").read_text())
    pl = plausibility(new)
    bad_pl = pl[~pl["ok"]]
    for _, r in bad_pl.iterrows():
        log(f"IMPLAUSIBLE {r['sample_id']} (batch {r['batch']}): {r['kpi']} = {r['value']:.4g} outside [{r['lo']}, {r['hi']}]")
    snaps = snapshots(out)
    (out / "snapshots").mkdir(exist_ok=True)
    if not snaps:
        p = out / "snapshots" / "kpis_1.csv"
        new.to_csv(p, index=False)
        log(f"no snapshot yet: wrote {p} ({len(new)} samples) as the first accepted snapshot")
        return 0
    last = snaps[-1]
    n = int(re.search(r"kpis_(\d+)", last.name).group(1))
    res = compare(pd.read_csv(last, dtype={"batch": str}), new, stats)
    ch = res["changes"]
    log(f"regress vs {last.name}: {len(res['common'])} common samples, {len(res['new_samples'])} new "
        f"({', '.join(res['new_samples']) or '-'}), {len(res['missing'])} missing; {len(ch)} behaviour change(s)")
    if len(res["baseline_driven"]):
        log(f"  ({len(res['baseline_driven'])} crack-KPI changes are baseline-driven: the crack threshold moved with the "
            "baseline membership; not counted until it is frozen)")
    if len(ch):
        summ = ch.groupby("kpi").agg(n=("sample_id", "size"), max_delta_scale=("delta_scale", "max")).sort_values("n", ascending=False)
        log(summ.to_string())
    status = 0 if len(ch) == 0 else 1
    if accept or (status == 0 and res["new_samples"]):
        p = out / "snapshots" / f"kpis_{n + 1}.csv"
        new.to_csv(p, index=False)
        log(f"{'accepted behaviour change; ' if status else ''}wrote snapshot {p.name}")
        status = 0
    elif status:
        log("FAIL: behaviour change. If intended, record why in DECISIONS.md and rerun with --accept.")
    return status
