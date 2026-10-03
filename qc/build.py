"""build-baseline and evaluate orchestration."""
from __future__ import annotations

import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from . import assign as assign_mod
from . import stats as stats_mod
from .analyze import SampleResult, _clean, analyze_sample, finalize_cracks
from .gates import GATE_KEYS
from .io import SampleFiles, discover, manifest_frame
from .kpis.registry import KPI_NAMES


def _worker(args):
    files, cfg, cache = args
    return analyze_sample(files, cfg, cache)


def analyze_all(samples: list[SampleFiles], cfg: dict, cache: Path, workers: int | None = None, log=print) -> list[SampleResult]:
    workers = workers or int(os.environ.get("QC_WORKERS", min(3, os.cpu_count() or 1)))
    results = []
    if workers <= 1 or len(samples) <= 1:
        for s in samples:
            t = time.time()
            results.append(analyze_sample(s, cfg, cache))
            log(f"  {s.sample_id} (batch {s.batch}) {time.time() - t:.0f}s")
        return results
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for s, r in zip(samples, ex.map(_worker, [(s, cfg, cache) for s in samples])):
            log(f"  {s.sample_id} (batch {s.batch}) done")
            results.append(r)
    return results


def results_frame(results: list[SampleResult]) -> pd.DataFrame:
    rows = []
    for r in results:
        row = dict(sample_id=r.sample_id, batch=str(r.batch))
        for k in [*KPI_NAMES, "crack_count", "hiZ_inclusion_count", "area_mm2", "si_n_sized"]:
            row[k] = r.kpis.get(k, np.nan)
        for k in KPI_NAMES:
            row[f"{k}_se"] = r.kpi_se.get(k, np.nan)
        for k in GATE_KEYS:
            row[k] = r.gates.get(k, np.nan)
        row["gate_flags"] = r.gates.get("gate_flags", "")
        row["flags"] = ";".join(r.info.get("flags", []))
        rows.append(row)
    return pd.DataFrame(rows)


def strips_frame(results: list[SampleResult]) -> pd.DataFrame:
    rows = []
    for r in results:
        for s in r.strips:
            rows.append(dict(sample_id=r.sample_id, batch=str(r.batch), **s))
    return pd.DataFrame(rows)


def crack_threshold_from(results: list[SampleResult]) -> float:
    base = [r.info["sato_p995"] for r in results if str(r.batch) == stats_mod.BASELINE_BATCH and r.info.get("sato_p995") is not None]
    if not base:
        base = [r.info["sato_p995"] for r in results if r.info.get("sato_p995") is not None]
    return float(np.median(base))


def fit_from_results(results: list[SampleResult], cfg: dict, crack_threshold: float | None = None) -> dict:
    """Everything that depends on the labelled table: baseline stats, LOO S, model, signatures."""
    df = results_frame(results)
    st = stats_mod.baseline_stats(df, cfg)
    st["crack_threshold"] = crack_threshold
    s_loo = stats_mod.loo_baseline_scores(df, cfg)
    model = assign_mod.build_model(df, st, cfg)
    return dict(df=df, stats=st, s_loo=s_loo, model=model)


def build_baseline(data: Path, out: Path, cfg: dict, log=print, workers=None) -> dict:
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    samples = discover(data)
    if not samples:
        raise SystemExit(f"no img_<id>_<det>.tif files found under {data}")
    man = manifest_frame(samples, cfg)
    man.to_csv(out / "manifest.csv", index=False)
    log(f"{len(samples)} samples: " + ", ".join(f"batch {b}: {n}" for b, n in man["batch"].value_counts().sort_index().items()))
    cache = out / "samples"
    results = analyze_all(samples, cfg, cache, workers, log)
    thr = crack_threshold_from(results)
    for r in results:
        if r.info.get("crack_threshold") != thr:
            finalize_cracks(r, cache / r.sample_id, thr, cfg)
    labelled = [r for r in results if str(r.batch) in ("1", "2", "3")]
    fit = fit_from_results(labelled, cfg, thr)
    df_all = results_frame(results)
    df_all.to_csv(out / "kpis.csv", index=False)
    strips_frame(results).to_csv(out / "strips.csv", index=False)
    st = fit["stats"]
    st["loo_S"] = fit["s_loo"]
    st["delta"] = fit["model"]["delta"]
    st["T"] = fit["model"]["T"]
    st["within_distance_percentiles"] = {
        str(q): float(np.percentile(fit["model"]["within_distances"], q)) for q in (50, 90, 97.5)
    } if fit["model"]["within_distances"] else {}
    (out / "baseline_stats.json").write_text(json.dumps(_clean(st), indent=1))
    (out / "model.json").write_text(json.dumps(_clean(fit["model"]), indent=1))
    (out / "signatures.json").write_text(json.dumps(_clean(fit["model"]["signatures"]), indent=1))
    from . import report

    figs = out / "figs"
    figs.mkdir(exist_ok=True)
    report.fig_signatures(fit["model"]["signatures"], figs / "signatures.png")
    report.fig_kpi_boxplots(df_all[df_all["batch"].isin(["1", "2", "3"])], figs / "kpi_boxplots.png")
    if st["n"] < 3:
        log(f"WARNING: only {st['n']} baseline (batch 3) image(s); z-scores rely on within-image strip SE and the 5 % floor.")
    log(f"baseline version {st['version']}  n={st['n']}  crack threshold {thr:.4f}  Δ={st['delta']} T={st['T']}")
    return fit


def load_baseline(out: Path) -> dict:
    out = Path(out)
    st = json.loads((out / "baseline_stats.json").read_text())
    model = json.loads((out / "model.json").read_text())
    df = pd.read_csv(out / "kpis.csv", dtype={"batch": str})
    from .analyze import _unclean

    return dict(stats=_unclean(st), model=_unclean(model), df=df, out=out)


def evaluate_files(files: SampleFiles, base: dict, cfg: dict, report_dir: Path | None = None, log=print) -> dict:
    t0 = time.time()
    st, model = base["stats"], base["model"]
    out = Path(base["out"])
    cache = out / "samples"
    # if the sample is part of the baseline table, score it against the baseline *without* it
    in_table = files.sample_id in set(base["df"]["sample_id"])
    held_out_note = None
    res = analyze_sample(files, cfg, cache, crack_threshold=st.get("crack_threshold"))
    if res.info.get("crack_threshold") != st.get("crack_threshold") and st.get("crack_threshold") is not None:
        finalize_cracks(res, cache / res.sample_id, st["crack_threshold"], cfg)
    x = dict(res.kpis)
    s_loo = list((st.get("loo_S") or {}).values())
    n_base_wo = int(((base["df"]["batch"] == "3") & (base["df"]["sample_id"] != files.sample_id)).sum())
    if in_table and n_base_wo == 0:
        held_out_note = ("sample is the only baseline image, so it is compared with itself: the verdict and assignment "
                         "are a pipeline check only, not evidence")
        in_table = False
    if in_table:
        df = base["df"]
        st = stats_mod.baseline_stats(df[df["batch"].isin(["1", "2", "3"])], cfg, exclude=files.sample_id)
        st["crack_threshold"] = base["stats"].get("crack_threshold")
        lab = df[(df["batch"].isin(["1", "2", "3"])) & (df["sample_id"] != files.sample_id)]
        model = assign_mod.build_model(lab, st, cfg)
        s_loo = list(stats_mod.loo_baseline_scores(lab, cfg).values())
        held_out_note = "sample is in the labelled set; it was held out of the baseline and the centroids for this evaluation"
    v = stats_mod.verdict(x, st, cfg, s_loo)
    g = stats_mod.gate_status(res.gates, st, cfg)
    a = assign_mod.assign(x, res.strips, model, st, cfg)
    elapsed = time.time() - t0
    verdict_json = report.build_verdict_json(res, st, model, v, g, a, cfg, held_out_note, elapsed)
    if report_dir:
        report.write_reports(verdict_json, res, st, model, base, Path(report_dir), cache / res.sample_id)
    log(f"{files.sample_id}: batch {a['assigned_batch']} p={a['probabilities'][a['assigned_batch']]:.2f} "
        f"{v['decision']} ({elapsed:.0f}s)")
    return verdict_json


from . import report  # noqa: E402  (circular-safe late import)
