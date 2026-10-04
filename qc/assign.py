"""Batch assignment: nearest shrunken centroid in baseline-sigma units, softmax with temperature,
uniform priors, strip stability, novelty score, and batch signatures (brief 5.3-5.4)."""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

from .kpis.registry import SHORTLIST, weight

EPS = 1e-3


def feature_vector(x: dict, stats: dict, cfg: dict, keys=SHORTLIST) -> np.ndarray:
    v = []
    for k in keys:
        s = stats["kpi"].get(k)
        val = x.get(k, np.nan)
        if s is None or val is None or not np.isfinite(val) or not np.isfinite(s["med"]):
            v.append(0.0)  # missing -> baseline median (no evidence either way)
        else:
            v.append((val - s["med"]) / s["scale"] * weight(k, cfg))
    return np.array(v, float)


def feature_matrix(df: pd.DataFrame, stats: dict, cfg: dict, keys=SHORTLIST) -> np.ndarray:
    return np.array([feature_vector(r, stats, cfg, keys) for r in df.to_dict("records")]).reshape(len(df), len(keys))


def fit_centroids(X: np.ndarray, y: np.ndarray, delta: float) -> dict:
    classes = sorted(set(y))
    raw = {c: np.median(X[y == c], axis=0) for c in classes}
    # grand centre = median of class centroids, so class sizes do not pull it (uniform priors)
    g = np.median(np.array(list(raw.values())), axis=0)
    out = {}
    for c, m in raw.items():
        d = m - g
        out[c] = g + np.sign(d) * np.maximum(np.abs(d) - delta, 0)
    return out


def distances(x: np.ndarray, C: dict) -> dict:
    return {c: float(np.sqrt(((x - m) ** 2).sum())) for c, m in C.items()}


def softmax_probs(d: dict, T: float) -> dict:
    keys = list(d)
    a = np.array([-(d[k] ** 2) / (2 * T) for k in keys])
    a -= a.max()
    p = np.exp(a)
    p /= p.sum()
    return {k: float(v) for k, v in zip(keys, p)}


def loio_predictions(X, y, delta, T, ids=None):
    """Leave-one-image-out probabilities with centroids refit without each image."""
    out = []
    for i in range(len(y)):
        mask = np.arange(len(y)) != i
        if len(set(y[mask])) == 0:
            out.append(None)
            continue
        C = fit_centroids(X[mask], y[mask], delta)
        d = distances(X[i], C)
        out.append((softmax_probs(d, T), d))
    return out


def nll(preds, y) -> float:
    vals = []
    for p, yi in zip(preds, y):
        if p is None:
            continue
        vals.append(-np.log(max(p[0].get(yi, 0.0), EPS)))
    return float(np.mean(vals)) if vals else np.nan


def choose_params(X, y, cfg) -> dict:
    """Δ and T by LOIO negative log-likelihood on the given labelled set."""
    best = None
    table = []
    if len(y) < 3 or len(set(y)) < 2:
        return dict(delta=float(cfg["shrinkage_grid"][0]), T=float(np.median(cfg["temperature_grid"])), table=[],
                    note="too few labelled images/classes to tune; using defaults")
    for delta, T in itertools.product(cfg["shrinkage_grid"], cfg["temperature_grid"]):
        v = nll(loio_predictions(X, y, delta, T), y)
        table.append(dict(delta=delta, T=T, nll=v))
        if best is None or v < best[2] - 1e-12:
            best = (delta, T, v)
    return dict(delta=float(best[0]), T=float(best[1]), nll=best[2], table=table)


def within_batch_distances(X, y, delta) -> list[float]:
    """LOIO distance of each labelled image to its own batch centroid (fit without it)."""
    out = []
    for i in range(len(y)):
        mask = np.arange(len(y)) != i
        if not np.any(y[mask] == y[i]):
            continue
        C = fit_centroids(X[mask], y[mask], delta)
        out.append(distances(X[i], C)[y[i]])
    return sorted(out)


def id_score(dmin: float, within: list[float]) -> float:
    if not within:
        return np.nan
    w = np.asarray(within)
    return float(1.0 - (np.sum(w <= dmin) / (len(w) + 1)))


def confidence_label(p_top: float, stability: float) -> str:
    if p_top >= 0.75 and stability >= 0.8:
        return "confident"
    if p_top >= 0.5:
        return "moderately confident"
    return "low confidence"


def assign(x: dict, strips: list[dict], model: dict, stats: dict, cfg: dict) -> dict:
    C = {c: np.asarray(v) for c, v in model["centroids"].items()}
    T = model["T"]
    xv = feature_vector(x, stats, cfg)
    d = distances(xv, C)
    p = softmax_probs(d, T)
    top = max(p, key=p.get)
    # stability: single strips and strip-bootstrap means
    calls = []
    sv = np.array([feature_vector(s, stats, cfg) for s in strips]) if strips else np.zeros((0, len(xv)))
    # Strip-level KPIs are biased relative to whole-image KPIs (narrower field: shorter correlation
    # lengths, particles cut at strip edges), so strip vectors are re-centred on the image vector:
    # only their *spread* is used for stability (see DECISIONS.md).
    if len(sv):
        sv = xv + (sv - sv.mean(axis=0))
    for v in sv:
        pp = softmax_probs(distances(v, C), T)
        calls.append(max(pp, key=pp.get))
    strip_agree = float(np.mean([c == top for c in calls])) if calls else np.nan
    rng = np.random.default_rng(cfg["seed"])
    boot_calls = []
    boot_pts = []
    if len(sv):
        idx = rng.integers(0, len(sv), size=(cfg["bootstrap_n"], len(sv)))
        means = sv[idx].mean(axis=1)
        for v in means:
            dd = distances(v, C)
            boot_calls.append(min(dd, key=dd.get))
        boot_pts = means
    stability = float(np.mean([c == top for c in boot_calls])) if boot_calls else np.nan
    dmin = min(d.values())
    ids = id_score(dmin, model.get("within_distances", []))
    within = model.get("within_distances", [])
    # With < 40 reference distances the empirical 97.5th percentile is just the maximum, which flags
    # ordinary members of a spread-out batch; then require 1.25x the largest distance ever seen.
    if not within:
        novelty = False
    elif len(within) >= 40:
        novelty = dmin > np.percentile(within, 97.5)
    else:
        novelty = dmin > 1.25 * max(within)
    sig = signature_match(x, stats, model.get("signatures", {}), cfg)
    beyond_range = bool(novelty)
    # far away, but in the assigned batch's own signature direction = a more extreme member, not a new way
    sm_top = sig["match"].get(top, np.nan)
    if novelty and np.isfinite(sm_top) and sm_top >= cfg["signature_min_effect"]:
        novelty = False
    if sig.get("new_way"):
        novelty = True
    label = confidence_label(p[top], stability if np.isfinite(stability) else 0.0)
    return dict(
        assigned_batch=top, probabilities=p, distances=d, stability=stability, strip_agreement=strip_agree,
        strip_calls=calls, confidence_label=label, novelty_flag=bool(novelty), in_distribution_score=ids,
        signature_match=sig["match"], signature_kpis_matched=sig["matched"], signature_kpis_missed=sig["missed"],
        signature_new_way=sig.get("new_way", False), beyond_labelled_range=beyond_range, priors="uniform", feature_vector=xv.tolist(),
        bootstrap_points=np.asarray(boot_pts).tolist() if len(boot_pts) else [],
        classes_available=sorted(C), T=T, delta=model["delta"],
    )


# ---------------------------------------------------------------- signatures

def batch_signatures(df: pd.DataFrame, stats: dict, cfg: dict) -> dict:
    """effect_bk = (median_b - baseline median)/scale for each non-baseline batch, with LOO sign stability."""
    from .stats import baseline_stats

    out = {}
    lab = df[df["batch"].astype(str).isin(["1", "2", "3"])]
    loo_stats = {sid: baseline_stats(lab, cfg, exclude=sid) for sid in lab["sample_id"]}
    for b in sorted(set(lab["batch"].astype(str)) - {"3"}):
        sub = lab[lab["batch"].astype(str) == b]
        effects = {}
        for k in SHORTLIST:
            s = stats["kpi"].get(k)
            if s is None or not np.isfinite(s["med"]):
                continue
            vals = sub[k].to_numpy(float)
            vals = vals[np.isfinite(vals)]
            if not len(vals):
                continue
            e = (np.median(vals) - s["med"]) / s["scale"]
            # LOO over all labelled images: drop each one, recompute baseline stats and batch median
            agree = []
            for sid in lab["sample_id"]:
                st2 = loo_stats[sid]
                s2 = st2["kpi"].get(k)
                v2 = sub[sub["sample_id"] != sid][k].to_numpy(float)
                v2 = v2[np.isfinite(v2)]
                if s2 is None or not len(v2) or not np.isfinite(s2["med"]):
                    agree.append(False)
                    continue
                e2 = (np.median(v2) - s2["med"]) / s2["scale"]
                agree.append(bool(np.sign(e2) == np.sign(e) and abs(e2) > 1))
            stab = float(np.mean(agree)) if agree else 0.0
            effects[k] = dict(effect=float(e), sign=int(np.sign(e)), stability=stab, weight=weight(k, cfg),
                              in_signature=bool(abs(e) >= cfg["signature_min_effect"] and stab >= cfg["signature_min_stability"]))
        out[b] = dict(n_images=int(len(sub)), effects=effects,
                      signature=[k for k, v in sorted(effects.items(), key=lambda kv: -abs(kv[1]["effect"])) if v["in_signature"]])
    # overlap between signatures
    if "1" in out and "2" in out:
        s1 = {(k, out["1"]["effects"][k]["sign"]) for k in out["1"]["signature"]}
        s2 = {(k, out["2"]["effects"][k]["sign"]) for k in out["2"]["signature"]}
        u = s1 | s2
        out["overlap_1_2"] = float(len(s1 & s2) / len(u)) if u else 0.0
    return out


def signature_match(x: dict, stats: dict, sigs: dict, cfg: dict) -> dict:
    from .stats import zscores

    z = zscores(x, stats)
    match, matched, missed = {}, {}, {}
    sig_keys = set()
    for b, s in sigs.items():
        if not isinstance(s, dict) or "signature" not in s:
            continue
        ks = s["signature"]
        sig_keys |= set(ks)
        if not ks:
            match[b] = np.nan
            matched[b], missed[b] = [], []
            continue
        num = den = 0.0
        matched[b], missed[b] = [], []
        for k in ks:
            e = s["effects"][k]
            zk = z.get(k, np.nan)
            if not np.isfinite(zk):
                continue
            num += e["weight"] * e["sign"] * zk
            den += e["weight"]
            (matched[b] if e["sign"] * zk > 1 else missed[b]).append(k)
        match[b] = float(num / den) if den else np.nan
    # "different in a new way": large |z| on shortlist KPIs that belong to no signature, and no signature matched
    off = [k for k in SHORTLIST if k not in sig_keys and np.isfinite(z.get(k, np.nan)) and abs(z[k]) > cfg["zones"]["reject"]]
    best = max([v for v in match.values() if np.isfinite(v)], default=np.nan)
    new_way = bool(off) and not (np.isfinite(best) and best > 1.5)
    return dict(match=match, matched=matched, missed=missed, new_way=new_way, off_signature_kpis=off)


def build_model(df: pd.DataFrame, stats: dict, cfg: dict) -> dict:
    lab = df[df["batch"].astype(str).isin(["1", "2", "3"])].reset_index(drop=True)
    X = feature_matrix(lab, stats, cfg)
    y = lab["batch"].astype(str).to_numpy()
    params = choose_params(X, y, cfg)
    C = fit_centroids(X, y, params["delta"]) if len(y) else {}
    within = within_batch_distances(X, y, params["delta"])
    sigs = batch_signatures(lab, stats, cfg)
    loio = loio_predictions(X, y, params["delta"], params["T"])
    acc = np.mean([p is not None and max(p[0], key=p[0].get) == yi for p, yi in zip(loio, y)]) if len(y) else np.nan
    return dict(centroids={c: v.tolist() for c, v in C.items()}, delta=params["delta"], T=params["T"],
                tuning=params, within_distances=within, signatures=sigs, n_labelled=int(len(y)),
                counts={c: int((y == c).sum()) for c in sorted(set(y))}, loio_accuracy_fixed_baseline=float(acc) if np.isfinite(acc) else None,
                features=SHORTLIST)
