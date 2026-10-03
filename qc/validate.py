"""Self-verification (brief section 6): nested LOIO, confound checks, sensitivity, synthetic truth."""
from __future__ import annotations

import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import assign as A
from . import stats as S
from .analyze import analyze_sample
from .gates import SHADOW_GATES
from .io import discover
from .kpis.registry import KPI_NAMES, SHORTLIST, rclass

R_KPIS = [k for k in KPI_NAMES if rclass(k) == "R" and k not in ("st_orientation_deg",)]


# ------------------------------------------------------------------ LOIO

def nested_loio(df: pd.DataFrame, strips: pd.DataFrame | None, cfg: dict) -> dict:
    lab = df[df["batch"].isin(["1", "2", "3"])].reset_index(drop=True)
    rows = []
    for i, r in lab.iterrows():
        train = lab.drop(index=i).reset_index(drop=True)
        st = S.baseline_stats(train, cfg)
        if st["n"] == 0:
            continue
        X = A.feature_matrix(train, st, cfg)
        y = train["batch"].to_numpy()
        params = A.choose_params(X, y, cfg)  # inner LOIO on the remaining images
        C = A.fit_centroids(X, y, params["delta"])
        model = dict(centroids={c: v.tolist() for c, v in C.items()}, T=params["T"], delta=params["delta"],
                     within_distances=A.within_batch_distances(X, y, params["delta"]),
                     signatures=A.batch_signatures(train, st, cfg))
        x = r.to_dict()
        st_rows = strips[strips["sample_id"] == r["sample_id"]].to_dict("records") if strips is not None else []
        a = A.assign(x, st_rows, model, st, cfg)
        s_loo = list(S.loo_baseline_scores(train, cfg).values())
        v = S.verdict(x, st, cfg, s_loo)
        p = {b: a["probabilities"].get(b, 0.0) for b in ("1", "2", "3")}
        top3 = [d["kpi"] for d in v["drivers"][:3]]
        rows.append(dict(sample_id=r["sample_id"], true=r["batch"], pred=a["assigned_batch"], p1=p["1"], p2=p["2"], p3=p["3"],
                         stability=a["stability"], strip_agreement=a["strip_agreement"],
                         strip_vs_truth=float(np.mean([c == r["batch"] for c in a["strip_calls"]])) if a["strip_calls"] else np.nan,
                         verdict=v["decision"], S=v["aggregate_score_S"], novelty=a["novelty_flag"], delta=params["delta"], T=params["T"],
                         top3=", ".join(top3), classes_trained="".join(sorted(C))))
    t = pd.DataFrame(rows)
    if t.empty:
        return dict(table=t)
    cm = pd.crosstab(t["true"], t["pred"], rownames=["true"], colnames=["pred"], dropna=False)
    nll = float(np.mean([-np.log(max(r[f"p{r['true']}"], A.EPS)) for _, r in t.iterrows()]))
    brier = float(np.mean([sum((r[f"p{b}"] - (r["true"] == b)) ** 2 for b in ("1", "2", "3")) for _, r in t.iterrows()]))
    base = t[t["true"] == "3"]
    other = t[t["true"] != "3"]
    return dict(
        table=t, confusion=cm, accuracy=float((t["true"] == t["pred"]).mean()), nll=nll, brier=brier,
        false_alarms=int((base["verdict"] != "ACCEPT").sum()), n_baseline=len(base),
        detections=int((other["verdict"] != "ACCEPT").sum()), n_other=len(other),
        strip_consistency=float(t["strip_agreement"].mean()),
    )


def shadow_classifier(df: pd.DataFrame, cfg: dict) -> dict:
    """Same nearest-centroid LOIO, but on acquisition gates only (confound check)."""
    lab = df[df["batch"].isin(["1", "2", "3"])].reset_index(drop=True)
    keys = [k for k in SHADOW_GATES if k in lab and np.isfinite(lab[k].to_numpy(float)).any()]
    preds = []
    for i in range(len(lab)):
        train = lab.drop(index=i)
        base = train[train["batch"] == "3"]
        if len(base) == 0:
            preds.append(None)
            continue
        med = base[keys].median()
        sc = np.maximum(1.4826 * (base[keys] - med).abs().median(), cfg["scale_rel_floor"] * med.abs()).replace(0, 1e-9)
        Xt = ((train[keys] - med) / sc).fillna(0).to_numpy()
        yt = train["batch"].to_numpy()
        C = A.fit_centroids(Xt, yt, 0.0)
        xi = ((lab.loc[i, keys].astype(float) - med) / sc).fillna(0).to_numpy()
        d = A.distances(xi, C)
        preds.append(min(d, key=d.get))
    y = lab["batch"].to_numpy()
    ok = [p == t for p, t in zip(preds, y) if p is not None]
    return dict(accuracy=float(np.mean(ok)) if ok else np.nan, n=len(ok), keys=keys, preds=preds)


def effect_sizes(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    lab = df[df["batch"].isin(["1", "2", "3"])]
    rows = []
    for a, b in (("1", "3"), ("2", "3"), ("1", "2")):
        A_ = lab[lab["batch"] == a]
        B_ = lab[lab["batch"] == b]
        if not len(A_) or not len(B_):
            continue
        for k in SHORTLIST:
            sa = S.robust_stats(A_[k].to_numpy(float), A_.get(f"{k}_se", pd.Series(dtype=float)).to_numpy(float), cfg["scale_rel_floor"])
            sb = S.robust_stats(B_[k].to_numpy(float), B_.get(f"{k}_se", pd.Series(dtype=float)).to_numpy(float), cfg["scale_rel_floor"])
            pooled = np.sqrt((sa["scale"] ** 2 + sb["scale"] ** 2) / 2)
            rows.append(dict(pair=f"{a} vs {b}", kpi=k, class_=rclass(k), median_a=sa["med"], median_b=sb["med"],
                             effect=abs(sa["med"] - sb["med"]) / pooled if pooled > 0 else np.nan))
    t = pd.DataFrame(rows)
    return t.sort_values(["pair", "effect"], ascending=[True, False]) if len(t) else t


# ------------------------------------------------------------------ sensitivity (re-segmentation)

def sensitivity(samples, df: pd.DataFrame, stats: dict, cfg: dict, log=print) -> pd.DataFrame:
    per_batch = cfg.get("validation_sensitivity_samples_per_batch", 1)
    chosen = []
    for b in ("1", "2", "3"):
        chosen += [s for s in samples if s.batch == b][:per_batch]
    rows = []
    for s in chosen:
        ref = df[df["sample_id"] == s.sample_id].iloc[0]
        variants = {"t1-5": (-5, 0), "t1+5": (5, 0), "t2-5": (0, -5), "t2+5": (0, 5)}
        for name, sh in variants.items():
            t = time.time()
            r = analyze_sample(s, cfg, None, crack_threshold=stats.get("crack_threshold"), t_shift=sh, save_cache=False)
            rows += _deltas(s.sample_id, name, r.kpis, ref, stats)
            log(f"  sensitivity {s.sample_id} {name} {time.time() - t:.0f}s")
        t = time.time()
        r = analyze_sample(s, cfg, None, crack_threshold=stats.get("crack_threshold"), downsample=2, save_cache=False)
        rows += _deltas(s.sample_id, "50nm/px (2x down)", r.kpis, ref, stats)
        log(f"  resolution {s.sample_id} {time.time() - t:.0f}s")
    return pd.DataFrame(rows)


def _deltas(sid, variant, kp, ref, stats):
    out = []
    for k in SHORTLIST + [k for k in R_KPIS if k not in SHORTLIST]:
        sc = stats["kpi"].get(k, {}).get("scale", np.nan)
        a, b = kp.get(k, np.nan), ref.get(k, np.nan)
        out.append(dict(sample_id=sid, variant=variant, kpi=k, class_=rclass(k),
                        delta_scale=float((a - b) / sc) if np.isfinite(sc) and sc > 0 else np.nan))
    return out


# ------------------------------------------------------------------ synthetic truth

def light_kpis(imgs: dict, cfg: dict) -> dict:
    """Segmentation + fractions + particles + dispersion + two-point only (no files, no Sato)."""
    from .gates import noise_sd
    from .kpis import correlation, dispersion, fractions, sizes
    from .kpis.context import Ctx
    from .segment import segment

    seg = segment(imgs["BSE"], imgs.get("ETD"), cfg)
    ctx = Ctx(px_um=0.025, bse=imgs["BSE"], bse_s=seg.bse_s, pore=seg.pore, si=seg.si, carbon=seg.carbon, ambiguous=seg.ambiguous,
              si_labels=seg.si_labels, t1=seg.t1, t2=seg.t2, noise_sd=noise_sd(imgs["BSE"], seg.carbon), cfg=cfg)
    out = fractions.fractions(ctx)
    k, lists = sizes.si_particles(ctx)
    out.update(k)
    out.update(dispersion.dispersion(ctx, lists.get("si_centroids_px")))
    out.update(correlation.two_point(ctx))
    return out


def synthetic_checks(cfg: dict, width: int = 3500, height: int = 1000, n_rep: int = 4, log=print) -> dict:
    from .synth import SynthParams, make_sample, perturb, truth_values, write_triple
    from .io import SampleFiles

    out = {"recovery": [], "clark_evans": {}, "se_check": {}, "perturb": []}
    tmp = Path(tempfile.mkdtemp(prefix="qc_synth_"))

    def run(imgs, sid):
        write_triple(tmp / sid, sid, imgs)
        f = discover(tmp / sid)[0]
        return analyze_sample(f, cfg, None, crack_threshold=None, save_cache=False)

    reps = []
    for i in range(n_rep):
        p = SynthParams(width=width, height=height, seed=500 + i)
        imgs, masks = make_sample(p)
        tv = truth_values(masks)
        r = run(imgs, f"rep{i}")
        reps.append(r)
        for k, tol, kind in (("si_frac_total", 0.005, "abs"), ("pore_frac_deep", 0.01, "abs"),
                             ("si_d50_aw_um", 0.10, "rel"), ("aniso_pore_chord_ratio", 0.10, "rel")):
            got = r.kpis[k]
            err = abs(got - tv[k]) if kind == "abs" else abs(got - tv[k]) / tv[k]
            out["recovery"].append(dict(rep=i, kpi=k, truth=tv[k], measured=got, error=err, tol=tol, kind=kind, ok=bool(err <= tol)))
        log(f"  synthetic rep {i} done")
    # SE check: 10 independent images with a Poisson number of particles (fill-to-target replicates
    # have almost no fraction variance by construction); fractions + two-point only, for speed.
    fr, se = [], []
    for i in range(10):
        imgs, _ = make_sample(SynthParams(width=width, height=height, seed=700 + i, poisson_count=True))
        k = light_kpis(imgs, cfg)
        fr.append(k["si_frac_total"])
        se.append(k["si_frac_se_imagerep"])
    fr, se = np.array(fr), np.array(se)
    out["se_check"] = dict(empirical_sd=float(fr.std(ddof=1)), mean_imagerep_se=float(se.mean()),
                           ratio=float(se.mean() / max(fr.std(ddof=1), 1e-12)))
    out["se_check"]["ok"] = bool(0.5 <= out["se_check"]["ratio"] <= 2.0)
    # Clark-Evans needs particles small relative to their spacing: at 7 % loading with 4 um discs the
    # hard-core exclusion alone pins R near 1 whatever the placement, so this check uses sparse 1 um
    # discs on full-size frames; the Poisson value is a mean over 3 seeds (SE of R ~ 0.05 per image).
    pois = []
    for sd in (901, 902, 903):
        imgs, _ = make_sample(SynthParams(width=7000, height=2000, seed=sd, si_d50_um=1.0, si_frac=0.015))
        pois.append(float(light_kpis(imgs, cfg)["si_clark_evans_R"]))
    out["clark_evans"]["poisson"] = float(np.mean(pois))
    out["clark_evans"]["poisson_each"] = pois
    imgs, _ = make_sample(SynthParams(width=7000, height=2000, seed=901, clustered=True, si_d50_um=1.0, si_frac=0.015))
    out["clark_evans"]["clustered"] = float(light_kpis(imgs, cfg)["si_clark_evans_R"])
    out["clark_evans"]["ok"] = bool(abs(out["clark_evans"]["poisson"] - 1) <= 0.1 and out["clark_evans"]["clustered"] < 0.8)
    # perturbation robustness: scale from the replicate set (strip SE / MAD as in the baseline)
    import pandas as _pd

    df = _pd.DataFrame([{"sample_id": f"rep{i}", "batch": "3", **r.kpis, **{f"{k}_se": v for k, v in r.kpi_se.items()}}
                        for i, r in enumerate(reps)])
    st = S.baseline_stats(df, cfg)
    base_imgs, _ = make_sample(SynthParams(width=width, height=height, seed=500))
    ref = reps[0].kpis
    variants = {"brightness+15%": dict(brightness=1.15), "contrast x0.85": dict(contrast=0.85), "gamma 1.2": dict(gamma=1.2),
                "noise SD 5": dict(noise=5)}
    for name, kw in variants.items():
        imgs = {d: perturb(a, **kw) for d, a in base_imgs.items()}
        r = run(imgs, f"pert{name.split()[0].replace('+', '').replace('%', '')}")
        worst = 0.0
        worst_k = None
        failing = []
        for k in R_KPIS:
            sc = st["kpi"].get(k, {}).get("scale", np.nan)
            if np.isfinite(sc) and sc > 0 and np.isfinite(r.kpis.get(k, np.nan)) and np.isfinite(ref.get(k, np.nan)):
                d = abs(r.kpis[k] - ref[k]) / sc
                if d >= 0.25:
                    failing.append(f"{k} {d:.2f}")
                if d > worst:
                    worst, worst_k = d, k
        v0 = S.verdict(ref, st, cfg)["decision"]
        v1 = S.verdict(r.kpis, st, cfg)["decision"]
        out["perturb"].append(dict(variant=name, worst_R_kpi=worst_k, worst_delta_scale=worst, ok=bool(worst < 0.25),
                                   R_kpis_over_0_25="; ".join(failing) or "none",
                                   verdict_ref=v0, verdict_perturbed=v1, verdict_same=v0 == v1))
        log(f"  perturbation {name} done")
    return out


# ------------------------------------------------------------------ writer

def write_validation_md(path: Path, loio: dict, shadow: dict, eff: pd.DataFrame, sens: pd.DataFrame | None, synth: dict | None,
                        counts: dict, stats: dict, model: dict):
    L = ["# VALIDATION", "", f"_Generated by `qc validate` — baseline version {stats.get('version')}; labelled images: "
         + ", ".join(f"batch {b}: {n}" for b, n in sorted(counts.items())) + "._", ""]
    n_base = counts.get("3", 0)
    if sum(counts.values()) < 6 or len(counts) < 3:
        L += ["> **Small or incomplete labelled set.** With these counts the LOIO numbers below are not a measure of accuracy "
              "on an unseen batch; they demonstrate the machinery. Add the remaining Drive images (`qc fetch`) and rerun.", ""]
    L += ["## 1. Leave-one-image-out (nested: Δ and T tuned inside each fold; baseline stats rebuilt without the held-out image)", ""]
    if loio.get("table") is not None and len(loio["table"]):
        L += [f"- Accuracy: **{loio['accuracy']:.2f}** ({(loio['table']['true'] == loio['table']['pred']).sum()}/{len(loio['table'])})",
              f"- Mean negative log-likelihood: {loio['nll']:.3f} (uniform guess = {np.log(3):.3f}); Brier score: {loio['brier']:.3f}",
              f"- **False alarms** (batch-3 images given INVESTIGATE/REJECT when held out): {loio['false_alarms']} / {loio['n_baseline']}",
              f"- **Detections** (batch-1/2 images given INVESTIGATE/REJECT): {loio['detections']} / {loio['n_other']}",
              f"- Strip consistency (share of strips agreeing with the parent image call): {loio['strip_consistency']:.2f} (target > 0.7)",
              f"- p-value floor: with n = {n_base} baseline images the smallest achievable rank p-value is {1 / (n_base + 1):.3f}; "
              "no significance is claimed beyond that.", "", "Confusion matrix (rows = true, columns = predicted):", "",
              loio["confusion"].to_markdown(), "", "Per-image results:", "",
              loio["table"].to_markdown(index=False, floatfmt=".3f"), ""]
    else:
        L += ["Not enough labelled images for LOIO.", ""]
    L += ["## 2. Which KPIs separate the batches (effect = |Δ median| / pooled robust scale)", ""]
    if len(eff):
        for pair, g in eff.groupby("pair"):
            L += [f"**{pair}**", "", g.drop(columns=["pair"]).head(8).to_markdown(index=False, floatfmt=".3f"), ""]
        L += ["![KPI box plots](out/figs/kpi_boxplots.png)", ""]
    else:
        L += ["Needs at least two labelled batches.", ""]
    sigs = model.get("signatures", {})
    L += ["## 3. Batch signatures", ""]
    for b in ("1", "2"):
        if b in sigs:
            L.append(f"- Batch {b} (n={sigs[b]['n_images']}): " + (", ".join(
                f"{k} {'+' if sigs[b]['effects'][k]['sign'] > 0 else '−'}{abs(sigs[b]['effects'][k]['effect']):.1f}σ"
                for k in sigs[b]["signature"]) or "no stable KPI with |effect| ≥ 1.5"))
    if "overlap_1_2" in sigs:
        L.append(f"- Overlap of batch-1 and batch-2 signatures: {sigs['overlap_1_2']:.2f} (1 = identical; high overlap → 1-vs-2 calls are low-confidence by construction)")
    L += ["", "## 4. Confound check — shadow classifier on acquisition gates only", ""]
    L.append(f"Gates used: {', '.join(shadow.get('keys', []))}. LOIO accuracy: **{shadow.get('accuracy', float('nan')):.2f}** (n={shadow.get('n', 0)}).")
    if loio.get("accuracy") is not None and np.isfinite(shadow.get("accuracy", np.nan)) and shadow["accuracy"] >= loio.get("accuracy", 1) - 0.1:
        L += ["", "> **WARNING: acquisition gates alone separate the batches about as well as the material KPIs.** The batches may "
              "differ by imaging session as much as by material; read material conclusions with that in mind."]
    L += [""]
    L += ["## 5. Threshold and resolution sensitivity (Δ in units of baseline scale)", ""]
    if sens is not None and len(sens):
        piv = sens.pivot_table(index=["kpi", "class_"], columns="variant", values="delta_scale", aggfunc=lambda v: np.nanmax(np.abs(v)))
        L += [piv.to_markdown(floatfmt=".2f"), ""]
        thr = sens[sens["variant"].str.startswith("t")]
        bad = thr[(thr["class_"] == "R") & (thr["delta_scale"].abs() > 0.5)]["kpi"].unique()
        L.append(f"R-class KPIs moving > 0.5 scale under ±5-level threshold shifts: {', '.join(bad) if len(bad) else 'none'} (target: none).")
    else:
        L.append("Not run (use `qc validate` without `--skip-sensitivity`).")
    L += ["", "## 6. Synthetic ground truth", ""]
    if synth:
        rec = pd.DataFrame(synth["recovery"])
        L += [rec.to_markdown(index=False, floatfmt=".4f"), "",
              f"- Clark–Evans R: Poisson placement {synth['clark_evans']['poisson']:.2f} (target 1 ± 0.1), clustered {synth['clark_evans']['clustered']:.2f} (target < 0.8) → {'PASS' if synth['clark_evans']['ok'] else 'FAIL'}",
              f"- ImageRep SE of Si fraction: mean {synth['se_check']['mean_imagerep_se']:.4f} vs empirical SD over replicates {synth['se_check']['empirical_sd']:.4f} (ratio {synth['se_check']['ratio']:.2f}; target within ×2) → {'PASS' if synth['se_check']['ok'] else 'FAIL'}",
              "", "Perturbation robustness (R-class KPIs, Δ in replicate-baseline scale units; target < 0.25, verdict unchanged):", "",
              pd.DataFrame(synth["perturb"]).to_markdown(index=False, floatfmt=".3f"), ""]
    else:
        L.append("Not run (use `qc validate` without `--skip-synthetic`).")
    path.write_text("\n".join(L) + "\n")
