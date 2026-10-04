"""Verdict JSON, Markdown + standalone HTML report, paste-ready DM text, figures."""
from __future__ import annotations

import base64
import html
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .analyze import _clean  # noqa: E402
from .kpis.registry import KPI_META, SHORTLIST, meaning, rclass, unit  # noqa: E402

CAVEATS = [
    "Bright particles are identified by backscatter Z-contrast only (no EDS); they are reported as Si-candidates and may include SiOx or Si-C composite.",
    "Pores appear open (not resin-infiltrated); the deep-pore fraction is a lower bound on porosity, not a porosity measurement.",
    "Graphite, carbon black and binder cannot be separated in these images; the 'carbon matrix' class contains all three.",
    "Coating thickness, surface roughness and collector delamination are not assessable: the fields lie entirely inside the coating.",
    "Anisotropy and transport indices are 2-D section quantities compared like-with-like against the baseline, not 3-D values; the through-plane direction is assumed vertical in the frame.",
    "Pixel size (25 nm/px) is taken from the TIFF resolution tags; vendor metadata is absent, so it is unverified.",
    "With N baseline images the smallest achievable rank p-value is 1/(N+1); verdicts are effect-size judgements against a small baseline.",
    "Image grey levels have been remapped after acquisition (comb histograms); any method relying on raw intensities would be confounded - this system uses BSE phase identity after smoothing and treats ETD/InLens levels as acquisition covariates.",
]
NOT_MEASURABLE = ["coating thickness", "surface roughness", "collector delamination", "binder/carbon-black distribution",
                  "Si vs SiOx identity", "true (total) porosity"]
SHORT_NAMES = {
    "si_frac_solid": "Si-candidate loading", "si_d50_aw_um": "Si-candidate D50", "si_d90_aw_um": "Si-candidate D90",
    "si_num_density_mm2": "Si-candidate number density", "si_solidity_med": "Si-candidate solidity",
    "si_quadrat_cv": "Si-candidate patchiness", "pore_frac_deep": "deep-pore fraction",
    "pore_chord_h_mean_um": "in-plane pore chord", "pore_chord_v_mean_um": "through-plane pore chord",
    "aniso_pore_chord_ratio": "pore anisotropy", "aniso_carbon_chord_ratio": "carbon/flake anisotropy",
    "s2_len_pore_h_px": "pore correlation length", "pore_percolating_frac_v": "through-plane pore percolation",
    "si_contact_pore_frac": "Si-candidate/pore contact", "si_contact_carbon_frac": "Si-candidate/carbon contact",
    "crack_density_um_per_mm2": "crack density",
}
BATCH_COLORS = {"1": "#2a6fdb", "2": "#d9822b", "3": "#3b9c5a", "unknown": "#888888"}


def short(k):
    return SHORT_NAMES.get(k, k)


def fmt(v, k=None):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "n/a"
    u = unit(k) if k else ""
    if u == "fraction":
        return f"{100 * v:.2f} %"
    a = abs(v)
    if a >= 1000:
        return f"{v:,.0f}"
    if a >= 10:
        return f"{v:.1f}"
    if a >= 1:
        return f"{v:.2f}"
    return f"{v:.3f}"


def fmt_u(v, k):
    u = unit(k)
    s = fmt(v, k)
    return s if u in ("fraction", "ratio", "") or s == "n/a" else f"{s} {u.replace('um', 'µm')}"


# ------------------------------------------------------------------ verdict JSON

def build_verdict_json(res, st, model, v, g, a, cfg, held_out_note, elapsed) -> dict:
    kp = {}
    z = v["z_all"]
    for k in KPI_META:
        if k in res.kpis:
            kp[k] = dict(value=res.kpis[k], se=res.kpi_se.get(k), ci95=res.kpi_ci.get(k), z=z.get(k), unit=unit(k),
                         robustness=rclass(k), shortlist=k in SHORTLIST)
    n_lab = model.get("n_labelled", 0)
    acc = model.get("loio_accuracy_fixed_baseline")
    calib = (f"T={model['T']} and shrinkage Δ={model['delta']} chosen by leave-one-image-out log-loss on {n_lab} labelled images"
             + (f"; LOIO accuracy {acc:.2f}" if acc is not None and np.isfinite(acc) else ""))
    missing_classes = sorted({"1", "2", "3"} - set(a["classes_available"]))
    out = {
        "sample_id": res.sample_id,
        "true_batch_if_known": res.batch if res.batch in ("1", "2", "3") else None,
        "baseline": {"batch": 3, "n_images": st.get("n"), "version": st.get("version")},
        "acquisition_gates": {"status": g["status"], "flags": g["flags"], "z": g["z"], "values": res.gates},
        "verdict_vs_baseline": {k: v[k] for k in ("decision", "reasons", "aggregate_score_S", "conformal_rank_p", "p_floor",
                                                   "n_kpis_beyond_2_5", "defect_alarms")} | {"drivers": v["drivers"]},
        "batch_assignment": {k: a[k] for k in ("assigned_batch", "probabilities", "distances", "stability", "strip_agreement",
                                               "confidence_label", "novelty_flag", "in_distribution_score", "signature_match",
                                               "signature_kpis_matched", "signature_kpis_missed", "beyond_labelled_range", "priors")}
        | {"calibration_note": calib, "classes_without_training_images": missing_classes},
        "kpis": kp,
        "image": {"shape_px": res.info["shape"], "px_nm": res.info["px_nm"], "flags": res.info["flags"],
                  "crop": res.info["crop"], "etd_is_se": res.info.get("etd_is_se")},
        "not_measurable": NOT_MEASURABLE,
        "caveats": CAVEATS,
        "runtime_s": round(elapsed, 1),
        "notes": [n for n in [held_out_note] if n],
        "_assign_internal": {"bootstrap_points": a["bootstrap_points"], "feature_vector": a["feature_vector"],
                             "strip_calls": a["strip_calls"]},
    }
    return _clean(out)


# ------------------------------------------------------------------ text

def signature_phrase(b: str, sigs: dict, n=3) -> str:
    s = sigs.get(b) or {}
    ks = s.get("signature", [])[:n]
    if not ks:
        return "no stable signature"
    arrows = [f"{short(k)} {'↑' if s['effects'][k]['sign'] > 0 else '↓'}" for k in ks]
    return ", ".join(arrows)


def lead_paragraph(vj: dict, model: dict) -> str:
    a = vj["batch_assignment"]
    v = vj["verdict_vs_baseline"]
    b = a["assigned_batch"]
    p = a["probabilities"]
    others = "; ".join(f"Batch {k}: {p[k]:.2f}" for k in sorted(p) if k != b) or "no other batch trained"
    stab = a.get("stability")
    stab_s = f"assignment stable in {100 * stab:.0f} % of strip-bootstrap resamples" if stab is not None else "stability n/a"
    s = f"**{vj['sample_id']} is most consistent with Batch {b}** (probability {p[b]:.2f}; {others}; {stab_s}; {a['confidence_label']})."
    sigs = model.get("signatures", {})
    parts = []
    for bb in sorted(a.get("signature_match", {})):
        m = a["signature_match"][bb]
        if m is None:
            continue
        parts.append(f"{m:+.1f} σ along Batch {bb}'s signature ({signature_phrase(bb, sigs)})")
    if parts:
        s += " It sits " + " and ".join(parts) + "."
    if a.get("novelty_flag"):
        s += (f" **Novelty: it does not resemble any known batch closely** (in-distribution score {a['in_distribution_score']:.2f}); "
              f"betting on Batch {b} as the nearest, treat as INVESTIGATE.")
    elif a.get("beyond_labelled_range"):
        s += (f" It lies further from the Batch {b} centroid than any labelled image, but in Batch {b}'s own signature "
              f"direction: a more extreme Batch-{b}-like sample rather than a new kind of variation.")
    if a.get("classes_without_training_images"):
        s += f" (No labelled images available for batch(es) {', '.join(a['classes_without_training_images'])}; they cannot be assigned.)"
    dec = v["decision"]
    n25 = v["n_kpis_beyond_2_5"]
    head = "within the baseline distribution" if dec == "ACCEPT" else "different from the baseline"
    s += f" **Relative to the Batch-3 baseline it is {head} — verdict {dec}**: {n25} of {len(SHORTLIST)} KPIs exceed 2.5 robust σ."
    drv = v["drivers"][:3]
    if drv:
        s += " Largest shifts: " + "; ".join(
            f"{short(d['kpi'])} {fmt_u(d['value'], d['kpi'])} vs {fmt_u(d['baseline_median'], d['kpi'])} ± {fmt(d['baseline_scale'], d['kpi'])} (z = {d['z']:+.1f})"
            for d in drv) + "."
        big = [d for d in drv if abs(d["z"]) > 2]
        if big:
            s += " Reads as: " + "; ".join(d["meaning"] for d in big[:2]) + "."
    nb = vj["baseline"]["n_images"] or 0
    s += (f" Baseline has {nb} image(s), so the smallest reportable rank p-value is {v['p_floor']:.3f}; "
          "this verdict is an effect-size judgement, not a significance test.")
    return s


def dm_block(vj: dict, model: dict) -> str:
    a = vj["batch_assignment"]
    v = vj["verdict_vs_baseline"]
    b = a["assigned_batch"]
    p = a["probabilities"]
    others = " | ".join(f"B{k} {p[k]:.2f}" for k in sorted(p) if k != b) or "no other batch trained"
    stab = a.get("stability")
    l1 = f"{vj['sample_id']} → Batch {b} (p={p[b]:.2f} | {others}) — {a['confidence_label']}, stability {100 * (stab or 0):.0f}%"
    rp = v.get("conformal_rank_p")
    rp_s = f"rank-p {rp:.3f}" + (" = floor" if rp is not None and abs(rp - v["p_floor"]) < 1e-9 else "") if rp is not None else "rank-p n/a"
    l2 = f"  vs baseline: {v['decision']} ({v['n_kpis_beyond_2_5']}/{len(SHORTLIST)} KPIs > 2.5σ; {rp_s})"
    drv = v["drivers"][:3]
    l3 = "  why: " + "; ".join(f"{short(d['kpi'])} {fmt(d['value'], d['kpi'])} vs {fmt(d['baseline_median'], d['kpi'])}±{fmt(d['baseline_scale'], d['kpi'])} (z {d['z']:+.1f})" for d in drv)
    big = [d for d in drv if abs(d["z"]) > 2]
    l4 = "  reads as: " + ("; ".join(d["meaning"].split(":")[0] for d in big[:2]) if big else "no material shift beyond 2σ")
    if a.get("novelty_flag"):
        l5 = f"  novelty: YES — far from all known batches (id-score {a['in_distribution_score']:.2f}); betting on nearest"
    else:
        ids = a.get("in_distribution_score")
        l5 = f"  novelty: no (in-distribution score {ids:.2f})" if ids is not None else "  novelty: n/a"
        if a.get("beyond_labelled_range"):
            l5 += f"; more extreme than any labelled Batch {b} image, same direction"
    return "\n".join([l1, l2, l3, l4, l5])


# ------------------------------------------------------------------ figures

def fig_zbars(vj: dict, path: Path):
    z = {k: vj["kpis"][k]["z"] for k in SHORTLIST if k in vj["kpis"] and vj["kpis"][k]["z"] is not None}
    ks = sorted(z, key=lambda k: abs(z[k]))
    fig, ax = plt.subplots(figsize=(7, 0.32 * len(ks) + 1.2))
    ax.axvspan(-2.5, 2.5, color="#e8f3ea", zorder=0)
    for lo, hi in ((2.5, 4), (-4, -2.5)):
        ax.axvspan(lo, hi, color="#fff2d6", zorder=0)
    lim = max(5, max((abs(v) for v in z.values()), default=0) + 0.5)
    for lo, hi in ((4, lim), (-lim, -4)):
        ax.axvspan(lo, hi, color="#fbe0de", zorder=0)
    cols = ["#c0392b" if abs(z[k]) > 4 else "#d68910" if abs(z[k]) > 2.5 else "#4a6b8a" for k in ks]
    ax.barh([f"{short(k)} [{rclass(k)}]" for k in ks], [z[k] for k in ks], color=cols)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlim(-lim, lim)
    ax.set_xlabel("robust z vs Batch-3 baseline (baseline σ units)")
    ax.set_title("Shortlist KPIs: shift from baseline")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def fig_pca(df, model, stats, cfg, vj, path: Path):
    from .assign import feature_matrix

    lab = df[df["batch"].isin(["1", "2", "3"])]
    if len(lab) < 2:
        return False
    X = feature_matrix(lab, stats, cfg)
    mu = X.mean(0)
    U, S, Vt = np.linalg.svd(X - mu, full_matrices=False)
    P = Vt[:2].T if Vt.shape[0] >= 2 else np.vstack([Vt[:1], np.zeros_like(Vt[:1])]).T
    fig, ax = plt.subplots(figsize=(6, 5))
    Y = (X - mu) @ P
    for b in ("1", "2", "3"):
        m = lab["batch"].to_numpy() == b
        if m.any():
            ax.scatter(Y[m, 0], Y[m, 1], c=BATCH_COLORS[b], label=f"Batch {b} ({m.sum()})", s=40, alpha=0.8, edgecolor="white")
    for b, c in model["centroids"].items():
        cy = (np.asarray(c) - mu) @ P
        ax.scatter(*cy, marker="X", s=160, c=BATCH_COLORS.get(b, "k"), edgecolor="k")
    bp = np.asarray(vj["_assign_internal"]["bootstrap_points"])
    if len(bp):
        by = (bp[:300] - mu) @ P
        ax.scatter(by[:, 0], by[:, 1], s=4, c="k", alpha=0.15, label="strip bootstrap")
    xv = (np.asarray(vj["_assign_internal"]["feature_vector"]) - mu) @ P
    ax.scatter(*xv, marker="*", s=320, c="#e6194b", edgecolor="k", label=vj["sample_id"], zorder=5)
    ev = S**2 / max((S**2).sum(), 1e-12)
    ax.set_xlabel(f"PC1 ({100 * ev[0]:.0f} %)")
    ax.set_ylabel(f"PC2 ({100 * ev[1]:.0f} %)" if len(ev) > 1 else "PC2")
    ax.set_title("Standardised shortlist KPIs (X = batch centroid)")
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return True


def _pooled_lists(base_out: Path, ids, key):
    vals = []
    for sid in ids:
        p = Path(base_out) / "cache" / sid / "lists.npz"
        if p.exists():
            with np.load(p) as z:
                if key in z.files:
                    vals.append(z[key])
    return np.concatenate(vals) if vals else np.array([])


def fig_distributions(res, stats, base_out: Path, path: Path):
    ids = [i for i in stats.get("ids", []) if i != res.sample_id]
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.4))
    e = res.lists.get("si_ecd_um", np.array([]))
    eb = _pooled_lists(base_out, ids, "si_ecd_um")
    bins = np.linspace(0, max(15, np.percentile(np.concatenate([e, eb, [1]]), 99)), 30)
    if len(eb):
        axs[0].hist(eb, bins=bins, weights=eb**2, density=True, alpha=0.5, color=BATCH_COLORS["3"], label="baseline (pooled)")
    if len(e):
        axs[0].hist(e, bins=bins, weights=e**2, density=True, histtype="step", lw=2, color="#e6194b", label="this sample")
    axs[0].set_xlabel("Si-candidate section ECD (µm), area-weighted")
    axs[0].legend(fontsize=8)
    for ax, key, lab in ((axs[1], "pore_chord_h_um", "horizontal (in-plane)"), (axs[2], "pore_chord_v_um", "vertical (through-plane)")):
        c = res.lists.get(key, np.array([]))
        cb = _pooled_lists(base_out, ids, key)
        bins = np.logspace(np.log10(0.025), np.log10(20), 40)
        if len(cb):
            ax.hist(cb, bins=bins, density=True, alpha=0.5, color=BATCH_COLORS["3"], label="baseline (pooled)")
        if len(c):
            ax.hist(c, bins=bins, density=True, histtype="step", lw=2, color="#e6194b", label="this sample")
        ax.set_xscale("log")
        ax.set_xlabel(f"pore chord length, {lab} (µm)")
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def fig_signatures(sigs: dict, path: Path):
    bs = [b for b in ("1", "2") if b in sigs]
    fig, axs = plt.subplots(1, max(len(bs), 1), figsize=(6 * max(len(bs), 1), 5), squeeze=False)
    for ax, b in zip(axs[0], bs):
        eff = sigs[b]["effects"]
        ks = sorted(eff, key=lambda k: abs(eff[k]["effect"]))
        vals = [eff[k]["effect"] for k in ks]
        bars = ax.barh([short(k) for k in ks], vals, color=BATCH_COLORS[b])
        for bar, k in zip(bars, ks):
            if not eff[k]["in_signature"]:
                bar.set_hatch("///")
                bar.set_alpha(0.45)
        ax.axvline(0, color="k", lw=0.8)
        for x in (-1.5, 1.5):
            ax.axvline(x, color="grey", ls=":", lw=0.8)
        ax.set_title(f"Batch {b} vs baseline (n={sigs[b]['n_images']}); filled = stable signature KPI")
        ax.set_xlabel("effect = (batch median − baseline median) / baseline scale")
    if not bs:
        axs[0][0].text(0.5, 0.5, "no labelled batch-1/2 images yet", ha="center", va="center")
        axs[0][0].axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def fig_kpi_boxplots(df, path: Path):
    ks = [k for k in SHORTLIST if k in df]
    n = len(ks)
    cols = 4
    rows = int(np.ceil(n / cols))
    fig, axs = plt.subplots(rows, cols, figsize=(14, 2.6 * rows), squeeze=False)
    for ax, k in zip(axs.ravel(), ks):
        data, labels = [], []
        for b in ("1", "2", "3"):
            v = df[df["batch"] == b][k].to_numpy(float)
            v = v[np.isfinite(v)]
            if len(v):
                data.append(v)
                labels.append(f"B{b}")
        if data:
            bp = ax.boxplot(data, tick_labels=labels, patch_artist=True, widths=0.6)
            for patch, lab in zip(bp["boxes"], labels):
                patch.set_facecolor(BATCH_COLORS[lab[1:]])
                patch.set_alpha(0.5)
            for i, v in enumerate(data):
                ax.scatter(np.full(len(v), i + 1) + np.linspace(-0.12, 0.12, len(v)), v, s=10, c="k", zorder=3)
        ax.set_title(f"{short(k)} [{unit(k)}]", fontsize=9)
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=90)
    plt.close(fig)


# ------------------------------------------------------------------ writers

def kpi_table_md(vj: dict, st: dict) -> str:
    lines = ["| KPI | value ± SE | baseline median ± scale | z | class | shortlist |", "|---|---|---|---|---|---|"]
    order = sorted(vj["kpis"], key=lambda k: (k not in SHORTLIST, list(KPI_META).index(k)))
    for k in order:
        d = vj["kpis"][k]
        s = st["kpi"].get(k, {})
        z = d.get("z")
        lines.append(f"| `{k}` | {fmt_u(d['value'], k)} ± {fmt(d['se'], k)} | {fmt(s.get('med'), k)} ± {fmt(s.get('scale'), k)} | "
                     f"{'n/a' if z is None else f'{z:+.2f}'} | {d['robustness']} | {'✓' if d['shortlist'] else ''} |")
    return "\n".join(lines)


def build_markdown(vj, st, model, figs: dict) -> str:
    a = vj["batch_assignment"]
    v = vj["verdict_vs_baseline"]
    g = vj["acquisition_gates"]
    md = [f"# QC report — {vj['sample_id']}", ""]
    if g["status"] == "caution":
        md += ["> **Acquisition drift or artefact suspected — compare with caution.** " + "; ".join(g["flags"]), ""]
    md += ["## 1. Verdict", "", lead_paragraph(vj, model), ""]
    md += ["## 2. Acquisition gates", "", f"Status: **{g['status']}**" + (f" — flags: {', '.join(g['flags'])}" if g["flags"] else ""), ""]
    zs = {k: val for k, val in (g.get("z") or {}).items() if val is not None}
    if zs:
        md += ["| gate | value | z |", "|---|---|---|"]
        for k in sorted(zs, key=lambda k: -abs(zs[k]))[:10]:
            md.append(f"| `{k}` | {fmt(g['values'].get(k))} | {zs[k]:+.1f} |")
        md.append("")
    md += ["## 3. KPIs", "", kpi_table_md(vj, st), ""]
    md += ["## 4. Drivers (shortlist, sorted by |z|)", ""]
    for d in v["drivers"][:8]:
        md.append(f"- **{short(d['kpi'])}** (`{d['kpi']}`, class {d['robustness']}): {fmt_u(d['value'], d['kpi'])} vs "
                  f"{fmt_u(d['baseline_median'], d['kpi'])} ± {fmt(d['baseline_scale'], d['kpi'])}, z = {d['z']:+.2f} ({d['direction']}) — {d['meaning']}")
    md += ["", "## 5. Batch assignment", ""]
    md.append("| batch | probability | distance | signature match (σ) |")
    md.append("|---|---|---|---|")
    for b in sorted(a["probabilities"]):
        sm = (a.get("signature_match") or {}).get(b)
        md.append(f"| {b} | {a['probabilities'][b]:.3f} | {a['distances'][b]:.2f} | {'—' if sm is None else f'{sm:+.2f}'} |")
    md += ["", f"- Assigned: **Batch {a['assigned_batch']}** — {a['confidence_label']}; strip-bootstrap stability "
           f"{fmt(a['stability'])}; single-strip agreement {fmt(a['strip_agreement'])}",
           f"- Novelty flag: **{'YES' if a['novelty_flag'] else 'no'}**; in-distribution score {fmt(a['in_distribution_score'])}",
           f"- Priors: {a['priors']}. {a['calibration_note']}"]
    for b, ks in (a.get("signature_kpis_matched") or {}).items():
        miss = (a.get("signature_kpis_missed") or {}).get(b, [])
        md.append(f"- Batch {b} signature — matched: {', '.join(map(short, ks)) or 'none'}; missed: {', '.join(map(short, miss)) or 'none'}")
    md += ["", "## 6. Figures", ""]
    for title, name in figs.items():
        md += [f"**{title}**", "", f"![{title}]({name})", ""]
    md += ["## 7. Not measurable from these images", "", ", ".join(vj["not_measurable"]) + ".", "", "## Caveats", ""]
    md += [f"{i + 1}. {c}" for i, c in enumerate(vj["caveats"])]
    if vj.get("notes"):
        md += ["", "Notes: " + "; ".join(vj["notes"])]
    md += ["", f"_Image flags: {', '.join(vj['image']['flags']) or 'none'}; runtime {vj['runtime_s']} s; baseline version {vj['baseline']['version']}._"]
    return "\n".join(md)


def md_to_html(md: str, figdir: Path, title: str) -> str:
    """Tiny Markdown subset -> HTML with images inlined as base64 (no external dependencies)."""
    import re

    def inline(s):
        s = html.escape(s)
        s = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s)
        s = re.sub(r"`(.+?)`", r"<code>\1</code>", s)
        s = re.sub(r"_(.+?)_$", r"<em>\1</em>", s)
        return s

    out, in_table, in_list, in_ol = [], False, False, False
    for line in md.splitlines():
        if line.startswith("|"):
            cells = [c.strip() for c in line.strip("|").split("|")]
            if all(set(c) <= set("-: ") for c in cells):
                continue
            if not in_table:
                out.append("<table>")
                in_table = True
                out.append("<tr>" + "".join(f"<th>{inline(c)}</th>" for c in cells) + "</tr>")
            else:
                out.append("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in cells) + "</tr>")
            continue
        if in_table:
            out.append("</table>")
            in_table = False
        m = re.match(r"!\[(.*?)\]\((.*?)\)", line)
        if m:
            p = figdir / m.group(2)
            if p.exists():
                b64 = base64.b64encode(p.read_bytes()).decode()
                out.append(f'<img alt="{html.escape(m.group(1))}" src="data:image/png;base64,{b64}">')
            continue
        if line.startswith("- "):
            if not in_list:
                out.append("<ul>")
                in_list = True
            out.append(f"<li>{inline(line[2:])}</li>")
            continue
        if in_list:
            out.append("</ul>")
            in_list = False
        m = re.match(r"(\d+)\. (.*)", line)
        if m:
            if not in_ol:
                out.append("<ol>")
                in_ol = True
            out.append(f"<li>{inline(m.group(2))}</li>")
            continue
        if in_ol:
            out.append("</ol>")
            in_ol = False
        if line.startswith("# "):
            out.append(f"<h1>{inline(line[2:])}</h1>")
        elif line.startswith("## "):
            out.append(f"<h2>{inline(line[3:])}</h2>")
        elif line.startswith("> "):
            out.append(f'<div class="banner">{inline(line[2:])}</div>')
        elif line.strip():
            out.append(f"<p>{inline(line)}</p>")
    if in_table:
        out.append("</table>")
    if in_list:
        out.append("</ul>")
    if in_ol:
        out.append("</ol>")
    css = """
:root{--bg:#ffffff;--fg:#1d2329;--muted:#5b6670;--line:#d9dee3;--accent:#1f5fa8;--warn:#fff3cd;--warnline:#e0a800}
body{background:var(--bg);color:var(--fg);font:15px/1.5 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;max-width:1100px;margin:24px auto;padding:0 16px}
h1{font-size:24px;border-bottom:2px solid var(--accent);padding-bottom:6px}h2{font-size:19px;margin-top:28px;color:var(--accent)}
table{border-collapse:collapse;font-size:13px;margin:8px 0;width:100%}th,td{border:1px solid var(--line);padding:4px 7px;text-align:left}
th{background:#f3f5f7}code{font-size:12px;background:#f3f5f7;padding:1px 4px;border-radius:3px}
img{max-width:100%;border:1px solid var(--line);margin:6px 0}.banner{background:var(--warn);border-left:4px solid var(--warnline);padding:8px 12px;margin:10px 0}
@media print{body{margin:0;max-width:none}h2{page-break-after:avoid}img{page-break-inside:avoid}}
"""
    return (f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>{html.escape(title)}</title><style>{css}</style></head><body>" + "\n".join(out) + "</body></html>")


def write_reports(vj, res, st, model, base, rdir: Path, cache_dir: Path):
    from .config import load_config

    cfg = load_config()
    rdir.mkdir(parents=True, exist_ok=True)
    figs = {}
    if (cache_dir / "overlay.png").exists():
        shutil.copy(cache_dir / "overlay.png", rdir / "overlay.png")
        figs["Segmentation overlay (blue = deep pore, orange = Si-candidate, red = crack skeleton; 4× downscaled)"] = "overlay.png"
    fig_zbars(vj, rdir / "zbars.png")
    figs["Robust z-scores vs baseline (green < 2.5σ, amber 2.5–4σ, red > 4σ)"] = "zbars.png"
    if fig_pca(base["df"], model, st, cfg, vj, rdir / "pca.png"):
        figs["Where the sample sits among labelled images (PCA of standardised KPIs)"] = "pca.png"
    fig_distributions(res, st, base["out"], rdir / "distributions.png")
    figs["Si-candidate size and pore chord distributions vs pooled baseline"] = "distributions.png"
    sig_png = Path(base["out"]) / "figs" / "signatures.png"
    if sig_png.exists():
        shutil.copy(sig_png, rdir / "signatures.png")
        figs["Batch signatures (what makes batches 1 and 2 different from baseline)"] = "signatures.png"
    md = build_markdown(vj, st, model, figs)
    (rdir / "report.md").write_text(md)
    (rdir / "report.html").write_text(md_to_html(md, rdir, f"QC {vj['sample_id']}"))
    pub = {k: v for k, v in vj.items() if not k.startswith("_")}
    (rdir / "verdict.json").write_text(json.dumps(pub, indent=1))
    (rdir / "dm.txt").write_text(dm_block(vj, model) + "\n")
