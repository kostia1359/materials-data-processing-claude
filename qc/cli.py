"""Command-line entry points:  python -m qc <command> ..."""
from __future__ import annotations

import csv
import json
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

import typer

from .config import load_config

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
app = typer.Typer(add_completion=False, help="Si-graphite anode SEM batch-QC")


@app.command("build-baseline")
def build_baseline(data: Path = typer.Option(..., help="data folder (Batch_1/, Batch_2/, Batch_3/ ...)"),
                   out: Path = typer.Option(Path("out")), config: Path = typer.Option(None),
                   workers: int = typer.Option(None, help="parallel processes (default min(3, cpus))")):
    """Analyse every sample, fix the crack threshold, fit baseline stats, centroids and signatures."""
    from .build import build_baseline as bb

    bb(data, out, load_config(config), workers=workers)


@app.command()
def evaluate(sample: Path = typer.Option(..., help="folder holding img_<id>_{BSE,ETD,Inlens}.tif (or one of the files)"),
             baseline: Path = typer.Option(Path("out")), report_dir: Path = typer.Option(None), config: Path = typer.Option(None)):
    """Verdict JSON + report.md + report.html + dm.txt for one sample."""
    from .build import evaluate_files, load_baseline
    from .io import discover

    cfg = load_config(config)
    folder = sample if sample.is_dir() else sample.parent
    found = discover(folder)
    if not sample.is_dir():
        sid = sample.name.split("_img_")[-1] if "_img_" in sample.name else sample.name
        found = [s for s in found if s.sample_id.split("img_")[-1] in sid] or found
    if not found:
        raise SystemExit(f"no sample found in {sample}")
    base = load_baseline(baseline)
    for s in found:
        rdir = report_dir or (baseline / "reports" / s.sample_id)
        vj = evaluate_files(s, base, cfg, rdir)
        typer.echo(f"report: {rdir / 'report.html'}")
        typer.echo(Path(rdir / "dm.txt").read_text())


@app.command()
def predict(folder: Path = typer.Option(...), baseline: Path = typer.Option(Path("out")), config: Path = typer.Option(None)):
    """Evaluate every sample in a test-drop folder -> predictions.json + dm.txt (paste-ready)."""
    from .build import evaluate_files, load_baseline
    from .io import discover

    cfg = load_config(config)
    base = load_baseline(baseline)
    samples = discover(folder)
    preds, blocks = [], []
    for s in samples:
        rdir = baseline / "reports" / s.sample_id
        vj = evaluate_files(s, base, cfg, rdir)
        a = vj["batch_assignment"]
        preds.append(dict(sample_id=s.sample_id, assigned_batch=a["assigned_batch"], probabilities=a["probabilities"],
                          confidence_label=a["confidence_label"], stability=a["stability"], novelty_flag=a["novelty_flag"],
                          in_distribution_score=a["in_distribution_score"], verdict=vj["verdict_vs_baseline"]["decision"],
                          report=str(rdir / "report.html")))
        blocks.append((rdir / "dm.txt").read_text().rstrip())
    meta = dict(created_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"), baseline_version=base["stats"].get("version"),
                folder=str(folder), n=len(preds))
    (baseline / "predictions.json").write_text(json.dumps(dict(meta=meta, predictions=preds), indent=1))
    (baseline / "dm.txt").write_text("\n\n".join(blocks) + "\n")
    typer.echo("\n\n".join(blocks))
    typer.echo(f"\nwrote {baseline / 'predictions.json'} and {baseline / 'dm.txt'}")


@app.command()
def score(predictions: Path = typer.Option(...), truth: Path = typer.Option(...),
          validation_md: Path = typer.Option(Path("VALIDATION.md"))):
    """Score pre-registered predictions against revealed labels (sample_id,batch with batch in 1,2,3,none)."""
    import numpy as np

    P = json.loads(predictions.read_text())
    with open(truth) as f:
        tr = {r["sample_id"].strip(): r["batch"].strip().lower() for r in csv.DictReader(f)}

    def norm(s):
        s = s.strip()
        return s if s.startswith("img_") else f"img_{s}"

    tr = {norm(k): v for k, v in tr.items()}
    rows, nll = [], []
    for p in P["predictions"]:
        t = tr.get(p["sample_id"])
        if t is None:
            continue
        ok = (p["assigned_batch"] == t)
        if t in ("1", "2", "3"):
            nll.append(-np.log(max(p["probabilities"].get(t, 0.0), 1e-3)))
        rows.append(f"| {p['sample_id']} | {t} | {p['assigned_batch']} | {p['probabilities'].get(t, float('nan')) if t != 'none' else float('nan'):.2f} "
                    f"| {p['confidence_label']} | {'yes' if p['novelty_flag'] else 'no'} | {'✓' if ok else ('novelty ✓' if t == 'none' and p['novelty_flag'] else '✗')} |")
    n_ok = sum(r.endswith("| ✓ |") or "novelty ✓" in r for r in rows)
    text = [f"\n## Test-drop score ({datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC; predictions made {P['meta']['created_utc']} with baseline {P['meta']['baseline_version']})", "",
            f"Correct: **{n_ok}/{len(rows)}**" + (f"; mean NLL on labelled truth {np.mean(nll):.3f}" if nll else ""), "",
            "| sample | truth | predicted | p(truth) | confidence | novelty | result |", "|---|---|---|---|---|---|---|", *rows, ""]
    with open(validation_md, "a") as f:
        f.write("\n".join(text))
    typer.echo("\n".join(text))
    typer.echo("Next: move the scored samples into their true Batch_N folders and re-run build-baseline + validate.")


@app.command()
def validate(data: Path = typer.Option(...), out: Path = typer.Option(Path("out")), config: Path = typer.Option(None),
             skip_sensitivity: bool = typer.Option(False), skip_synthetic: bool = typer.Option(False),
             validation_md: Path = typer.Option(Path("VALIDATION.md"))):
    """LOIO + confound + sensitivity + synthetic checks -> VALIDATION.md."""
    import pandas as pd

    from . import validate as V
    from .build import build_baseline as bb
    from .build import load_baseline
    from .io import discover

    cfg = load_config(config)
    if not (out / "baseline_stats.json").exists():
        bb(data, out, cfg)
    base = load_baseline(out)
    df = base["df"]
    strips = pd.read_csv(out / "strips.csv", dtype={"batch": str})
    t = time.time()
    loio = V.nested_loio(df, strips, cfg)
    typer.echo(f"LOIO done ({time.time() - t:.0f}s): accuracy {loio.get('accuracy')}")
    shadow = V.shadow_classifier(df, cfg)
    eff = V.effect_sizes(df, cfg)
    sens = None if skip_sensitivity else V.sensitivity(discover(data), df, base["stats"], cfg, log=typer.echo)
    synth = None if skip_synthetic else V.synthetic_checks(cfg, log=typer.echo)
    counts = df[df["batch"].isin(["1", "2", "3"])]["batch"].value_counts().to_dict()
    V.write_validation_md(validation_md, loio, shadow, eff, sens, synth, counts, base["stats"], base["model"])
    if loio.get("table") is not None and len(loio["table"]):
        loio["table"].to_csv(out / "loio.csv", index=False)
    if sens is not None:
        sens.to_csv(out / "sensitivity.csv", index=False)
    typer.echo(f"wrote {validation_md}")


@app.command()
def regress(data: Path = typer.Option(Path("data")), out: Path = typer.Option(Path("out")), accept: bool = typer.Option(False),
            config: Path = typer.Option(None), workers: int = typer.Option(None)):
    """Recompute KPIs (cached stages reused) and compare with the last accepted snapshot."""
    from .build import build_baseline as bb
    from .regress import run

    bb(data, out, load_config(config), workers=workers)
    raise typer.Exit(run(out, accept, log=typer.echo))


@app.command()
def quicklook(data: Path = typer.Option(...), out: Path = typer.Option(Path("out/quicklook")), config: Path = typer.Option(None)):
    """Downscaled 3-channel composite + histograms per sample, and a facts table (DATA_FACTS input)."""
    from .quicklook import run

    run(data, out, load_config(config))


@app.command()
def fetch(manifest: Path = typer.Option(Path("data/drive_manifest.csv")), data: Path = typer.Option(Path("data")),
          verify_only: bool = typer.Option(False),
          only: list[str] = typer.Option(None, help="sample ids to fetch (repeatable); default all")):
    """Download the Drive files listed in the manifest (GOOGLE_API_KEY via googleapis.com, else gdown/drive.google.com)."""
    from .download import run

    run(manifest, data, verify_only, only=only or None)


@app.command()
def synth(out: Path = typer.Option(Path("synthdata")), per_batch: int = typer.Option(4), width: int = typer.Option(3500),
          height: int = typer.Option(1000)):
    """Write a SYNTHETIC 3-batch dataset (for pipeline tests and demos only)."""
    from .synth import make_dataset

    make_dataset(out, per_batch, width, height)
    typer.echo(f"wrote synthetic dataset to {out}")


def main():
    app()


if __name__ == "__main__":
    main()
