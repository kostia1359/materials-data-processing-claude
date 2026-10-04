# anode-qc — Si–graphite anode SEM batch QC

Interpretable QC for SEM cross-sections of a silicon–graphite Li-ion anode. For each sample (three co-registered
images `img_<id>_BSE.tif`, `img_<id>_ETD.tif`, `img_<id>_Inlens.tif`) it:

1. extracts ~40 material KPIs with physical units and per-sample uncertainty, plus a 16-KPI shortlist,
2. compares against the **batch-3 baseline** → ACCEPT / INVESTIGATE / REJECT with the KPIs driving it,
3. **always bets on a batch (1, 2 or 3)** with a calibrated probability, stability, novelty flag and a plain-language
   explanation through batch *signatures* ("what makes batch 2 different from baseline"),
4. validates itself (nested leave-one-image-out, synthetic ground truth, confound and sensitivity checks → `VALIDATION.md`).

No neural networks: every number is traceable to a mask on the image. Phase identity comes from BSE only; ETD
confirms pore floors; InLens/ETD intensities are acquisition covariates ("gates").

## Quick start

```bash
pip install -e .            # or: pip install numpy scipy scikit-image tifffile imagecodecs pandas scikit-learn matplotlib pyyaml typer tabulate pytest
python -m qc fetch                                  # download the Drive folder in data/drive_manifest.csv (gdown; GOOGLE_API_KEY route optional); resumable, size-verified
python -m qc quicklook --data data                  # composites + facts table per sample (re-verify DATA_FACTS.md)
python -m qc build-baseline --data data --out out   # KPIs, baseline stats, centroids, signatures + simulation layer (--no-sim to skip)
python -m qc regress --data data --out out          # Section 0.3 loop: recompute (cached) and compare with the last KPI snapshot
python -m qc validate --data data --out out         # -> VALIDATION.md
python -m qc evaluate --sample ./new/ --baseline out  # one sample -> out/reports/<id>/report.html (+ .md, verdict.json, dm.txt)
```

Data layout: `data/Batch_1/`, `data/Batch_2/`, `data/Batch_3/` (any folder whose name contains `batch N`);
samples in `data/test/` or `data/batch_unknown/` are labelled `unknown` and never used for fitting.
A `Batch_N_` prefix in front of `img_` is ignored; an `_SE.tif` is accepted in place of `_ETD.tif` (flagged).

## Test-drop workflow (organiser drops a folder, truth arrives next morning)

```bash
python -m qc predict --folder ./test_drop --baseline out     # -> out/predictions.json + out/dm.txt (paste into the DM)
# next morning, truth.csv = "sample_id,batch" with batch in {1,2,3,none}
python -m qc score --predictions out/predictions.json --truth truth.csv   # appended to VALIDATION.md
# then move the scored samples into their true Batch_N folders and re-run build-baseline + validate
```

`predictions.json` stores the baseline version hash and a UTC timestamp, so the morning score is a pre-registered test.

## Development loop (brief §0.3)

New triples are added to `data/` one or two at a time; `qc regress` recomputes only the cache stages whose code or
config changed (`out/cache/<id>/{seg.npz,sato.npy,result.json,sim.json}` keyed by file SHA-1 + stage version) and
fails if any KPI of an already-processed sample moves by ≥ 0.1 × its baseline scale (or ≥ 1 % for fractions). An
intended change is accepted with `--accept` after writing the reason in `DECISIONS.md`; snapshots live in
`out/snapshots/kpis_<n>.csv` (committed). It also flags KPIs outside the plausible ranges of §0.3.2.

## Simulation layer (brief §10, `qc/sim/`)

Additive, never feeds the classifier: an evidence-fused phase map (BSE + ETD, InLens only gates smoothing) with an
explicit *uncertain* class → nested bound maps L_solid ⊆ L_mid ⊆ L_pore → ×4 block-majority maps (audited) →
two-conductivity finite-volume Laplace solves (pore 1, carbon D_c = 0.05, Si 0; refuses non-spanning masks),
electronic network, Si swelling scenarios (f_A 1.6 / 2.43), closed-form indices and a PyBaMM composite
graphite–Si DFN (Chen2020_composite) under both porosity conventions. Everything is reported as a ratio to the
baseline median with A/B/C trust grades and L_solid…L_pore bounds; outputs `out/sim.csv`, `out/sim_bounds.csv`,
`out/sim_baseline.json` and a `simulation` block in every `verdict.json`.

## Outputs

| file | content |
|---|---|
| `out/manifest.csv` | discovered samples, batch, paths, size, pixel size |
| `out/kpis.csv`, `out/strips.csv` | one row per sample / per (sample, strip): KPIs, `<kpi>_se`, gates |
| `out/baseline_stats.json` | per-KPI `n, med, mad, within_se, scale`, crack threshold, Δ, T, LOO scores, version |
| `out/model.json`, `out/signatures.json` | shrunken centroids, tuning table, within-batch distances; batch signatures |
| `out/figs/` | `signatures.png`, `kpi_boxplots.png` |
| `out/reports/<id>/` | `report.html` (standalone, figures embedded), `report.md`, `verdict.json`, `dm.txt`, figures |

An example report for the starter sample is in `examples/img_0grcilhi/` (open `report.html`). With only one
baseline image it is a self-comparison and demonstrates layout only.

## How it works (one paragraph per stage)

- **Hygiene** (`qc/io.py`): drop edge columns where R≠G (green edge artefact), crop un-sectioned edge bands,
  read 25 nm/px from TIFF tags (resample if different).
- **Segmentation** (`qc/segment.py`): Gaussian σ 1.5 px on BSE, 3-class multi-Otsu → deep pore (hysteresis +
  ETD-dark confirmation) / carbon matrix / Si-candidate (opening + interior check). Watershed splits touching particles.
- **KPIs** (`qc/kpis/`): fractions, area-weighted Si D10/50/90, shape, dispersion (quadrat CV, Clark–Evans,
  agglomerates), pore chords and local thickness, anisotropy (chord ratios, structure tensor), two-point statistics
  (correlation lengths, integral range → ImageRep SE), percolation/Euler, interfaces/contacts, cracks (Sato ridges),
  inclusions, vertical gradients. Each KPI has unit, robustness class R/M/W and a battery meaning
  (`qc/kpis/registry.py`). Uncertainty: 5 vertical strips → SE and block bootstrap.
- **Baseline & verdict** (`qc/stats.py`): robust z with `scale = max(MAD, median strip SE, 5 % of median)`;
  REJECT/INVESTIGATE zones; conformal rank p with its floor 1/(n+1) always stated.
- **Assignment** (`qc/assign.py`): nearest shrunken centroid on the 16 shortlist KPIs in baseline-σ units,
  softmax with temperature; Δ and T by leave-one-image-out log-loss; uniform priors; strip-bootstrap stability;
  in-distribution score and novelty flag; signatures explain the bet.
- **Validation** (`qc/validate.py`, `qc/synth.py`, `tests/`): see `VALIDATION.md`.

Parameters live in `config.yaml`; every judgment call is in `DECISIONS.md`; verified data facts in `DATA_FACTS.md`.

## Runtime

KPI pipeline ≈ 50–70 s per 7000 × ~2000 sample on 4 CPU cores (Sato ridge filter ≈ 25 s of it); simulation layer
≈ 60–100 s more (11 sparse solves on the ×4 map + PyBaMM ≈ 15 s). A full-resolution solve was measured at 354 s /
7.8 GB, hence the audited ×4 path. `build-baseline` runs 3 samples in
parallel (`--workers`, or `QC_WORKERS`); results are cached in `out/samples/<id>/` and reused when files are unchanged.
`pytest` (14 tests, synthetic ground truth) takes ~3 min.

## What this cannot measure

Coating thickness, surface roughness, collector delamination (no collector or surface in frame), binder/carbon-black
distribution (same Z as graphite), Si vs SiOx identity (no EDS), true porosity (open pores → lower bound only).
These caveats are printed in every report.
