# DECISIONS

One line per judgment call: **what** — why — alternative considered. Post-hoc changes are marked **[post-hoc]**.

## Data and I/O
- **Sample ID = `img_<id>`, any prefix before `img_` ignored** — 7 of 31 Drive samples carry a `Batch_N_` prefix and the rest do not; the brief's regex would miss them. Alternative: rename files on disk (rejected: the judged drop may use either form).
- **`SE` accepted as the ETD channel and flagged (`etd_is_SE_detector`)** — 4 Drive samples (`rxax5ozo`, `vc2whyaq`, `x77cy643`, `utfgcjfa`) ship `_SE.tif` instead of `_ETD.tif`. Both are secondary-electron images; ETD is only used to confirm pore floors and as a gate, so the role is the same. Alternative: treat those samples as BSE+InLens only (loses pore confirmation).
- **Batch label from the first `batch\s*_?\d` in the path; `test/` or `batch_unknown/` → unknown.**
- **RGB-disagreement columns/rows are cropped only at the frame edge** (observed: 2-px green column at the right edge, fraction-of-rows disagreement > 50 %). Interior disagreement is flagged, not cropped. Alternative: drop every disagreeing column (would split the frame).
- **All masking is rectangular cropping** (bad edge columns, edge bands, manual masks), so the valid area is the whole working array and every per-area normalisation is exact.
- **Edge-band rule: row flag smoothed over 9 rows, band must touch the frame edge and be ≤ 10 % of height.** Not triggered on the starter sample; the brief reports ~90 nodular rows on another sample. Override via `manual_masks` in `config.yaml`.
- **Pixel size read from X/YResolution (25.0005 nm/px on the starter triple)**; resampling only if > 2 % off 25 nm.
- **Raw data are not committed** (93 × ~18 MB). `data/drive_manifest.csv` (made from a Drive listing on 2026-10-03) + `qc fetch` reproduce the folder; `qc fetch --verify-only` checks sizes.

## Segmentation
- **3-class multi-Otsu on the σ = 1.5 px smoothed BSE** (brief). Starter sample: t1 = 38.2, t2 = 81.6 (brief expected ~[40, 82]).
- **Pore = hysteresis-grown BSE core AND ETD darker than the 10th percentile of ETD in the carbon class.** On the starter sample this keeps 16.1 % of 17.1 % grown pixels; 1.0 % are "ambiguous" (BSE-dark, ETD-bright).
- **Si-candidate = opening(disk 1) of BSE ≥ t2, then the per-component interior check (median of 1-px-eroded interior ≥ t2 + 10).** Si fraction drops from 7.1 % (raw t2) to 5.3 % of the frame after opening + interior check (rejects pore-edge halos), which is the intended behaviour.
- **Watershed splitting done per connected component bounding box** (identical result to a global watershed, ~10× faster); components < 200 px are not split.
- **Shading diagnostic**: the starter sample's per-band t2 drifts from 84 to 68 (bottom band darker; BSE row-mean gradient 18 %), so `shading_flag` is set. KPIs still use global thresholds (brief); the flag is reported as a gate.

## KPIs
- **Strips = 5 equal-width strips of the cropped frame** (1399–1400 px), not exactly 1400 px, because the frame is 6998 px after removing the 2 bad columns.
- **Si contact ring taken 2–4 px from the particle, not the 1-px ring** — the smoothed BSE always passes through the carbon grey band between a bright particle and a dark pore, so the 1-px ring classified 100 % carbon on the real sample (Si–pore contact = 0.000). With the 2–4 px ring: 1.3 % pore contact. Alternative: unsmoothed BSE ring (noisier).
- **Local thickness by successive openings (two EDTs per radius, radii 1–20 then coarser) on the ×2-downsampled pore mask; D50 interpolated within the quantisation bin.** porespy not required.
- **Integral range A₂ = sum of the autocovariance inside the first zero crossing of its radial average** **[post-hoc, from synthetic truth]** — summing the connected positive lobe let small positive noise percolate across the frame and over-estimated the SE ~7×.
- **`_px` two-point lengths are expressed in 25-nm reference pixels** (scaled by px_nm/25) so they survive a magnification change; names kept from the brief.
- **Structure tensor computed on the 3-phase map (pore 0 / carbon 1 / Si 2), not on raw BSE** **[post-hoc, from synthetic perturbation test]** — on intensity, `st_coherence` moved 0.76 scale under gamma 1.2 (an R-class KPI must move < 0.25). σ scaled with pixel size.
- **Clark–Evans uses Donnelly's edge-corrected expectation** **[post-hoc, deviation from Appendix B]** — the uncorrected ratio read 1.26 for Poisson truth in a 25 µm-tall frame. The synthetic test uses 1 µm particles at 1.5 % (mean of 3 seeds), because 4 µm particles at 7 % are close to hard-core packing and R ≈ 1 regardless of clustering.
- **Crack threshold = median over baseline images of each image's 99.5th percentile Sato response** (brief); Sato response is cached as float16 so the threshold can be fixed after all baseline images are seen. With 1 baseline image it is that image's own percentile.
- **Curtaining exclusion only applied when the Fourier peak in 50–80° is ≥ 3× the median angular power** — the starter InLens shows only a 1.6× peak at 51.5° (no visible curtaining); excluding a direction there would delete real cracks.
- **`si_agglomerate_frac` kept as specified even though it reads 0.90 on the starter sample** — heavy-tailed component sizes put most area above 3× the median; it is compared like-for-like and not on the shortlist.
- **SE floor for scale uses `within_se` = median per-image strip SE** (brief). On half-size synthetic frames this dominates the scale (strip SE of Si fraction ≈ 2.5 pp) — real frames are 2× wider so strips hold more particles.

## Statistics and assignment
- **Grand centre for shrinkage = median of class centroids** (not of all images), so batch sizes cannot pull it — consistent with uniform priors.
- **Missing KPI value → feature set to the baseline median (0 after standardisation)** — "no evidence either way".
- **Strip vectors re-centred on the image vector for stability/bootstrap** — strip KPIs are biased relative to whole-image KPIs (narrower field, particles cut at strip edges); using raw strip means gave 20 % stability for a p = 1.00 synthetic call. Only the strip spread is used.
- **Inner tuning keeps the outer fold's baseline stats fixed** (one level of nesting for stats, two for Δ/T).
- **id_score = 1 − ECDF of LOIO within-batch distances; novelty when min distance > their 97.5th percentile, or when the sample has |z| > 4 on KPIs outside every signature and matches no signature (> 1.5 σ).**
- **Signature stability = share of leave-one-image-out refits (baseline stats rebuilt too) with the same sign and |effect| > 1; signature = stability ≥ 0.8 and |effect| ≥ 1.5.**
- **Defect rule uses event counts (cracks, inclusions) scaled by area against the baseline maximum count + 3 Poisson SD.**
- **When the evaluated sample is in the labelled table it is held out** of baseline stats, centroids and signatures for that evaluation; if it is the only baseline image, the report says the comparison is self-referential.

## Robustness classes (frozen a priori, not changed)
- **`si_contact_pore_frac` stays class R as in the brief, although synthetic perturbation moves it 0.35–1.0 scale** (its mask edge moves with t1). Changing the class would change its weight post hoc; it is listed here as the first candidate for demotion to M once real baseline scales are known.
- **Pore chord KPIs (R) move 1–3 scale under ±5-level t1 shifts on synthetic data** because the synthetic baseline scale is tiny (homogeneous replicates → 5 % floor). Re-check on the real baseline in VALIDATION.md before trusting them as R.

## Not done (time-boxed)
- `tau_index_v` (2-D diffusion solve), FastAPI/React front end — Phase 6 stretch items.
- The physics-simulation layer from research report 2 (brief Section 10) — not in the brief text supplied, so not built.

## Added after the synthetic 3-batch demonstration
- **Novelty by distance: `min_b d_b` > 97.5th percentile of LOIO within-batch distances only when ≥ 40 reference distances exist; otherwise > 1.25 × the largest one** **[post-hoc]** — with ~11 distances the 97.5th percentile equals the maximum, and a synthetic batch-2 sample matching batch 2's signature at +5.8 σ was flagged "novel". `id_score` is still reported as in the brief. A distance-based novelty is also withdrawn when the sample matches the assigned batch's signature at ≥ 1.5 σ (a more extreme member of that batch); the report then says "beyond the labelled range, same direction".
- **Minimum-size cut-offs (Si 100/500 px², pore 20 px², inclusions 50 px², crack length 80 px / width 6 px, quadrat 700 px) are in 25-nm reference pixels and rescale with pixel size** **[post-hoc, from the 2× resolution check]** — at 50 nm/px the unscaled cut-offs changed Si number density by ~6 scale units.
