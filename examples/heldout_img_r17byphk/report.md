# QC report — img_r17byphk

> **Acquisition drift or artefact suspected — compare with caution.** sat255_bse z=+1579.7; sat0_inlens z=+8.1; shading_flag

## 1. Verdict

**img_r17byphk is most consistent with Batch 3** (probability 0.61; Batch 1: 0.16; Batch 2: 0.23; assignment stable in 99 % of strip-bootstrap resamples; moderately confident). **Relative to the Batch-3 baseline it is within the baseline distribution — verdict ACCEPT**: 0 of 16 KPIs exceed the investigate zone (8 robust σ; reject at 12.8). Largest shifts: Si-candidate number density 34,615 mm^-2 vs 26,433 mm^-2 ± 5,008 (z = +1.6); Si-candidate/pore contact 2.63 % vs 1.71 % ± 0.70 % (z = +1.3); pore correlation length 39.7 px vs 33.3 px ± 5.00 (z = +1.3). Baseline has 17 image(s), so the smallest reportable rank p-value is 0.056; this verdict is an effect-size judgement, not a significance test.

## 2. Acquisition gates

Status: **caution** — flags: sat255_bse z=+1579.7, sat0_inlens z=+8.1, shading_flag

| gate | value | z |
|---|---|---|
| `sat255_bse` | 0.000 | +1579.7 |
| `sat0_inlens` | 0.009 | +8.1 |
| `sat255_etd` | 0.000 | +3.4 |
| `sat0_bse` | 0.024 | +2.6 |
| `sat0_etd` | 0.037 | +2.5 |
| `etd_c_med` | 82.0 | +2.0 |
| `n_grey_levels_bse` | 209.0 | +1.6 |
| `bse_row_gradient_pct` | 8.47 | +1.6 |
| `t1` | 37.7 | -0.8 |
| `noise_sd_bse` | 10.7 | +0.7 |

## 3. KPIs

| KPI | value ± SE | baseline median ± scale | z | class | shortlist |
|---|---|---|---|---|---|
| `si_frac_solid` | 5.79 % ± 1.02 % | 6.71 % ± 0.84 % | -1.09 | M | ✓ |
| `si_d50_aw_um` | 2.90 µm ± 0.510 | 3.19 ± 0.334 | -0.85 | R | ✓ |
| `si_d90_aw_um` | 6.35 µm ± 0.620 | 5.82 ± 0.776 | +0.69 | R | ✓ |
| `si_num_density_mm2` | 34,615 mm^-2 ± 8,531 | 26,433 ± 5,008 | +1.63 | M | ✓ |
| `si_solidity_med` | 0.856 ± 0.014 | 0.884 ± 0.044 | -0.63 | R | ✓ |
| `si_quadrat_cv` | 0.897 ± 0.072 | 0.749 ± 0.180 | +0.82 | M | ✓ |
| `pore_frac_deep` | 13.39 % ± 0.69 % | 12.11 % ± 1.14 % | +1.12 | W | ✓ |
| `pore_chord_h_mean_um` | 0.702 µm ± 0.042 | 0.644 ± 0.075 | +0.78 | R | ✓ |
| `pore_chord_v_mean_um` | 0.572 µm ± 0.029 | 0.540 ± 0.053 | +0.61 | R | ✓ |
| `aniso_pore_chord_ratio` | 1.23 ± 0.032 | 1.18 ± 0.059 | +0.85 | R | ✓ |
| `aniso_carbon_chord_ratio` | 1.17 ± 0.047 | 1.20 ± 0.060 | -0.53 | R | ✓ |
| `s2_len_pore_h_px` | 39.7 px ± 3.85 | 33.3 ± 5.00 | +1.29 | R | ✓ |
| `pore_percolating_frac_v` | 0.00 % ± 0.00 % | 0.00 % ± 0.50 % | +0.00 | M | ✓ |
| `si_contact_pore_frac` | 2.63 % ± 0.44 % | 1.71 % ± 0.70 % | +1.31 | R | ✓ |
| `si_contact_carbon_frac` | 97.37 % ± 0.44 % | 98.29 % ± 4.91 % | -0.19 | R | ✓ |
| `crack_density_um_per_mm2` | 1,497 µm/mm^2 ± 757.6 | 2,541 ± 1,707 | -0.61 | M | ✓ |
| `si_frac_total` | 5.02 % ± 0.92 % | 5.89 % ± 0.87 % | -1.00 | M |  |
| `si_bse_contrast` | 1.27 ± 0.034 | 1.73 ± 0.234 | -1.95 | M |  |
| `si_d10_aw_um` | 1.23 µm ± 0.178 | 1.31 ± 0.158 | -0.52 | M |  |
| `si_span` | 1.76 ± 0.141 | 1.41 ± 0.216 | +1.61 | M |  |
| `si_aspect_med` | 1.80 ± 0.092 | 1.70 ± 0.085 | +1.18 | R |  |
| `si_circularity_med` | 0.541 ± 0.027 | 0.604 ± 0.030 | -2.07 | R |  |
| `si_interior_texture` | 0.672 ± 0.021 | 0.733 ± 0.075 | -0.81 | M |  |
| `si_clark_evans_R` | 0.723 ± 0.055 | 0.754 ± 0.062 | -0.49 | M |  |
| `si_agglomerate_frac` | 92.12 % ± 2.14 % | 87.85 % ± 4.39 % | +0.97 | M |  |
| `pore_frac_slope_per_level` | 0.003 1/level ± n/a | 0.004 ± 0.001 | -0.39 | W |  |
| `ambiguous_frac` | 0.76 % ± 0.06 % | 0.62 % ± 0.50 % | +0.28 | W |  |
| `pore_chord_h_p90_um` | 1.68 µm ± 0.094 | 1.55 ± 0.185 | +0.67 | M |  |
| `pore_lt_d50_um` | 0.799 µm ± 0.039 | 0.694 ± 0.056 | +1.87 | M |  |
| `st_coherence` | 0.146 ± 0.017 | 0.122 ± 0.028 | +0.84 | R |  |
| `st_orientation_deg` | -2.92 deg ± 1.98 | 2.73 ± 3.42 | -1.65 | R |  |
| `s2_len_pore_v_px` | 23.1 px ± 0.959 | 20.5 ± 3.40 | +0.76 | R |  |
| `s2_len_si_px` | 62.9 px ± 6.95 | 66.7 ± 7.69 | -0.49 | R |  |
| `s2_integral_range_pore_px2` | 8,801 px^2 ± 822.8 | 6,172 ± 1,730 | +1.52 | M |  |
| `pore_percolating_frac_h` | 0.00 % ± 0.00 % | 0.00 % ± 0.50 % | +0.00 | M |  |
| `pore_euler_density_mm2` | 140,217 mm^-2 ± 17,761 | 141,234 ± 23,120 | -0.04 | M |  |
| `interface_pore_solid_um_per_um2` | 0.665 µm^-1 ± 0.009 | 0.651 ± 0.033 | +0.44 | M |  |
| `largest_void_ecd_um` | 10.1 µm ± 0.331 | 9.94 ± 3.33 | +0.04 | M |  |
| `hiZ_inclusion_count_mm2` | 0.000 mm^-2 ± 0.000 | 0.000 ± 0.000 | +0.00 | M |  |
| `vertical_pore_slope` | 0.251 pp/10µm ± 0.712 | 0.127 ± 0.580 | +0.21 | W |  |
| `vertical_si_slope` | -0.900 pp/10µm ± 0.880 | -0.164 ± 0.623 | -1.18 | W |  |
| `si_frac_se_imagerep` | 0.86 % ± 0.25 % | 1.01 % ± 0.50 % | -0.29 | R |  |
| `pore_frac_se_imagerep` | 0.84 % ± 0.10 % | 0.70 % ± 0.50 % | +0.27 | R |  |

## 4. Drivers (shortlist, sorted by |z|)

- **Si-candidate number density** (`si_num_density_mm2`, class M): 34,615 mm^-2 vs 26,433 mm^-2 ± 5,008, z = +1.63 (higher) — more Si-candidate particles per area (finer or more loaded)
- **Si-candidate/pore contact** (`si_contact_pore_frac`, class R): 2.63 % vs 1.71 % ± 0.70 %, z = +1.31 (higher) — more Si-candidate surface facing pores: more room to swell, but more SEI
- **pore correlation length** (`s2_len_pore_h_px`, class R): 39.7 px vs 33.3 px ± 5.00, z = +1.29 (higher) — pore network stretched in-plane (longer correlation length): flatter, more aligned pores; consistent with heavier calendering
- **deep-pore fraction** (`pore_frac_deep`, class W): 13.39 % vs 12.11 % ± 1.14 %, z = +1.12 (higher) — more deep porosity (lower bound): better electrolyte transport, lower energy density
- **Si-candidate loading** (`si_frac_solid`, class M): 5.79 % vs 6.71 % ± 0.84 %, z = -1.09 (lower) — less Si-candidate in the solid: lower capacity, less swelling
- **pore anisotropy** (`aniso_pore_chord_ratio`, class R): 1.23 vs 1.18 ± 0.059, z = +0.85 (higher) — pores more aligned parallel to the collector: higher through-plane tortuosity, slower fast charge; consistent with heavier calendering
- **Si-candidate D50** (`si_d50_aw_um`, class R): 2.90 µm vs 3.19 µm ± 0.334, z = -0.85 (lower) — finer Si-candidate: more surface/SEI, less fracture
- **Si-candidate patchiness** (`si_quadrat_cv`, class M): 0.897 vs 0.749 ± 0.180, z = +0.82 (higher) — Si-candidate less evenly distributed (mixing/dispersion issue): local swelling hot-spots

## 5. Batch assignment

| batch | probability | distance | signature match (σ) |
|---|---|---|---|
| 1 | 0.156 | 4.53 | — |
| 2 | 0.229 | 4.18 | — |
| 3 | 0.615 | 3.10 | — |

- Assigned: **Batch 3** — moderately confident; strip-bootstrap stability 0.991; single-strip agreement 0.800
- Novelty flag: **no**; in-distribution score 0.645
- Priors: uniform. T=4.0 and shrinkage Δ=0.0 chosen by leave-one-image-out log-loss on 30 labelled images; LOIO accuracy 0.47
- Batch 1 signature — matched: none; missed: none
- Batch 2 signature — matched: none; missed: none

## 6. Figures

**Segmentation overlay (blue = deep pore, orange = Si-candidate, red = crack skeleton; 4× downscaled)**

![Segmentation overlay (blue = deep pore, orange = Si-candidate, red = crack skeleton; 4× downscaled)](overlay.png)

**Robust z-scores vs baseline (green < 8σ, amber 8–12.8σ, red > 12.8σ)**

![Robust z-scores vs baseline (green < 8σ, amber 8–12.8σ, red > 12.8σ)](zbars.png)

**Where the sample sits among labelled images (PCA of standardised KPIs)**

![Where the sample sits among labelled images (PCA of standardised KPIs)](pca.png)

**Si-candidate size and pore chord distributions vs pooled baseline**

![Si-candidate size and pore chord distributions vs pooled baseline](distributions.png)

**Batch signatures (what makes batches 1 and 2 different from baseline)**

![Batch signatures (what makes batches 1 and 2 different from baseline)](signatures.png)

**Simulation: transport index vs assumed D_c, with L_solid/L_pore bounds, over the baseline band**

![Simulation: transport index vs assumed D_c, with L_solid/L_pore bounds, over the baseline band](sim_sweep.png)

**Simulation: Si-candidate swelling before/after and per-particle constraint**

![Simulation: Si-candidate swelling before/after and per-particle constraint](swelling.png)

## 6b. Simulation layer — relative performance indices (ratios to baseline, frozen assumptions)

| index | grade | value | baseline median | ratio | bound interval (L_solid … L_pore) |
|---|---|---|---|---|---|
| Through-plane electrolyte transport (2-conductivity D_eff) (higher = easier ionic transport) | B | 0.059 | 0.059 | 1.01× | 0.054 … 0.064 |
| Through-plane tortuosity index (higher = more tortuous) | B | 2.18 | 2.00 | 1.09× |  |
| In-plane / through-plane transport anisotropy (higher = more aligned (calendering)) | B | 1.32 | 1.21 | 1.09× |  |
| Limiting current proxy (higher = better fast charge) | B | n/a | n/a | 1.01× |  |
| Plating-risk ionic-resistance ratio (higher = more plating risk) | B | n/a | n/a | 0.99× |  |
| Electronic network conductance (Si insulating) (higher = better wired) | B | 0.241 | 0.298 | 0.81× |  |
| Swelling: particles with enough pore buffer (full lithiation) (higher = better) | B | 0.013 | 0.036 | 0.37× | 0.013 … 0.048 |
| Swelling: constraint index (area-weighted) (higher = more matrix displaced) | B | 0.818 | 0.812 | 1.01× |  |
| Swelling: 2-D thickness-change bound (higher = more electrode swelling) | B | 0.046 | 0.056 | 0.82× |  |
| First-cycle-loss proxy (exposed interface) (higher = more SEI) | B | 0.495 | 0.485 | 1.02× |  |
| Cycle-life index (heuristic) (higher = worse) | C | n/a | n/a | 1.00× |  |
| Energy density ratio [offset_baseline_to_0.30] | A | ε = 0.313 |  | 0.952× |  |
| PyBaMM 1C / 3C CC-charge capacity [offset_baseline_to_0.30] | B | 0.695 / 0.302 |  |  |  |
| PyBaMM min negative surface Δφ at separator, 3C (V; < 0 = plating indicator) [offset_baseline_to_0.30] | B | -0.107 |  |  |  |
| Energy density ratio [deep_pore_as_is] | A | ε = 0.134 |  | 0.956× |  |
| PyBaMM 1C / 3C CC-charge capacity [deep_pore_as_is] | B | 0.630 / 0.142 |  |  |  |
| PyBaMM min negative surface Δφ at separator, 3C (V; < 0 = plating indicator) [deep_pore_as_is] | B | -0.362 |  |  |  |

**Marked unstable** (sample ranking changes between conventions, brief 10.4): Q_CC_3C_over_Q_C10 rank: porosity offset vs as-is; Q_CC_1C_over_Q_C10 rank: porosity offset vs as-is. Read those rows as direction-free.

Checks: bounds ordered True, per-strip ordered True, flux imbalance 2.5e-11, downsampling audit ok True; uncertain pixels 5.4 %; assumptions hash f5d1bb3ceb. Grades: A = arithmetic on measured quantities; B = direction supported, level set by an assumption (ratios only); C = heuristic.

## 7. Not measurable from these images

coating thickness, surface roughness, collector delamination, binder/carbon-black distribution, Si vs SiOx identity, true (total) porosity.

## Caveats

1. Bright particles are identified by backscatter Z-contrast only (no EDS); they are reported as Si-candidates and may include SiOx or Si-C composite.
2. Pores appear open (not resin-infiltrated); the deep-pore fraction is a lower bound on porosity, not a porosity measurement.
3. Graphite, carbon black and binder cannot be separated in these images; the 'carbon matrix' class contains all three.
4. Coating thickness, surface roughness and collector delamination are not assessable: the fields lie entirely inside the coating.
5. Anisotropy and transport indices are 2-D section quantities compared like-with-like against the baseline, not 3-D values; the through-plane direction is assumed vertical in the frame.
6. Pixel size (25 nm/px) is taken from the TIFF resolution tags; vendor metadata is absent, so it is unverified.
7. With N baseline images the smallest achievable rank p-value is 1/(N+1); verdicts are effect-size judgements against a small baseline.
8. Image grey levels have been remapped after acquisition (comb histograms); any method relying on raw intensities would be confounded - this system uses BSE phase identity after smoothing and treats ETD/InLens levels as acquisition covariates.
9. Transport and mechanics indices are 2-D effective-medium quantities computed with a fixed matrix diffusivity D_c = 0.05 and fixed moduli; they are ratios to the baseline under identical assumptions, not electrode tortuosity, conductivity or stress values.
10. The deep-pore phase does not percolate in 2-D, so no pore-only tortuosity is reported; the two-conductivity index's absolute level is set by D_c and only its ratios are meaningful.
11. Cell-level outputs come from a PyBaMM composite graphite-Si model with a frozen LG-M50-type parameter set and are relative rate-capability and plating-indicator shifts, not predictions of the real cell.

_Image flags: none; runtime 8.1 s; baseline version a381ba39._