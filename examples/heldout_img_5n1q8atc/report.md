# QC report — img_5n1q8atc

> **Acquisition drift or artefact suspected — compare with caution.** sat255_etd z=+13.9; sat0_inlens z=+15.3; bse_si_med z=-4.0

## 1. Verdict

**img_5n1q8atc is most consistent with Batch 1** (probability 0.92; Batch 2: 0.07; Batch 3: 0.01; assignment stable in 100 % of strip-bootstrap resamples; confident). **Novelty: it does not resemble any known batch closely** (in-distribution score 0.06); betting on Batch 1 as the nearest, treat as INVESTIGATE. **Relative to the Batch-3 baseline it is different from the baseline — verdict REJECT**: 3 of 16 KPIs exceed the investigate zone (8 robust σ; reject at 12.8). Largest shifts: Si-candidate number density 202,516 mm^-2 vs 26,433 mm^-2 ± 5,008 (z = +35.2); Si-candidate loading 20.25 % vs 6.71 % ± 0.84 % (z = +16.1); Si-candidate/pore contact 7.40 % vs 1.71 % ± 0.70 % (z = +8.1). Reads as: more Si-candidate particles per area (finer or more loaded); more Si-candidate in the solid: higher capacity but more swelling, first-cycle loss and N/P shift. Baseline has 17 image(s), so the smallest reportable rank p-value is 0.056; this verdict is an effect-size judgement, not a significance test.

## 2. Acquisition gates

Status: **caution** — flags: sat255_etd z=+13.9, sat0_inlens z=+15.3, bse_si_med z=-4.0

| gate | value | z |
|---|---|---|
| `sat0_inlens` | 0.017 | +15.3 |
| `sat255_etd` | 0.002 | +13.9 |
| `bse_si_med` | 91.0 | -4.0 |
| `t2` | 71.5 | -2.8 |
| `height_px` | 2,316 | +2.5 |
| `sat0_bse` | 0.022 | +2.3 |
| `etd_si_med` | 113.0 | -2.2 |
| `sat0_etd` | 0.033 | +2.1 |
| `etd_c_med` | 80.0 | +1.6 |
| `noise_sd_bse` | 11.6 | +1.5 |

## 3. KPIs

| KPI | value ± SE | baseline median ± scale | z | class | shortlist |
|---|---|---|---|---|---|
| `si_frac_solid` | 20.25 % ± 2.29 % | 6.71 % ± 0.84 % | +16.10 | M | ✓ |
| `si_d50_aw_um` | 4.17 µm ± 0.167 | 3.19 ± 0.334 | +2.93 | R | ✓ |
| `si_d90_aw_um` | 8.43 µm ± 0.653 | 5.82 ± 0.776 | +3.37 | R | ✓ |
| `si_num_density_mm2` | 202,516 mm^-2 ± 14,199 | 26,433 ± 5,008 | +35.16 | M | ✓ |
| `si_solidity_med` | 0.732 ± 0.012 | 0.884 ± 0.044 | -3.45 | R | ✓ |
| `si_quadrat_cv` | 0.582 ± 0.039 | 0.749 ± 0.180 | -0.93 | M | ✓ |
| `pore_frac_deep` | 8.25 % ± 0.48 % | 12.11 % ± 1.14 % | -3.38 | W | ✓ |
| `pore_chord_h_mean_um` | 0.458 µm ± 0.016 | 0.644 ± 0.075 | -2.48 | R | ✓ |
| `pore_chord_v_mean_um` | 0.437 µm ± 0.015 | 0.540 ± 0.053 | -1.92 | R | ✓ |
| `aniso_pore_chord_ratio` | 1.05 ± 0.006 | 1.18 ± 0.059 | -2.22 | R | ✓ |
| `aniso_carbon_chord_ratio` | 1.05 ± 0.019 | 1.20 ± 0.060 | -2.44 | R | ✓ |
| `s2_len_pore_h_px` | 28.9 px ± 1.48 | 33.3 ± 5.00 | -0.89 | R | ✓ |
| `pore_percolating_frac_v` | 0.00 % ± 0.00 % | 0.00 % ± 0.50 % | +0.00 | M | ✓ |
| `si_contact_pore_frac` | 7.40 % ± 0.08 % | 1.71 % ± 0.70 % | +8.09 | R | ✓ |
| `si_contact_carbon_frac` | 92.60 % ± 0.08 % | 98.29 % ± 4.91 % | -1.16 | R | ✓ |
| `crack_density_um_per_mm2` | 204.8 µm/mm^2 ± 204.8 | 2,541 ± 1,707 | -1.37 | M | ✓ |
| `si_frac_total` | 18.58 % ± 2.18 % | 5.89 % ± 0.87 % | +14.58 | M |  |
| `si_bse_contrast` | 0.804 ± 0.024 | 1.73 ± 0.234 | -3.95 | M |  |
| `si_d10_aw_um` | 0.991 µm ± 0.086 | 1.31 ± 0.158 | -2.02 | M |  |
| `si_span` | 1.79 ± 0.153 | 1.41 ± 0.216 | +1.72 | M |  |
| `si_aspect_med` | 1.83 ± 0.071 | 1.70 ± 0.085 | +1.50 | R |  |
| `si_circularity_med` | 0.310 ± 0.008 | 0.604 ± 0.030 | -9.73 | R |  |
| `si_interior_texture` | 0.534 ± 0.019 | 0.733 ± 0.075 | -2.66 | M |  |
| `si_clark_evans_R` | 0.850 ± 0.022 | 0.754 ± 0.062 | +1.54 | M |  |
| `si_agglomerate_frac` | 93.33 % ± 1.11 % | 87.85 % ± 4.39 % | +1.25 | M |  |
| `pore_frac_slope_per_level` | 0.003 1/level ± n/a | 0.004 ± 0.001 | -1.50 | W |  |
| `ambiguous_frac` | 0.43 % ± 0.03 % | 0.62 % ± 0.50 % | -0.38 | W |  |
| `pore_chord_h_p90_um` | 1.15 µm ± 0.056 | 1.55 ± 0.185 | -2.16 | M |  |
| `pore_lt_d50_um` | 0.660 µm ± 0.031 | 0.694 ± 0.056 | -0.60 | M |  |
| `st_coherence` | 0.038 ± 0.004 | 0.122 ± 0.028 | -2.97 | R |  |
| `st_orientation_deg` | 13.5 deg ± 5.14 | 2.73 ± 3.42 | +3.14 | R |  |
| `s2_len_pore_v_px` | 22.3 px ± 0.895 | 20.5 ± 3.40 | +0.53 | R |  |
| `s2_len_si_px` | 77.4 px ± 3.31 | 66.7 ± 7.69 | +1.40 | R |  |
| `s2_integral_range_pore_px2` | 6,719 px^2 ± 810.8 | 6,172 ± 1,730 | +0.32 | M |  |
| `pore_percolating_frac_h` | 0.00 % ± 0.00 % | 0.00 % ± 0.50 % | +0.00 | M |  |
| `pore_euler_density_mm2` | 293,904 mm^-2 ± 14,720 | 141,234 ± 23,120 | +6.60 | M |  |
| `interface_pore_solid_um_per_um2` | 0.576 µm^-1 ± 0.018 | 0.651 ± 0.033 | -2.28 | M |  |
| `largest_void_ecd_um` | 6.46 µm ± 0.206 | 9.94 ± 3.33 | -1.04 | M |  |
| `hiZ_inclusion_count_mm2` | 0.000 mm^-2 ± 0.000 | 0.000 ± 0.000 | +0.00 | M |  |
| `vertical_pore_slope` | -0.135 pp/10µm ± 0.284 | 0.127 ± 0.580 | -0.45 | W |  |
| `vertical_si_slope` | 0.108 pp/10µm ± 1.07 | -0.164 ± 0.623 | +0.44 | W |  |
| `si_frac_se_imagerep` | 2.29 % ± 0.42 % | 1.01 % ± 0.50 % | +2.55 | R |  |
| `pore_frac_se_imagerep` | 0.56 % ± 0.06 % | 0.70 % ± 0.50 % | -0.28 | R |  |

## 4. Drivers (shortlist, sorted by |z|)

- **Si-candidate number density** (`si_num_density_mm2`, class M): 202,516 mm^-2 vs 26,433 mm^-2 ± 5,008, z = +35.16 (higher) — more Si-candidate particles per area (finer or more loaded)
- **Si-candidate loading** (`si_frac_solid`, class M): 20.25 % vs 6.71 % ± 0.84 %, z = +16.10 (higher) — more Si-candidate in the solid: higher capacity but more swelling, first-cycle loss and N/P shift
- **Si-candidate/pore contact** (`si_contact_pore_frac`, class R): 7.40 % vs 1.71 % ± 0.70 %, z = +8.09 (higher) — more Si-candidate surface facing pores: more room to swell, but more SEI
- **Si-candidate solidity** (`si_solidity_med`, class R): 0.732 vs 0.884 ± 0.044, z = -3.45 (lower) — more irregular/porous Si-candidate shapes (grade change?)
- **deep-pore fraction** (`pore_frac_deep`, class W): 8.25 % vs 12.11 % ± 1.14 %, z = -3.38 (lower) — less deep porosity (lower bound): denser coating, slower ionic transport; consistent with heavier calendering
- **Si-candidate D90** (`si_d90_aw_um`, class R): 8.43 µm vs 5.82 µm ± 0.776, z = +3.37 (higher) — more coarse Si-candidate particles in the tail: fracture/swelling hot-spots
- **Si-candidate D50** (`si_d50_aw_um`, class R): 4.17 µm vs 3.19 µm ± 0.334, z = +2.93 (higher) — coarser Si-candidate (area-weighted D50): higher fracture and swelling risk, more particle isolation
- **in-plane pore chord** (`pore_chord_h_mean_um`, class R): 0.458 µm vs 0.644 µm ± 0.075, z = -2.48 (lower) — narrower in-plane pores; consistent with heavier calendering

## 5. Batch assignment

| batch | probability | distance | signature match (σ) |
|---|---|---|---|
| 1 | 0.924 | 28.56 | — |
| 2 | 0.070 | 28.92 | — |
| 3 | 0.006 | 29.26 | — |

- Assigned: **Batch 1** — confident; strip-bootstrap stability 1.00; single-strip agreement 1.00
- Novelty flag: **YES**; in-distribution score 0.065
- Priors: uniform. T=4.0 and shrinkage Δ=0.0 chosen by leave-one-image-out log-loss on 30 labelled images; LOIO accuracy 0.33
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
| Through-plane electrolyte transport (2-conductivity D_eff) (higher = easier ionic transport) | B | 0.029 | 0.059 | 0.49× | 0.027 … 0.031 |
| Through-plane tortuosity index (higher = more tortuous) | B | 2.65 | 2.00 | 1.32× |  |
| In-plane / through-plane transport anisotropy (higher = more aligned (calendering)) | B | 0.958 | 1.21 | 0.79× |  |
| Limiting current proxy (higher = better fast charge) | B | n/a | n/a | 0.49× |  |
| Plating-risk ionic-resistance ratio (higher = more plating risk) | B | n/a | n/a | 2.03× |  |
| Electronic network conductance (Si insulating) (higher = better wired) | B | 0.166 | 0.298 | 0.56× |  |
| Swelling: particles with enough pore buffer (full lithiation) (higher = better) | B | 0.191 | 0.036 | 5.32× | 0.184 … 0.225 |
| Swelling: constraint index (area-weighted) (higher = more matrix displaced) | B | 0.861 | 0.812 | 1.06× |  |
| Swelling: 2-D thickness-change bound (higher = more electrode swelling) | B | 0.151 | 0.056 | 2.67× |  |
| First-cycle-loss proxy (exposed interface) (higher = more SEI) | B | 0.421 | 0.485 | 0.87× |  |
| Cycle-life index (heuristic) (higher = worse) | C | n/a | n/a | 2.06× |  |
| Energy density ratio [offset_baseline_to_0.30] | A | ε = 0.261 |  | 1.52× |  |
| PyBaMM 1C / 3C CC-charge capacity [offset_baseline_to_0.30] | B | 0.700 / 0.302 |  |  |  |
| PyBaMM min negative surface Δφ at separator, 3C (V; < 0 = plating indicator) [offset_baseline_to_0.30] | B | -0.223 |  |  |  |
| Energy density ratio [deep_pore_as_is] | A | ε = 0.083 |  | 1.51× |  |
| PyBaMM 1C / 3C CC-charge capacity [deep_pore_as_is] | B | 0.388 / 0.066 |  |  |  |
| PyBaMM min negative surface Δφ at separator, 3C (V; < 0 = plating indicator) [deep_pore_as_is] | B | -0.431 |  |  |  |

**Marked unstable** (sample ranking changes between conventions, brief 10.4): Q_CC_3C_over_Q_C10 rank: porosity offset vs as-is; Q_CC_1C_over_Q_C10 rank: porosity offset vs as-is. Read those rows as direction-free.

Checks: bounds ordered True, per-strip ordered True, flux imbalance 3.3e-11, downsampling audit ok True; uncertain pixels 3.6 %; assumptions hash f5d1bb3ceb. Grades: A = arithmetic on measured quantities; B = direction supported, level set by an assumption (ratios only); C = heuristic.

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

_Image flags: none; runtime 60.8 s; baseline version a381ba39._