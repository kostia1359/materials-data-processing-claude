# QC report — img_0grcilhi

> **Acquisition drift or artefact suspected — compare with caution.** sat0_inlens z=+10.0; shading_flag

## 1. Verdict

**img_0grcilhi is most consistent with Batch 3** (probability 0.58; Batch 1: 0.11; Batch 2: 0.31; assignment stable in 100 % of strip-bootstrap resamples; moderately confident). **Relative to the Batch-3 baseline it is within the baseline distribution — verdict ACCEPT**: 0 of 16 KPIs exceed the investigate zone (8 robust σ; reject at 12.8). Largest shifts: deep-pore fraction 16.14 % vs 12.06 % ± 0.75 % (z = +5.5); in-plane pore chord 0.898 µm vs 0.630 µm ± 0.056 (z = +4.8); through-plane pore chord 0.738 µm vs 0.535 µm ± 0.043 (z = +4.7). Reads as: more deep porosity (lower bound): better electrolyte transport, lower energy density; wider in-plane pores. Baseline has 16 image(s), so the smallest reportable rank p-value is 0.059; this verdict is an effect-size judgement, not a significance test.

## 2. Acquisition gates

Status: **caution** — flags: sat0_inlens z=+10.0, shading_flag

| gate | value | z |
|---|---|---|
| `sat0_inlens` | 0.011 | +10.0 |
| `bse_row_gradient_pct` | 18.0 | +3.8 |
| `height_px` | 1,904 | -1.6 |
| `sat0_bse` | 0.016 | +1.5 |
| `blur_bse` | 0.221 | +1.5 |
| `sat0_etd` | 0.025 | +1.4 |
| `etd_si_med` | 119.0 | -1.3 |
| `etd_c_med` | 71.0 | -0.9 |
| `t1` | 38.2 | -0.7 |
| `reg_shift_etd_px` | 0.200 | -0.7 |

## 3. KPIs

| KPI | value ± SE | baseline median ± scale | z | class | shortlist |
|---|---|---|---|---|---|
| `si_frac_solid` | 6.29 % ± 1.02 % | 6.74 % ± 0.85 % | -0.53 | M | ✓ |
| `si_d50_aw_um` | 3.35 µm ± 0.349 | 3.18 ± 0.347 | +0.49 | R | ✓ |
| `si_d90_aw_um` | 6.09 µm ± 0.537 | 5.77 ± 0.818 | +0.38 | R | ✓ |
| `si_num_density_mm2` | 23,055 mm^-2 ± 1,817 | 26,599 ± 4,836 | -0.73 | M | ✓ |
| `si_solidity_med` | 0.882 ± 0.019 | 0.885 ± 0.044 | -0.05 | R | ✓ |
| `si_quadrat_cv` | 0.857 ± 0.113 | 0.738 ± 0.175 | +0.67 | M | ✓ |
| `pore_frac_deep` | 16.14 % ± 0.99 % | 12.06 % ± 0.75 % | +5.48 | W | ✓ |
| `pore_chord_h_mean_um` | 0.898 µm ± 0.063 | 0.630 ± 0.056 | +4.81 | R | ✓ |
| `pore_chord_v_mean_um` | 0.738 µm ± 0.065 | 0.535 ± 0.043 | +4.71 | R | ✓ |
| `aniso_pore_chord_ratio` | 1.22 ± 0.021 | 1.18 ± 0.059 | +0.71 | R | ✓ |
| `aniso_carbon_chord_ratio` | 1.20 ± 0.012 | 1.20 ± 0.060 | -0.08 | R | ✓ |
| `s2_len_pore_h_px` | 53.2 px ± 5.16 | 32.6 ± 4.69 | +4.40 | R | ✓ |
| `pore_percolating_frac_v` | 0.00 % ± 0.00 % | 0.00 % ± 0.50 % | +0.00 | M | ✓ |
| `si_contact_pore_frac` | 1.27 % ± 0.17 % | 1.73 % ± 0.79 % | -0.58 | R | ✓ |
| `si_contact_carbon_frac` | 98.73 % ± 0.17 % | 98.27 % ± 4.91 % | +0.09 | R | ✓ |
| `crack_density_um_per_mm2` | 297.2 µm/mm^2 ± 297.1 | 2,696 ± 1,711 | -1.40 | M | ✓ |
| `si_frac_total` | 5.28 % ± 0.86 % | 5.92 % ± 0.76 % | -0.85 | M |  |
| `si_bse_contrast` | 1.57 ± 0.047 | 1.73 ± 0.215 | -0.75 | M |  |
| `si_d10_aw_um` | 1.35 µm ± 0.146 | 1.29 ± 0.144 | +0.41 | M |  |
| `si_span` | 1.41 ± 0.136 | 1.47 ± 0.220 | -0.24 | M |  |
| `si_aspect_med` | 1.76 ± 0.060 | 1.70 ± 0.085 | +0.78 | R |  |
| `si_circularity_med` | 0.599 ± 0.034 | 0.605 ± 0.031 | -0.21 | R |  |
| `si_interior_texture` | 0.783 ± 0.028 | 0.725 ± 0.066 | +0.87 | M |  |
| `si_clark_evans_R` | 0.833 ± 0.066 | 0.749 ± 0.059 | +1.41 | M |  |
| `si_agglomerate_frac` | 89.84 % ± 3.26 % | 87.59 % ± 4.38 % | +0.51 | M |  |
| `pore_frac_slope_per_level` | 0.004 1/level ± n/a | 0.004 ± 0.001 | +0.64 | W |  |
| `ambiguous_frac` | 0.97 % ± 0.05 % | 0.62 % ± 0.50 % | +0.71 | W |  |
| `pore_chord_h_p90_um` | 2.23 µm ± 0.235 | 1.51 ± 0.148 | +4.81 | M |  |
| `pore_lt_d50_um` | 1.09 µm ± 0.151 | 0.683 ± 0.064 | +6.42 | M |  |
| `st_coherence` | 0.133 ± 0.009 | 0.121 ± 0.031 | +0.36 | R |  |
| `st_orientation_deg` | 9.73 deg ± 4.52 | 2.46 ± 3.34 | +2.17 | R |  |
| `s2_len_pore_v_px` | 29.8 px ± 3.25 | 20.5 ± 3.39 | +2.77 | R |  |
| `s2_len_si_px` | 66.2 px ± 6.01 | 67.4 ± 8.17 | -0.15 | R |  |
| `s2_integral_range_pore_px2` | 15,828 px^2 ± 2,949 | 6,132 ± 1,534 | +6.32 | M |  |
| `pore_percolating_frac_h` | 0.00 % ± 10.24 % | 0.00 % ± 0.50 % | +0.00 | M |  |
| `pore_euler_density_mm2` | 78,530 mm^-2 ± 16,186 | 141,240 ± 20,592 | -3.05 | M |  |
| `interface_pore_solid_um_per_um2` | 0.624 µm^-1 ± 0.021 | 0.652 ± 0.033 | -0.86 | M |  |
| `largest_void_ecd_um` | 18.7 µm ± 1.35 | 9.38 ± 3.00 | +3.11 | M |  |
| `hiZ_inclusion_count_mm2` | 0.000 mm^-2 ± 0.000 | 0.000 ± 0.000 | +0.00 | M |  |
| `vertical_pore_slope` | 1.58 pp/10µm ± 1.01 | 0.083 ± 0.541 | +2.76 | W |  |
| `vertical_si_slope` | -1.69 pp/10µm ± 0.604 | -0.068 ± 0.632 | -2.57 | W |  |
| `si_frac_se_imagerep` | 1.07 % ± 0.31 % | 1.01 % ± 0.50 % | +0.14 | R |  |
| `pore_frac_se_imagerep` | 1.27 % ± 0.33 % | 0.70 % ± 0.50 % | +1.14 | R |  |

## 4. Drivers (shortlist, sorted by |z|)

- **deep-pore fraction** (`pore_frac_deep`, class W): 16.14 % vs 12.06 % ± 0.75 %, z = +5.48 (higher) — more deep porosity (lower bound): better electrolyte transport, lower energy density
- **in-plane pore chord** (`pore_chord_h_mean_um`, class R): 0.898 µm vs 0.630 µm ± 0.056, z = +4.81 (higher) — wider in-plane pores
- **through-plane pore chord** (`pore_chord_v_mean_um`, class R): 0.738 µm vs 0.535 µm ± 0.043, z = +4.71 (higher) — taller through-plane pores: easier through-plane transport
- **pore correlation length** (`s2_len_pore_h_px`, class R): 53.2 px vs 32.6 px ± 4.69, z = +4.40 (higher) — pore network stretched in-plane (longer correlation length): flatter, more aligned pores; consistent with heavier calendering
- **crack density** (`crack_density_um_per_mm2`, class M): 297.2 µm/mm^2 vs 2,696 µm/mm^2 ± 1,711, z = -1.40 (lower) — fewer cracks
- **Si-candidate number density** (`si_num_density_mm2`, class M): 23,055 mm^-2 vs 26,599 mm^-2 ± 4,836, z = -0.73 (lower) — fewer Si-candidate particles per area (coarser or less loaded)
- **pore anisotropy** (`aniso_pore_chord_ratio`, class R): 1.22 vs 1.18 ± 0.059, z = +0.71 (higher) — pores more aligned parallel to the collector: higher through-plane tortuosity, slower fast charge; consistent with heavier calendering
- **Si-candidate patchiness** (`si_quadrat_cv`, class M): 0.857 vs 0.738 ± 0.175, z = +0.67 (higher) — Si-candidate less evenly distributed (mixing/dispersion issue): local swelling hot-spots

## 5. Batch assignment

| batch | probability | distance | signature match (σ) |
|---|---|---|---|
| 1 | 0.111 | 9.93 | — |
| 2 | 0.313 | 9.06 | — |
| 3 | 0.576 | 8.50 | — |

- Assigned: **Batch 3** — moderately confident; strip-bootstrap stability 1.00; single-strip agreement 1.00
- Novelty flag: **no**; in-distribution score 0.129
- Priors: uniform. T=8.0 and shrinkage Δ=0.0 chosen by leave-one-image-out log-loss on 30 labelled images; LOIO accuracy 0.37
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
| Through-plane electrolyte transport (2-conductivity D_eff) (higher = easier ionic transport) | B | 0.067 | 0.058 | 1.14× | 0.059 … 0.074 |
| Through-plane tortuosity index (higher = more tortuous) | B | 2.52 | 2.00 | 1.26× |  |
| In-plane / through-plane transport anisotropy (higher = more aligned (calendering)) | B | 1.42 | 1.21 | 1.17× |  |
| Limiting current proxy (higher = better fast charge) | B | n/a | n/a | 1.14× |  |
| Plating-risk ionic-resistance ratio (higher = more plating risk) | B | n/a | n/a | 0.88× |  |
| Electronic network conductance (Si insulating) (higher = better wired) | B | 0.142 | 0.302 | 0.47× |  |
| Swelling: particles with enough pore buffer (full lithiation) (higher = better) | B | 0.036 | 0.035 | 1.03× | 0.036 … 0.051 |
| Swelling: constraint index (area-weighted) (higher = more matrix displaced) | B | 0.801 | 0.813 | 0.99× |  |
| Swelling: 2-D thickness-change bound (higher = more electrode swelling) | B | 0.046 | 0.057 | 0.82× |  |
| First-cycle-loss proxy (exposed interface) (higher = more SEI) | B | 0.491 | 0.485 | 1.01× |  |
| Cycle-life index (heuristic) (higher = worse) | C | n/a | n/a | 0.96× |  |
| Energy density ratio [offset_baseline_to_0.30] | A | ε = 0.341 |  | 0.928× |  |
| PyBaMM 1C / 3C CC-charge capacity [offset_baseline_to_0.30] | B | 0.689 / 0.292 |  |  |  |
| PyBaMM min negative surface Δφ at separator, 3C (V; < 0 = plating indicator) [offset_baseline_to_0.30] | B | -0.111 |  |  |  |
| Energy density ratio [deep_pore_as_is] | A | ε = 0.161 |  | 0.940× |  |
| PyBaMM 1C / 3C CC-charge capacity [deep_pore_as_is] | B | 0.639 / 0.151 |  |  |  |
| PyBaMM min negative surface Δφ at separator, 3C (V; < 0 = plating indicator) [deep_pore_as_is] | B | -0.355 |  |  |  |

**Marked unstable** (sample ranking changes between conventions, brief 10.4): Q_CC_3C_over_Q_C10 rank: porosity offset vs as-is; Q_CC_1C_over_Q_C10 rank: porosity offset vs as-is. Read those rows as direction-free.

Checks: bounds ordered True, per-strip ordered True, flux imbalance 2.0e-11, downsampling audit ok True; uncertain pixels 7.2 %; assumptions hash f5d1bb3ceb. Grades: A = arithmetic on measured quantities; B = direction supported, level set by an assumption (ratios only); C = heuristic.

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

_Image flags: dropped_rgb_mismatch_cols:2; runtime 66.3 s; baseline version 64639794._