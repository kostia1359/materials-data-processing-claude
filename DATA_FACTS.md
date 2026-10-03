# DATA_FACTS

Status on 2026-10-03: **only the starter triple (`Batch_3_img_0grcilhi`, batch 3) has been verified pixel by pixel**.
The remaining 90 files could not be downloaded into the build container (`drive.google.com` is blocked by its
network policy). The inventory below comes from a Drive listing; run `python -m qc fetch` then
`python -m qc quicklook --data data` to re-verify Section 1 of the brief on every file
(`out/quicklook/facts.md` is the per-sample table this file should be updated from).

## Inventory (Drive listing, 2026-10-03) — `data/drive_manifest.csv`

| | Batch 1 | Batch 2 | Batch 3 (baseline) | Total |
|---|---|---|---|---|
| samples (triples) | 7 | 7 | 17 | **31** |
| files | 21 | 21 | 51 | 93 (+ a stray `.DS_Store`) |

- File size 12.8–23.2 MB each; all `.tif`.
- **Naming is not uniform**: 7 samples carry a `Batch_N_` prefix (`Batch_1_img_4ih2ggld_BSE.tif`), the rest do not (`img_ffwubibz_BSE.tif`).
  The loader strips any prefix before `img_`.
- **4 samples have an `SE` image instead of `ETD`**: batch 2 `rxax5ozo`; batch 3 `vc2whyaq`, `x77cy643`, `utfgcjfa`.
  Treated as the ETD role and flagged `etd_is_SE_detector` (see DECISIONS.md). If SE-vs-ETD acquisition
  differences show up in the gates, these 4 must be checked separately in the shadow-classifier analysis.

## Verified on the starter triple (`img_0grcilhi`, batch 3)

| Fact (brief §1) | Brief | Measured | Consequence |
|---|---|---|---|
| Container | 8-bit RGB, R=G=B, LZW, written by tifffile | ✔ 8-bit RGB, `Software=tifffile.py`; R=B everywhere, **G≠R on the last 2 columns** in all three detectors | 2 right-edge columns cropped (6998 px wide) |
| Size | 7000 px wide, height varies | ✔ 7000 × 1904 | per-area normalisation everywhere |
| Pixel size | 25.0 nm from X/YResolution | ✔ 25.0005 nm/px (126997216/125 px per inch) | no resampling |
| Data bar | none | ✔ none | — |
| Co-registration | ≤ 0.1 px | ✔ BSE↔ETD 0.2 px, BSE↔InLens 0.1 px (phase correlation on gradient magnitude, 1024² centre crop) | per-pixel triplets valid |
| Comb histograms | 95–211 levels | BSE **177**, ETD 224, InLens 226 distinct levels; BSE comb step 1 inside 20–200 | smoothing before thresholding; no 256-bin features |
| BSE phase levels | graphite ~56–59, Si ~108–115, pore ~11–18 | carbon median **58**, Si-candidate median **115**, Otsu t1/t2 **38.2 / 81.6** | consistent with brief |
| ETD / InLens | not reproducible between samples | ETD carbon 71 / Si 119; InLens carbon 71 / Si 156 | gates only |
| Vertical shading | one sample darkens 17 % top→bottom | **this sample: BSE row mean 59.4 → 49.4 (gradient 18 %)**; per-band t2 drifts 84 → 68 | `shading_flag` set; logged as covariate |
| Edge band | ~90 nodular rows on one sample | none detected here | rule kept; `manual_masks` override |
| Curtaining | 60–70° stripes in ETD/InLens | InLens angular Fourier peak only 1.6× median (at 51.5°): **no clear curtaining in this sample** | crack-direction exclusion only when peak ≥ 3× |
| Microstructure | flakes 10–30 µm, blocky bright 2–12 µm, open pores | ✔ (overlay `examples/img_0grcilhi/overlay.png`); Si-candidate area-weighted D50 3.3 µm, D90 6.1 µm; 113 sized particles | — |
| Phases | 3 separable classes | deep pore 16.1 %, Si-candidate 5.3 % of frame (6.3 % of solid), ambiguous 1.0 % | exactly 3 classes |
| Collector / free surface | none in frame | ✔ none visible | thickness/roughness/delamination not measurable |

Nothing found on the starter triple contradicts the brief. Re-verify per batch once all files are present.
