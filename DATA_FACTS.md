# DATA_FACTS

All **31 samples / 93 TIFFs** from the Drive folder `data_all` were downloaded on 2026-10-04 and verified
byte-exact against `data/drive_manifest.csv`. Per-sample tables: `out/manifest.csv`, `out/kpis.csv` (gates
included) and `out/quicklook/facts.csv`. Composites and histograms for every sample are in `out/quicklook/*.png`.

## Download log

| route | files | bytes | status |
|---|---|---|---|
| gdown by file id (`qc fetch`, drive.google.com) | 93 | 1.73 GB | 93/93 present, size = `size_bytes`, every file opens with tifffile |

A Drive re-list on 2026-10-04 found no files added since the 2026-10-03 manifest. `data/` (the ingested set)
was filled by hard links from the download folder in the Section 0.3 order: starter → Batch_3 → Batch_1 → Batch_2.
The manifest and the Drive listing agree on all 31 samples; no extra or missing files.

## Inventory and naming

| batch | samples | width (px) | height (px) | px size (nm) | `Batch_N_` prefix | `_SE.tif` instead of `_ETD.tif` |
|---|---|---|---|---|---|---|
| 1 | 7 | 6960–7000 | 1780–2316 | 25.000 | 2 | 0 |
| 2 | 7 | 7000 | 2048–2272 | 25.000 | 2 | 1 (`rxax5ozo`) |
| 3 (baseline) | 17 | 6996–7000 | 1612–2272 | 25.000 | 3 | 3 (`utfgcjfa`, `vc2whyaq`, `x77cy643`) |

- **Width is not always 7000 px** (6960 and 6996 occur); every KPI is per area, so this is harmless.
- Pixel size is 24.9992–25.0005 nm/px from X/YResolution on every file, so nothing is resampled.
- No filename prefix disagrees with its folder.

## Brief §1 facts, re-verified on all 31

| fact | result |
|---|---|
| 8-bit RGB, R=G=B, LZW, tifffile | ✔ all. G≠R edge columns (1–3 px) on 13 samples, dropped automatically |
| Co-registration | ✔ BSE↔ETD ≤ 0.32 px, BSE↔InLens ≤ 0.91 px (flag at 2 px never fires) |
| Comb histograms | ✔ 97–218 distinct BSE levels; comb step 1–2 inside 20–200 |
| BSE phase levels | ✔ carbon median 51–63, Si-candidate 108–116 in baseline; **two Batch_1 samples have a darker bright phase (≈ 86–90, see below)** |
| No data bar / no collector / open pores | ✔ visually on all quicklooks and overlays |
| Edge bands | rule never fired; checked by eye on the 2088-px samples (`9luzk4jm` bottom is cleanly sectioned); `kbdh4tri` has ~40 rough top rows (≈ 0.1 pp extra Si, logged in DECISIONS.md) |
| Curtaining | not detected: the InLens 50–80° Fourier peak is 1.1–2.5× the median angular power on all 31 images, below the 3× presence test, so no crack direction is excluded |

## Per-batch acquisition statistics (medians) — **session effects are real**

| gate | Batch 1 | Batch 2 | Batch 3 |
|---|---|---|---|
| BSE noise SD (grey levels) | 11.10 | 10.74 | 9.92 |
| BSE noise range | 10.66–11.66 | 10.15–11.22 | **7.32–10.99** |
| BSE saturated at 0 (fraction) | 0.022 | 0.024 | 0.006 |
| InLens carbon median | 110 | 103 | 75 |
| ETD carbon median | 80 | 78 | 73 |
| BSE row gradient top→bottom (%) | −10.5 | −6.3 | +1.7 |
| Otsu t1 / t2 | 36.3 / 85.5 | 36.7 / 82.2 | 40.5 / 83.0 |
| BSE carbon / Si-candidate median | 58 / 115 | 55 / 112 | 59 / 115 |

**The batches are almost perfectly separated by acquisition noise.** Only 6 images (1 / 3 / 2) fall in the noise
range shared by all three batches. Any KPI that responds to noise will "separate batches" for the wrong reason.
This is exactly what the crack KPI did until it was noise-normalised (DECISIONS.md, VALIDATION.md §4).

## ETD vs SE label check (brief §0.2)

The three `_SE` baseline samples sit inside the ETD-labelled baseline range on every channel statistic:
ETD carbon median 73/78/73 (ETD range 70–85), Si-candidate 128/129/128 (111–131), saturation 0.013–0.016 (0–0.038),
deep-pore fraction 0.120–0.122 (0.103–0.177), ambiguous fraction 0.006–0.007 (0.003–0.010). There is no evidence of
a different detector setting. The label is kept as covariate `etd_is_se`.

## Material observations

- **Batch_1 `4ih2ggld` and `5n1q8atc`** carry a separate bright BSE mode at grey ≈ 86–90 (baseline Si-candidate ≈ 108–116) with ~3× the
  baseline bright-phase mass: Si-candidate 18.6 % / 20.2 % of solids, `si_bse_contrast` 0.73 / 0.80 (baseline 1.66 ± 0.23).
  This is a lower-Z bright phase, consistent with SiOx or Si–C composite (no EDS, so it can't be confirmed). The other five Batch_1 samples look
  baseline-like on Si loading.
- The deep-pore phase never spans the frame (2-D percolation), on any sample (VALIDATION.md §8).
