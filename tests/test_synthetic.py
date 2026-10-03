"""Synthetic ground-truth recovery (brief 6.1). Smaller frames than the real 7000 px images keep this fast."""
import numpy as np
import pytest

from qc.config import load_config
from qc.synth import SynthParams, make_sample, perturb, truth_values
from qc.validate import light_kpis

CFG = load_config()


@pytest.fixture(scope="module")
def sample():
    imgs, masks = make_sample(SynthParams(width=3500, height=1000, seed=11))
    return imgs, truth_values(masks), light_kpis(imgs, CFG)


def test_fractions(sample):
    _, t, k = sample
    assert abs(k["si_frac_total"] - t["si_frac_total"]) <= 0.005
    assert abs(k["pore_frac_deep"] - t["pore_frac_deep"]) <= 0.01


def test_si_d50(sample):
    _, t, k = sample
    assert abs(k["si_d50_aw_um"] - t["si_d50_aw_um"]) / t["si_d50_aw_um"] <= 0.10


def test_anisotropy():
    from qc.kpis.pores import chords
    from qc.segment import segment

    imgs, masks = make_sample(SynthParams(width=3500, height=1000, seed=12))
    t = truth_values(masks)
    seg = segment(imgs["BSE"], imgs["ETD"], CFG)
    ratio = chords(seg.pore, 1).mean() / chords(seg.pore, 0).mean()
    assert abs(ratio - t["aniso_pore_chord_ratio"]) / t["aniso_pore_chord_ratio"] <= 0.10


def test_clark_evans_clustered_vs_poisson():
    pois = np.mean([light_kpis(make_sample(SynthParams(width=7000, height=2000, seed=s, si_d50_um=1.0, si_frac=0.015))[0], CFG)["si_clark_evans_R"]
                    for s in (901, 902, 903)])
    clus = light_kpis(make_sample(SynthParams(width=7000, height=2000, seed=901, clustered=True, si_d50_um=1.0, si_frac=0.015))[0], CFG)["si_clark_evans_R"]
    assert abs(pois - 1) <= 0.1
    assert clus < 0.8


def test_imagerep_se_matches_replicates():
    fr, se = [], []
    for i in range(10):
        imgs, _ = make_sample(SynthParams(width=3500, height=1000, seed=700 + i, poisson_count=True))
        k = light_kpis(imgs, CFG)
        fr.append(k["si_frac_total"])
        se.append(k["si_frac_se_imagerep"])
    ratio = np.mean(se) / np.std(fr, ddof=1)
    assert 0.5 <= ratio <= 2.0


@pytest.mark.parametrize("kw", [dict(brightness=1.15), dict(contrast=0.85), dict(gamma=1.2), dict(noise=5)])
def test_si_kpis_robust_to_intensity_changes(sample, kw):
    imgs, _, k0 = sample
    k1 = light_kpis({d: perturb(a, **kw) for d, a in imgs.items()}, CFG)
    # Si fraction and size: the Z-contrast step is large, so these must barely move
    assert abs(k1["si_frac_total"] - k0["si_frac_total"]) < 0.004
    assert abs(k1["si_d50_aw_um"] - k0["si_d50_aw_um"]) / k0["si_d50_aw_um"] < 0.05
