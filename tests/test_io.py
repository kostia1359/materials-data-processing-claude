from pathlib import Path

import numpy as np

from qc.config import load_config
from qc.io import FILE_RE, discover, load_sample
from qc.synth import SynthParams, make_sample, write_triple


def test_filename_regex_variants():
    for name, sid, det in [
        ("img_abc123_BSE.tif", "abc123", "BSE"),
        ("Batch_3_img_0grcilhi_Inlens.tif", "0grcilhi", "Inlens"),
        ("img_rxax5ozo_SE.tif", "rxax5ozo", "SE"),
        ("IMG_X1_etd.TIFF", "X1", "etd"),
    ]:
        m = FILE_RE.match(name)
        assert m and m.group("id") == sid and m.group("det") == det


def test_discover_and_load(tmp_path: Path):
    imgs, _ = make_sample(SynthParams(width=1400, height=500, seed=1))
    # simulate the green edge column artefact on BSE
    write_triple(tmp_path / "Batch_2", "t1", imgs)
    se_dir = tmp_path / "Batch_3"
    write_triple(se_dir, "t2", imgs)
    (se_dir / "img_t2_ETD.tif").rename(se_dir / "img_t2_SE.tif")
    found = {s.sample_id: s for s in discover(tmp_path)}
    assert found["img_t1"].batch == "2" and found["img_t1"].complete
    assert found["img_t2"].etd_is_se and found["img_t2"].complete
    cfg = load_config()
    L = load_sample(found["img_t1"], cfg)
    assert L.bse.shape == (500, 1400)
    assert abs(L.px_nm - 25.0) < 0.01 and not L.px_nm_assumed
    assert L.bse.dtype == np.uint8
