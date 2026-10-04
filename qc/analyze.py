"""Per-sample orchestration: load -> segment -> gates -> KPIs (image + strips) -> cache.

Crack KPIs need an absolute Sato threshold that is fixed from the baseline images, so the Sato
response is cached per sample and crack KPIs are (re)finalised once the threshold is known.
"""
from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import gates as gates_mod
from . import segment as seg_mod
from . import tiles
from .io import SampleFiles, load_sample
from .kpis.anisotropy import structure_tensor_maps
from .kpis.context import Ctx
from .kpis.defects import crack_mask, curtaining_angle, sato_response
from .kpis.pores import local_thickness_map, pore_threshold_slope

CACHE_VERSION = 5
QC_DIR = Path(__file__).resolve().parent

# Source files each cached stage depends on: a stage is recomputed only when one of these (or the
# config) changes, so `qc regress` over 31 samples re-runs just the stages whose code changed.
STAGE_SOURCES = {
    "seg": ["io.py", "segment.py"],
    "sato": ["io.py", "segment.py", "kpis/defects.py"],
    "kpis": ["io.py", "segment.py", "gates.py", "tiles.py", "analyze.py", "kpis/*.py"],
}


# config keys used only after the per-image KPIs (statistics, decisions, validation) — changing them
# must not invalidate the image-level cache
_POST_KPI_KEYS = {"sim", "crack_threshold_frozen", "zones", "weights", "shrinkage_grid", "temperature_grid", "gate_caution_z",
                  "signature_min_effect", "signature_min_stability", "scale_rel_floor", "scale_abs_floor",
                  "validation_sensitivity_samples_per_batch", "manual_masks"}


def _cfg_core(cfg: dict) -> dict:
    return {k: v for k, v in cfg.items() if k not in _POST_KPI_KEYS}


def stage_version(stage: str, cfg: dict) -> str:
    h = hashlib.sha1(str(CACHE_VERSION).encode())
    for pat in STAGE_SOURCES[stage]:
        for f in sorted(QC_DIR.glob(pat)):
            h.update(f.read_bytes())
    h.update(json.dumps(_cfg_core(cfg), sort_keys=True, default=str).encode())
    return h.hexdigest()[:10]


def pipeline_version(cfg: dict) -> str:
    return stage_version("kpis", cfg)


@dataclass
class SampleResult:
    sample_id: str
    batch: str
    kpis: dict
    kpi_se: dict
    kpi_ci: dict
    strips: list
    gates: dict
    info: dict
    lists: dict = field(default_factory=dict)

    def to_json(self) -> dict:
        return dict(sample_id=self.sample_id, batch=self.batch, kpis=self.kpis, kpi_se=self.kpi_se,
                    kpi_ci=self.kpi_ci, strips=self.strips, gates=self.gates, info=self.info)


def _clean(o):
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def _unclean(o):
    if isinstance(o, dict):
        return {k: _unclean(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_unclean(v) for v in o]
    return np.nan if o is None else o


def build_ctx(sample, seg, cfg, crack_threshold=None, with_cracks=True, crack_resp=None) -> Ctx:
    g_noise = gates_mod.noise_sd(sample.bse, seg.carbon)
    ctx = Ctx(px_um=sample.px_um, bse=sample.bse, bse_s=seg.bse_s, pore=seg.pore, si=seg.si, carbon=seg.carbon,
              ambiguous=seg.ambiguous, si_labels=seg.si_labels, t1=seg.t1, t2=seg.t2, noise_sd=g_noise, cfg=cfg)
    ctx.lt_map = local_thickness_map(seg.pore, sample.px_um)
    # structure tensor on the 3-phase map (not raw intensity) so it is invariant to LUT/gamma/noise,
    # with sigma scaled to keep the same physical size if the pixel size differs from 25 nm
    phase_map = seg.carbon.astype(np.float32) + 2.0 * seg.si.astype(np.float32)
    st_sigma = cfg["structure_tensor_sigma_px"] * cfg["px_nm_target"] / sample.px_nm
    ctx.st = structure_tensor_maps(phase_map, st_sigma)
    if with_cracks:
        ctx.crack_resp = crack_resp if crack_resp is not None else sato_response(seg.bse_s, cfg["crack_sato_sigmas"])
        ctx.crack_threshold = crack_threshold
        src = sample.inlens if sample.inlens is not None else sample.bse
        lo, hi = cfg["curtaining_search_deg"]
        ang, strength = curtaining_angle(src, lo, hi)
        ctx.extras["curtain_strength"] = strength
        ctx.extras["curtain_measured_deg"] = ang
        # only exclude a direction when curtaining is actually visible as a spectral peak
        ctx.curtain_deg = ang if strength >= cfg["curtaining_min_strength"] else None
    return ctx


def _seg_to_npz(seg, path: Path, key: str):
    np.savez_compressed(path, key=key, pore=np.packbits(seg.pore), si=np.packbits(seg.si), amb=np.packbits(seg.ambiguous),
                        si_labels=seg.si_labels, shape=np.array(seg.pore.shape), t1=seg.t1, t2=seg.t2,
                        band=np.array(seg.band_thresholds, float), shading=seg.shading_flag,
                        etd_cut=np.nan if seg.etd_pore_cut is None else seg.etd_pore_cut)


def _seg_from_npz(path: Path, key: str, sample, cfg):
    if not path.exists():
        return None
    z = np.load(path)
    if str(z["key"]) != key:
        return None
    shape = tuple(z["shape"])
    n = shape[0] * shape[1]

    def unpack(a):
        return np.unpackbits(a)[:n].reshape(shape).astype(bool)

    pore, si, amb = unpack(z["pore"]), unpack(z["si"]), unpack(z["amb"])
    sig = cfg["smoothing_sigma_px"]
    bse_s = seg_mod.smooth(sample.bse, sig)
    etd_s = seg_mod.smooth(sample.etd, sig) if sample.etd is not None else None
    cut = float(z["etd_cut"])
    return seg_mod.Segmentation(bse_s, etd_s, float(z["t1"]), float(z["t2"]), pore, amb, si, ~pore & ~si,
                                z["si_labels"], z["band"].tolist(), bool(z["shading"]), None if np.isnan(cut) else cut)


def analyze_sample(files: SampleFiles, cfg: dict, cache_dir: Path | None = None, crack_threshold: float | None = None,
                   force: bool = False, t_shift=(0.0, 0.0), downsample: int = 1, save_cache: bool = True) -> SampleResult:
    """Cache layout out/cache/<sample_id>/: seg.npz, sato.npy(+sato.key), result.json, lists.npz, overlay.png.
    Each stage is keyed by (file SHA-1, stage code+config version)."""
    sdir = Path(cache_dir) / files.sample_id if cache_dir else None
    plain = t_shift == (0.0, 0.0) and downsample == 1
    use_cache = sdir is not None and plain
    fp = files.fingerprint()
    man = (cfg.get("manual_masks") or {}).get(files.sample_id)
    if man:  # a manual mask invalidates only its own sample's cache
        fp = fp + ":" + hashlib.sha1(json.dumps(man, sort_keys=True).encode()).hexdigest()[:6]
    kv = stage_version("kpis", cfg)
    if use_cache and not force:
        cached = load_cached(sdir, fp, kv)
        if cached is not None:
            if crack_threshold is not None and cached.info.get("crack_threshold") != crack_threshold:
                finalize_cracks(cached, sdir, crack_threshold, cfg)
            return cached
    t0 = time.time()
    sample = load_sample(files, cfg)
    if downsample > 1:
        from skimage.transform import downscale_local_mean

        for name in ("bse", "etd", "inlens"):
            a = getattr(sample, name)
            if a is not None:
                setattr(sample, name, np.clip(downscale_local_mean(a.astype(np.float32), (downsample, downsample)), 0, 255).astype(np.uint8))
        sample.px_nm *= downsample
    seg_key = f"{fp}:{stage_version('seg', cfg)}"
    seg = _seg_from_npz(sdir / "seg.npz", seg_key, sample, cfg) if use_cache and not force else None
    stages = {"seg": "cached" if seg is not None else "computed"}
    if seg is None:
        seg = seg_mod.segment(sample.bse, sample.etd, cfg, t_shift=t_shift, px_nm=sample.px_nm)
    t_seg = time.time() - t0
    g = gates_mod.compute_gates(sample, seg, cfg)
    g["etd_is_se"] = float(files.etd_is_se)
    sato_key = f"{fp}:{stage_version('sato', cfg)}"
    resp = None
    if use_cache and not force and (sdir / "sato.key").exists() and (sdir / "sato.key").read_text() == sato_key:
        resp = np.load(sdir / "sato.npy")
    stages["sato"] = "cached" if resp is not None else "computed"
    ctx = build_ctx(sample, seg, cfg, crack_threshold, crack_resp=resp)
    resp = ctx.crack_resp
    kpis, lists = tiles.compute_kpis(ctx)
    kpis["pore_frac_slope_per_level"] = pore_threshold_slope(seg.bse_s, seg.etd_s, seg.t1, seg.t2, cfg)
    strips = tiles.strip_kpis(ctx, cfg["strips"])
    keys = [k for k in kpis if isinstance(kpis[k], (int, float, np.floating, np.integer))]
    se = tiles.strip_se(strips, keys)
    ci = tiles.bootstrap_ci(strips, keys, cfg["bootstrap_n"], cfg["seed"])
    info = dict(
        fingerprint=fp, kpi_version=kv, cache_version=CACHE_VERSION, stages=stages, flags=sample.flags, crop=sample.crop,
        shape=list(sample.bse.shape), raw_shape=list(sample.raw_shape), px_nm=sample.px_nm, px_nm_assumed=sample.px_nm_assumed,
        etd_is_se=files.etd_is_se, has_batch_prefix=files.has_batch_prefix, complete=files.complete, t_shift=list(t_shift),
        downsample=downsample, crack_threshold=crack_threshold,
        sato_p995=float(np.percentile(resp[::2, ::2].astype(np.float32), cfg["crack_threshold_default_quantile"])) if resp is not None else None,
        curtain_deg=ctx.curtain_deg, curtain_measured_deg=ctx.extras.get("curtain_measured_deg"),
        curtain_strength=ctx.extras.get("curtain_strength"),
        etd_pore_cut=seg.etd_pore_cut, band_thresholds=seg.band_thresholds,
        seconds=dict(segment=round(t_seg, 1), total=round(time.time() - t0, 1)),
        paths={k: str(v) for k, v in files.paths.items()},
    )
    res = SampleResult(files.sample_id, files.batch, kpis, se, ci, strips, g, info, lists)
    if use_cache and save_cache:
        sdir.mkdir(parents=True, exist_ok=True)
        if stages["seg"] == "computed":
            _seg_to_npz(seg, sdir / "seg.npz", seg_key)
        if stages["sato"] == "computed":
            np.save(sdir / "sato.npy", resp)
            (sdir / "sato.key").write_text(sato_key)
        save_overlay(sample.bse, seg, ctx.extras.get("crack_skeleton"), sdir / "overlay.png")
        save_cached(res, sdir)
    return res


def save_overlay(bse, seg, crack_skel, path: Path, factor: int = 4):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rgb = seg_mod.overlay_rgb(bse, seg, factor)
    if crack_skel is not None:
        from scipy import ndimage as ndi

        ck = ndi.binary_dilation(crack_skel, iterations=2)[::factor, ::factor]
        rgb[ck] = [230, 20, 60]
    plt.imsave(path, rgb)


def save_cached(res: SampleResult, sdir: Path):
    with open(sdir / "result.json", "w") as f:
        json.dump(_clean(res.to_json()), f)
    np.savez_compressed(sdir / "lists.npz", **{k: np.asarray(v, np.float32) for k, v in res.lists.items()})


def load_cached(sdir: Path, fingerprint: str | None = None, kpi_version: str | None = None) -> SampleResult | None:
    p = sdir / "result.json"
    if not p.exists():
        return None
    d = _unclean(json.loads(p.read_text()))
    if d["info"].get("cache_version") != CACHE_VERSION:
        return None
    if fingerprint and d["info"].get("fingerprint") != fingerprint:
        return None
    if kpi_version and d["info"].get("kpi_version") != kpi_version:
        return None
    lists = {}
    if (sdir / "lists.npz").exists():
        with np.load(sdir / "lists.npz") as z:
            lists = {k: z[k] for k in z.files}
    return SampleResult(d["sample_id"], d["batch"], d["kpis"], d["kpi_se"], d["kpi_ci"], d["strips"], d["gates"], d["info"], lists)


def finalize_cracks(res: SampleResult, sdir: Path, thr: float, cfg: dict):
    """Recompute crack KPIs (image + strips) from the cached Sato response at threshold `thr`."""
    resp = np.load(sdir / "sato.npy")
    area_mm2 = res.kpis["area_mm2"]
    px_um = res.info["px_nm"] / 1000
    curtain = res.info.get("curtain_deg")
    kept, n, length = crack_mask(resp, thr, cfg, curtain, 25.0 / res.info["px_nm"])
    res.kpis["crack_density_um_per_mm2"] = length * px_um / area_mm2
    res.kpis["crack_count"] = n
    for s, sl in zip(res.strips, tiles.strip_slices(resp.shape[1], len(res.strips))):
        _, n_s, len_s = crack_mask(resp[:, sl], thr, cfg, curtain, 25.0 / res.info["px_nm"])
        s["crack_density_um_per_mm2"] = len_s * px_um / s["area_mm2"]
        s["crack_count"] = n_s
    keys = ["crack_density_um_per_mm2", "crack_count"]
    res.kpi_se.update(tiles.strip_se(res.strips, keys))
    res.kpi_ci.update(tiles.bootstrap_ci(res.strips, keys, cfg["bootstrap_n"], cfg["seed"]))
    res.info["crack_threshold"] = thr
    save_cached(res, sdir)
    ov = sdir / "overlay.png"
    if ov.exists():
        import matplotlib.pyplot as plt
        from scipy import ndimage as ndi

        rgb = (plt.imread(ov)[..., :3] * 255).astype(np.uint8)
        f = int(round(resp.shape[1] / rgb.shape[1])) or 1
        ck = ndi.binary_dilation(kept, iterations=2)[::f, ::f][: rgb.shape[0], : rgb.shape[1]]
        rgb[: ck.shape[0], : ck.shape[1]][ck] = [230, 20, 60]
        plt.imsave(ov, rgb)
