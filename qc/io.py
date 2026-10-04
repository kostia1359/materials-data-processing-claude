"""Discovery, manifest, loading and image hygiene.

Every downstream step works on a *rectangular* valid region: bad columns/rows (channel
disagreement), detected edge bands and manual masks are all cropped away here, so the
valid-area mask is the full array of what is returned.
"""
from __future__ import annotations

import hashlib
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

# Files in the Drive folder sometimes carry a "Batch_N_" prefix and four samples have an
# "SE" image where the others have "ETD" (both are secondary-electron detectors).
FILE_RE = re.compile(
    r"^(?:Batch_(?P<pb>\d+)_)?img_(?P<id>[a-z0-9]+)_(?P<det>BSE|ETD|SE|Inlens)\.tiff?$", re.IGNORECASE
)
BATCH_RE = re.compile(r"batch[\s_-]*(\d)", re.IGNORECASE)
DETECTORS = ("BSE", "ETD", "Inlens")


@dataclass
class SampleFiles:
    sample_id: str
    batch: str  # "1" | "2" | "3" | "unknown"
    paths: dict = field(default_factory=dict)  # det -> Path ; det in BSE/ETD/Inlens
    etd_is_se: bool = False
    has_batch_prefix: bool = False

    @property
    def etd_label_in_filename(self) -> str:
        return "SE" if self.etd_is_se else ("ETD" if "ETD" in self.paths else "")

    @property
    def complete(self) -> bool:
        return all(d in self.paths for d in DETECTORS)

    def fingerprint(self) -> str:
        """SHA-1 over the file contents of all channels (stable across copies, renames and mtimes)."""
        h = hashlib.sha1()
        for d in sorted(self.paths):
            h.update(d.encode())
            h.update(file_sha1(self.paths[d]).encode())
        return h.hexdigest()[:12]


_SHA_CACHE: dict = {}


def file_sha1(path) -> str:
    p = Path(path)
    st = p.stat()
    key = (str(p.resolve()), st.st_size, st.st_mtime_ns)
    if key not in _SHA_CACHE:
        h = hashlib.sha1()
        with open(p, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 22), b""):
                h.update(chunk)
        _SHA_CACHE[key] = h.hexdigest()
    return _SHA_CACHE[key]


def _norm_det(det: str) -> tuple[str, bool]:
    d = det.lower()
    if d == "bse":
        return "BSE", False
    if d == "inlens":
        return "Inlens", False
    if d == "etd":
        return "ETD", False
    return "ETD", True  # SE used in place of ETD


def batch_from_path(path: Path, root: Path) -> str:
    """Batch label = first `batch N` in the parent path (closest folder wins)."""
    try:
        rel = path.parent.relative_to(root)
        parts = [root.name, *rel.parts]
    except ValueError:
        parts = list(path.parent.parts)
    for part in reversed(parts):
        m = BATCH_RE.search(part)
        if m:
            return m.group(1)
    return "unknown"


def discover(data_dir: str | Path) -> list[SampleFiles]:
    root = Path(data_dir)
    samples: dict[str, SampleFiles] = {}
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        m = FILE_RE.match(p.name)
        if not m:
            continue
        sid = f"img_{m.group('id').lower()}"
        det, is_se = _norm_det(m.group("det"))
        batch = batch_from_path(p, root)
        if any(part.lower() in ("test", "batch_unknown", "unknown") for part in p.relative_to(root).parts[:-1]):
            batch = "unknown"
        pb = m.group("pb")
        if pb is not None and batch != "unknown" and pb != batch:
            raise SystemExit(f"{p}: filename prefix says batch {pb} but the folder says batch {batch}; "
                             "the folder is authoritative - fix the file location before continuing")
        s = samples.setdefault(sid, SampleFiles(sid, batch))
        s.has_batch_prefix = s.has_batch_prefix or pb is not None
        if det in s.paths:
            warnings.warn(f"{sid}: duplicate {det} file {p.name}; keeping {s.paths[det].name}")
            continue
        s.paths[det] = p
        s.etd_is_se = s.etd_is_se or is_se
        if s.batch != batch:
            warnings.warn(f"{sid}: files in different batch folders ({s.batch} vs {batch})")
    out = []
    for s in samples.values():
        if "BSE" not in s.paths:
            warnings.warn(f"{s.sample_id}: no BSE image, skipped")
            continue
        if not s.complete:
            warnings.warn(f"{s.sample_id}: incomplete triple {sorted(s.paths)}")
        out.append(s)
    return sorted(out, key=lambda s: (s.batch, s.sample_id))


def read_px_nm(page) -> tuple[float | None, float | None]:
    """Pixel size in nm from X/YResolution tags (None when absent)."""
    out = []
    for tag in ("XResolution", "YResolution"):
        if tag not in page.tags:
            out.append(None)
            continue
        num, den = page.tags[tag].value
        res = num / den if den else 0
        unit = page.tags["ResolutionUnit"].value if "ResolutionUnit" in page.tags else 2
        unit = int(unit)
        if res <= 1 or unit == 1:
            out.append(None)
            continue
        per_m = res / 0.0254 if unit == 2 else res / 0.01
        out.append(1e9 / per_m)
    return out[0], out[1]


@dataclass
class Channel:
    img: np.ndarray  # uint8 2-D
    px_nm: float
    px_nm_assumed: bool
    bad_cols: np.ndarray
    bad_rows: np.ndarray
    n_levels: int


def load_channel(path: str | Path, px_default: float = 25.0) -> Channel:
    with tifffile.TiffFile(path) as tf:
        page = tf.pages[0]
        a = page.asarray()
        px_x, px_y = read_px_nm(page)
    if a.ndim == 3:
        a = a[..., :3]
        disagree = (a[..., 0] != a[..., 1]) | (a[..., 0] != a[..., 2])
        bad_cols = np.where(disagree.mean(axis=0) > 0.5)[0]
        bad_rows = np.where(disagree.mean(axis=1) > 0.5)[0]
        g = a[..., 0]
    else:
        g = a
        bad_cols = bad_rows = np.array([], dtype=int)
    if g.dtype != np.uint8:
        g = np.clip(np.round(g.astype(float) / g.max() * 255), 0, 255).astype(np.uint8)
    assumed = px_x is None
    px = px_x if px_x is not None else px_default
    if px_y is not None and px_x is not None and abs(px_x - px_y) / px_x > 0.02:
        warnings.warn(f"{path}: anisotropic pixels {px_x:.2f}x{px_y:.2f} nm; using X")
    return Channel(g, float(px), assumed, bad_cols, bad_rows, int(len(np.unique(g))))


@dataclass
class LoadedSample:
    files: SampleFiles
    bse: np.ndarray
    etd: np.ndarray | None
    inlens: np.ndarray | None
    px_nm: float
    px_nm_assumed: bool
    resampled: bool
    raw_shape: tuple
    crop: dict  # top, bottom, left, right rows/cols removed
    n_levels: dict
    flags: list

    @property
    def px_um(self) -> float:
        return self.px_nm / 1000.0

    @property
    def area_um2(self) -> float:
        return self.bse.size * self.px_um**2


def _runs_from_edge(flag: np.ndarray, max_len: int) -> tuple[int, int]:
    """Length of contiguous True runs from the start and the end, each capped at max_len."""
    top = 0
    while top < len(flag) and flag[top]:
        top += 1
    bot = 0
    while bot < len(flag) and flag[len(flag) - 1 - bot]:
        bot += 1
    return (top if top <= max_len else 0), (bot if bot <= max_len else 0)


def detect_edge_bands(inlens: np.ndarray | None, etd: np.ndarray | None, cfg: dict) -> tuple[int, int]:
    """Un-sectioned nodular band at top/bottom: InLens dark AND ETD locally rough (brief 4.1.3)."""
    if inlens is None or etd is None:
        return 0, 0
    h = inlens.shape[0]
    max_len = int(cfg["edge_band_max_frac"] * h)
    row_med = np.median(inlens, axis=1).astype(float)
    glob_med = float(np.median(inlens))
    e = etd.astype(np.float32)
    local_var = ndi.uniform_filter(e * e, 5) - ndi.uniform_filter(e, 5) ** 2
    row_var = np.median(local_var, axis=1)
    glob_var = float(np.median(local_var))
    flag = (row_med < cfg["edge_band_inlens_ratio"] * glob_med) & (row_var > cfg["edge_band_etd_var_ratio"] * glob_var)
    # smooth the row flag over 9 rows so single noisy rows do not break a band
    flag = ndi.uniform_filter1d(flag.astype(float), 9) > 0.5
    return _runs_from_edge(flag, max_len)


def load_sample(files: SampleFiles, cfg: dict) -> LoadedSample:
    from skimage.transform import rescale

    flags: list[str] = []
    chans = {d: load_channel(p, cfg["px_nm_default"]) for d, p in files.paths.items()}
    bse = chans["BSE"]
    shapes = {d: c.img.shape for d, c in chans.items()}
    if len(set(shapes.values())) > 1:
        flags.append(f"channel_shape_mismatch:{shapes}")
    h = min(s[0] for s in shapes.values())
    w = min(s[1] for s in shapes.values())
    bad_cols = np.unique(np.concatenate([c.bad_cols for c in chans.values()]))
    bad_rows = np.unique(np.concatenate([c.bad_rows for c in chans.values()]))
    # bad columns/rows are only ever cropped from the frame edges (the observed artefact is a
    # 2-px column at the right edge); interior disagreement is flagged but kept.
    left = right = top = bottom = 0
    cols = set(bad_cols.tolist())
    while left in cols:
        left += 1
    while (w - 1 - right) in cols:
        right += 1
    rows = set(bad_rows.tolist())
    while top in rows:
        top += 1
    while (h - 1 - bottom) in rows:
        bottom += 1
    interior_bad = len(cols) - left - right
    if interior_bad > 0:
        flags.append(f"interior_channel_disagreement_cols:{interior_bad}")
    if left or right:
        flags.append(f"dropped_rgb_mismatch_cols:{left + right}")

    def crop(a):
        return a[top : h - bottom, left : w - right]

    imgs = {d: crop(c.img[:h, :w]) for d, c in chans.items()}
    px_nm = bse.px_nm
    assumed = bse.px_nm_assumed
    pxs = {d: c.px_nm for d, c in chans.items()}
    if max(pxs.values()) - min(pxs.values()) > cfg["resample_tolerance"] * px_nm:
        flags.append(f"px_size_differs_between_channels:{pxs}")
    resampled = False
    target = cfg["px_nm_target"]
    if abs(px_nm - target) / target > cfg["resample_tolerance"]:
        scale = px_nm / target
        for d in imgs:
            imgs[d] = np.clip(
                np.round(rescale(imgs[d].astype(np.float32), scale, order=1, anti_aliasing=scale < 1, preserve_range=True)),
                0, 255,
            ).astype(np.uint8)
        flags.append(f"resampled_from_{px_nm:.2f}nm")
        px_nm = target
        resampled = True
    if assumed:
        flags.append("px_nm_assumed")
    # edge bands of un-sectioned material
    band_top, band_bot = detect_edge_bands(imgs.get("Inlens"), imgs.get("ETD"), cfg)
    man = (cfg.get("manual_masks") or {}).get(files.sample_id)
    m_l = m_r = 0
    if man:
        band_top = int(man.get("top", band_top))
        band_bot = int(man.get("bottom", band_bot))
        m_l, m_r = int(man.get("left", 0)), int(man.get("right", 0))
        flags.append("manual_mask")
    if band_top or band_bot:
        flags.append(f"edge_band_masked:top={band_top},bottom={band_bot}")
    hh, ww = imgs["BSE"].shape
    for d in imgs:
        imgs[d] = np.ascontiguousarray(imgs[d][band_top : hh - band_bot, m_l : ww - m_r])
    if files.etd_is_se:
        flags.append("etd_is_SE_detector")
    return LoadedSample(
        files=files,
        bse=imgs["BSE"],
        etd=imgs.get("ETD"),
        inlens=imgs.get("Inlens"),
        px_nm=px_nm,
        px_nm_assumed=assumed,
        resampled=resampled,
        raw_shape=(h, w),
        crop=dict(top=top + band_top, bottom=bottom + band_bot, left=left + m_l, right=right + m_r,
                  band_top=band_top, band_bottom=band_bot),
        n_levels={d: c.n_levels for d, c in chans.items()},
        flags=flags,
    )


def manifest_frame(samples: list[SampleFiles], cfg: dict) -> pd.DataFrame:
    rows = []
    for s in samples:
        with tifffile.TiffFile(s.paths["BSE"]) as tf:
            page = tf.pages[0]
            px, _ = read_px_nm(page)
            shape = page.shape
        rows.append(dict(
            sample_id=s.sample_id, batch=s.batch,
            path_bse=str(s.paths.get("BSE", "")), path_etd=str(s.paths.get("ETD", "")),
            path_inlens=str(s.paths.get("Inlens", "")), etd_is_se=s.etd_is_se,
            etd_label_in_filename=s.etd_label_in_filename, has_batch_prefix=s.has_batch_prefix,
            sha1_bse=file_sha1(s.paths["BSE"]),
            width=shape[1], height=shape[0], px_nm=px if px else cfg["px_nm_default"],
            px_nm_assumed=px is None, complete=s.complete,
        ))
    return pd.DataFrame(rows)
