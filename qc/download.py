"""drive_manifest.csv -> data/ : resume-safe, size-verified downloads of the Drive folder.

Routes, tried in order:
  1. Drive API v3 media endpoint (www.googleapis.com) when GOOGLE_API_KEY is set — works for files
     shared "Anyone with the link" and only needs googleapis.com egress.
  2. gdown by file id (drive.google.com / drive.usercontent.google.com).
  3. plain HTTPS to drive.usercontent.google.com.
Files are written to <data>/<folder>/<title>.part, size-checked, opened with tifffile, then renamed.
"""
from __future__ import annotations

import csv
import os
import shutil
import time
import urllib.request
from pathlib import Path

API_URL = "https://www.googleapis.com/drive/v3/files/{id}?alt=media&supportsAllDrives=true&key={key}"
USERCONTENT_URL = "https://drive.usercontent.google.com/download?id={id}&export=download&confirm=t"


def _http(url: str, dest: Path):
    with urllib.request.urlopen(url, timeout=300) as resp, open(dest, "wb") as f:
        shutil.copyfileobj(resp, f)


def _gdown(file_id: str, dest: Path):
    import gdown

    out = gdown.download(id=file_id, output=str(dest), quiet=True)
    if out is None:
        raise RuntimeError("gdown returned None")


def _routes():
    key = os.environ.get("GOOGLE_API_KEY")
    if key:
        yield "drive-api", lambda fid, d: _http(API_URL.format(id=fid, key=key), d)
    yield "gdown", _gdown
    yield "usercontent", lambda fid, d: _http(USERCONTENT_URL.format(id=fid), d)


def _valid_tiff(p: Path) -> bool:
    try:
        import tifffile

        with tifffile.TiffFile(p) as tf:
            return len(tf.pages) >= 1
    except Exception:
        return False


def run(manifest: Path, data: Path, verify_only: bool = False, only: list[str] | None = None, retries: int = 2) -> list[dict]:
    rows = list(csv.DictReader(open(manifest)))
    if only:
        rows = [r for r in rows if r["sample_id"] in only or r["sample_id"].removeprefix("img_") in only]
    log = []
    for r in rows:
        dest = Path(data) / r["folder"] / r["title"]
        size = int(r["size_bytes"])
        if dest.exists() and dest.stat().st_size == size:
            log.append(dict(file=r["title"], bytes=size, status="present"))
            continue
        if verify_only:
            log.append(dict(file=r["title"], bytes=dest.stat().st_size if dest.exists() else 0, status="missing"))
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_name(dest.name + ".part")
        status, errors = "failed", []
        for name, fn in _routes():
            for attempt in range(retries + 1):
                try:
                    fn(r["file_id"], tmp)
                    got = tmp.stat().st_size
                    if got != size:
                        raise RuntimeError(f"size {got} != {size} (HTML error page?)")
                    if not _valid_tiff(tmp):
                        raise RuntimeError("not a readable TIFF")
                    tmp.rename(dest)
                    status = f"downloaded ({name})"
                    break
                except Exception as e:  # noqa: BLE001 - every failure is logged and the next route tried
                    errors.append(f"{name}: {str(e)[:120]}")
                    time.sleep(2**attempt)
            if status != "failed":
                break
        if tmp.exists() and status == "failed":
            tmp.unlink()
        log.append(dict(file=r["title"], bytes=size, status=status if status != "failed" else "failed: " + " | ".join(errors[-3:])))
        print(f"{r['title']}: {log[-1]['status']}")
    ok = sum(l["status"] in ("present",) or l["status"].startswith("downloaded") for l in log)
    print(f"{ok}/{len(log)} files present and verified")
    return log
