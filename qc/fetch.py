"""Download the Drive data folder from data/drive_manifest.csv (public share link required).

Run this where drive.google.com is reachable. Files land in <data>/<folder>/<title> and sizes are
verified against the manifest.
"""
from __future__ import annotations

import csv
import shutil
import urllib.request
from pathlib import Path

URL = "https://drive.usercontent.google.com/download?id={id}&export=download&confirm=t"


def run(manifest: Path, data: Path, verify_only: bool = False):
    rows = list(csv.DictReader(open(manifest)))
    ok = bad = 0
    for r in rows:
        dest = Path(data) / r["folder"] / r["title"]
        size = int(r["size_bytes"])
        if dest.exists() and dest.stat().st_size == size:
            ok += 1
            continue
        if verify_only:
            print(f"missing/size mismatch: {dest}")
            bad += 1
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_suffix(dest.suffix + ".part")
        print(f"downloading {r['title']} ({size / 1e6:.1f} MB)")
        with urllib.request.urlopen(URL.format(id=r["file_id"])) as resp, open(tmp, "wb") as f:
            shutil.copyfileobj(resp, f)
        if tmp.stat().st_size != size:
            print(f"  size mismatch ({tmp.stat().st_size} != {size}); left as {tmp}")
            bad += 1
            continue
        tmp.rename(dest)
        ok += 1
    print(f"{ok} files present and verified, {bad} problems")
