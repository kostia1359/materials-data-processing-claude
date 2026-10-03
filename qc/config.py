from __future__ import annotations

from pathlib import Path

import yaml

DEFAULT_PATH = Path(__file__).resolve().parent.parent / "config.yaml"


def load_config(path: str | Path | None = None) -> dict:
    with open(DEFAULT_PATH) as f:
        cfg = yaml.safe_load(f)
    if path and Path(path).resolve() != DEFAULT_PATH:
        with open(path) as f:
            cfg.update(yaml.safe_load(f) or {})
    return cfg
