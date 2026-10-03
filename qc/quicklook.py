"""Quicklook PNGs (R=BSE, G=ETD, B=InLens composite + histograms) and a per-sample facts table."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .gates import comb_step, registration_shift  # noqa: E402
from .io import discover, load_sample  # noqa: E402


def run(data: Path, out: Path, cfg: dict):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for s in discover(data):
        L = load_sample(s, cfg)
        chans = {"BSE": L.bse, "ETD": L.etd, "Inlens": L.inlens}
        f = 8
        comp = np.stack([(c if c is not None else L.bse)[::f, ::f] for c in chans.values()], -1)
        fig, axs = plt.subplots(2, 1, figsize=(12, 6), gridspec_kw={"height_ratios": [2, 1]})
        axs[0].imshow(comp)
        axs[0].set_title(f"{s.sample_id} (batch {s.batch}) — R=BSE G={'SE' if s.etd_is_se else 'ETD'} B=InLens, 1/{f} scale")
        axs[0].axis("off")
        for (name, c), col in zip(chans.items(), ("r", "g", "b")):
            if c is not None:
                axs[1].hist(c.ravel()[::7], bins=256, range=(0, 256), histtype="step", color=col, label=name)
        axs[1].set_yscale("log")
        axs[1].legend()
        axs[1].set_xlabel("grey level (note comb histograms = post-acquisition LUT)")
        fig.tight_layout()
        fig.savefig(out / f"{s.sample_id}.png", dpi=90)
        plt.close(fig)
        r = dict(sample_id=s.sample_id, batch=s.batch, etd_is_se=s.etd_is_se, width=L.raw_shape[1], height=L.raw_shape[0],
                 px_nm=round(L.px_nm, 3), px_nm_assumed=L.px_nm_assumed, flags=";".join(L.flags))
        for name, c in chans.items():
            if c is None:
                continue
            r[f"levels_{name}"] = L.n_levels.get(name)
            r[f"comb_{name}"] = comb_step(c)
            r[f"median_{name}"] = float(np.median(c))
        r["reg_etd_px"] = registration_shift(L.bse, L.etd) if L.etd is not None else np.nan
        r["reg_inlens_px"] = registration_shift(L.bse, L.inlens) if L.inlens is not None else np.nan
        rows.append(r)
        print(f"  {s.sample_id}")
    df = pd.DataFrame(rows)
    df.to_csv(out / "facts.csv", index=False)
    (out / "facts.md").write_text(df.to_markdown(index=False, floatfmt=".2f") + "\n\nPer-batch medians:\n\n"
                                  + df.groupby("batch")[[c for c in df if c.startswith(("median_", "levels_", "height"))]].median().to_markdown(floatfmt=".1f") + "\n")
    print(f"wrote {out / 'facts.csv'}")
