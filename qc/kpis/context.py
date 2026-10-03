from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np


@dataclass
class Ctx:
    """Masks and precomputed maps for one region (whole image or a vertical strip)."""

    px_um: float
    bse: np.ndarray  # uint8 raw BSE
    bse_s: np.ndarray  # float32 smoothed BSE
    pore: np.ndarray
    si: np.ndarray
    carbon: np.ndarray
    ambiguous: np.ndarray
    si_labels: np.ndarray
    t1: float
    t2: float
    noise_sd: float
    cfg: dict
    lt_map: np.ndarray | None = None  # pore local thickness (um), full resolution
    st: tuple | None = None  # structure tensor (Arr, Arc, Acc)
    crack_resp: np.ndarray | None = None  # Sato ridge response (float16)
    crack_threshold: float | None = None
    curtain_deg: float | None = None
    extras: dict = field(default_factory=dict)

    def px_len(self, n_ref: float) -> float:
        """A length given in 25-nm reference pixels, converted to this image's pixels."""
        return n_ref * 0.025 / self.px_um

    def px_area(self, n_ref: float) -> float:
        """An area given in 25-nm reference pixels², converted to this image's pixels²."""
        return n_ref * (0.025 / self.px_um) ** 2

    @property
    def shape(self):
        return self.bse.shape

    @property
    def area_px(self) -> int:
        return int(self.bse.size)

    @property
    def area_um2(self) -> float:
        return self.area_px * self.px_um**2

    @property
    def area_mm2(self) -> float:
        return self.area_um2 * 1e-6

    def sub(self, cols: slice) -> "Ctx":
        def c(a):
            return None if a is None else a[:, cols]

        return replace(
            self,
            bse=self.bse[:, cols], bse_s=self.bse_s[:, cols], pore=self.pore[:, cols], si=self.si[:, cols],
            carbon=self.carbon[:, cols], ambiguous=self.ambiguous[:, cols], si_labels=self.si_labels[:, cols],
            lt_map=c(self.lt_map), st=None if self.st is None else tuple(a[:, cols] for a in self.st),
            crack_resp=c(self.crack_resp), extras={},
        )
