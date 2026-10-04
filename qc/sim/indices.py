"""Closed-form relative indices (brief 10.3, indices.py). Each carries a trust grade:
A = arithmetic on measured quantities; B = direction supported, level set by an assumption (ratios only);
C = heuristic."""
from __future__ import annotations

import numpy as np

F = 96485.0

GRADES = {
    "energy_density_mAh_cm3": "A",
    "D_eff_rel_TP": "B", "D_eff_rel_IP": "B", "aniso_ratio": "B", "tau_p": "B", "N_M": "B",
    "t_diffusion_s": "B", "i_lim_A_m2": "B", "plating_risk_index": "B", "icl_index": "B",
    "sigma_eff_rel_TP": "B", "carbon_spanning_frac": "B", "exposed_si_frac": "B",
    "pore_closure_frac": "B", "constraint_index": "B", "constraint_index_aw": "B", "buffer_sufficiency_frac": "B",
    "si_touch_frac": "B", "largest_merged_cluster_ecd_um": "B", "si_area_in_clusters_gt_5um": "B", "dH_H_bound": "B",
    "Q_CC_1C_over_Q_C10": "B", "Q_CC_3C_over_Q_C10": "B", "min_neg_surface_dphi_sep_3C_V": "B", "min_electrolyte_conc_3C": "B",
    "cli_index": "C",
}
# direction in which a higher value is *worse* for the battery (+1) or better (-1); used for wording only
WORSE_IF_HIGHER = {"tau_p": 1, "N_M": 1, "t_diffusion_s": 1, "plating_risk_index": 1, "icl_index": 1, "constraint_index": 1,
                   "constraint_index_aw": 1, "si_touch_frac": 1, "largest_merged_cluster_ecd_um": 1, "si_area_in_clusters_gt_5um": 1,
                   "dH_H_bound": 1, "cli_index": 1, "pore_closure_frac": 1,
                   "D_eff_rel_TP": -1, "i_lim_A_m2": -1, "sigma_eff_rel_TP": -1, "buffer_sufficiency_frac": -1, "energy_density_mAh_cm3": -1,
                   "Q_CC_1C_over_Q_C10": -1, "Q_CC_3C_over_Q_C10": -1, "min_neg_surface_dphi_sep_3C_V": -1, "min_electrolyte_conc_3C": -1}


def transport_proxies(D_eff_rel: float, cfg: dict) -> dict:
    el = cfg["sim"]["electrolyte"]
    L = cfg["sim"]["electrode"]["L_um"] * 1e-6
    D = el["D_bulk_m2s"] * D_eff_rel
    if not np.isfinite(D) or D <= 0:
        return dict(t_diffusion_s=np.nan, i_lim_A_m2=np.nan)
    return dict(t_diffusion_s=L**2 / D, i_lim_A_m2=2 * F * el["c0_molm3"] * D / ((1 - el["t_plus"]) * L))


def energy_density(eps: float, f_si: float, cfg: dict) -> float:
    """Volumetric negative-electrode capacity (mAh/cm3) from *volume* fractions:
    Q = (1 - eps) * sum_i phi_i rho_i q_i, phi over solids with CBD inactive."""
    en = cfg["sim"]["energy"]
    cbd = cfg["sim"]["electrode"]["cbd_frac_of_solids"]
    phi_si = f_si
    phi_gr = max(1 - f_si - cbd, 0.0)
    return float((1 - eps) * (phi_si * en["rho_si"] * en["q_si_pract"] + phi_gr * en["rho_gr"] * en["q_gr"]))


def icl_index(labels, px_um: float) -> float:
    """First-cycle-loss proxy: exposed interface per area, Si-pore interface weighted x3 (SEI on Si grows
    and re-forms with every cycle) plus carbon-pore interface. Absolute level arbitrary - ratio only."""
    from scipy import ndimage as ndi

    from .laplace import CARBON, PORE, SI

    four = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], bool)
    pore = labels == PORE
    ring = ndi.binary_dilation(pore, structure=four) & ~pore
    si_if = (ring & (labels == SI)).sum()
    c_if = (ring & (labels == CARBON)).sum()
    return float((3 * si_if + c_if) / labels.size / px_um)


def plating_risk_index(D_eff_rel: float, D_eff_rel_base: float) -> float:
    """Ionic-resistance ratio vs baseline (R_ion ∝ L / (kappa * D_eff_rel), same L and kappa)."""
    return float(D_eff_rel_base / D_eff_rel) if D_eff_rel and np.isfinite(D_eff_rel) and D_eff_rel > 0 else np.nan


def cli_index(x: dict, base: dict, cfg: dict) -> float:
    """Cycle-life HEURISTIC (grade C): weighted sum of sample/baseline ratios of 'worse' drivers; 1 = baseline."""
    w = cfg["sim"]["cli_weights"]
    terms = {
        "f_si": (x.get("si_frac_solid"), base.get("si_frac_solid")),
        "d50": (x.get("si_d50_aw_um"), base.get("si_d50_aw_um")),
        "buffer": (1 - (x.get("buffer_sufficiency_frac") or 0), 1 - (base.get("buffer_sufficiency_frac") or 0)),
        "agglomerate": (x.get("si_agglomerate_frac"), base.get("si_agglomerate_frac")),
        "contact": (1 - (x.get("si_contact_carbon_frac") or 0), 1 - (base.get("si_contact_carbon_frac") or 0)),
    }
    tot = wsum = 0.0
    for k, (a, b) in terms.items():
        if a is None or b is None or not np.isfinite(a) or not np.isfinite(b) or b <= 0:
            continue
        tot += w[k] * a / b
        wsum += w[k]
    return float(tot / wsum) if wsum else np.nan


def porosity(pore_frac_deep: float, base_median: float | None, convention: str) -> float:
    """offset_baseline_to_0.30: shift so the baseline median deep-pore fraction maps to 0.30 (the deep-pore
    fraction is a lower bound); deep_pore_as_is: use the lower bound directly."""
    if convention.startswith("offset"):
        return float(np.clip(0.30 + (pore_frac_deep - (base_median if base_median is not None else pore_frac_deep)), 0.05, 0.8))
    return float(np.clip(pore_frac_deep, 0.03, 0.8))
