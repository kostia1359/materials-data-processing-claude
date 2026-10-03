"""KPI metadata: unit, robustness class, one-line battery meaning, shortlist membership."""
from __future__ import annotations

# name: (unit, robustness class, meaning when HIGHER than baseline, meaning when LOWER, family)
KPI_META = {
    "si_frac_solid": ("fraction", "M", "more Si-candidate in the solid: higher capacity but more swelling, first-cycle loss and N/P shift",
                      "less Si-candidate in the solid: lower capacity, less swelling", "Si loading"),
    "si_frac_total": ("fraction", "M", "more Si-candidate per image area", "less Si-candidate per image area", "Si loading"),
    "si_d50_aw_um": ("um", "R", "coarser Si-candidate (area-weighted D50): higher fracture and swelling risk, more particle isolation",
                     "finer Si-candidate: more surface/SEI, less fracture", "Si size"),
    "si_d90_aw_um": ("um", "R", "more coarse Si-candidate particles in the tail: fracture/swelling hot-spots",
                     "fewer coarse Si-candidate particles", "Si size"),
    "si_d10_aw_um": ("um", "M", "fine fraction coarser", "more fine Si-candidate", "Si size"),
    "si_span": ("ratio", "M", "broader Si-candidate size distribution", "narrower Si-candidate size distribution", "Si size"),
    "si_num_density_mm2": ("mm^-2", "M", "more Si-candidate particles per area (finer or more loaded)",
                           "fewer Si-candidate particles per area (coarser or less loaded)", "Si size"),
    "si_solidity_med": ("ratio", "R", "more compact/blocky Si-candidate shapes", "more irregular/porous Si-candidate shapes (grade change?)", "Si shape"),
    "si_aspect_med": ("ratio", "R", "more elongated Si-candidate particles", "more equiaxed Si-candidate particles", "Si shape"),
    "si_circularity_med": ("ratio", "R", "rounder Si-candidate outlines", "more angular Si-candidate outlines", "Si shape"),
    "si_interior_texture": ("ratio", "M", "more internal contrast in Si-candidate (composite/porous grade?)", "smoother Si-candidate interiors", "Si shape"),
    "si_quadrat_cv": ("ratio", "M", "Si-candidate less evenly distributed (mixing/dispersion issue): local swelling hot-spots",
                      "Si-candidate more evenly distributed", "Si dispersion"),
    "si_clark_evans_R": ("ratio", "M", "Si-candidate more regularly spaced", "Si-candidate more clustered", "Si dispersion"),
    "si_agglomerate_frac": ("fraction", "M", "more Si-candidate in agglomerates: they swell together, lose contact and debond neighbours",
                            "less agglomerated Si-candidate", "Si dispersion"),
    "pore_frac_deep": ("fraction", "W", "more deep porosity (lower bound): better electrolyte transport, lower energy density",
                       "less deep porosity (lower bound): denser coating, slower ionic transport; consistent with heavier calendering", "Porosity"),
    "pore_frac_slope_per_level": ("1/level", "W", "pore fraction more threshold-sensitive", "pore fraction less threshold-sensitive", "Porosity"),
    "ambiguous_frac": ("fraction", "W", "more BSE-dark/ETD-bright pixels (shallow gaps)", "fewer ambiguous pixels", "Porosity"),
    "pore_chord_h_mean_um": ("um", "R", "wider in-plane pores", "narrower in-plane pores; consistent with heavier calendering", "Pore size"),
    "pore_chord_v_mean_um": ("um", "R", "taller through-plane pores: easier through-plane transport",
                             "flatter through-plane pores: more tortuous transport; consistent with heavier calendering", "Pore size"),
    "pore_chord_h_p90_um": ("um", "M", "more large in-plane pores", "fewer large in-plane pores", "Pore size"),
    "pore_lt_d50_um": ("um", "M", "thicker pores (local thickness): better wetting and fast-charge paths",
                       "thinner pores: harder wetting, slower fast charge", "Pore size"),
    "aniso_pore_chord_ratio": ("ratio", "R", "pores more aligned parallel to the collector: higher through-plane tortuosity, slower fast charge; consistent with heavier calendering",
                               "pores less aligned: lower through-plane tortuosity", "Anisotropy"),
    "aniso_carbon_chord_ratio": ("ratio", "R", "graphite/carbon more aligned parallel to the collector: higher through-plane tortuosity; consistent with heavier calendering",
                                 "carbon matrix less aligned", "Anisotropy"),
    "st_coherence": ("ratio", "R", "more strongly oriented texture", "less oriented texture", "Anisotropy"),
    "st_orientation_deg": ("deg", "R", "texture tilted counter-clockwise", "texture tilted clockwise", "Anisotropy"),
    "s2_len_pore_h_px": ("px", "R", "pore network stretched in-plane (longer correlation length): flatter, more aligned pores; consistent with heavier calendering",
                         "shorter in-plane pore correlation length: finer, less stretched pore network", "Two-point"),
    "s2_len_pore_v_px": ("px", "R", "longer through-plane pore correlation length", "shorter through-plane pore correlation length", "Two-point"),
    "s2_len_si_px": ("px", "R", "Si-candidate features larger/more widely correlated", "Si-candidate features smaller", "Two-point"),
    "s2_integral_range_pore_px2": ("px^2", "M", "coarser pore network (larger RVE)", "finer pore network", "Two-point"),
    "pore_percolating_frac_v": ("fraction", "M", "more through-plane connected deep pores: better fast-charge transport",
                                "fewer through-plane connected deep pores: transport bottleneck, plating risk", "Connectivity"),
    "pore_percolating_frac_h": ("fraction", "M", "more in-plane connected deep pores", "fewer in-plane connected deep pores", "Connectivity"),
    "pore_euler_density_mm2": ("mm^-2", "M", "more isolated pores (less connected)", "more loops/connected pores", "Connectivity"),
    "interface_pore_solid_um_per_um2": ("um^-1", "M", "more exposed pore/solid interface: more SEI, larger first-cycle loss",
                                        "less exposed interface", "Interfaces"),
    "si_contact_pore_frac": ("fraction", "R", "more Si-candidate surface facing pores: more room to swell, but more SEI",
                             "less buffer space around Si-candidate: swelling stresses the matrix", "Interfaces"),
    "si_contact_carbon_frac": ("fraction", "R", "Si-candidate better wired to the carbon matrix",
                               "Si-candidate less wired to carbon: risk of electronic isolation", "Interfaces"),
    "crack_density_um_per_mm2": ("um/mm^2", "M", "more cracks: impedance rise, crack cascade in Si", "fewer cracks", "Defects"),
    "largest_void_ecd_um": ("um", "M", "larger largest void", "smaller largest void", "Defects"),
    "hiZ_inclusion_count_mm2": ("mm^-2", "M", "more high-Z inclusions (contaminants): separator puncture risk", "fewer high-Z inclusions", "Defects"),
    "vertical_pore_slope": ("pp/10um", "W", "porosity increases down the frame (diagnostic)", "porosity decreases down the frame (diagnostic)", "Vertical"),
    "vertical_si_slope": ("pp/10um", "W", "Si-candidate increases down the frame (diagnostic)", "Si-candidate decreases down the frame (diagnostic)", "Vertical"),
    "si_frac_se_imagerep": ("fraction", "R", "larger representativeness error on Si fraction", "smaller representativeness error", "Two-point"),
    "pore_frac_se_imagerep": ("fraction", "R", "larger representativeness error on pore fraction", "smaller representativeness error", "Two-point"),
}

SHORTLIST = [
    "si_frac_solid", "si_d50_aw_um", "si_d90_aw_um", "si_num_density_mm2", "si_solidity_med", "si_quadrat_cv",
    "pore_frac_deep", "pore_chord_h_mean_um", "pore_chord_v_mean_um", "aniso_pore_chord_ratio",
    "aniso_carbon_chord_ratio", "s2_len_pore_h_px", "pore_percolating_frac_v", "si_contact_pore_frac",
    "si_contact_carbon_frac", "crack_density_um_per_mm2",
]
DEFECT_COUNTS = {"crack_density_um_per_mm2": "crack_count", "hiZ_inclusion_count_mm2": "hiZ_inclusion_count"}
KPI_NAMES = list(KPI_META)


def unit(k):
    return KPI_META.get(k, ("", "M", "", "", ""))[0]


def rclass(k):
    return KPI_META.get(k, ("", "M", "", "", ""))[1]


def meaning(k, direction: str) -> str:
    m = KPI_META.get(k)
    if not m:
        return ""
    return m[2] if direction == "higher" else m[3]


def weight(k, cfg) -> float:
    return float(cfg["weights"][rclass(k)])
