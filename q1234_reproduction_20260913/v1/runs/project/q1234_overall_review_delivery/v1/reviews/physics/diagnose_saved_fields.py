"""Read frozen arrays; do algebra and spatial quadrature only. Never solve PDEs.

Run from the repository root with Python, NumPy, SciPy and openpyxl.
Writes only beside this script. All reported physical interpretations are conditional.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform
import subprocess

import numpy as np
import scipy
from scipy.fft import dct
from scipy.special import exp1
from openpyxl import load_workbook

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[3]
FILES = {
    "q4_saved_field": "q4_complete_delivery/v1/output/main.npz",
    "q1_saved_field": "q1_complete_delivery/q1_delivery/output/q1_unrounded.npz",
    "environment": "q4_complete_delivery/v1/inputs/attachment1.xlsx",
    "radius_input": "q4_complete_delivery/v1/inputs/attachment2.xlsx",
    "q4_model": "q4_complete_delivery/v1/第四问模型与文献采用.md",
    "q4_source": "q4_complete_delivery/v1/source/q4_spectral.py",
    "problem_text": "readable/problem/fulltext.md",
}


def array_workbook(path):
    wb = load_workbook(path, read_only=True, data_only=True)
    rows = list(wb.worksheets[0].values)
    wb.close()
    return np.asarray(rows[1:], float)


def reconstructed_c_integral(c, gauss_order):
    """Independent Gauss integration of density from the primitive polynomial.

    Uses NumPy Chebyshev evaluation, not the production quadrature or PDE operator.
    C interpolation follows the frozen definition: interpolate F(C), then invert it.
    """
    primitive = c * np.exp(-0.30 / c) - 0.30 * exp1(0.30 / c)
    coeff = dct(primitive[::-1], type=1) / (len(c) - 1)
    coeff[[0, -1]] *= 0.5
    points, weights = np.polynomial.legendre.leggauss(gauss_order)
    target = np.polynomial.chebyshev.chebval(points, coeff)
    cc = np.interp((points + 1) / 2, (1 - np.cos(np.linspace(0, np.pi, len(c)))) / 2, c)
    for _ in range(50):
        step = (cc * np.exp(-0.30 / cc) - 0.30 * exp1(0.30 / cc) - target) / np.exp(-0.30 / cc)
        while np.any(cc - step <= 0):
            bad = cc - step <= 0
            step[bad] *= 0.5
        cc -= step
        if np.max(np.abs(step)) < 1e-13:
            break
    else:
        raise RuntimeError("Primitive inverse did not converge")
    assert np.max(np.abs(cc * np.exp(-0.30 / cc) - 0.30 * exp1(0.30 / cc) - target)) < 1e-12
    return float(np.dot(weights, (760 + 90 * cc) / (1 + cc)) / 2)


def main():
    base_commit = "anonymous-code-snapshot"
    provenance = {
        key: {"path": value, "sha256": hashlib.sha256((ROOT / value).read_bytes()).hexdigest()}
        for key, value in FILES.items()
    }
    data = array_workbook(ROOT / FILES["environment"])
    future = data[data[:, 0] >= 10800, 1:].mean(axis=0)
    q4 = np.load(ROOT / FILES["q4_saved_field"])
    t, radii, fields, stored_w = (q4[k] for k in ("time_s", "radius_m", "full_TC", "quadrature_weights"))
    rho_d0 = (760 + 90 * 2.55) / 3.55  # initial matching convention, not a measurement
    q_implied = (760 + 90 * fields[:, :, 1]) / (1 + fields[:, :, 1])
    ratios = (radii / 0.02) ** 2 * (q_implied @ stored_w) / rho_d0
    rho_d_material = rho_d0 * (0.02 / radii) ** 2
    rows = []
    for i in sorted({int(np.argmin(ratios)), int(np.searchsorted(t, 7 * 3600)), len(t) - 1}):
        quad = {str(n): reconstructed_c_integral(fields[i, :, 1], n) for n in (192, 384)}
        ratio_independent = (radii[i] / 0.02) ** 2 * quad["384"] / rho_d0
        rows.append({
            "time_s": float(t[i]), "radius_m": float(radii[i]),
            "implied_dry_mass_ratio_nodal": float(ratios[i]),
            "implied_dry_mass_ratio_gauss384": float(ratio_independent),
            "gauss192_384_rho_d_integral_difference_kg_m3": quad["192"] - quad["384"],
            "length_for_global_mass_only_cm": float(25 / ratio_independent),
            "physical_rho_d_minmax_if_empirical_rho_real": [float(q_implied[i].min()), float(q_implied[i].max())],
            "affine_material_rho_d_kg_m3": float(rho_d_material[i]),
            "qualification": "Changing length alone preserves these saved C and R only diagnostically; not a recomputed model. Local consistency does not follow from global correction.",
        })
    energy = []
    for hour in (6, 12, 24, 48):
        i = int(np.searchsorted(t, hour * 3600))
        ts, cs = fields[i, -1]
        j = rho_d_material[i] * 8e-7 * (cs - future[1])
        energy.append({
            "time_s": float(t[i]), "surface_T_degC": float(ts), "surface_C": float(cs),
            "convective_heat_input_on_saved_temperature_W_m2": float(25 * (future[0] - ts)),
            "conditional_material_water_flux_kg_m2_s": float(j),
            "conditional_latent_demand_at_2_45_MJ_kg_W_m2": float(2.45e6 * j),
        })
    q1 = np.load(ROOT / FILES["q1_saved_field"])
    volume = np.pi * 0.02**2 * 0.25
    mass_loss = 820 / 3.55 * volume * (2.55 - q1["average_C"][-1])
    sensible = 820 * 2600 * volume * (q1["average_temperature_degC"][-1] - 28)
    qconv = float(np.trapezoid(2 * np.pi * 0.02 * 0.25 * 25 * (q1["environment"][:, 0] - q1["temperature_degC"][:, -1]), q1["time_s"]))
    scale = (0.02 / radii[-1]) ** 2
    result = {
        "base_commit": base_commit,
        "scope": "Algebra and saved-field quadrature only; no PDE integration, no new temperature or drying-time solution, no calibration.",
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "sources": provenance,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "q4_future_environment_T_C": future.tolist(),
        "conditional_air_state_if_W_and_101_325kPa": {
            "vapor_partial_pressure_kPa": float(101.325 * future[1] / (0.621945 + future[1])),
            "relative_humidity_using_FAO_eq11": float((101.325 * future[1] / (0.621945 + future[1])) / (0.6108 * np.exp(17.27 * future[0] / (future[0] + 237.3)))),
            "qualification": "Only if the workbook kg/kg means water per dry air, with stated pressure. This is not a verified interpretation of the problem input.",
        },
        "q4_density_contradiction_diagnostics": rows,
        "q4_conditional_latent_flux": energy,
        "q4_boundary_coefficient_scaling": {
            "terminal_R_over_R0": float(radii[-1] / 0.02),
            "rho_d_and_per_area_K_C_multiplier_at_fixed_effective_hm": float(scale),
            "total_mass_rate_multiplier_at_fixed_surface_C_difference": float(0.02 / radii[-1]),
            "hypothetical_effective_hm_if_K_C_constant_m_s": float(8e-7 / scale),
            "qualification": "Hypothetical hm is a sensitivity contract only; not a fitted or recommended replacement parameter.",
        },
        "q1_conditional_energy_readback": {
            "initial_dry_mass_kg": float(820 / 3.55 * volume),
            "water_loss_kg": float(mass_loss), "latent_J": float(2.45e6 * mass_loss),
            "baseline_sensible_J": float(sensible), "baseline_convective_integral_J": qconv,
            "latent_over_baseline_sensible": float(2.45e6 * mass_loss / sensible),
            "qualification": "rho=820 interpreted as initial wet density only; all effective moisture loss assumed evaporative; latent value is a scale diagnostic.",
        },
        "assertions": {
            "density_gauss192_384_close": all(abs(x["gauss192_384_rho_d_integral_difference_kg_m3"]) < 1e-8 for x in rows),
            "q4_required_length_exceeds_initial_at_7h_and_end": all(x["length_for_global_mass_only_cm"] > 25 for x in rows),
            "rho_d_shape_not_spatially_uniform_at_7h": bool(np.ptp(q_implied[int(np.searchsorted(t, 7*3600))]) > 1),
        },
    }
    assert all(result["assertions"].values()), result["assertions"]
    (OUT / "saved_field_physics_diagnostics.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with (OUT / "q4_implied_dry_mass_curve.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time_s", "radius_m", "implied_dry_mass_ratio_if_empirical_rho_real", "length_cm_for_global_mass_only"])
        writer.writerows(zip(t, radii, ratios, 25 / ratios))
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
