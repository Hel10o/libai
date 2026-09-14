"""Cold-start sealed-boundary pure-shrinkage check over observed 72 hours."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import solve_ivp


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    sys.path.insert(0, str(a.root / "q4_complete_delivery/v1/source"))
    from geometry_solver import GeometryFV, Inputs
    inputs = Inputs(a.root / "q4_complete_delivery/v1/inputs")
    op = GeometryFV(40, 0, "q4", inputs)
    op.side[:] = 0.0
    op.end[:] = 0.0
    initial = np.tile([28., 2.55], op.m)
    times = np.linspace(0, 259200, 145)
    solution = solve_ivp(lambda t, y: op.rhs(t, y, inputs.ambient), (0, times[-1]), initial,
                         method="BDF", jac=lambda t, y: op.rhs(t, y, inputs.ambient, True),
                         rtol=2e-11, atol=1e-13, t_eval=times, max_step=1800)
    assert solution.success
    change = float(np.max(abs(solution.y - initial[:, None])))
    ratios = np.array([inputs.radius(t) / .02 for t in times])
    rho0 = (760 + 90 * 2.55) / 3.55
    dry_mass = rho0 / ratios**2 * ratios**2
    water_mass = dry_mass * solution.y[1::2].mean(axis=0)
    dry_error = float(np.max(abs(dry_mass / dry_mass[0] - 1)))
    water_error = float(np.max(abs(water_mass / water_mass[0] - 1)))
    assert change < 1e-12 and max(dry_error, water_error) < 1e-14
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps({"scope": "zero heat/moisture exchange, uniform initial state, prescribed observed radius",
        "end_h": 72, "temperature_moisture_max_change": change,
        "relative_dry_mass_change": dry_error, "relative_water_mass_change": water_error,
        "pass": True}, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
