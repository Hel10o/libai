"""Read frozen geometry fields; targeted continuations, never full cold starts.

Uses the frozen FV operator, so this is not an independent discretization.
All generated files remain next to this script. No original data is modified.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import platform
import subprocess
import time
import numpy as np
import scipy
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = ROOT / "q4_complete_delivery/v1/source/geometry_solver_executed.py"
GEOM = ROOT / "q4_complete_delivery/v1/validation/geometry"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_operator():
    spec = importlib.util.spec_from_file_location("frozen_geometry", SOURCE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    if (HERE / "endpoint_audit.json").exists():
        raise FileExistsError("Refusing to overwrite completed endpoint audit")
    mod = load_operator()
    inputs = mod.Inputs()
    ambient = lambda t: inputs.ambient(t, True)
    report = {
        "git_commit_read": "anonymous-code-snapshot",
        "scope": "Targeted continuations from previously saved 2D states. Shared frozen FV operator. No full PDE cold start or continuum error bound.",
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "source_sha256": sha(SOURCE), "cases": {}, "sources": {},
    }
    for case, nr, nz, official in (("q23", 80, 256, 206906.76), ("q4", 80, 128, 183931.56)):
        source_npz = GEOM / f"{case}_{nr}x{nz}_integral.npz"
        source_json = source_npz.with_suffix(".json")
        report["sources"][str(source_npz.relative_to(ROOT))] = sha(source_npz)
        report["sources"][str(source_json.relative_to(ROOT))] = sha(source_json)
        saved = np.load(source_npz)
        op = mod.GeometryFV(nr, nz, case, inputs, flux="integral")
        assert saved["event_states"].shape[1:] == op.shape + (2,)
        if case == "q23":
            initial_t = float(saved["event_times_s"][0])
            initial = saved["event_states"][0].ravel()
            evaluation_times = np.unique(np.r_[saved["event_times_s"], official])
            max_step = 0.5
        else:
            prior = np.flatnonzero(saved["snapshot_times_s"] < official)[-1]
            initial_t = float(saved["snapshot_times_s"][prior])
            initial = saved["snapshots"][prior].ravel()
            evaluation_times = np.array([official, *saved["event_times_s"]])
            evaluation_times = np.unique(evaluation_times[evaluation_times >= initial_t])
            max_step = 30.0
            # The existing radius is a constant 0.012 m across this short window.
            checks = np.r_[initial_t, evaluation_times[-1], inputs.rad[(inputs.rad[:, 0] >= initial_t) & (inputs.rad[:, 0] <= evaluation_times[-1]), 0]]
            assert np.ptp([inputs.radius(t) for t in checks]) == 0
        started = time.perf_counter()
        print(f"START {case}: saved {initial_t:.9f} s -> {evaluation_times[-1]:.9f} s", flush=True)
        solution = solve_ivp(lambda t, y: op.rhs(t, y, ambient),
            (initial_t, float(evaluation_times[-1])), initial, method="BDF",
            jac=lambda t, y: op.rhs(t, y, ambient, True),
            rtol=2e-11, atol=np.tile([1e-10, 1e-13], op.m),
            max_step=max_step, t_eval=evaluation_times)
        assert solution.success, solution.message
        states = solution.y.T.reshape((-1,) + op.shape + (2,))
        official_idx = np.flatnonzero(solution.t == official)[0]
        final = states[official_idx]
        c = final[:, :, 1]
        ix = np.unravel_index(np.argmax(c), c.shape)
        item = {
            "initial_saved_time_s": initial_t, "execution_s": official,
            "rtol": 2e-11, "atol_T": 1e-10, "atol_C": 1e-13,
            "max_step_s": max_step, "nfev": solution.nfev, "nlu": solution.nlu,
            "elapsed_s": time.perf_counter() - started,
            "direct_FV_Cmax_at_execution": float(c[ix]),
            "threshold_margin_0p15_minus_Cmax": float(.15-c[ix]),
            "wettest_r_z_m": [float(op.r[ix[1]] * (inputs.radius(official) / .02 if case == "q4" else 1)), float(op.z[ix[0]])],
            "all_nodes_strictly_pass": bool(np.all(c < .15)),
            "radial_max_positive_difference": float(np.max(np.diff(c, axis=1))),
            "axial_max_positive_difference": float(np.max(np.diff(c, axis=0))),
            "upstream_error": "Unchanged saved coarse-grid history; local tight solve does not eliminate upstream spatial/temporal errors.",
        }
        if case in ("q23", "q4"):
            repeated = []
            for t, y in zip(saved["event_times_s"], saved["event_states"]):
                idx = np.flatnonzero(solution.t == t)[0]
                repeated.append({"time_s": float(t), "Cmax_new": float(np.max(states[idx, :, :, 1])),
                    "Cmax_saved": float(np.max(y[:, :, 1])),
                    "max_abs_C_difference": float(np.max(np.abs(states[idx, :, :, 1]-y[:, :, 1]))),
                    "max_abs_T_difference": float(np.max(np.abs(states[idx, :, :, 0]-y[:, :, 0])))})
            item["continuation_crosscheck_with_frozen_events"] = repeated
        np.savez_compressed(HERE / f"{case}_direct_execution_field.npz", time_s=solution.t, states=states,
            reference_r_m=op.r, z_m=op.z, source_initial_time_s=initial_t)
        report["cases"][case] = item
        print(json.dumps(item, ensure_ascii=False), flush=True)

    # A deterministic observed nonuniform Q2 field demonstrates mixed-sign
    # cross-temperature effects. This refutes naive cooperative-system claims,
    # not the observed direction of the actual 1D/2D event shift.
    saved = np.load(GEOM / "q23_80x256_integral.npz")
    ix = int(np.flatnonzero(saved["snapshot_times_s"] == 10800)[0])
    op = mod.GeometryFV(80, 256, "q23", inputs, flux="integral")
    matrix = op.rhs(10800, saved["snapshots"][ix].ravel(), inputs.ambient, True).tocoo()
    select = (matrix.row % 2 == 1) & (matrix.col % 2 == 0)
    values = matrix.data[select]
    report["comparison_principle_diagnostic"] = {
        "state": "saved Q23 80x256 field at 10800 s", "block": "d(Cdot)/d(T)",
        "minimum": float(values.min()), "maximum": float(values.max()),
        "negative_entries": int(np.sum(values < 0)), "positive_entries": int(np.sum(values > 0)),
        "meaning": "The coupled semidiscrete operator is not cooperative in componentwise (T,C) order. A scalar D>0 comparison argument alone does not establish 1D >= 2D moisture for the coupled model.",
    }
    (HERE / "endpoint_audit.json").write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print("COMPLETE", flush=True)


if __name__ == "__main__":
    main()
