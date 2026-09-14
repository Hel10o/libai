"""Refined Q3 reference solver used for the frozen refinement delivery.

Main features relative to the historical reference solver:
- future-environment scenario is explicit (mean / nominal / time_mean / last);
- strict execution endpoint is derived automatically from the event root,
  numerical budget and output quantum; no historical endpoint is hard-coded;
- the execution endpoint is actually evaluated on the same continuous solution;
- all official 60 s / 0.1 cm output values come from one unrounded trajectory;
- output refuses overwrite if either JSON or NPZ already exists.

The spatial method remains the audited global Chebyshev collocation in x=(r/R)^2
with primitive-gradient moisture flux and algebraic Robin surface conditions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import tempfile
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from q3_reference import Reference
from q3_solver import Environment, load_input


class RefinedEnvironment(Environment):
    def __init__(self, data, scenario="mean", kind="linear"):
        if scenario == "time_mean":
            # Initialize interpolation machinery without accepting an unknown scenario.
            super().__init__(data, scenario="mean", kind=kind)
            mask = self.data[:, 0] >= 10800.0
            tail = self.data[mask]
            duration = tail[-1, 0] - tail[0, 0]
            self.future = np.trapezoid(tail[:, 1:], tail[:, 0], axis=0) / duration
            self.scenario = scenario
        else:
            super().__init__(data, scenario=scenario, kind=kind)


def _atomic_npz(path: Path, **arrays) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp.npz", dir=path.parent)
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        np.savez_compressed(tmp, **arrays)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def _atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        tmp.write_text(text, encoding="utf-8")
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def derive_execution_time(root_s: float, budget_s: float, quantum_s: float) -> float:
    if budget_s < 0 or quantum_s <= 0:
        raise ValueError("budget_s must be >= 0 and quantum_s must be > 0")
    # Tiny subtraction avoids an accidental extra quantum from binary representation
    # when the mathematical quotient is already an integer.
    return quantum_s * math.ceil((root_s + budget_s) / quantum_s - 1e-13)


def run(
    *,
    n: int,
    out: Path,
    input_file: Path,
    scenario: str = "mean",
    method: str = "Radau",
    rtol: float = 2e-11,
    atol_c: float = 2e-13,
    max_step: float = 120.0,
    budget_s: float = 0.03,
    quantum_s: float = 0.36,
    interpolation: str = "linear",
):
    json_path = out.with_suffix(".json")
    npz_path = out.with_suffix(".npz")
    if json_path.exists() or npz_path.exists():
        raise FileExistsError(f"Refusing overwrite: {json_path} / {npz_path}")

    data = load_input(input_file)
    env = RefinedEnvironment(data, scenario=scenario, kind=interpolation)
    op = Reference(n)
    y = np.tile([28.0, 2.55], n)

    # Segment at all observed boundary knots; after 4 h, use wider blocks only.
    cuts = np.unique(np.r_[data[:, 0], np.arange(21600.0, 360001.0, 21600.0)])
    times = [0.0]
    profiles = [np.tile([28.0, 2.55], (21, 1))]
    snapshot_times = [0.0]
    snapshot_full = [np.tile([28.0, 2.55], (n + 1, 1))]
    event_full = []
    event_times = []
    extrema_log = []
    next_output = 60.0

    tic = time.perf_counter()
    nsteps = nfev = nlu = 0
    critical = None
    execution = None
    execution_full = None

    for a, b in zip(cuts[:-1], cuts[1:]):
        seg_env = env.segment(a, b)
        fun = lambda t, yy: op.evaluate(t, yy, seg_env)

        def event(t, yy):
            return float(np.max(yy[1::2]) - 0.15)

        event.terminal = False
        event.direction = -1

        sol = solve_ivp(
            fun,
            (a, b),
            y,
            method=method,
            jac=lambda t, yy: op.evaluate(t, yy, seg_env, True),
            rtol=rtol,
            atol=np.tile([2e-11, atol_c], n),
            max_step=min(max_step, 30.0) if a < 14400.0 else max_step,
            dense_output=True,
            events=event,
        )
        if not sol.success:
            raise RuntimeError(sol.message)
        nsteps += len(sol.t) - 1
        nfev += sol.nfev
        nlu += sol.nlu

        segment_end = b
        if len(sol.t_events[0]):
            root = float(sol.t_events[0][0])
            execution = derive_execution_time(root, budget_s, quantum_s)
            if execution > b:
                raise RuntimeError("Derived execution endpoint crosses unresolved integration segment")
            segment_end = execution

            y_root = sol.sol(root)
            root_full = op.surface(root, y_root.reshape(n, 2), seg_env)
            y_minus = sol.sol(root - 0.01)
            y_plus = sol.sol(root + 0.01)
            minus_full = op.surface(root - 0.01, y_minus.reshape(n, 2), seg_env)
            plus_full = op.surface(root + 0.01, y_plus.reshape(n, 2), seg_env)
            y_exec = sol.sol(execution)
            execution_full = op.surface(execution, y_exec.reshape(n, 2), seg_env)

            critical = {
                "critical_s": root,
                "critical_h": root / 3600.0,
                "t_minus_s": root - 0.01,
                "max_minus": float(np.max(minus_full[:, 1])),
                "critical_max": float(np.max(root_full[:, 1])),
                "t_plus_s": root + 0.01,
                "max_plus": float(np.max(plus_full[:, 1])),
                "critical_slope_kgkg_s": float(fun(root, y_root)[1::2][np.argmax(y_root[1::2])]),
                "execution_s": float(execution),
                "execution_h": float(execution / 3600.0),
                "execution_max": float(np.max(execution_full[:, 1])),
                "execution_argmax_x": float(op.extrema(execution_full)["x_max"]),
                "strict_execution_pass": bool(np.max(execution_full[:, 1]) < 0.15),
                "budget_s": budget_s,
                "quantum_s": quantum_s,
            }
            if not critical["strict_execution_pass"]:
                raise AssertionError("Derived execution endpoint is not strictly qualified")
            event_times = [root - 0.01, root, root + 0.01, execution]
            event_full = [minus_full, root_full, plus_full, execution_full]

        out_times = np.arange(next_output, segment_end + 1e-8, 60.0)
        # For the event segment the execution endpoint is not an integer minute.
        if critical is not None:
            out_times = out_times[out_times < execution - 1e-9]
        for t in out_times:
            full = op.surface(t, sol.sol(t).reshape(n, 2), seg_env)
            times.append(float(t))
            profiles.append(op.profiles(full))
        if len(out_times):
            next_output = 60.0 * (np.floor(out_times[-1] / 60.0) + 1.0)

        if b in (1800.0, 3600.0, 5400.0, 7200.0, 9000.0, 10800.0, 14400.0, 21600.0) or b % 21600.0 == 0 or critical is not None:
            snap_t = segment_end if critical is not None else b
            full = op.surface(snap_t, sol.sol(snap_t).reshape(n, 2), seg_env)
            snapshot_times.append(float(snap_t))
            snapshot_full.append(full)
            extrema_log.append({"t_s": float(snap_t), **op.extrema(full)})

        y = sol.sol(segment_end)
        if critical is not None:
            times.append(float(execution))
            profiles.append(op.profiles(execution_full))
            break

    if critical is None:
        raise RuntimeError("Threshold event not reached")

    stats = {
        "model": "global Chebyshev collocation in x=(r/R)^2 with primitive-gradient moisture flux",
        "n": n,
        "method": method,
        "rtol": rtol,
        "atol_C": atol_c,
        "max_step_s": max_step,
        "scenario": scenario,
        "interpolation": interpolation,
        "future": env.future.tolist(),
        "event": critical,
        "elapsed_s": time.perf_counter() - tic,
        "nsteps": nsteps,
        "nfev": nfev,
        "nlu": nlu,
        "min_scalar_boundary_derivative": op.min_boundary_derivative,
        "max_boundary_residual": op.max_bc_residual,
        "max_boundary_newton": op.max_newton,
        "boundary_fallbacks": op.boundary_fallbacks,
        "extrema": extrema_log,
        "input_sha256": hashlib.sha256(input_file.read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }

    arrays = {
        "time_s": np.asarray(times),
        "radius_cm": np.linspace(0.0, 2.0, 21),
        "profile_TC": np.asarray(profiles),
        "x": op.x,
        "snapshot_time_s": np.asarray(snapshot_times),
        "snapshot_full_TC": np.asarray(snapshot_full),
        "event_time_s": np.asarray(event_times),
        "event_full_TC": np.asarray(event_full),
    }
    _atomic_npz(npz_path, **arrays)
    _atomic_text(json_path, json.dumps(stats, ensure_ascii=False, indent=2) + "\n")
    return stats


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--n", type=int, default=80)
    p.add_argument("--scenario", choices=["mean", "nominal", "time_mean", "last"], default="mean")
    p.add_argument("--method", default="Radau")
    p.add_argument("--rtol", type=float, default=2e-11)
    p.add_argument("--atol-c", type=float, default=2e-13)
    p.add_argument("--max-step", type=float, default=120.0)
    p.add_argument("--budget-s", type=float, default=0.03)
    p.add_argument("--quantum-s", type=float, default=0.36)
    p.add_argument("--interpolation", choices=["linear", "pchip"], default="linear")
    a = p.parse_args()
    result = run(
        n=a.n,
        out=a.out,
        input_file=a.input,
        scenario=a.scenario,
        method=a.method,
        rtol=a.rtol,
        atol_c=a.atol_c,
        max_step=a.max_step,
        budget_s=a.budget_s,
        quantum_s=a.quantum_s,
        interpolation=a.interpolation,
    )
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()
