"""Conditional future-environment probes from frozen 4 h states.

Reuses the audited effective-model operators; not an independent physical model
or a replacement for any official result. Q4 radius remains the observed input.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
import time
import numpy as np
import scipy
from scipy.integrate import solve_ivp

REPO = Path(__file__).resolve().parents[3]
DELIVERY = REPO / 'q4_complete_delivery/v1'
sys.path.insert(0, str(DELIVERY / 'source'))
sys.path.insert(0, str(DELIVERY / 'q123_closeout/source/legacy_q3'))
import q4_spectral
import q3_reference
from q4_common import Radius, Environment, load_input, read_numeric

SCENARIOS = {'baseline': (0., 0.), 'temperature_minus_1K': (-1., 0.),
             'temperature_plus_1K': (1., 0.), 'boundary_minus_0p005': (0., -.005),
             'boundary_plus_0p005': (0., .005)}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(case, n, scenario, out):
    stem = out / f'{case}_{scenario}_n{n}'
    if stem.with_suffix('.json').exists() or stem.with_suffix('.npz').exists():
        raise FileExistsError(stem)
    out.mkdir(parents=True, exist_ok=True)
    env_path = DELIVERY / 'inputs/attachment1.xlsx'
    future = Environment(load_input(env_path)).future
    boundary = future + np.array(SCENARIOS[scenario])
    env = lambda t: boundary
    if case == 'q3':
        original = DELIVERY / 'q123_closeout/output/q23_unified.npz'
        z = np.load(original)
        idx = int(np.flatnonzero(z['snapshot_time_s'] == 14400.)[0])
        full = z['snapshot_full_TC'][idx]
        op = q3_reference.Reference(n)
        official, official_root = 206906.76, 206906.38993232802
        module = q3_reference
    else:
        original = DELIVERY / 'output/main.npz'
        z = np.load(original)
        idx = int(np.flatnonzero(z['time_s'] == 14400.)[0])
        full = z['full_TC'][idx]
        radius = Radius(read_numeric(DELIVERY / 'inputs/attachment2.xlsx', 2)[1])
        op = q4_spectral.Reference(n, radius)
        official, official_root = 183931.56, 183931.21171548913
        module = q4_spectral
    old_n = len(z['x']) - 1
    # Use only exact nested Chebyshev nodes, avoiding an interpolation error.
    assert old_n % n == 0
    start = full[::old_n // n].copy()
    assert np.max(abs(z['x'][::old_n // n] - op.x)) < 1e-14
    y = start[:-1].ravel()
    saved_t, saved_full = [14400.], [op.surface(14400., y.reshape(n, 2), env)]
    critical = None
    official_full = None
    event_full = None
    nfev = nlu = 0
    tic = time.perf_counter()
    # Q4 never extends the observed radius beyond 72 h in these probes.
    cuts = np.unique(np.r_[np.arange(14400., 259201., 1800.), official])
    for a, b in zip(cuts[:-1], cuts[1:]):
        def event(t, state):
            return float(np.max(state[1::2]) - .15)
        event.direction = -1
        event.terminal = False
        sol = solve_ivp(lambda t, u: op.evaluate(t, u, env), (a, b), y,
                        method='Radau', jac=lambda t, u: op.evaluate(t, u, env, True),
                        rtol=2e-9, atol=np.tile([2e-10, 2e-12], n),
                        max_step=300., dense_output=True, events=event)
        if not sol.success:
            raise RuntimeError(sol.message)
        nfev += sol.nfev
        nlu += sol.nlu
        if critical is None and len(sol.t_events[0]):
            critical = float(sol.t_events[0][0])
            event_full = op.surface(critical, sol.sol(critical).reshape(n, 2), env)
        y = sol.y[:, -1]
        state = op.surface(b, y.reshape(n, 2), env)
        saved_t.append(float(b))
        saved_full.append(state)
        if b == official:
            official_full = state.copy()
        if critical is not None and official_full is not None:
            break
    result = {'case': case, 'scenario': scenario, 'n': n,
              'scope': 'Actual coupled 1D continuation from frozen 4 h; observed Q4 R(t) fixed across environment probes; not a physical confidence interval.',
              'baseline_future_T_b': future.tolist(), 'future_T_b': boundary.tolist(),
              'perturbation_T_b': list(SCENARIOS[scenario]),
              'perturbation_status': 'Analyst-selected stress size, not a measured uncertainty or a confidence interval.',
              'start_s': 14400., 'source_state_node_count': old_n + 1,
              'initial_state_projection': 'Exact nested Chebyshev nodes; same source state for all scenarios.',
              'source_npz': original.relative_to(REPO).as_posix(), 'source_npz_sha256': sha(original),
              'operator_source': Path(module.__file__).relative_to(REPO).as_posix(),
              'operator_sha256': sha(Path(module.__file__)),
              'environment_input_sha256': sha(env_path), 'script_sha256': sha(Path(__file__)),
              'method': 'Radau', 'rtol': 2e-9, 'atol_TC': [2e-10, 2e-12], 'max_step_s': 300.,
              'critical_s': critical, 'critical_h': critical / 3600 if critical else None,
              'difference_from_frozen_official_root_s': critical - official_root if critical else None,
              'polynomial_extrema_at_event': op.extrema(event_full) if event_full is not None else None,
              'official_execution_s': official,
              'Cmax_at_official_execution': op.extrema(official_full)['max_C'] if official_full is not None else None,
              'strict_at_official_execution_in_this_probe': bool(op.extrema(official_full)['max_C'] < .15) if official_full is not None else None,
              'max_boundary_residual': op.max_bc_residual,
              'nfev': nfev, 'nlu': nlu, 'elapsed_s': time.perf_counter() - tic,
              'runtime': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__}}
    assert critical is not None and official_full is not None, 'No event within observed radius span; report requires explicit unresolved state.'
    states = np.asarray(saved_full)
    assert np.isfinite(states).all() and states[:, :, 1].min() > 0
    np.savez_compressed(stem.with_suffix('.npz'), time_s=np.asarray(saved_t), full_TC=states,
                        x=op.x, source_initial_full_TC=start, event_full_TC=event_full,
                        official_full_TC=official_full)
    stem.with_suffix('.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--case', choices=['q3', 'q4'], required=True)
    p.add_argument('--n', type=int, default=40)
    p.add_argument('--scenario', choices=list(SCENARIOS) + ['all'], default='all')
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    for scenario in SCENARIOS if a.scenario == 'all' else [a.scenario]:
        run(a.case, a.n, scenario, a.out)
