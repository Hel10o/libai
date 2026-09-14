"""Independent postprocessing of 14 future-environment probe result pairs.

Does not import the probe or solver operators and does not run any PDE.
Reconstructs primitive polynomials with least-squares Chebyshev interpolation
(the solver used a DCT) and independently checks Robin residuals and maxima.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import numpy as np
from numpy.polynomial import Chebyshev
from scipy.optimize import brentq
from scipy.special import exp1
from openpyxl import load_workbook

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BASE = ROOT / "q1234_overall_review_delivery/v1"
DATA = BASE / "evidence/environment"
SCRIPT = BASE / "source/environment_probe.py"
OUT = HERE / "environment_evidence_audit.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_xlsx(path, columns):
    book = load_workbook(path, read_only=True, data_only=True)
    rows = [r[:columns] for r in list(book.active.values)[1:] if r[0] is not None]
    book.close()
    return np.asarray(rows, float)


def primitive(c, b):
    return c*np.exp(-b/c)-b*exp1(b/c)


def check_field(case, x, state, time, ambient, radius_data):
    n = len(x)-1
    b = .45 if case == "q3" else .30
    t, c = state.T
    u = primitive(c, b)
    # Fit around a constant reference to reduce cancellation in derivatives.
    p = Chebyshev.fit(x, u-u[-1], n, domain=[0, 1])
    dt = Chebyshev.fit(x, t-t[-1], n, domain=[0, 1]).deriv()(1)
    du = p.deriv()(1)
    roots = p.deriv().roots()
    roots = roots[np.abs(roots.imag) < 1e-9].real
    roots = roots[(roots > 1e-12) & (roots < 1-1e-12)]
    candidates = np.r_[0, roots, 1]
    values = p(candidates)+u[-1]
    index = np.argmax(values)
    maximum = brentq(lambda value: primitive(value, b)-values[index], 1e-7, 10., xtol=5e-15)
    radius = .02 if case == "q3" else np.interp(time, radius_data[:, 0], radius_data[:, 1])*.01
    if case == "q3":
        k = .21+.38*c[-1]/(1+c[-1]); prefactor=.0024
    else:
        k = .12+.20*c[-1]/(1+c[-1]); prefactor=.00042
    thermal = 2/radius*k*dt+25*(t[-1]-ambient[0])
    mass = 2/radius*prefactor*np.exp(-3850/(t[-1]+273.15))*du+8e-7*(c[-1]-ambient[1])
    return {"max_C": float(maximum), "x_at_max": float(candidates[index]),
            "stationary_points": int(len(roots)), "node_max_C": float(c.max()),
            "heat_boundary_residual_W_m2": float(thermal),
            "moisture_boundary_residual_m_s": float(mass), "radius_m": float(radius)}


def main():
    if OUT.exists():
        raise FileExistsError(OUT)
    env_path = ROOT / "q4_complete_delivery/v1/inputs/attachment1.xlsx"
    rad_path = ROOT / "q4_complete_delivery/v1/inputs/attachment2.xlsx"
    env = read_xlsx(env_path, 3)
    rad = read_xlsx(rad_path, 2)
    future = env[env[:, 0]>=10800, 1:].mean(axis=0)
    paths = sorted(DATA.glob("*.json"))
    assert len(paths) == 14 and len(list(DATA.glob("*.npz"))) == 14
    scenarios = {"baseline": [0.,0.], "temperature_minus_1K": [-1.,0.],
                 "temperature_plus_1K": [1.,0.], "boundary_minus_0p005": [0.,-.005],
                 "boundary_plus_0p005": [0.,.005]}
    rows=[]; hash_records={}; originals={}; starts={}
    for path in paths:
        info=json.loads(path.read_text(encoding="utf-8"))
        raw=np.load(path.with_suffix(".npz"))
        case,n=info["case"],info["n"]
        original=ROOT/info["source_npz"]
        for file in (path,path.with_suffix(".npz"), original, ROOT/info["operator_source"], SCRIPT,env_path,rad_path):
            hash_records[str(file.relative_to(ROOT))]=sha(file)
        assert sha(SCRIPT)==info["script_sha256"]
        assert sha(original)==info["source_npz_sha256"]
        assert sha(ROOT/info["operator_source"])==info["operator_sha256"]
        assert sha(env_path)==info["environment_input_sha256"]
        assert np.array_equal(future,np.asarray(info["baseline_future_T_b"]))
        ambient=future+np.asarray(scenarios[info["scenario"]])
        assert np.array_equal(ambient,np.asarray(info["future_T_b"]))
        if original not in originals:
            source=np.load(original)
            key="snapshot_time_s" if case=="q3" else "time_s"
            fieldkey="snapshot_full_TC" if case=="q3" else "full_TC"
            idx=int(np.flatnonzero(source[key]==14400)[0])
            originals[original]=(source["x"],source[fieldkey][idx])
        old_x,old_field=originals[original]
        step=(len(old_x)-1)//n
        expected=old_field[::step]
        assert expected.shape==(n+1,2)
        assert np.array_equal(expected,raw["source_initial_full_TC"])
        assert np.max(np.abs(old_x[::step]-raw["x"]))<1e-14
        assert np.array_equal(raw["full_TC"][0,:-1],expected[:-1])
        pair=(case,n)
        if pair in starts:assert np.array_equal(starts[pair],raw["source_initial_full_TC"])
        starts[pair]=raw["source_initial_full_TC"].copy()
        assert raw["time_s"][0]==14400 and np.all(np.diff(raw["time_s"])>0)
        ix=int(np.flatnonzero(raw["time_s"]==info["official_execution_s"])[0])
        assert np.array_equal(raw["full_TC"][ix],raw["official_full_TC"])
        assert 14400 < info["critical_s"] <= raw["time_s"][-1] <= 259200
        event=check_field(case,raw["x"],raw["event_full_TC"],info["critical_s"],ambient,rad)
        official=check_field(case,raw["x"],raw["official_full_TC"],info["official_execution_s"],ambient,rad)
        initial=check_field(case,raw["x"],raw["full_TC"][0],14400,ambient,rad)
        assert abs(event["max_C"]-.15)<1e-10
        assert abs(event["max_C"]-info["polynomial_extrema_at_event"]["max_C"])<1e-10
        assert abs(official["max_C"]-info["Cmax_at_official_execution"])<1e-10
        assert (official["max_C"]<.15)==info["strict_at_official_execution_in_this_probe"]
        for diagnostic in (event,official,initial):
            assert abs(diagnostic["heat_boundary_residual_W_m2"])<1e-6
            assert abs(diagnostic["moisture_boundary_residual_m_s"])<1e-14
        rows.append({"name":path.stem,"case":case,"scenario":info["scenario"],"n":n,
            "critical_s":info["critical_s"],"critical_h":info["critical_h"],
            "event_independent":event,"official_independent":official,
            "initial_independent":initial,
            "surface_TC_reset_from_source_at_4h":(raw["full_TC"][0,-1]-expected[-1]).tolist(),
            "source_initial_interior_max_difference":float(np.max(np.abs(raw["full_TC"][0,:-1]-expected[:-1]))),
            "baseline_error_vs_frozen_root_s":info["difference_from_frozen_official_root_s"] if info["scenario"]=="baseline" else None})
    by={(r["case"],r["scenario"],r["n"]):r for r in rows}
    for r in rows:
        base=by.get((r["case"],"baseline",r["n"]))
        r["difference_from_same_n_baseline_s"]=r["critical_s"]-base["critical_s"]
        r["difference_from_same_n_baseline_h"]=r["difference_from_same_n_baseline_s"]/3600
    refinement={}
    for case,n in (("q3",80),("q4",60)):
        base40=by[case,"baseline",40];basefine=by[case,"baseline",n]
        low40=by[case,"temperature_minus_1K",40];lowfine=by[case,"temperature_minus_1K",n]
        shift40=low40["critical_s"]-base40["critical_s"]
        shiftfine=lowfine["critical_s"]-basefine["critical_s"]
        refinement[case]={"coarse_n":40,"fine_n":n,
            "baseline_refinement_change_s":basefine["critical_s"]-base40["critical_s"],
            "minus1K_refinement_change_s":lowfine["critical_s"]-low40["critical_s"],
            "minus1K_effect_coarse_s":shift40,"minus1K_effect_fine_s":shiftfine,
            "effect_refinement_change_s":shiftfine-shift40,"effect_fine_h":shiftfine/3600,
            "fine_baseline_error_vs_frozen_root_s":basefine["baseline_error_vs_frozen_root_s"]}
    report={"status":"PASS","git_commit_read":"anonymous-code-snapshot",
        "scope":"Independent artifact and polynomial/boundary postprocessing only. No PDE solve, no independent physical calibration.",
        "pairs_checked":14,"scenarios_n40":10,"refinement_runs":4,
        "records":rows,"refinement":refinement,"sha256":hash_records,
        "limitations":["Nested node restriction is exact at retained nodes, not an error-free lower-degree representation.",
            "Only baseline and minus-1 K were spatially refined; no temporal refinement was added in this audit.",
            "Events were detected at interior nodes; independent whole-radial primitive-polynomial extrema validate the saved event states, not continuum maxima.",
            "4 h future-platform changes are step interventions. Surface is algebraically reset; bulk state is unchanged.",
            "Q4 radius is held to observed R(t) under each changed environment, so results are conditional, not feedback predictions."]}
    OUT.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":"PASS","pairs_checked":14,"refinement":refinement},ensure_ascii=False,indent=2))
    for r in rows: print(f'{r["name"]}: {r["critical_h"]:.10f} h; delta same-grid baseline {r["difference_from_same_n_baseline_s"]:.9f} s; official Cmax {r["official_independent"]["max_C"]:.12f}')


if __name__=="__main__":main()
