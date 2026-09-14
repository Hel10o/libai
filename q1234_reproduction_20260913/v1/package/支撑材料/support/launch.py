#!/usr/bin/env python3
"""Anonymous source snapshot: safe smoke by default, explicit numerical runs.

Run from the paper directory: python support/launch.py smoke
Use `list` to see numerical tasks. A `run` requires a new output directory.
This launcher never treats syntax checks as numerical validation.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PROJECT = HERE / "project"

TASKS = {
    "q1-main": "Q1 10240-cell main solution; includes independent flux quadrature",
    "q1-full": "Q1 complete solver/validation/refinement/workbook/report chain",
    "q2-sensitivity": "Seven 320-cell FV runs: baseline, hm/hT +/-20%, PCHIP, frozen properties",
    "q23-main": "Current Q2/Q3 unified trajectory, 80-node Radau, one-second records",
    "q3-main": "Current Q3 80-node trajectory, threshold budget and strict execution",
    "q4-main": "Current Q4 120-node shrinking-radius main solution and balance audit",
}

SAFE_HELP = [
    "q1_complete_delivery/q1_delivery/run_all.py",
    *["q1_complete_delivery/q1_delivery/source/" + p for p in
      ("q1_solver.py", "q1_validate.py", "q1_refine.py", "q1_excel.py", "q1_report.py")],
    "q2_refinement_delivery/runtime/source/q2_core.py",
    "q2_refinement_delivery/runtime/source/q2_spectral.py",
    "q2_refinement_delivery/source/prepare_reused.py",
    "q2_refinement_delivery/source/q2_mechanism.py",
    "q4_complete_delivery/v1/q123_closeout/source/solve_unified_q23.py",
    "q4_complete_delivery/v1/q123_closeout/source/audit_unified_q23.py",
    *["q3_refinement_delivery/source/" + p for p in
      ("q3_refined_reference.py", "future_scenarios_from_checkpoint.py",
       "export_workbook.py", "validate_xlsx.py", "verify_frozen.py")],
    *["q4_complete_delivery/v1/source/" + p for p in
      ("run_q4.py", "geometry_solver.py", "geometry_solver_executed.py",
       "geometry_checks.py", "check_numerics.py", "finalize_q4.py", "package_workbooks.py")],
    "q1234_overall_review_delivery/v1/source/environment_probe.py",
]


def child_environment():
    env = os.environ.copy()
    env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONUTF8="1",
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    return env


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def smoke(log=None, node=None):
    """Parse every Python file and import only through reviewed --help paths."""
    manifest = json.loads((HERE / "snapshot_manifest.json").read_text(encoding="utf-8"))
    checks = []
    for entry in manifest["files"]:
        p = PROJECT / entry["path"]
        assert p.is_file() and sha(p) == entry["snapshot_sha256"], entry["path"]
        if p.suffix == ".py":
            ast.parse(p.read_text(encoding="utf-8-sig"), filename=entry["path"])
    ast.parse(Path(__file__).read_text(encoding="utf-8"), filename="support/launch.py")
    for name in SAFE_HELP:
        p = subprocess.run([sys.executable, str(PROJECT / name), "--help"],
                           cwd=PROJECT, env=child_environment(), capture_output=True,
                           text=True, encoding="utf-8", errors="replace", timeout=45)
        # Store relative identity and status only; do not leak runtime home paths.
        checks.append({"path": name, "check": "--help", "returncode": p.returncode})
        if p.returncode:
            raise RuntimeError(f"Safe entry check failed: {name}; return code {p.returncode}")
    if node:
        name = "q4_complete_delivery/v1/source/export_workbooks.mjs"
        p = subprocess.run([node, "--check", str(PROJECT / name)], cwd=PROJECT,
                           capture_output=True, text=True, timeout=30)
        checks.append({"path": name, "check": "node --check", "returncode": p.returncode})
        if p.returncode:
            raise RuntimeError("JavaScript syntax check failed")
    versions = {"python": platform.python_version()}
    for package in ("numpy", "scipy", "matplotlib", "openpyxl"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not installed"
    result = {"status": "PASS", "scope": "hash/AST/safe --help only; no PDE or workbook export",
              "snapshot_files_checked": len(manifest["files"]),
              "python_files_parsed": 1 + sum(e["path"].endswith(".py") for e in manifest["files"]),
              "versions": versions, "entry_checks": checks}
    if log:
        log = Path(log)
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=True, indent=2))
    return result


def command(task, out):
    q1 = PROJECT / "q1_complete_delivery/q1_delivery"
    q23 = PROJECT / "q4_complete_delivery/v1/q123_closeout"
    q3 = PROJECT / "q3_refinement_delivery"
    q4 = PROJECT / "q4_complete_delivery/v1"
    if task == "q1-full":
        return [str(q1 / "run_all.py"), "--out", str(out)]
    if task == "q1-main":
        return [str(q1 / "source/q1_solver.py"), "--input", str(q1 / "input/附件1.xlsx"),
                "--out", str(out), "--quadrature"]
    if task == "q23-main":
        return [str(q23 / "source/solve_unified_q23.py"), "--input", str(q23 / "inputs/attachment1.xlsx"),
                "--out", str(out / "q23_unified"), "--n", "80", "--method", "Radau",
                "--rtol", "2e-13", "--max-step", "30", "--step", "1", "--end", "206906.76"]
    if task == "q3-main":
        return [str(q3 / "source/q3_refined_reference.py"), "--input", str(q3 / "inputs/attachment1.xlsx"),
                "--out", str(out / "solution"), "--n", "80", "--scenario", "mean",
                "--method", "Radau", "--rtol", "2e-11", "--atol-c", "2e-13",
                "--max-step", "120", "--budget-s", "0.03", "--quantum-s", "0.36"]
    if task == "q4-main":
        return [str(q4 / "source/run_q4.py"), "--out", str(out / "main"), "--n", "120",
                "--kind", "linear", "--method", "Radau", "--rtol", "2e-11", "--max-step", "90", "--audit"]
    raise ValueError(task)


def run(task, out):
    out = Path(out).resolve()
    if out.exists():
        raise FileExistsError("Choose a new output directory; existing results are never overwritten.")
    if out == PROJECT or PROJECT in out.parents:
        raise ValueError("Numerical output must be outside the immutable project snapshot.")
    out.mkdir(parents=True)
    print(TASKS[task], flush=True)
    if task == "q2-sensitivity":
        for key, value in child_environment().items():
            os.environ[key] = value
        sys.path.insert(0, str(PROJECT / "q2_refinement_delivery/runtime/source"))
        from q2_core import Environment, load_environment
        from q2_tests import run_sensitivity
        data, _ = load_environment(PROJECT / "q2_refinement_delivery/runtime/inputs/附件1.xlsx")
        run_sensitivity(Environment(data), out)
        return
    subprocess.run([sys.executable, *command(task, out)], cwd=PROJECT,
                   env=child_environment(), check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action")
    s = sub.add_parser("smoke", help="safe default: no numerical solving")
    s.add_argument("--log", type=Path)
    s.add_argument("--node", help="optional Node executable for JavaScript syntax only")
    sub.add_parser("list", help="list explicit, potentially expensive numerical tasks")
    r = sub.add_parser("run", help="explicit numerical integration, never the default")
    r.add_argument("task", choices=TASKS)
    r.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.action == "run":
        run(args.task, args.out)
    elif args.action == "list":
        for name, description in TASKS.items():
            print(f"{name}: {description}")
    else:
        smoke(getattr(args, "log", None), getattr(args, "node", None))


if __name__ == "__main__":
    main()
