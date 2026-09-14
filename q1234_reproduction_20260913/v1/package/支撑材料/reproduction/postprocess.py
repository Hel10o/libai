"""Complete dependent calculations, export and compare without external data."""
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from reproduce import execute, save, sha


def extra_tasks(package, root):
    q4 = root / "q4_complete_delivery/v1"
    q3 = root / "q3_refinement_delivery"
    review = root / "q1234_overall_review_delivery/v1"
    tasks = {}
    def add(name, cmd, products):
        tasks[name] = ([str(v) for v in cmd], list(products))
    destination = root / "q2_refinement_delivery/reproduced_sensitivity"
    add("q2-sensitivity", [package / "reproduction/supplement.py", "--root", root,
                            "--out", destination, "--case", "q2-sensitivity"],
        [destination / "sensitivity.json", destination / "sensitivity_baseline.npz"])
    for case, refined_n in [("q3", 80), ("q4", 60)]:
        for n, scenarios in [(40, ["all"]), (refined_n, ["baseline", "temperature_minus_1K"])]:
            for scenario in scenarios:
                names = ["baseline", "temperature_minus_1K", "temperature_plus_1K", "boundary_minus_0p005", "boundary_plus_0p005"] if scenario == "all" else [scenario]
                add(f"environment-{case}-{n}-{scenario}",
                    [review / "source/environment_probe.py", "--case", case, "--n", n,
                     "--scenario", scenario, "--out", review / "evidence/environment"],
                    [review / "evidence/environment" / f"{case}_{s}_n{n}{ext}" for s in names for ext in (".json", ".npz")])
    add("q3-time-mean", [q3 / "source/q3_refined_reference.py", "--input", q3 / "inputs/attachment1.xlsx",
                         "--out", q3 / "validation/time_mean", "--n", 80, "--scenario", "time_mean",
                         "--rtol", "2e-11", "--max-step", 120],
        [q3 / "validation" / ("time_mean" + ext) for ext in (".json", ".npz")])
    add("q4-fv640-extra-tight", [q4 / "source/geometry_solver.py", "--case", "q4", "--nr", 640,
                                "--nz", 0, "--flux", "integral", "--rtol", "2e-11", "--max-step", 120,
                                "--input-dir", q4 / "inputs", "--out-dir", q4 / "validation/geometry"],
        [q4 / "validation/geometry" / ("q4_640x0_integral_rtol2e-11" + ext) for ext in (".json", ".npz")])
    add("geometry-algebra-checks", [q4 / "source/geometry_checks.py", "--out-dir", q4 / "validation/algebra"],
        [q4 / "validation/algebra" / s for s in ("jacobian_flux_checks.json", "dimension_reduction_check.json")])
    add("material-motion-72h", [package / "reproduction/material_motion.py", "--root", root,
                              "--out", q4 / "validation/material_motion.json"],
        [q4 / "validation/material_motion.json"])
    add("geometry-endpoint-continuations", [review / "reviews/geometry/audit_geometry_endpoints.py"],
        [review / "reviews/geometry/endpoint_audit.json", review / "reviews/geometry/q23_direct_execution_field.npz",
         review / "reviews/geometry/q4_direct_execution_field.npz"])
    return tasks


def finish(package, out, root, signature):
    q4 = root / "q4_complete_delivery/v1"
    # All budget terms are regenerated; finalize refuses a state not actually evaluated.
    execute("q4-finalize", [str(q4 / "source/finalize_q4.py"), "--base", str(q4)],
            [q4 / "output/end_event.json", q4 / "output/table6_unrounded.csv"], root, out, signature)
    tasks = extra_tasks(package, root)
    errors = []
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = {pool.submit(execute, name, cmd, products, root, out, signature): name
                   for name, (cmd, products) in tasks.items()}
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as error:
                errors.append(str(error))
                print("FAIL " + str(error), flush=True)
    if errors:
        raise RuntimeError("\n".join(errors))
    execute("q4-numeric-comparison", [str(q4 / "source/check_numerics.py"), "--base", str(q4)],
            [q4 / "validation/numeric_comparison.json"], root, out, signature)
    physics = root / "q1234_overall_review_delivery/v1/reviews/physics"
    # This diagnostic computes quantities from newly generated fields, not historical arrays.
    execute("physics-diagnostics", [str(physics / "diagnose_saved_fields.py")],
            [physics / "saved_field_physics_diagnostics.json"], root, out, signature)
    from reproduction.verify import verify
    result = verify(package, root, out)
    save(out / "comparison.json", result)
    if not result["pass"]:
        raise AssertionError("Paper comparison failed; inspect comparison.json")
    from reproduction.portable_workbooks import export
    export(root, out / "results")
    from reproduction.verify import verify_exports
    book_record = out / "workbook_verification.json"
    if book_record.exists():
        prior = json.loads(book_record.read_text(encoding="utf-8"))
        assert prior["pass"] and prior.get("export_manifest_sha256") == sha(out / "results/export_manifest.json")
    else:
        save(book_record, verify_exports(root, out))
    from reproduction.figures import build
    build(root, out)
    print("COMPLETE: paper comparisons and portable workbooks passed", flush=True)
