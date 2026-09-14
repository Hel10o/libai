"""Recompute the paper from packaged source and inputs in a new directory.

Usage: python -X utf8 reproduce.py --out ../reproduced --jobs 3
No saved answer is used as a numerical initial state. The results/reference
directories are read only by the final comparison, never by PDE solvers.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def prepare(out):
    source = HERE / "support/project"
    manifest = {p.relative_to(source).as_posix(): sha(p) for p in source.rglob("*") if p.is_file()}
    record = out / "source_manifest.json"
    if record.exists():
        if json.loads(record.read_text(encoding="utf-8")) != manifest:
            raise RuntimeError("Packaged source/input changed; choose a new output directory")
        for name, digest in manifest.items():
            if sha(out / "project" / name) != digest:
                raise RuntimeError("Working source/input changed: " + name)
    else:
        if out.exists() and any(out.iterdir()):
            raise FileExistsError("Choose an empty output directory")
        out.mkdir(parents=True, exist_ok=True)
        shutil.copytree(source, out / "project")
        save(record, manifest)
    # The preserved executed geometry snapshot uses the original attachment
    # layout. Create byte-identical aliases from packaged inputs, never read
    # the author's external working directory.
    aliases = out / "project/A题/附件"
    aliases.mkdir(parents=True, exist_ok=True)
    for original, target in [("attachment1.xlsx", "附件1.xlsx"), ("attachment2.xlsx", "附件2.xlsx")]:
        src = source / "q4_complete_delivery/v1/inputs" / original
        dst = aliases / target
        if dst.exists():
            if sha(dst) != sha(src):
                raise RuntimeError("Input alias differs from packaged attachment")
        else:
            shutil.copyfile(src, dst)
    return out / "project", manifest


def commands(root):
    q1 = root / "q1_complete_delivery/q1_delivery"
    q3 = root / "q3_refinement_delivery"
    q4 = root / "q4_complete_delivery/v1"
    q23 = q4 / "q123_closeout"
    geom = q4 / "validation/geometry"
    tasks = {}
    def add(name, args, products):
        tasks[name] = ([str(a) for a in args], [Path(p) for p in products])
    def geometry(name, case, nr, nz, flux="integral", rtol="2e-9", step="240"):
        label = f"{case}_{nr}x{nz}" + ("_integral" if flux == "integral" else "")
        if float(rtol) != 2e-9:
            label += "_rtol" + format(float(rtol), "g")
        add(name, [q4 / "source/geometry_solver.py", "--case", case, "--nr", nr, "--nz", nz,
                   "--flux", flux, "--rtol", rtol, "--max-step", step,
                   "--input-dir", q4 / "inputs", "--out-dir", geom],
            [geom / (label + ext) for ext in (".json", ".npz")])
    # Start the two longest independent calculations before the small jobs.
    geometry("geometry-q23-2d", "q23", 80, 256)
    geometry("geometry-q4-2d", "q4", 80, 128)
    add("q1-main", [q1 / "source/q1_solver.py", "--input", q1 / "input/附件1.xlsx",
                    "--out", q1 / "output", "--n", 10240, "--quadrature"],
        [q1 / "output" / x for x in ("q1_unrounded.npz", "run.json", "table_temperature.csv", "table_moisture.csv")])
    add("q23-main", [q23 / "source/solve_unified_q23.py", "--input", q23 / "inputs/attachment1.xlsx",
                     "--out", q23 / "output/q23_unified", "--n", 80, "--method", "Radau",
                     "--rtol", "2e-13", "--max-step", 30, "--step", 1, "--end", "206906.76"],
        [q23 / "output" / ("q23_unified" + ext) for ext in (".json", ".npz")])
    add("q3-main", [q3 / "source/q3_refined_reference.py", "--input", q3 / "inputs/attachment1.xlsx",
                    "--out", q3 / "output/solution", "--n", 80, "--scenario", "mean", "--method", "Radau",
                    "--rtol", "2e-11", "--atol-c", "2e-13", "--max-step", 120,
                    "--budget-s", "0.03", "--quantum-s", "0.36"],
        [q3 / "output" / ("solution" + ext) for ext in (".json", ".npz")])
    for name, stem, n, kind, method, rtol, step, audit in [
        ("q4-main", "output/main", 120, "linear", "Radau", "2e-11", 90, True),
        ("q4-spectral80", "validation/spectral80", 80, "linear", "Radau", "2e-10", 180, False),
        ("q4-spectral160", "validation/spectral160", 160, "linear", "Radau", "2e-11", 90, False),
        ("q4-time120-bdf", "validation/time120_bdf", 120, "linear", "BDF", "1e-11", 45, False),
        ("q4-fixed80", "validation/fixed80", 80, "fixed", "Radau", "2e-10", 180, False),
        ("q4-pchip80", "validation/pchip80", 80, "pchip", "Radau", "2e-10", 180, False),
    ]:
        add(name, [q4 / "source/run_q4.py", "--out", q4 / stem, "--n", n, "--kind", kind,
                   "--method", method, "--rtol", rtol, "--max-step", step] + (["--audit"] if audit else []),
            [q4 / (stem + ext) for ext in (".json", ".npz")])
    geometry("q4-fv320", "q4", 320, 0)
    geometry("q4-fv320-tight", "q4", 320, 0, rtol="2e-10", step="120")
    geometry("q4-fv640-tight", "q4", 640, 0, rtol="2e-10", step="120")
    geometry("geometry-q1-1d", "q1", 80, 0, "harmonic")
    geometry("geometry-q1-2d", "q1", 80, 128, "harmonic")
    geometry("geometry-q23-1d", "q23", 80, 0)
    geometry("geometry-q4-1d", "q4", 80, 0)
    return tasks


def execute(name, args, products, root, out, source_signature):
    record = out / "logs" / (name + ".json")
    command = [sys.executable, "-X", "utf8", "-B", *args]
    identity = {"command": [a.replace(str(out), "<OUTPUT>") for a in command],
                "source_signature": source_signature}
    if record.exists():
        old = json.loads(record.read_text(encoding="utf-8"))
        if old.get("returncode") == 0 and all(old.get(k) == v for k, v in identity.items()):
            if all(p.is_file() and sha(p) == old["outputs"].get(p.relative_to(out).as_posix()) for p in products):
                print("VERIFIED RESUME " + name, flush=True)
                return old
        raise RuntimeError("Incomplete or changed run: " + name + "; use a new output directory")
    env = os.environ.copy()
    env.update(PYTHONUTF8="1", PYTHONDONTWRITEBYTECODE="1", OPENBLAS_NUM_THREADS="1",
               OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", MPLBACKEND="Agg")
    record.parent.mkdir(parents=True, exist_ok=True)
    start = time.time()
    print("START " + name, flush=True)
    with (record.with_suffix(".log")).open("wb") as log:
        process = subprocess.run(command, cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT)
    result = dict(identity, returncode=process.returncode, elapsed_s=time.time() - start,
                  outputs={p.relative_to(out).as_posix(): sha(p) for p in products if p.is_file()})
    save(record, result)
    if process.returncode or len(result["outputs"]) != len(products):
        raise RuntimeError("Run failed: " + name + "; see " + str(record.with_suffix(".log")))
    print("PASS " + name + " (" + str(round(result["elapsed_s"], 1)) + " s)", flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--stage", choices=["solve", "postprocess", "all"], default="all")
    args = parser.parse_args()
    out = args.out.resolve()
    if out == HERE or HERE in out.parents:
        raise ValueError("Output must be outside the submission package")
    root, manifest = prepare(out)
    signature = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()
    if args.stage in ("solve", "all"):
        tasks = commands(root)
        errors = []
        with ThreadPoolExecutor(max_workers=max(1, min(args.jobs, 4))) as pool:
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
    if args.stage in ("postprocess", "all"):
        from reproduction.postprocess import finish
        finish(HERE, out, root, signature)


if __name__ == "__main__":
    main()
