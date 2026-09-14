"""Build a compact submission archive only after actual reproduction passes."""
from pathlib import Path
import hashlib
import json
import shutil
import zipfile

BASE = Path(__file__).resolve().parent
PACKAGE = BASE / "package/支撑材料"
RUN = BASE / "runs"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def main():
    comparison = load(RUN / "comparison.json")
    workbooks = load(RUN / "workbook_verification.json")
    assert comparison["pass"] and workbooks["pass"]
    logs = {p.stem: load(p) for p in sorted((RUN / "logs").glob("*.json"))}
    assert len(logs) >= 28 and all(v["returncode"] == 0 for v in logs.values())
    for log in logs.values():
        for name, digest in log["outputs"].items():
            assert sha(RUN / name) == digest, name
    summary = {
        "status": "PASS", "scope": "Cold-start package-only numerical reproduction and portable workbook read-back",
        "paper_numeric_cells": comparison["paper_numeric_cells"],
        "paper_tables": {k: {"pass": v["pass"], "numeric_cells": v["numeric_cells"]} for k, v in comparison["tables"].items()},
        "claims": comparison["claims"], "full_outputs": comparison["full_outputs"],
        "geometry_table": comparison["geometry_table"], "environment_table": comparison["environment_table"],
        "geometry_endpoint_reproduction": comparison["geometry_endpoint_reproduction"],
        "exported_workbook_cells_read_back": workbooks["total_cells"],
        "completed_tasks": {name: {"returncode": v["returncode"], "elapsed_s": v["elapsed_s"]} for name, v in logs.items()},
        "comparison_sha256": sha(RUN / "comparison.json"),
        "workbook_verification_sha256": sha(RUN / "workbook_verification.json"),
        "solver_input_source_sha256": sha(RUN / "source_manifest.json"),
        "note": "Numerical reproduction does not establish a physical calibration or a continuous two-dimensional strict guarantee."
    }
    (PACKAGE / "复现验收结果.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    listing = sorted({p.relative_to(PACKAGE).as_posix() for p in PACKAGE.rglob("*")
                      if p.is_file() and "__pycache__" not in p.parts} | {"MANIFEST.json"})
    (PACKAGE / "支撑材料文件列表.md").write_text("# 支撑材料文件列表\n\n完整复现入口：`reproduce.py`。运行方法见 `复现说明.md`，实际验收见 `复现验收结果.json`。文件哈希见 `MANIFEST.json`（清单不对自身求哈希）。\n\n" + "\n".join("- `" + name + "`" for name in listing) + "\n", encoding="utf-8")
    sources = {p.relative_to(PACKAGE).as_posix(): sha(p) for p in sorted(PACKAGE.rglob("*"))
               if p.is_file() and "__pycache__" not in p.parts and p.name != "MANIFEST.json"}
    (PACKAGE / "MANIFEST.json").write_text(json.dumps(sources, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    target = BASE / "out/A题_支撑材料_复现验收版.zip"
    target.parent.mkdir(exist_ok=True)
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for name in sources:
            z.write(PACKAGE / name, "支撑材料/" + name)
        z.write(PACKAGE / "MANIFEST.json", "支撑材料/MANIFEST.json")
    with zipfile.ZipFile(target) as z:
        assert z.testzip() is None
        for name, digest in sources.items():
            assert hashlib.sha256(z.read("支撑材料/" + name)).hexdigest() == digest
    assert target.stat().st_size < 20_000_000, "Compact submission target exceeded"
    target.with_suffix(".zip.sha256").write_text(sha(target) + "  " + target.name + "\n", encoding="utf-8")
    report = {"zip": target.relative_to(BASE).as_posix(), "bytes": target.stat().st_size,
              "sha256": sha(target), "members": len(sources) + 1, "all_members_read_back": True}
    (BASE / "evidence/verified_package.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False))
    paper = BASE / "paper_source"
    # The paper source has the identical runnable support tree and extra TeX files.
    for source in PACKAGE.rglob("*"):
        if source.is_file() and "__pycache__" not in source.parts:
            dest = paper / source.relative_to(PACKAGE)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, dest)
    source_files = {p.relative_to(paper).as_posix(): sha(p) for p in sorted(paper.rglob("*"))
                    if p.is_file() and "build" not in p.relative_to(paper).parts and "__pycache__" not in p.parts}
    source_zip = BASE / "out/四问论文_复现验收版_源稿.zip"
    with zipfile.ZipFile(source_zip, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for name in source_files:
            z.write(paper / name, name)
        z.writestr("SOURCE_MANIFEST.json", json.dumps(source_files, ensure_ascii=False, indent=2))
    with zipfile.ZipFile(source_zip) as z:
        assert z.testzip() is None
        for name, digest in source_files.items():
            assert hashlib.sha256(z.read(name)).hexdigest() == digest
    source_zip.with_suffix(".zip.sha256").write_text(sha(source_zip) + "  " + source_zip.name + "\n", encoding="utf-8")
    (BASE / "evidence/paper_source_package.json").write_text(json.dumps({"bytes": source_zip.stat().st_size,
        "sha256": sha(source_zip), "members": len(source_files) + 1, "files": source_files}, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
