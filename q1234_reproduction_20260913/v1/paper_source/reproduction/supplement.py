"""Additional paper scenarios using only the newly calculated workspace."""
import argparse
import sys
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--case", choices=["q2-sensitivity"], required=True)
    a = p.parse_args()
    path = a.root / "q2_refinement_delivery/runtime/source"
    sys.path.insert(0, str(path))
    from q2_core import Environment, load_environment
    from q2_tests import run_sensitivity
    data, _ = load_environment(a.root / "q2_refinement_delivery/runtime/inputs/附件1.xlsx")
    run_sensitivity(Environment(data), a.out)


if __name__ == "__main__":
    main()
