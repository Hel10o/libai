"""Read saved same-grid 1D/2D fields and independently integrate stored weights."""
from pathlib import Path
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BASE = ROOT / "q4_complete_delivery/v1/validation/geometry"
target = HERE / "saved_geometry_scope.json"
if target.exists():
    raise FileExistsError(target)
results = []
for two_name, one_name, when in (("q1_80x128", "q1_80x0", 1800),
    ("q23_80x256_integral", "q23_80x0_integral", 10800),
    ("q4_80x128_integral", "q4_80x0_integral", 21600)):
    items = []
    sources = {}
    for name in (two_name, one_name):
        path = BASE / (name + ".npz")
        sources[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        field = np.load(path)
        ix = int(np.flatnonzero(field["snapshot_times_s"] == when)[0])
        state = field["snapshots"][ix].reshape(-1, 2)
        weights = field["weights"]
        assert len(state) == len(weights) and np.all(weights > 0)
        # Common shrinkage Jacobian cancels from the normalized average.
        items.append((state * weights[:, None]).sum(axis=0) / weights.sum())
    two, one = items
    item = {"two_dimensional": two_name, "paired_one_dimensional": one_name,
        "saved_time_s": when, "mean_T_C_2d": two.tolist(), "mean_T_C_1d": one.tolist(),
        "mean_T_C_2d_minus_1d": (two-one).tolist(),
        "normalized_water_loss_2d_relative_increase": float((2.55-two[1])/(2.55-one[1])-1),
        "source_sha256": sources,
        "interpretation": "Postprocessing saved states, not new PDE runs. Volume means equal dry-mass means only for the current spatially uniform dry density assumption. Thermal mean is arithmetic volume mean, not enthalpy."}
    results.append(item)
    print(json.dumps(item, ensure_ascii=False), flush=True)
target.write_text(json.dumps(results, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
