"""Regenerate seven numerical figures from the new calculations.

The original publication PDFs remain in figures/. These review figures use
portable fonts and expose the same physical quantities, not pixel identity.
Five schematic/workflow diagrams have no numerical solver output.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from reproduction.verify import at, data, load
from reproduce import sha, save


def build(root, out):
    dest = out / "figures"
    dest.mkdir(parents=True, exist_ok=True)
    q1, q2, q3, q4 = data(root)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "figure.constrained_layout.use": True, "pdf.fonttype": 42})
    outputs = []
    def finish(fig, name, sources):
        for ext in ("png", "pdf"):
            path = dest / (name + "." + ext)
            fig.savefig(path, dpi=170)
            outputs.append({"path": path.relative_to(out).as_posix(), "sha256": sha(path), "sources": sources})
        plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for time in [100, 300, 600, 900, 1200, 1500, 1800]:
        idx = at(q1["time_s"], time)
        for ax, key in zip(axes, ["temperature_degC", "moisture_dry_basis"]):
            ax.plot(q1["radius_cm"], q1[key][idx], label=f"{time} s")
    for ax, ylabel in zip(axes, ["Temperature / deg C", "Dry-basis moisture / (kg/kg)"]):
        ax.set(xlabel="Radius / cm", ylabel=ylabel)
        ax.grid(alpha=.2)
    axes[0].legend(ncols=2, fontsize=8)
    finish(fig, "q1_profiles", ["q1_unrounded.npz"])
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    sel = (q2["time_s"] <= 10800) & (q2["time_s"] % 30 == 0)
    for ax, channel, title in zip(axes, [0, 1], ["Temperature / deg C", "Dry-basis moisture / (kg/kg)"]):
        mesh = ax.pcolormesh(q2["radius_cm"], q2["time_s"][sel] / 3600,
                             q2["profile_TC"][sel, :, channel], shading="auto", cmap="viridis")
        ax.set(xlabel="Radius / cm", ylabel="Time / h", title=title)
        fig.colorbar(mesh, ax=ax)
    finish(fig, "q2_fields", ["q23_unified.npz, first 3h"])
    sel = q2["time_s"] <= 10800
    fields = q2["profile_TC"][sel]
    times = q2["time_s"][sel] / 3600
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    mechanism = {}
    for ax, i, location in zip(axes, [0, -1], ["Center", "Surface"]):
        T, C = fields[:, i, 0], fields[:, i, 1]
        H = 3850 * (1 / 301.15 - 1 / (T + 273.15))
        M = .45 * (1 / 2.55 - 1 / C)
        ax.plot(times, H, label="Temperature term H")
        ax.plot(times, M, label="Moisture term M")
        ax.plot(times, H + M, label="ln(D/D0)", lw=2)
        ax.set(xlabel="Time / h", ylabel="Log ratio", title=location)
        ax.grid(alpha=.2)
        ratio = np.exp(H + M)
        mechanism[location] = {"ratio_end": float(ratio[-1]), "peak_time_h": float(times[np.argmax(ratio)]),
                               "peak_ratio": float(ratio.max())}
    axes[0].legend(fontsize=8)
    finish(fig, "q2_mechanism", ["q23_unified.npz; analytical log decomposition"])
    save(out / "mechanism.json", mechanism)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax in axes:
        ax.plot(q2["time_s"] / 3600, q2["profile_TC"][:, 0, 1], label="Q3 center")
        ax.plot(q4["time_s"] / 3600, q4["max_C"], label="Q4 radial maximum")
        ax.axhline(.15, ls="--", color="black", lw=.8)
        ax.set(xlabel="Time / h", ylabel="Dry-basis moisture / (kg/kg)")
        ax.grid(alpha=.2)
    axes[0].legend()
    axes[1].set(xlim=(45, 59), ylim=(.145, .19), title="Endpoint region")
    finish(fig, "q3_q4_drying", ["q23_unified.npz", "main.npz"])
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot(q4["time_s"] / 3600, q4["radius_m"] * 100)
    axes[0].set(xlabel="Time / h", ylabel="Observed radius / cm")
    fixed = load(root / "q4_complete_delivery/v1/validation/fixed80.json")["event"]["critical_h"]
    shrinking = load(root / "q4_complete_delivery/v1/output/end_event.json")["critical_h"]
    axes[1].bar(["Fixed radius", "Observed shrinkage"], [fixed, shrinking], color=["#D55E00", "#0072B2"])
    axes[1].set(ylabel="Critical time / h")
    for i, value in enumerate([fixed, shrinking]):
        axes[1].text(i, value + 2, f"{value:.4f}", ha="center")
    finish(fig, "q4_shrinkage", ["main.npz", "fixed80.json", "end_event.json"])
    fig, axes = plt.subplots(2, 3, figsize=(12, 6))
    geometry = root / "q4_complete_delivery/v1/validation/geometry"
    for col, (case, pair, time) in enumerate([
        ("Q1, 1800 s", ["q1_80x128", "q1_80x0"], 1800),
        ("Q2, 3 h", ["q23_80x256_integral", "q23_80x0_integral"], 10800),
        ("Q4, 6 h", ["q4_80x128_integral", "q4_80x0_integral"], 21600)]):
        values = [np.load(geometry / (p + ".npz")) for p in pair]
        fields = [z["snapshots"][at(z["snapshot_times_s"], time), ..., 1] for z in values]
        radius = values[0]["reference_r_m"] * values[0]["radius_m"][at(values[0]["time_s"], time)] / .02 * 100
        axes[0, col].plot(radius, fields[0][0] - fields[1][0])
        axes[0, col].set(title=case, xlabel="Radius / cm", ylabel="C(2D)-C(1D) / (kg/kg)")
        means = [float(f.ravel() @ z["weights"] / z["weights"].sum()) for f, z in zip(fields, values)]
        axes[1, col].bar(["2D", "1D"], means, color=["#0072B2", "#D55E00"])
        axes[1, col].set(ylabel="Whole-cylinder mean / (kg/kg)")
        for i, value in enumerate(means):
            axes[1, col].text(i, value, f"{value:.6f}", ha="center", va="bottom", fontsize=8)
    finish(fig, "dimension_reduction", ["six newly generated matched 1D/2D fields"])
    fig, ax = plt.subplots(figsize=(9, 4.5))
    envdir = root / "q1234_overall_review_delivery/v1/evidence/environment"
    scenarios = ["temperature_minus_1K", "temperature_plus_1K", "boundary_minus_0p005", "boundary_plus_0p005"]
    for case, offset, marker in [("q3", -.12, "o"), ("q4", .12, "s")]:
        base = load(envdir / f"{case}_baseline_n40.json")["critical_s"]
        deltas = [(load(envdir / f"{case}_{s}_n40.json")["critical_s"] - base) / 3600 for s in scenarios]
        ax.scatter(deltas, np.arange(4) + offset, marker=marker, label=case.upper())
        for i, value in enumerate(deltas):
            ax.annotate(f"{value:+.6f}", (value, i + offset), xytext=(6, 0), textcoords="offset points", fontsize=8)
    ax.axvline(0, color="gray", lw=.8)
    ax.set_yticks(range(4), ["Temperature -1 K", "Temperature +1 K", "Boundary b -0.005", "Boundary b +0.005"])
    ax.set(xlabel="Change of critical time / h", xlim=(-2, 2.45))
    ax.legend()
    finish(fig, "environment_sensitivity", ["ten newly generated n40 environmental scenarios"])
    save(out / "figure_manifest.json", {"numerical_figures": 7, "layout": "Portable numerical review plots; not pixel-identical publication layout", "files": outputs})
