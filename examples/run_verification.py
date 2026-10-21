"""Produce figures, animation, and summary from the C++ demo output.

Runs the compiled demo (build/ignition_demo.exe), reads
build/ignition_history.csv, and writes docs/results/:
- ignition_runaway.png/.svg (cover): fuel/temperature/rate histories with
  the adiabatic limit and energy audit panel
- ignition_evolution.gif: 41 frames of the computed ignition runaway
- energy_audit.csv, run_summary.json

The chemistry itself is computed by the native C++ sources; this script
only builds, runs, and plots them.
"""

from __future__ import annotations

import json
import pathlib
import platform
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.animation as mpl_animation
import matplotlib.pyplot as plt
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
RESULTS = ROOT / "docs" / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

BLUE, RED, GREEN, GRAY = "#1f4e79", "#b23a48", "#0b6e4f", "#555555"
plt.rcParams.update(
    {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 10,
    }
)

CMAKE = None  # resolved portably at first use via shutil.which("cmake")
CTEST = None


def _find_cmake():
    """Locate cmake/ctest portably: PATH first, known fallbacks second."""
    global CMAKE, CTEST
    import shutil

    found = shutil.which("cmake")
    if not found:
        fallback = (r"C:\Program Files (x86)\Microsoft Visual Studio\2022"
                    r"\BuildTools\Common7\IDE\CommonExtensions\Microsoft"
                    r"\CMake\CMake\bin\cmake.exe")
        if pathlib.Path(fallback).exists():
            found = fallback
    if not found:
        raise RuntimeError("cmake not found on PATH or at the documented "
                           "fallback; install CMake to regenerate results")
    CMAKE = found
    ctest = shutil.which("ctest")
    if not ctest:
        ctest = found.replace("cmake.exe", "ctest.exe")
    CTEST = ctest


BUILD_DIR = "build-cmake"


def build_and_run():
    if CMAKE is None:
        _find_cmake()
    build = ROOT / BUILD_DIR
    subprocess.run([CMAKE, "-S", str(ROOT), "-B", str(build)], check=True,
                   capture_output=True)
    subprocess.run([CMAKE, "--build", str(build), "--config", "Release",
                    "--parallel"], check=True, capture_output=True)
    subprocess.run([CTEST, "--test-dir", str(build), "-C", "Release",
                    "--output-on-failure"], check=True)
    demo = build / "Release" / "ignition_demo.exe"
    if not demo.exists():
        demo = build / "ignition_demo"
    hist_dir = demo.parent
    # The demo writes build/ignition_history.csv relative to its cwd.
    (hist_dir / "build").mkdir(exist_ok=True)
    subprocess.run([str(demo)], check=True, cwd=str(hist_dir))
    out = hist_dir / "build" / "ignition_history.csv"
    # Commit the raw C++ history alongside the derived results.
    raw = np.genfromtxt(out, delimiter=",", names=True)
    (RESULTS / "ignition_history.csv").write_text(
        out.read_text(encoding="utf-8"), encoding="utf-8")
    return raw


def make_cover(data):
    t = data["time_s"]
    fig, ax = plt.subplots(1, 3, figsize=(13.0, 4.0))
    ax[0].plot(t * 1e3, data["Y_fuel"], color=BLUE, lw=1.8)
    ax[0].set_xlabel("t [ms]")
    ax[0].set_ylabel("fuel mass fraction Y [-]")
    ax[0].set_title("Fuel consumption (thermal runaway)")

    ax[1].plot(t * 1e3, data["T_K"], color=RED, lw=1.8)
    ax[1].axhline(1233.333, color=GRAY, ls="--", lw=1.2,
                  label="adiabatic limit T₀ + qY₀/cₚ")
    ax[1].set_xlabel("t [ms]")
    ax[1].set_ylabel("T [K]")
    ax[1].set_title("Temperature rise")
    ax[1].legend(frameon=False, fontsize=8)

    ax[2].semilogy(t * 1e3, data["k_1_s"], color=GREEN, lw=1.8)
    ax[2].set_xlabel("t [ms]")
    ax[2].set_ylabel("Arrhenius rate k(T) [1/s]")
    ax[2].set_title("Rate feedback (Arrhenius)")
    fig.suptitle("Ignition Zero — computed adiabatic runaway of the reduced "
                 "reacting system", fontsize=13)
    fig.tight_layout()
    fig.savefig(RESULTS / "ignition_runaway.png", dpi=150, bbox_inches="tight")
    fig.savefig(RESULTS / "ignition_runaway.svg", dpi=180, bbox_inches="tight")
    plt.close(fig)


def make_energy_audit(data):
    t = data["time_s"]
    released = 2.0e6 * (0.08 - data["Y_fuel"])
    sensible = 1200.0 * (data["T_K"] - 1100.0)
    residual = 2.0e6 * data["Y_fuel"]
    rows = ["time_s,released_J_per_kg,sensible_J_per_kg,residual_chemical_J_per_kg"]
    for ti, r, s, c in zip(t, released, sensible, residual):
        rows.append(f"{ti:.9e},{r:.9e},{s:.9e},{c:.9e}")
    (RESULTS / "energy_audit.csv").write_text("\n".join(rows), encoding="utf-8")

    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    ax.plot(t * 1e3, released, color=RED, lw=1.8, label="released q(Y₀−Y)")
    ax.plot(t * 1e3, sensible, "--", color=BLUE, lw=1.6, label="sensible cₚ(T−T₀)")
    ax.plot(t * 1e3, residual, ":", color=GREEN, lw=1.8, label="stored qY")
    ax.set_xlabel("t [ms]")
    ax.set_ylabel("energy [J/kg]")
    ax.set_title("Energy audit: released = sensible (J/kg balance)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS / "energy_audit.png", dpi=150, bbox_inches="tight")
    fig.savefig(RESULTS / "energy_audit.svg", dpi=180, bbox_inches="tight")
    plt.close(fig)
    return float(np.max(np.abs(released - sensible)))


def make_gif(data, frames=41):
    t = data["time_s"]
    idx = np.linspace(0, len(t) - 1, frames).astype(int)
    fig, ax = plt.subplots(1, 2, figsize=(10.6, 4.2), dpi=140)
    (line_y,) = ax[0].plot([], [], color=BLUE, lw=1.8)
    ax[0].set_xlim(0.0, t[-1] * 1e3)
    ax[0].set_ylim(0.0, 0.09)
    ax[0].set_xlabel("t [ms]")
    ax[0].set_ylabel("Y [-]")
    title_y = ax[0].set_title("Fuel at t = 0.00 ms")
    (line_t,) = ax[1].plot([], [], color=RED, lw=1.8)
    ax[1].set_xlim(0.0, t[-1] * 1e3)
    ax[1].set_ylim(1050.0, 1300.0)
    ax[1].set_xlabel("t [ms]")
    ax[1].set_ylabel("T [K]")
    title_t = ax[1].set_title("Temperature at t = 0.00 ms")

    writer = mpl_animation.PillowWriter(fps=8)
    with writer.saving(fig, RESULTS / "ignition_evolution.gif", dpi=140):
        for j, i in enumerate(idx):
            line_y.set_data(t[: i + 1] * 1e3, data["Y_fuel"][: i + 1])
            line_t.set_data(t[: i + 1] * 1e3, data["T_K"][: i + 1])
            title_y.set_text(f"Fuel at t = {t[i]*1e3:.2f} ms")
            title_t.set_text(f"Temperature at t = {t[i]*1e3:.2f} ms")
            writer.grab_frame()
    plt.close(fig)
    return frames


def main():
    data = build_and_run()
    make_cover(data)
    imbalance = make_energy_audit(data)
    frames = make_gif(data)
    summary = {
        "runtime": {
            "python": platform.python_version(),
            "build": "CMake + CTest (Release), MSVC 19.44 toolchain",
            "system": platform.platform(),
        },
        "case": {
            "A_1_s": 1.0e6, "Ea_J_per_mol": 8.0e4, "q_J_per_kg": 2.0e6,
            "cp_J_per_kgK": 1200.0, "T0_K": 1100.0, "Y0": 0.08,
            "dt_s": 1.0e-5, "steps": 12000,
        },
        "results": {
            "T_final_K": float(data["T_K"][-1]),
            "adiabatic_limit_K": 1100.0 + 2.0e6 * 0.08 / 1200.0,
            "Y_final": float(data["Y_fuel"][-1]),
            "energy_audit_max_imbalance_J_per_kg": imbalance,
        },
        "animation": {"file": "ignition_evolution.gif", "frames": frames},
    }
    (RESULTS / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
