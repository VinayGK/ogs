#!/usr/bin/env python3
"""Plot the MS33 dry-density sweep against the CIEMAT/Villar reference fit.

The script reads the reduced Model-I history, overlays the Villar Eq. (7)
curve, and writes both a CSV export and the comparison plot.
"""

from __future__ import annotations

import csv
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mplconfig-"))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.mkdtemp(prefix="xdg-cache-"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
SUMMARY_JSON = ROOT / "ms33_model_i_summary.json"
OUT_CSV = ROOT / "ms33_vs_ciemat_swelling_pressure_vs_dry_density.csv"
OUT_PNG = ROOT / "ms33_vs_ciemat_swelling_pressure_vs_dry_density.png"


def ciemat_villar_eq7_ps_mpa(qd_g_cm3: np.ndarray) -> np.ndarray:
    """Compute the Villar Eq. (7) reference curve in MPa."""
    return np.exp(6.77 * qd_g_cm3 - 9.07)


def load_ms33_points() -> list[tuple[float, float]]:
    """Load dry-density and final-stress pairs from the reduced Model-I summary."""
    payload = json.loads(SUMMARY_JSON.read_text())
    rows = []
    for case_data in payload["cases"].values():
        rho_d_kg_m3 = float(case_data["rho_dry"])
        qd_g_cm3 = rho_d_kg_m3 / 1000.0
        # In the reduced Model-I scaffold this is the isotropic confined response proxy.
        ps_mpa = float(case_data["final"]["mean_total_stress_MPa"])
        rows.append((qd_g_cm3, ps_mpa))
    return sorted(rows, key=lambda item: item[0])


def write_overlay_csv(ms33_points: list[tuple[float, float]]) -> None:
    """Export the reference curve and simulation points as a flat CSV table."""
    qd_ref = np.array([1.40, 1.60, 1.80], dtype=float)
    ps_ref = ciemat_villar_eq7_ps_mpa(qd_ref)

    qd_curve = np.linspace(1.35, 1.85, 101)
    ps_curve = ciemat_villar_eq7_ps_mpa(qd_curve)

    with OUT_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "series",
                "kind",
                "dry_density_g_cm3",
                "dry_density_kg_m3",
                "swelling_pressure_MPa",
                "basis",
            ]
        )
        for qd, ps in zip(qd_curve, ps_curve, strict=True):
            writer.writerow(
                [
                    "CIEMAT/Villar Eq(7)",
                    "curve",
                    f"{qd:.6f}",
                    f"{qd * 1000.0:.3f}",
                    f"{ps:.9f}",
                    "Lloret et al. 2007 Eq.(7): Ps=exp(6.77*qd-9.07)",
                ]
            )
        for qd, ps in zip(qd_ref, ps_ref, strict=True):
            writer.writerow(
                [
                    "CIEMAT/Villar Eq(7)",
                    "reference_point",
                    f"{qd:.6f}",
                    f"{qd * 1000.0:.3f}",
                    f"{ps:.9f}",
                    "reference densities 1.40/1.60/1.80 g/cm^3",
                ]
            )
        for qd, ps in ms33_points:
            writer.writerow(
                [
                    "OGS MS33 reduced Model-I",
                    "simulation_point",
                    f"{qd:.6f}",
                    f"{qd * 1000.0:.3f}",
                    f"{ps:.9f}",
                    "final mean total stress at 120 d (constant-volume scaffold)",
                ]
            )


def plot_overlay(ms33_points: list[tuple[float, float]]) -> None:
    """Render the dry-density comparison plot used in the note."""
    qd_curve = np.linspace(1.35, 1.85, 101)
    ps_curve = ciemat_villar_eq7_ps_mpa(qd_curve)
    qd_ref = np.array([1.40, 1.60, 1.80], dtype=float)
    ps_ref = ciemat_villar_eq7_ps_mpa(qd_ref)

    x_ms33 = np.array([row[0] for row in ms33_points], dtype=float)
    y_ms33 = np.array([row[1] for row in ms33_points], dtype=float)

    fig, ax = plt.subplots(figsize=(7.6, 4.8))

    ax.plot(
        qd_curve,
        ps_curve,
        color="black",
        linewidth=2.0,
        label="CIEMAT/Villar Eq. (7) fit",
    )
    ax.scatter(
        qd_ref,
        ps_ref,
        marker="o",
        facecolors="white",
        edgecolors="black",
        s=44,
        zorder=3,
        label="CIEMAT reference dry densities",
    )
    ax.plot(
        x_ms33,
        y_ms33,
        color="#d62728",
        linewidth=1.8,
        marker="s",
        markersize=6,
        label="OGS reduced MS33 sweep (120 d)",
    )

    for x, y in zip(x_ms33, y_ms33, strict=True):
        ax.annotate(
            f"{int(round(1000.0 * x))}",
            (x, y),
            textcoords="offset points",
            xytext=(0, 7),
            ha="center",
            fontsize=8.5,
            color="#d62728",
        )

    ax.set_title("Swelling Pressure vs Dry Density: CIEMAT/Villar vs OGS MS33")
    ax.set_xlabel("Dry density, $\\rho_d$ (g/cm$^3$)")
    ax.set_ylabel("Swelling pressure, $P_s$ (MPa)")
    ax.set_xlim(1.35, 1.85)
    ax.set_ylim(0.0, max(float(ps_ref.max()), float(y_ms33.max())) * 1.08)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", frameon=True)

    # Secondary x-axis in kg/m^3 for direct use with project files.
    secax = ax.secondary_xaxis("top", functions=(lambda x: 1000.0 * x, lambda x: x / 1000.0))
    secax.set_xlabel("Dry density, $\\rho_d$ (kg/m$^3$)")

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=220)
    plt.close(fig)


def main() -> None:
    """Load reduced results, refresh the overlay CSV, and write the figure."""
    if not SUMMARY_JSON.exists():
        raise FileNotFoundError(
            f"Missing {SUMMARY_JSON.name}. Run summarize_ms33_model_i.py first."
        )

    ms33_points = load_ms33_points()
    write_overlay_csv(ms33_points)
    plot_overlay(ms33_points)

    print(OUT_CSV)
    print(OUT_PNG)


if __name__ == "__main__":
    main()
