#!/usr/bin/env python3
"""Plot consolidated native-vs-MFront calibrated benchmark comparisons.

Inputs:
- calibrated_benchmark_pairwise_deltas.csv produced by
  run_calibrated_beacon_anchors_comparison.py

Outputs:
- calibrated_beacon_native_mfront_comparison.png
- calibrated_anchors_native_mfront_parity.png

These figures are meant for direct insertion into the transition note.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTPUT_ROOT = Path(__file__).resolve().parent / "_outputs" / "calibrated_native_mfront_comparison"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pairwise-csv",
        type=Path,
        default=OUTPUT_ROOT / "calibrated_benchmark_pairwise_deltas.csv",
    )
    parser.add_argument(
        "--beacon-png",
        type=Path,
        default=OUTPUT_ROOT / "calibrated_beacon_native_mfront_comparison.png",
    )
    parser.add_argument(
        "--anchors-png",
        type=Path,
        default=OUTPUT_ROOT / "calibrated_anchors_native_mfront_parity.png",
    )
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def to_float(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def plot_beacon(rows: list[dict[str, str]], out_png: Path) -> None:
    case_order = ["1a01", "1b"]
    by_case = {r["case_id"]: r for r in rows if r["benchmark_group"] == "BEACON_report"}
    cases = [c for c in case_order if c in by_case]
    if not cases:
        raise RuntimeError("No BEACON rows found in pairwise CSV.")

    native_ps = np.array([to_float(by_case[c]["native_final_swelling_pressure_mpa"]) for c in cases])
    mfront_ps = np.array([to_float(by_case[c]["mfront_final_swelling_pressure_mpa"]) for c in cases])
    native_mult = np.array([to_float(by_case[c]["native_multiplier"]) for c in cases])
    mfront_mult = np.array([to_float(by_case[c]["mfront_multiplier"]) for c in cases])
    delta_ax = np.array([to_float(by_case[c]["delta_axial_sigma_kpa_native_minus_mfront"]) for c in cases])
    delta_rad = np.array([to_float(by_case[c]["delta_radial_sigma_kpa_native_minus_mfront"]) for c in cases])

    x = np.arange(len(cases))
    w = 0.34

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0))

    ax = axes[0, 0]
    ax.bar(x - w / 2, native_ps, width=w, label="Native", color="#1f77b4")
    ax.bar(x + w / 2, mfront_ps, width=w, label="MFront", color="#d62728")
    ax.set_xticks(x, cases)
    ax.set_ylabel("Final swelling pressure (MPa)")
    ax.set_title("BEACON calibrated replay: final swelling pressure")
    ax.grid(axis="y", alpha=0.3)
    ax.legend(frameon=False)

    ax = axes[0, 1]
    ax.bar(x - w / 2, native_mult, width=w, label="Native", color="#1f77b4")
    ax.bar(x + w / 2, mfront_mult, width=w, label="MFront", color="#d62728")
    ax.set_yscale("log")
    ax.set_xticks(x, cases)
    ax.set_ylabel("Effective vdW multiplier (-)")
    ax.set_title("Calibrated multipliers (log scale)")
    ax.grid(axis="y", which="both", alpha=0.3)

    ax = axes[1, 0]
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.bar(x - w / 2, delta_ax, width=w, label=r"$\Delta \sigma_{ax}$", color="#9467bd")
    ax.bar(x + w / 2, delta_rad, width=w, label=r"$\Delta \sigma_{rad}$", color="#2ca02c")
    ax.set_xticks(x, cases)
    ax.set_ylabel("Native - MFront stress (kPa)")
    ax.set_title("Boundary stress deltas")
    ax.grid(axis="y", alpha=0.3)
    ax.legend(frameon=False)

    ax = axes[1, 1]
    ratio = np.array([to_float(by_case[c]["native_to_mfront_multiplier_ratio"]) for c in cases])
    ax.plot(cases, ratio, marker="o", linestyle="-", color="#ff7f0e")
    ax.axhline(1.0, color="black", linewidth=1.0, linestyle="--")
    ax.set_ylabel("Native / MFront multiplier ratio")
    ax.set_title("Multiplier ratio by BEACON case")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def plot_anchors(rows: list[dict[str, str]], out_png: Path) -> None:
    anchors = [r for r in rows if r["benchmark_group"] == "ANCHORS_MS33_ModelI"]
    if not anchors:
        raise RuntimeError("No ANCHORS rows found in pairwise CSV.")

    anchors = sorted(anchors, key=lambda r: to_float(r["dry_density_kg_m3"]))
    dd = np.array([to_float(r["dry_density_kg_m3"]) for r in anchors])
    ratio = np.array([to_float(r["native_to_mfront_multiplier_ratio"]) for r in anchors])
    delta_ps_kpa = np.array(
        [1e3 * to_float(r["delta_swelling_pressure_mpa_native_minus_mfront"]) for r in anchors]
    )

    fig, axes = plt.subplots(2, 1, figsize=(9.5, 7.0), sharex=True)

    ax = axes[0]
    ax.plot(dd, ratio, marker="o", markersize=3.8, linewidth=1.5, color="#1f77b4")
    ax.axhline(1.0, color="black", linewidth=1.0, linestyle="--")
    ax.set_ylabel("Native / MFront multiplier ratio")
    ax.set_title("ANCHORS dense calibration parity")
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.plot(dd, delta_ps_kpa, marker="s", markersize=3.6, linewidth=1.3, color="#d62728")
    ax.set_xlabel(r"Dry density $\rho_d$ (kg/m$^3$)")
    ax.set_ylabel(r"$\Delta p_{cal}$ (kPa), native - MFront")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    rows = read_rows(args.pairwise_csv)
    plot_beacon(rows, args.beacon_png)
    plot_anchors(rows, args.anchors_png)
    print(args.beacon_png)
    print(args.anchors_png)


if __name__ == "__main__":
    main()
