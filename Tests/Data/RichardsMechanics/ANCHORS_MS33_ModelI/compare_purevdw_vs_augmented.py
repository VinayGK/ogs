#!/usr/bin/env python3
"""Compare pure-vdW and augmented-vdW Villar dry-density calibrations.

Reads CSV files produced by:
  run_villar_dense_dd_native_purevdw_calibration.py
  run_villar_dense_dd_native_augmented_calibration.py

and produces multi-panel comparison figures.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from pathlib import Path

import matplotlib
import numpy as np

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mplconfig-"))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.mkdtemp(prefix="xdg-cache-"))
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parent

HAMAKER_LIT = 5.1e-21   # J — for unit labelling
SPECIFIC_SURFACE = 523.0  # m²/g (code convention; NOT converted to m²/kg)
RHO_SOLID = 2780.0        # kg/m³
# lambda unit note: lambda_code = lambda_physical × 1000 (because Sa uses m²/g).
# lambda_code=1e-6 → physical 1 nm decay. lambda_code=1e-9 → physical 1 pm (sub-atomic).


def villar_curve(dd_kg_m3: np.ndarray) -> np.ndarray:
    return np.exp(6.77 * dd_kg_m3 / 1000.0 - 9.07)


def read_csv(path: Path) -> dict[str, np.ndarray]:
    with path.open() as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    if not rows:
        raise ValueError(f"Empty CSV: {path}")
    out: dict[str, list] = {k: [] for k in rows[0]}
    for row in rows:
        for k, v in row.items():
            try:
                out[k].append(float(v))
            except (ValueError, TypeError):
                out[k].append(v)
    return {k: np.array(v) for k, v in out.items()}


def find_augmented_csvs() -> list[Path]:
    return sorted(ROOT.glob("villar_dense_dd_native_augmented_*_calibration.csv"))


def find_purevdw_csv() -> Path | None:
    p = ROOT / "villar_dense_dd_native_purevdw_calibration.csv"
    return p if p.exists() else None


def xi0_from_n_l0(n_l0: float, lam: float, n_s: float) -> float:
    return n_l0 / (lam * n_s * RHO_SOLID * SPECIFIC_SURFACE)


def plot_swelling_pressure(
    ax,
    dd_dense: np.ndarray,
    vdw_data: dict | None,
    aug_datasets: list[dict],
) -> None:
    dd_fine = np.linspace(dd_dense.min(), dd_dense.max(), 200)
    ax.plot(dd_fine, villar_curve(dd_fine), "k-", lw=2, label="Villar target")
    if vdw_data is not None:
        ax.plot(
            vdw_data["dry_density_kg_m3"],
            vdw_data["calibrated_MPa"],
            "o--",
            color="#1f77b4",
            ms=5,
            lw=1.6,
            label="Pure vdW (calibrated $m_{vdW}$)",
        )
    colors = ["#d62728", "#ff7f0e", "#2ca02c", "#9467bd"]
    for i, aug in enumerate(aug_datasets):
        lam = float(aug["lambda_m"][0])
        c = colors[i % len(colors)]
        ax.plot(
            aug["dry_density_kg_m3"],
            aug["calibrated_MPa"],
            "s--",
            color=c,
            ms=4,
            lw=1.4,
            label=f"Augmented (λ={lam:.0e} m)",
        )
    ax.set_xlabel("Dry density (kg/m³)")
    ax.set_ylabel("Swelling pressure (MPa)")
    ax.set_yscale("log")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(which="both", alpha=0.25)
    ax.set_title("Villar swelling pressure vs dry density")


def plot_multiplier(ax, vdw_data: dict | None) -> None:
    if vdw_data is None:
        ax.text(0.5, 0.5, "No pure-vdW data", ha="center", va="center",
                transform=ax.transAxes)
        return
    ax.semilogy(
        vdw_data["dry_density_kg_m3"],
        vdw_data["vdw_multiplier"],
        "^-",
        color="#1f77b4",
        ms=5,
        lw=1.8,
        label="$m_{vdW}$ (pure)",
    )
    ax.axhline(1.0, color="grey", lw=1, ls="--", label="$m_{vdW}=1$ (literature A)")
    ax.set_xlabel("Dry density (kg/m³)")
    ax.set_ylabel("$m_{vdW}$ = $A_{eff}/A_{lit}$")
    ax.set_title("Required vdW multiplier (pure model)")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(which="both", alpha=0.25)


def plot_K_curves(ax, aug_datasets: list[dict]) -> None:
    if not aug_datasets:
        ax.text(0.5, 0.5, "No augmented data", ha="center", va="center",
                transform=ax.transAxes)
        return
    colors = ["#d62728", "#ff7f0e", "#2ca02c", "#9467bd"]
    for i, aug in enumerate(aug_datasets):
        lam = float(aug["lambda_m"][0])
        c = colors[i % len(colors)]
        ax.semilogy(
            aug["dry_density_kg_m3"],
            aug["K_calibrated_J_kg"],
            "v-",
            color=c,
            ms=4,
            lw=1.6,
            label=f"λ={lam:.0e} m",
        )
    ax.set_xlabel("Dry density (kg/m³)")
    ax.set_ylabel("K  (J/kg)")
    ax.set_title("Augmented prefactor K vs dry density")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(which="both", alpha=0.25)


def plot_xi0(ax, aug_datasets: list[dict]) -> None:
    """Plot initial dimensionless film thickness xi0 = h0/lambda."""
    if not aug_datasets:
        ax.text(0.5, 0.5, "No augmented data", ha="center", va="center",
                transform=ax.transAxes)
        return
    colors = ["#d62728", "#ff7f0e", "#2ca02c", "#9467bd"]
    for i, aug in enumerate(aug_datasets):
        lam = float(aug["lambda_m"][0])
        c = colors[i % len(colors)]
        xi0 = aug.get("xi0", None)
        if xi0 is None:
            n_l0 = aug["n_l0"]
            n_s = aug["n_s_ref"]
            xi0 = n_l0 / (lam * n_s * RHO_SOLID * SPECIFIC_SURFACE)
        ax.plot(
            aug["dry_density_kg_m3"],
            xi0,
            "D-",
            color=c,
            ms=4,
            lw=1.4,
            label=f"λ={lam:.0e} m",
        )
    ax.axhline(1.0, color="k", lw=1, ls=":", label="ξ₀ = 1 (thin-film limit)")
    ax.set_xlabel("Dry density (kg/m³)")
    ax.set_ylabel("ξ₀ = $n_{l0}$ / (λ·$n_S$·$ρ_{SR}$·$S_a$)")
    ax.set_title("Initial film-thickness parameter ξ₀")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.25)


def plot_relative_errors(ax, vdw_data: dict | None, aug_datasets: list[dict]) -> None:
    if vdw_data is not None:
        rel_err = (
            100 * np.abs(vdw_data["delta_MPa"]) / np.maximum(vdw_data["target_villar_MPa"], 1e-12)
        )
        ax.plot(
            vdw_data["dry_density_kg_m3"],
            rel_err,
            "o--",
            color="#1f77b4",
            ms=5,
            lw=1.4,
            label="Pure vdW",
        )
    colors = ["#d62728", "#ff7f0e", "#2ca02c", "#9467bd"]
    for i, aug in enumerate(aug_datasets):
        lam = float(aug["lambda_m"][0])
        c = colors[i % len(colors)]
        rel_err = (
            100 * np.abs(aug["delta_MPa"]) / np.maximum(aug["target_villar_MPa"], 1e-12)
        )
        ax.plot(
            aug["dry_density_kg_m3"],
            rel_err,
            "s--",
            color=c,
            ms=4,
            lw=1.2,
            label=f"Augmented λ={lam:.0e}",
        )
    ax.axhline(2.0, color="grey", lw=1, ls="--", label="2% tol")
    ax.set_xlabel("Dry density (kg/m³)")
    ax.set_ylabel("Relative error (%)")
    ax.set_title("Calibration relative error")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.25)


def main() -> None:
    ap = argparse.ArgumentParser(description="Compare pure-vdW vs augmented calibrations")
    ap.add_argument("--out-dir", type=Path, default=ROOT)
    ap.add_argument("--purevdw-csv", type=Path, default=None)
    ap.add_argument("--augmented-csv", type=Path, nargs="*", default=None)
    args = ap.parse_args()

    purevdw_csv = args.purevdw_csv or find_purevdw_csv()
    aug_csvs = args.augmented_csv if args.augmented_csv is not None else find_augmented_csvs()

    if purevdw_csv is None and not aug_csvs:
        print("No calibration CSVs found — run the calibration scripts first.")
        return

    vdw_data = read_csv(purevdw_csv) if purevdw_csv and purevdw_csv.exists() else None
    aug_datasets = [read_csv(p) for p in aug_csvs if p.exists()]

    if vdw_data is None and not aug_datasets:
        print("No readable CSV data found.")
        return

    if vdw_data is not None:
        dd_arr = vdw_data["dry_density_kg_m3"]
    else:
        dd_arr = aug_datasets[0]["dry_density_kg_m3"]

    if vdw_data:
        print(f"Pure vdW: {len(vdw_data['dry_density_kg_m3'])} points")
        print(f"  m_vdW range: {vdw_data['vdw_multiplier'].min():.2e} – "
              f"{vdw_data['vdw_multiplier'].max():.2e}")
        mean_re = 100 * np.mean(np.abs(vdw_data["delta_MPa"]) /
                                np.maximum(vdw_data["target_villar_MPa"], 1e-12))
        print(f"  Mean relative error: {mean_re:.2f}%")

    for aug in aug_datasets:
        lam = float(aug["lambda_m"][0])
        mean_re = 100 * np.mean(np.abs(aug["delta_MPa"]) /
                                np.maximum(aug["target_villar_MPa"], 1e-12))
        print(f"Augmented λ={lam:.0e}: {len(aug['dry_density_kg_m3'])} points, "
              f"mean err={mean_re:.2f}%, "
              f"K={aug['K_calibrated_J_kg'].min():.2e}–{aug['K_calibrated_J_kg'].max():.2e} J/kg")

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Figure 1: 4-panel overview ---
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle("DSM_native: pure vdW vs augmented vdW — Villar calibration", fontsize=12)

    plot_swelling_pressure(axes[0, 0], dd_arr, vdw_data, aug_datasets)
    plot_multiplier(axes[0, 1], vdw_data)
    plot_K_curves(axes[1, 0], aug_datasets)
    plot_xi0(axes[1, 1], aug_datasets)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_dir / "compare_purevdw_vs_augmented_overview.png", dpi=220)
    plt.close(fig)
    print(f"\nWrote: {out_dir / 'compare_purevdw_vs_augmented_overview.png'}")

    # --- Figure 2: swelling pressure + relative errors, 2 panels ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("DSM_native calibration accuracy — pure vdW vs augmented", fontsize=11)

    plot_swelling_pressure(axes[0], dd_arr, vdw_data, aug_datasets)
    plot_relative_errors(axes[1], vdw_data, aug_datasets)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_dir / "compare_purevdw_vs_augmented_accuracy.png", dpi=220)
    plt.close(fig)
    print(f"Wrote: {out_dir / 'compare_purevdw_vs_augmented_accuracy.png'}")

    # --- Figure 3: direct K vs m_vdW comparison ---
    if vdw_data is not None and aug_datasets:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.set_title("Calibration parameter magnitude comparison")

        ax2 = ax.twinx()

        lines1 = ax.semilogy(
            vdw_data["dry_density_kg_m3"],
            vdw_data["vdw_multiplier"],
            "^-",
            color="#1f77b4",
            ms=5,
            lw=1.8,
            label="$m_{vdW}$ (pure, left axis)",
        )

        colors = ["#d62728", "#ff7f0e", "#2ca02c", "#9467bd"]
        lines2 = []
        for i, aug in enumerate(aug_datasets):
            lam = float(aug["lambda_m"][0])
            c = colors[i % len(colors)]
            line, = ax2.semilogy(
                aug["dry_density_kg_m3"],
                aug["K_calibrated_J_kg"],
                "v--",
                color=c,
                ms=4,
                lw=1.4,
                label=f"K (augmented λ={lam:.0e}, right axis)",
            )
            lines2.append(line)

        ax.set_xlabel("Dry density (kg/m³)")
        ax.set_ylabel("$m_{vdW}$ (dimensionless)", color="#1f77b4")
        ax2.set_ylabel("K (J/kg)")
        ax.tick_params(axis="y", colors="#1f77b4")

        all_lines = lines1 + lines2
        all_labels = [l.get_label() for l in all_lines]
        ax.legend(all_lines, all_labels, frameon=False, fontsize=8)
        ax.grid(which="both", alpha=0.25)
        fig.tight_layout()
        fig.savefig(out_dir / "compare_purevdw_vs_augmented_param.png", dpi=220)
        plt.close(fig)
        print(f"Wrote: {out_dir / 'compare_purevdw_vs_augmented_param.png'}")


if __name__ == "__main__":
    main()
