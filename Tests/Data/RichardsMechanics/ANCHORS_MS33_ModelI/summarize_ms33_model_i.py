#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
from pathlib import Path

import vtk
from vtk.util.numpy_support import vtk_to_numpy


ROOT = Path(__file__).resolve().parent
RHO_SOLID = 2780.0
CASES = {
    "dd1400": {
        "rho_dry": 1400.0,
        "specific_weight_kN_m3": 14.0,
        "phi0": 0.49640287769784175,
        "k0": 1.2264011098493141e-20,
    },
    "dd1600": {
        "rho_dry": 1600.0,
        "specific_weight_kN_m3": 16.0,
        "phi0": 0.4244604316546763,
        "k0": 5.870260900123441e-21,
    },
    "dd1800": {
        "rho_dry": 1800.0,
        "specific_weight_kN_m3": 18.0,
        "phi0": 0.3525179856115108,
        "k0": 2.6569572325678573e-21,
    },
}


def read_grid(path: Path) -> vtk.vtkUnstructuredGrid:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    return reader.GetOutput()


def get_array(grid: vtk.vtkUnstructuredGrid, name: str):
    arr = grid.GetPointData().GetArray(name)
    if arr is None:
        arr = grid.GetCellData().GetArray(name)
    if arr is None:
        raise KeyError(name)
    return vtk_to_numpy(arr)


def clean_scalar(value: float, tol: float = 1e-12) -> float:
    return 0.0 if abs(value) < tol else value


def reduce_snapshot(path: Path):
    grid = read_grid(path)
    pressure = get_array(grid, "pressure")
    saturation = get_array(grid, "saturation")
    porosity = get_array(grid, "porosity")
    dry_density = get_array(grid, "dry_density_solid")
    sigma = get_array(grid, "sigma")
    swelling = get_array(grid, "swelling_stress")

    mean_total_stress = float((-sigma[:, 0] - sigma[:, 1] - sigma[:, 2]).mean() / 3.0)
    mean_swelling_stress = float((-swelling[:, 0] - swelling[:, 1] - swelling[:, 2]).mean() / 3.0)

    return {
        "pressure_Pa_avg": clean_scalar(float(pressure.mean())),
        "pressure_Pa_min": clean_scalar(float(pressure.min())),
        "pressure_Pa_max": clean_scalar(float(pressure.max())),
        "suction_MPa_avg": clean_scalar(float(-pressure.mean() / 1e6)),
        "saturation_avg": clean_scalar(float(saturation.mean())),
        "porosity_avg": clean_scalar(float(porosity.mean())),
        "porosity_spread": clean_scalar(float(porosity.max() - porosity.min())),
        "dry_density_avg_kg_m3": clean_scalar(float(dry_density.mean())),
        "dry_density_spread_kg_m3": clean_scalar(float(dry_density.max() - dry_density.min())),
        "mean_total_stress_MPa": clean_scalar(mean_total_stress / 1e6),
        "mean_swelling_stress_MPa": clean_scalar(mean_swelling_stress / 1e6),
    }


def main() -> None:
    history_rows = []
    summary = {"rho_solid_kg_m3": RHO_SOLID, "cases": {}}

    for case, meta in CASES.items():
        snapshots = sorted(
            ROOT.glob(f"ms33_model_i_{case}_ts_*_t_*.vtu"),
            key=lambda path: int(path.stem.split("_ts_")[1].split("_t_")[0]),
        )
        if not snapshots:
            raise FileNotFoundError(case)

        case_history = []
        for path in snapshots:
            stem = path.stem
            timestep = int(stem.split("_ts_")[1].split("_t_")[0])
            time_s = float(stem.split("_t_")[1])
            snap = reduce_snapshot(path)
            snap["timestep"] = timestep
            snap["time_s"] = time_s
            snap["time_days"] = time_s / 86400.0
            case_history.append(snap)
            history_rows.append(
                {
                    "case": case,
                    "specific_weight_kN_m3": meta["specific_weight_kN_m3"],
                    "rho_dry_kg_m3": meta["rho_dry"],
                    "time_days": snap["time_days"],
                    "suction_MPa": snap["suction_MPa_avg"],
                    "mean_total_stress_MPa": snap["mean_total_stress_MPa"],
                    "mean_swelling_stress_MPa": snap["mean_swelling_stress_MPa"],
                    "saturation": snap["saturation_avg"],
                    "porosity": snap["porosity_avg"],
                }
            )

        summary["cases"][case] = {
            **meta,
            "history": case_history,
            "final": case_history[-1],
            "initial": case_history[0],
        }

    (ROOT / "ms33_model_i_summary.json").write_text(json.dumps(summary, indent=2))
    with (ROOT / "ms33_model_i_history.csv").open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "case",
                "specific_weight_kN_m3",
                "rho_dry_kg_m3",
                "time_days",
                "suction_MPa",
                "mean_total_stress_MPa",
                "mean_swelling_stress_MPa",
                "saturation",
                "porosity",
            ],
        )
        writer.writeheader()
        writer.writerows(history_rows)


if __name__ == "__main__":
    main()
