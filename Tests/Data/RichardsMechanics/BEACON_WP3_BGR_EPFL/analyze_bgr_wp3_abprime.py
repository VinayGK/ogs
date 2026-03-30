#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


TARGETS = {
    "experiment_swelling_pressure_mpa_range": [3.12, 3.55],
    "bgr_boundary_value_swelling_pressure_mpa": 2.71,
    "initial_void_ratio_range": [0.83, 0.85],
    "initial_total_suction_mpa_range": [90.0, 110.0],
    "water_pressure_boundary_kpa": 20.0,
}


def load(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    points = np.array([grid.GetPoint(i) for i in range(grid.GetNumberOfPoints())])
    data = {
        grid.GetPointData().GetArrayName(i): vtk_to_numpy(grid.GetPointData().GetArray(i))
        for i in range(grid.GetPointData().GetNumberOfArrays())
    }
    return points, data


def boundary_mean(points: np.ndarray, values: np.ndarray, axis: int, side: str) -> np.ndarray:
    target = points[:, axis].max() if side == "max" else points[:, axis].min()
    mask = np.isclose(points[:, axis], target, atol=1e-10)
    return np.asarray(values[mask]).mean(axis=0)


def summarise_single(path: Path) -> dict:
    points, data = load(path)
    sigma_top = boundary_mean(points, data["sigma"], axis=1, side="max")
    sigma_side = boundary_mean(points, data["sigma"], axis=0, side="max")
    summary = {
        "path": str(path),
        "axial_swelling_pressure_mpa": abs(float(sigma_top[1])) / 1.0e6,
        "radial_stress_mpa": abs(float(sigma_side[0])) / 1.0e6,
        "mean_saturation": float(np.asarray(data["saturation"]).mean()),
    }
    for field in ("porosity", "transport_porosity", "dry_density_solid", "swelling_stress"):
        if field in data:
            summary[f"mean_{field}"] = float(np.asarray(data[field]).mean())
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--mfront", type=Path, required=True)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()

    result = {
        "targets": TARGETS,
        "native": summarise_single(args.native),
        "mfront": summarise_single(args.mfront),
    }

    text = json.dumps(result, indent=2, sort_keys=True)
    if args.json_out:
        args.json_out.write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
