#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


REPORT_TARGETS = {
    "1a01": {
        "axial_swelling_pressure_kpa": 604.0,
        "radial_swelling_pressure_kpa": 994.0,
        "stage1_dry_density": 1655.0,
        "post_mortem_density_profile": [1466.0, 1454.0, 1427.0, 1353.0],
    },
    "1b": {
        "dry_density_after_volume_adjustment": 1520.0,
        "swelling_pressure_note": "nonzero and stabilized after about 500 days",
    },
}


def load_point_data(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray], int]:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()

    points = np.array([grid.GetPoint(i) for i in range(grid.GetNumberOfPoints())])
    point_data = {
        grid.GetPointData().GetArrayName(i): vtk_to_numpy(grid.GetPointData().GetArray(i))
        for i in range(grid.GetPointData().GetNumberOfArrays())
    }
    return points, point_data, int(grid.GetNumberOfCells())


def boundary_mean(points: np.ndarray, values: np.ndarray, axis: int, atol: float = 1e-8) -> np.ndarray:
    target = points[:, axis].max()
    mask = np.isclose(points[:, axis], target, atol=atol)
    return np.asarray(values[mask]).mean(axis=0)


def density_profile(points: np.ndarray, density: np.ndarray, n_bins: int) -> list[dict[str, float]]:
    y = points[:, 1]
    ymin = float(y.min())
    ymax = float(y.max())
    edges = np.linspace(ymin, ymax, n_bins + 1)
    out = []
    for i in range(n_bins):
        lo = edges[i]
        hi = edges[i + 1]
        if i == n_bins - 1:
            mask = (y >= lo) & (y <= hi)
        else:
            mask = (y >= lo) & (y < hi)
        center_mm = 1000.0 * 0.5 * (lo + hi)
        out.append(
            {
                "center_mm": float(center_mm),
                "mean_kg_m3": float(np.asarray(density[mask]).mean()),
                "count": int(mask.sum()),
            }
        )
    return out


def max_abs_diff(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.abs(np.asarray(a) - np.asarray(b)).max())


def summarise(case: str, native_path: Path, bridge_path: Path) -> dict:
    native_points, native, native_n_cells = load_point_data(native_path)
    bridge_points, bridge, bridge_n_cells = load_point_data(bridge_path)

    if not np.allclose(native_points, bridge_points):
        raise RuntimeError("Native and bridge outputs are not on the same point set.")
    if native_n_cells != bridge_n_cells:
        raise RuntimeError("Native and bridge outputs are not on the same cell set.")

    bins = 4 if case == "1a01" else 10
    native_sigma_top = boundary_mean(native_points, native["sigma"], axis=1)
    native_sigma_side = boundary_mean(native_points, native["sigma"], axis=0)
    bridge_sigma_top = boundary_mean(bridge_points, bridge["sigma"], axis=1)
    bridge_sigma_side = boundary_mean(bridge_points, bridge["sigma"], axis=0)

    native_sw_top = boundary_mean(native_points, native["swelling_stress"], axis=1)
    native_sw_side = boundary_mean(native_points, native["swelling_stress"], axis=0)
    bridge_sw_top = boundary_mean(bridge_points, bridge["swelling_stress"], axis=1)
    bridge_sw_side = boundary_mean(bridge_points, bridge["swelling_stress"], axis=0)

    fields_to_compare = ["pressure", "saturation", "dry_density_solid", "sigma", "swelling_stress"]
    field_diffs = {field: max_abs_diff(native[field], bridge[field]) for field in fields_to_compare}

    return {
        "case": case,
        "report_targets": REPORT_TARGETS[case],
        "mesh": {
            "n_points": int(native_points.shape[0]),
            "n_cells": native_n_cells,
        },
        "native": {
            "axial_swelling_pressure_kpa": abs(float(native_sigma_top[1])) / 1000.0,
            "radial_swelling_pressure_kpa": abs(float(native_sigma_side[0])) / 1000.0,
            "top_swelling_stress_kpa": abs(float(native_sw_top[1])) / 1000.0,
            "side_swelling_stress_kpa": abs(float(native_sw_side[0])) / 1000.0,
            "dry_density_mean_kg_m3": float(np.asarray(native["dry_density_solid"]).mean()),
            "dry_density_profile_kg_m3": density_profile(native_points, native["dry_density_solid"], bins),
        },
        "bridge": {
            "axial_swelling_pressure_kpa": abs(float(bridge_sigma_top[1])) / 1000.0,
            "radial_swelling_pressure_kpa": abs(float(bridge_sigma_side[0])) / 1000.0,
            "top_swelling_stress_kpa": abs(float(bridge_sw_top[1])) / 1000.0,
            "side_swelling_stress_kpa": abs(float(bridge_sw_side[0])) / 1000.0,
            "dry_density_mean_kg_m3": float(np.asarray(bridge["dry_density_solid"]).mean()),
            "dry_density_profile_kg_m3": density_profile(bridge_points, bridge["dry_density_solid"], bins),
        },
        "native_vs_bridge_max_abs_diff": field_diffs,
    }

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=["1a01", "1b"], required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--bridge", type=Path, required=True)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()

    result = summarise(args.case, args.native, args.bridge)

    text = json.dumps(result, indent=2, sort_keys=True)
    if args.json_out:
        args.json_out.write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
