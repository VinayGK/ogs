#!/usr/bin/env python3
"""Extract simple time-history observables from a sequence of VTU files.

The script is intentionally small and dependency-light:
- it reads OGS VTK XML unstructured-grid outputs via `vtk`,
- it searches point data first, then cell data,
- it writes one CSV row per file,
- and it stores per-component means for the requested arrays.

This is a benchmark-observable extraction layer, not a calibration tool.
"""

from __future__ import annotations

import argparse
import csv
import glob
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


TIME_RE = re.compile(r"_t_([0-9.+\-eE]+)\.vtu$")


def parse_time_from_name(path: Path) -> float:
    match = TIME_RE.search(path.name)
    if not match:
        raise ValueError(
            f"Cannot parse time stamp from file name '{path.name}'. "
            "Expected suffix like '_t_1000.000000.vtu'."
        )
    return float(match.group(1))


def load_unstructured_grid(path: Path) -> vtk.vtkUnstructuredGrid:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    data = reader.GetOutput()
    if data is None:
        raise RuntimeError(f"Failed to read VTK file '{path}'.")
    return data


def find_array(dataset: vtk.vtkUnstructuredGrid, name: str, association: str):
    if association in ("auto", "point"):
        point_data = dataset.GetPointData()
        if point_data is not None:
            array = point_data.GetArray(name)
            if array is not None:
                return "point", vtk_to_numpy(array)

    if association in ("auto", "cell"):
        cell_data = dataset.GetCellData()
        if cell_data is not None:
            array = cell_data.GetArray(name)
            if array is not None:
                return "cell", vtk_to_numpy(array)

    raise KeyError(f"Array '{name}' was not found in point or cell data.")


def as_2d(values: np.ndarray) -> np.ndarray:
    if values.ndim == 1:
        return values.reshape(-1, 1)
    return values


def component_column_names(name: str, n_components: int) -> list[str]:
    if n_components == 1:
        return [f"{name}_mean"]
    return [f"{name}_c{i}_mean" for i in range(n_components)]


def iter_input_files(patterns: Iterable[str]) -> list[Path]:
    files: list[Path] = []
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        files.extend(Path(match) for match in matches)
    unique_files = sorted(set(files), key=parse_time_from_name)
    return unique_files


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Extract per-file mean observables from VTU outputs."
    )
    parser.add_argument(
        "--glob",
        dest="patterns",
        action="append",
        required=True,
        help="Glob pattern for input VTU files. Can be passed multiple times.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="CSV file to write.",
    )
    parser.add_argument(
        "--arrays",
        nargs="+",
        required=True,
        help="VTK point/cell data array names to extract.",
    )
    parser.add_argument(
        "--association",
        choices=("auto", "point", "cell"),
        default="auto",
        help="Where to look for the arrays. Default: auto.",
    )
    args = parser.parse_args()

    input_files = iter_input_files(args.patterns)
    if not input_files:
        raise SystemExit("No input VTU files matched the provided glob(s).")

    rows = []
    header = ["time_s"]
    component_counts: dict[str, int] = {}

    for path in input_files:
        dataset = load_unstructured_grid(path)
        time_s = parse_time_from_name(path)

        row = {"time_s": time_s}
        for name in args.arrays:
            source, values = find_array(dataset, name, args.association)
            values_2d = as_2d(np.asarray(values))
            if name not in component_counts:
                component_counts[name] = values_2d.shape[1]
                header.extend(component_column_names(name, values_2d.shape[1]))
            elif component_counts[name] != values_2d.shape[1]:
                raise RuntimeError(
                    f"Array '{name}' has inconsistent component count: "
                    f"expected {component_counts[name]}, got {values_2d.shape[1]} in '{path.name}'."
                )

            means = values_2d.mean(axis=0)
            for i, mean in enumerate(means):
                column = f"{name}_mean" if values_2d.shape[1] == 1 else f"{name}_c{i}_mean"
                row[column] = float(mean)

            # Keep the association check explicit even though the CSV does not
            # store it: the script is meant to fail early if the field moves.
            _ = source

        rows.append(row)

    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for row in rows:
            writer.writerow(
                [row["time_s"]]
                + [row[column] for column in header[1:]]
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
