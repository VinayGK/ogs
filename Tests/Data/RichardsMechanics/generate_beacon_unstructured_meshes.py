#!/usr/bin/env python3
"""Generate the small unstructured BEACON meshes used by the report runs.

The script rescales the source GML footprint into the 1a01 and 1b benchmark
domains and writes both the VTU meshes and their matching geometry files.
"""

from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import vtk


BASE_DIR = Path(__file__).resolve().parent


def _read_rectangle(gml_path: Path) -> tuple[float, float]:
    """Extract rectangle width and height from a simple GML geometry."""
    root = ET.parse(gml_path).getroot()
    xs = [float(point.attrib["x"]) for point in root.findall(".//point")]
    ys = [float(point.attrib["y"]) for point in root.findall(".//point")]
    return max(xs) - min(xs), max(ys) - min(ys)


def _build_mesh(width: float, height: float, *, nx: int, ny: int, seed: int) -> vtk.vtkUnstructuredGrid:
    """Create a lightly perturbed triangulation of a rectangular domain."""
    rng = np.random.default_rng(seed)
    points = vtk.vtkPoints()

    for j in range(ny):
        y = height * j / (ny - 1)
        for i in range(nx):
            x = width * i / (nx - 1)
            if 0 < i < nx - 1 and 0 < j < ny - 1:
                x += rng.uniform(-0.18, 0.18) * width / (nx - 1)
                y += rng.uniform(-0.18, 0.18) * height / (ny - 1)
            points.InsertNextPoint(float(x), float(y), 0.0)

    poly = vtk.vtkPolyData()
    poly.SetPoints(points)

    delaunay = vtk.vtkDelaunay2D()
    delaunay.SetInputData(poly)
    delaunay.SetTolerance(1e-8)
    delaunay.Update()

    triangulation = delaunay.GetOutput()
    append = vtk.vtkAppendFilter()
    append.AddInputData(triangulation)
    append.Update()
    return append.GetOutput()


def _write_mesh(mesh: vtk.vtkUnstructuredGrid, path: Path) -> None:
    """Persist a VTU mesh in ASCII form for reproducible diffs."""
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(mesh)
    writer.SetDataModeToAscii()
    if writer.Write() != 1:
        raise RuntimeError(f"failed to write {path}")


def main() -> None:
    """Generate all benchmark meshes and report their sizes."""
    cases = [
        ("beacon_1a01.gml", "beacon_1a01_domain_unstructured_162e.vtu", 10, 10, 101),
        ("beacon_1b.gml", "beacon_1b_domain_unstructured_162e.vtu", 10, 10, 202),
    ]

    for gml_name, mesh_name, nx, ny, seed in cases:
        width, height = _read_rectangle(BASE_DIR / gml_name)
        mesh = _build_mesh(width, height, nx=nx, ny=ny, seed=seed)
        out = BASE_DIR / mesh_name
        _write_mesh(mesh, out)
        print(f"{mesh_name}: {mesh.GetNumberOfPoints()} points, {mesh.GetNumberOfCells()} cells")


if __name__ == "__main__":
    main()
