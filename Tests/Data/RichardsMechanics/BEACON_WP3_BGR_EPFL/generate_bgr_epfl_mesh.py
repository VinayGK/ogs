#!/usr/bin/env python3

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE_MESH = HERE.parent / "beacon_1a01_domain_unstructured_162e.vtu"
TARGET_MESH = HERE / "bgr_wp3_epfl_domain_162e.vtu"
TARGET_MESH_STRUCTURED = HERE / "bgr_wp3_epfl_domain_2e.vtu"
TARGET_GML = HERE / "bgr_wp3_epfl_geometry.gml"

SOURCE_WIDTH = 0.025
SOURCE_HEIGHT = 0.02
TARGET_WIDTH = 0.0175
TARGET_HEIGHT = 0.0125


def scale_mesh() -> None:
    tree = ET.parse(SOURCE_MESH)
    root = tree.getroot()
    data_array = root.find(".//Points/DataArray")
    if data_array is None or data_array.text is None:
        raise RuntimeError("Could not find point coordinates in source VTU.")

    values = list(map(float, data_array.text.split()))
    sx = TARGET_WIDTH / SOURCE_WIDTH
    sy = TARGET_HEIGHT / SOURCE_HEIGHT

    for i in range(0, len(values), 3):
        values[i] *= sx
        values[i + 1] *= sy

    data_array.text = " ".join(f"{v:.16g}" for v in values)
    tree.write(TARGET_MESH, encoding="ISO-8859-1", xml_declaration=True)


def write_gml() -> None:
    text = f"""<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>
<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>bgr_wp3_epfl_geometry</name>
    <points>
        <point id="0" x="0.0" y="0.0" z="0.0" name="origin"/>
        <point id="1" x="0.0" y="{TARGET_HEIGHT}" z="0.0"/>
        <point id="2" x="{TARGET_WIDTH}" y="0.0" z="0.0"/>
        <point id="3" x="{TARGET_WIDTH}" y="{TARGET_HEIGHT}" z="0.0"/>
    </points>
    <polylines>
        <polyline id="0" name="left"><pnt>0</pnt><pnt>1</pnt></polyline>
        <polyline id="1" name="right"><pnt>2</pnt><pnt>3</pnt></polyline>
        <polyline id="2" name="bottom"><pnt>0</pnt><pnt>2</pnt></polyline>
        <polyline id="3" name="top"><pnt>1</pnt><pnt>3</pnt></polyline>
    </polylines>
</OpenGeoSysGLI>
"""
    TARGET_GML.write_text(text)


def write_structured_column_mesh() -> None:
    ny = 2
    points = []
    for j in range(ny + 1):
        y = TARGET_HEIGHT * j / ny
        points.append(f"0 {y:.16g} 0")
        points.append(f"{TARGET_WIDTH:.16g} {y:.16g} 0")

    connectivity = []
    offsets = []
    for j in range(ny):
        base = 2 * j
        connectivity.extend([base, base + 1, base + 3, base + 2])
        offsets.append(4 * (j + 1))

    material_ids = " ".join("0" for _ in range(ny))
    connectivity_text = " ".join(str(v) for v in connectivity)
    offsets_text = " ".join(str(v) for v in offsets)
    types_text = " ".join("9" for _ in range(ny))
    points_text = "\n                    ".join(points)

    text = f"""<?xml version="1.0" encoding="ISO-8859-1"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
    <UnstructuredGrid>
        <Piece NumberOfPoints="{2 * (ny + 1)}" NumberOfCells="{ny}">
            <PointData/>
            <CellData Scalars="MaterialIDs">
                <DataArray type="Int32" Name="MaterialIDs" format="ascii">{material_ids}</DataArray>
            </CellData>
            <Points>
                <DataArray type="Float64" NumberOfComponents="3" format="ascii">
                    {points_text}
                </DataArray>
            </Points>
            <Cells>
                <DataArray type="Int64" Name="connectivity" format="ascii">{connectivity_text}</DataArray>
                <DataArray type="Int64" Name="offsets" format="ascii">{offsets_text}</DataArray>
                <DataArray type="UInt8" Name="types" format="ascii">{types_text}</DataArray>
            </Cells>
        </Piece>
    </UnstructuredGrid>
</VTKFile>
"""
    TARGET_MESH_STRUCTURED.write_text(text)


def main() -> None:
    scale_mesh()
    write_structured_column_mesh()
    write_gml()
    print(TARGET_MESH)
    print(TARGET_MESH_STRUCTURED)
    print(TARGET_GML)


if __name__ == "__main__":
    main()
