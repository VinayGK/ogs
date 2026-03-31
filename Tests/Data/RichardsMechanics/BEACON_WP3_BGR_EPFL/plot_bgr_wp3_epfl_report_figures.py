#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mplconfig-"))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.mkdtemp(prefix="xdg-cache-"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import vtk
from PIL import Image, ImageDraw, ImageFont
from vtk.util.numpy_support import vtk_to_numpy


HERE = Path(__file__).resolve().parent
MFRONT_DEFAULT = Path("/Users/vinaykumar/git/build/release-mfront-tpm/bin/ogs")
NATIVE_DEFAULT = MFRONT_DEFAULT
NATIVE_PROJECT = HERE / "bgr_wp3_p2_1_abprime_native.prj"
MFRONT_PROJECT = HERE / "bgr_wp3_p2_1_abprime_mfront.prj"
EXPERIMENT_SWELLING_PRESSURE_RANGE_MPA = (3.12, 3.55)


@dataclass
class Series:
    label: str
    times_s: np.ndarray
    times_min: np.ndarray
    axial_swelling_pressure_mpa: np.ndarray
    vertical_stress_mpa: np.ndarray
    void_ratio: np.ndarray
    saturation: np.ndarray


def make_output_times(t_end_s: float, num_times: int) -> list[float]:
    raw = np.geomspace(60.0, t_end_s, num=num_times)
    times = np.unique(np.round(raw, 6))
    times = np.concatenate(([0.0], times))
    times[-1] = t_end_s
    return [float(t) for t in times]


def write_project_copy(source: Path, target: Path, prefix: str, output_times: list[float]) -> None:
    tree = ET.parse(source)
    root = tree.getroot()

    output = root.find("./time_loop/output")
    if output is None:
        raise RuntimeError(f"No output block found in {source}")

    prefix_node = output.find("prefix")
    if prefix_node is None:
        raise RuntimeError(f"No output prefix found in {source}")
    prefix_node.text = prefix

    fixed_times = output.find("fixed_output_times")
    if fixed_times is None:
        fixed_times = ET.SubElement(output, "fixed_output_times")
    fixed_times.text = " ".join(f"{t:.6f}" for t in output_times[1:])

    tree.write(target, encoding="ISO-8859-1", xml_declaration=True)


def run_ogs(ogs: Path, project: Path, workdir: Path) -> None:
    subprocess.run(
        [str(ogs), str(project)],
        cwd=workdir,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )


def load_vtu(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
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


def time_from_filename(path: Path) -> float:
    stem = path.stem
    marker = "_t_"
    if marker not in stem:
        raise RuntimeError(f"Could not parse time from {path.name}")
    return float(stem.split(marker, 1)[1])


def collect_series(run_dir: Path, prefix: str, label: str) -> Series:
    paths = sorted(run_dir.glob(f"{prefix}_t_*.vtu"), key=time_from_filename)
    if not paths:
        raise RuntimeError(f"No VTU outputs found for {prefix} in {run_dir}")

    times_s = []
    axial_swelling_pressure_mpa = []
    vertical_stress_mpa = []
    void_ratio = []
    saturation = []

    for path in paths:
        points, data = load_vtu(path)
        sigma_top = boundary_mean(points, data["sigma"], axis=1, side="max")
        phi = np.asarray(data["porosity"]).mean()
        times_s.append(time_from_filename(path))
        axial_swelling_pressure_mpa.append(abs(float(sigma_top[1])) / 1.0e6)
        vertical_stress_mpa.append(abs(float(sigma_top[1])) / 1.0e6)
        void_ratio.append(float(phi / (1.0 - phi)))
        saturation.append(float(np.asarray(data["saturation"]).mean()))

    times_s_arr = np.asarray(times_s)
    order = np.argsort(times_s_arr)
    times_s_arr = times_s_arr[order]

    return Series(
        label=label,
        times_s=times_s_arr,
        times_min=times_s_arr / 60.0,
        axial_swelling_pressure_mpa=np.asarray(axial_swelling_pressure_mpa)[order],
        vertical_stress_mpa=np.asarray(vertical_stress_mpa)[order],
        void_ratio=np.asarray(void_ratio)[order],
        saturation=np.asarray(saturation)[order],
    )


def write_series_csv(path: Path, series: list[Series]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "series",
                "time_s",
                "time_min",
                "axial_swelling_pressure_mpa",
                "vertical_stress_mpa",
                "void_ratio",
                "mean_saturation",
            ]
        )
        for s in series:
            for row in zip(
                s.times_s,
                s.times_min,
                s.axial_swelling_pressure_mpa,
                s.vertical_stress_mpa,
                s.void_ratio,
                s.saturation,
                strict=True,
                ):
                    writer.writerow([s.label, *row])


def read_series_csv(path: Path) -> list[Series]:
    buckets: dict[str, dict[str, list[float]]] = {}
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            label = row["series"]
            if label not in buckets:
                buckets[label] = {
                    "time_s": [],
                    "time_min": [],
                    "axial_swelling_pressure_mpa": [],
                    "vertical_stress_mpa": [],
                    "void_ratio": [],
                    "mean_saturation": [],
                }
            buckets[label]["time_s"].append(float(row["time_s"]))
            buckets[label]["time_min"].append(float(row["time_min"]))
            buckets[label]["axial_swelling_pressure_mpa"].append(float(row["axial_swelling_pressure_mpa"]))
            buckets[label]["vertical_stress_mpa"].append(float(row["vertical_stress_mpa"]))
            buckets[label]["void_ratio"].append(float(row["void_ratio"]))
            buckets[label]["mean_saturation"].append(float(row["mean_saturation"]))

    series = []
    for label, values in buckets.items():
        series.append(
            Series(
                label=label,
                times_s=np.asarray(values["time_s"]),
                times_min=np.asarray(values["time_min"]),
                axial_swelling_pressure_mpa=np.asarray(values["axial_swelling_pressure_mpa"]),
                vertical_stress_mpa=np.asarray(values["vertical_stress_mpa"]),
                void_ratio=np.asarray(values["void_ratio"]),
                saturation=np.asarray(values["mean_saturation"]),
            )
        )
    return series


def plot_swelling_pressure(series: list[Series], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.2))

    ax.axhspan(
        EXPERIMENT_SWELLING_PRESSURE_RANGE_MPA[0],
        EXPERIMENT_SWELLING_PRESSURE_RANGE_MPA[1],
        color="#d9d9d9",
        alpha=0.5,
        label="EPFL report range",
    )

    colors = {"Native OGS": "#d62728", "MFront OGS": "#1f77b4"}
    for s in series:
        ax.plot(
            s.times_min,
            s.axial_swelling_pressure_mpa,
            marker="o",
            markersize=3.2,
            linewidth=1.9,
            color=colors[s.label],
            label=s.label,
        )

    ax.set_title("EPFL P2-1 / A-B': Vertical Swelling Pressure With Time")
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Vertical swelling pressure (MPa)")
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=True, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_void_ratio(series: list[Series], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.6))

    colors = {"Native OGS": "#d62728", "MFront OGS": "#1f77b4"}
    for s in series:
        ax.plot(
            s.vertical_stress_mpa,
            s.void_ratio,
            marker="o",
            markersize=3.4,
            linewidth=1.9,
            color=colors[s.label],
            label=s.label,
        )
        ax.scatter(
            [s.vertical_stress_mpa[0], s.vertical_stress_mpa[-1]],
            [s.void_ratio[0], s.void_ratio[-1]],
            color=colors[s.label],
            s=24,
            zorder=3,
        )

    ax.axvspan(
        EXPERIMENT_SWELLING_PRESSURE_RANGE_MPA[0],
        EXPERIMENT_SWELLING_PRESSURE_RANGE_MPA[1],
        color="#d9d9d9",
        alpha=0.35,
        label="EPFL B' stress range",
    )

    ax.set_xscale("log")
    ax.set_title("EPFL P2-1 / A-B': Void Ratio - Vertical Stress Path")
    ax.set_xlabel("Total vertical stress (MPa)")
    ax.set_ylabel("Void ratio (-)")
    ax.set_xlim(left=1.0e-2)
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(frameon=True, loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_composite(swelling_path: Path, void_ratio_path: Path, out_path: Path) -> None:
    img1 = Image.open(swelling_path)
    img2 = Image.open(void_ratio_path)
    width = max(img1.width, img2.width)
    margin = 24
    title_h = 44
    canvas = Image.new(
        "RGB",
        (width + 2 * margin, title_h + img1.height + img2.height + 3 * margin),
        "white",
    )
    draw = ImageDraw.Draw(canvas)
    draw.text(
        (margin, 12),
        "EPFL report-style figures for reduced BGR P2-1 / A-B'",
        fill="black",
        font=ImageFont.load_default(),
    )
    canvas.paste(img1, (margin, title_h))
    canvas.paste(img2, (margin, title_h + img1.height + margin))
    canvas.save(out_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--native-ogs", type=Path, default=NATIVE_DEFAULT)
    parser.add_argument("--mfront-ogs", type=Path, default=MFRONT_DEFAULT)
    parser.add_argument("--num-times", type=int, default=48)
    parser.add_argument("--t-end-s", type=float, default=2.0e5)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    parser.add_argument("--from-csv", type=Path)
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    series_csv = output_dir / "bgr_wp3_p2_1_abprime_report_style_series.csv"
    swelling_png = output_dir / "bgr_wp3_p2_1_abprime_fig4_5b_like.png"
    void_ratio_png = output_dir / "bgr_wp3_p2_1_abprime_fig4_6_like.png"
    composite_png = output_dir / "bgr_wp3_p2_1_abprime_report_style_panel.png"

    if args.from_csv:
        print(f"Loading series from {args.from_csv}", flush=True)
        series = read_series_csv(args.from_csv)
    else:
        print("Running native and MFront EPFL replays", flush=True)
        output_times = make_output_times(args.t_end_s, args.num_times)

        with tempfile.TemporaryDirectory(prefix="bgr-epfl-native-") as native_tmp, tempfile.TemporaryDirectory(
            prefix="bgr-epfl-mfront-"
        ) as mfront_tmp:
            native_run = Path(native_tmp)
            mfront_run = Path(mfront_tmp)

            for src in (HERE / "bgr_wp3_epfl_domain_2e.vtu", HERE / "bgr_wp3_epfl_geometry.gml"):
                shutil.copy2(src, native_run / src.name)
                shutil.copy2(src, mfront_run / src.name)

            native_prefix = "bgr_wp3_p2_1_abprime_native_dense"
            mfront_prefix = "bgr_wp3_p2_1_abprime_mfront_dense"
            native_project = native_run / NATIVE_PROJECT.name
            mfront_project = mfront_run / MFRONT_PROJECT.name

            write_project_copy(NATIVE_PROJECT, native_project, native_prefix, output_times)
            write_project_copy(MFRONT_PROJECT, mfront_project, mfront_prefix, output_times)

            run_ogs(args.native_ogs, native_project, native_run)
            run_ogs(args.mfront_ogs, mfront_project, mfront_run)

            native_series = collect_series(native_run, native_prefix, "Native OGS")
            mfront_series = collect_series(mfront_run, mfront_prefix, "MFront OGS")
            series = [native_series, mfront_series]
            write_series_csv(series_csv, series)

    print(f"Writing {swelling_png.name}", flush=True)
    plot_swelling_pressure(series, swelling_png)
    print(f"Writing {void_ratio_png.name}", flush=True)
    plot_void_ratio(series, void_ratio_png)
    print(f"Writing {composite_png.name}", flush=True)
    plot_composite(swelling_png, void_ratio_png, composite_png)

    print(series_csv if series_csv.exists() else args.from_csv)
    print(swelling_png)
    print(void_ratio_png)
    print(composite_png)


if __name__ == "__main__":
    main()
