#!/usr/bin/env python3
"""Run the native notebook-style dense dry-density sweep.

The script builds per-density projects, runs OGS, and records the calibrated
native-branch response along with a comparison curve.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mplconfig-"))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.mkdtemp(prefix="xdg-cache-"))

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
RICHARDS_DATA_ROOT = ROOT.parent

RHO_SOLID = 2780.0  # kg/m^3
RHO_SOLID_REFERENCE = 2780.0  # kg/m^3 (aligned with MFront bridge run)
RHO_LR_REF = 1000.0  # kg/m^3

MS33_K0_REF = 5.6e-21
MS33_PHI_REF = 0.42

HAMAKER_LITERATURE = 5.1e-21  # J
SPECIFIC_SURFACE_MASS = 523.0  # m^2/g

TOTAL_SUCTION_MPA = 100.0
MACRO_SUCTION_MPA = 1.0
MICRO_SUCTION_MPA = TOTAL_SUCTION_MPA - MACRO_SUCTION_MPA

MASS_EXCHANGE_COEFFICIENT = 1e-13
MICRO_SWELLING_SLOPE = 0.1
PRESSURE_IC_PA = -MACRO_SUCTION_MPA * 1e6
TIME_END_S = 120 * 86400
PRESSURE_TOLERANCE_PA = 1e-12

AREA_FACTOR_TULLER = 1.0
PORE_AREA_SHAPE_FACTOR_TULLER = 0.8584073464102069
CHARACTERISTIC_PORE_SIZE = 1e-5
SURFACE_TENSION = 0.0715
MICRO_LIQUID_DENSITY_REFERENCE = 1e-6
MICRO_LIQUID_DENSITY_A = 1e-16
MICRO_LIQUID_DENSITY_B = 1.0

NATIVE_NOTEBOOK_SOURCE = Path("/Users/vinaykumar/git/ogs-native-dsm-transition")
DEFAULT_NATIVE_OGS = Path("/Users/vinaykumar/git/build/release-native-beacon/bin/ogs")
DEFAULT_MFRONT_CALIBRATION_CSV = ROOT / "villar_dense_dd_calibration.csv"


@dataclass(frozen=True)
class Case:
    dry_density: float

    @property
    def dry_density_g_cm3(self) -> float:
        return self.dry_density / 1000.0

    @property
    def phi0(self) -> float:
        return 1.0 - self.dry_density / RHO_SOLID

    @property
    def phi_micro_assumed(self) -> float:
        return 0.99 * self.phi0

    @property
    def phi_macro_assumed(self) -> float:
        return 0.01 * self.phi0

    @property
    def intrinsic_permeability(self) -> float:
        return (
            MS33_K0_REF
            * ((1.0 - MS33_PHI_REF) ** 2 / MS33_PHI_REF**3)
            * (self.phi0**3 / (1.0 - self.phi0) ** 2)
        )

    @property
    def villar_target_swelling_mpa(self) -> float:
        return math.exp(6.77 * self.dry_density_g_cm3 - 9.07)


def git_short_hash(repo: Path) -> str:
    """Return a short Git hash for provenance tracking."""
    try:
        return (
            subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"], cwd=repo, text=True
            )
            .strip()
        )
    except Exception:
        return ""


def read_grid(path: Path) -> vtk.vtkUnstructuredGrid:
    """Read a VTU snapshot into a VTK unstructured grid."""
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    return reader.GetOutput()


def get_array(grid: vtk.vtkUnstructuredGrid, name: str):
    """Fetch a field from point or cell data."""
    arr = grid.GetPointData().GetArray(name)
    if arr is None:
        arr = grid.GetCellData().GetArray(name)
    if arr is None:
        raise KeyError(name)
    return vtk_to_numpy(arr)


def mean_total_stress_mpa(vtu_path: Path) -> float:
    """Compute the mean isotropic stress in MPa from the final VTU."""
    grid = read_grid(vtu_path)
    sigma = get_array(grid, "sigma")
    p_mean = float((-sigma[:, 0] - sigma[:, 1] - sigma[:, 2]).mean() / 3.0)
    return p_mean / 1e6


def extract_last_vtu(prefix: str) -> Path:
    """Select the final VTU produced by a case run."""
    candidates = sorted(ROOT.glob(f"{prefix}_ts_*_t_*.vtu"))
    if not candidates:
        raise FileNotFoundError(f"No VTU outputs found for prefix {prefix}")

    def key(path: Path):
        stem = path.stem
        match = re.search(r"_ts_(\d+)_t_([-0-9eE+.]+)$", stem)
        if not match:
            return (-1, -1.0)
        return (int(match.group(1)), float(match.group(2)))

    return sorted(candidates, key=key)[-1]


def cleanup_runtime(prefix: str, project_path: Path) -> None:
    """Remove transient outputs from a calibration trial."""
    for pattern in (f"{prefix}.pvd", f"{prefix}_ts_*_t_*.vtu"):
        for p in ROOT.glob(pattern):
            p.unlink(missing_ok=True)
    project_path.unlink(missing_ok=True)


def n_l0_from_micro_suction(phi0: float, hamaker_eff: float) -> float:
    """Derive the initial micro-scale water content from the suction split."""
    # Target micro suction potential: mu = p / rho.
    mu_abs = MICRO_SUCTION_MPA * 1e6 / RHO_LR_REF
    n_s = 1.0 - phi0
    prefactor = (
        abs(hamaker_eff)
        * (SPECIFIC_SURFACE_MASS * n_s * RHO_SOLID_REFERENCE) ** 3
        / (6.0 * math.pi)
    )
    n_l0 = (prefactor / mu_abs) ** (1.0 / 3.0)
    return max(1e-12, n_l0)


def write_native_notebook_project(case: Case, project_path: Path) -> dict:
    """Create a temporary native notebook project for one density case."""
    prefix = project_path.stem
    hamaker_eff = HAMAKER_LITERATURE
    n_l0 = n_l0_from_micro_suction(case.phi0, hamaker_eff)
    n_s = 1.0 - case.phi0

    xml = f"""<?xml version='1.0' encoding='ISO-8859-1'?>
<OpenGeoSysProject>
    <meshes>
        <mesh axially_symmetric="true">../square_1x1_quad_1e0.vtu</mesh>
        <mesh axially_symmetric="true">../square_1x1_quad_1e0_left.vtu</mesh>
        <mesh axially_symmetric="true">../square_1x1_quad_1e0_right.vtu</mesh>
        <mesh axially_symmetric="true">../square_1x1_quad_1e0_top.vtu</mesh>
        <mesh axially_symmetric="true">../square_1x1_quad_1e0_bottom.vtu</mesh>
    </meshes>
    <processes>
        <process>
            <name>RM</name>
            <type>RICHARDS_MECHANICS</type>
            <integration_order>2</integration_order>
            <jacobian_assembler><type>Analytical</type></jacobian_assembler>
            <micro_porosity>
                <mass_exchange_coefficient>{MASS_EXCHANGE_COEFFICIENT:.16g}</mass_exchange_coefficient>
                <nonlinear_solver>
                    <maximum_iterations>100</maximum_iterations>
                    <residuum_tolerance>1e-8</residuum_tolerance>
                    <increment_tolerance>1e-20</increment_tolerance>
                </nonlinear_solver>
            </micro_porosity>
            <vk_potential_exchange>
                <enabled>true</enabled>
                <mode>full_potential</mode>
                <pressure_tolerance>{PRESSURE_TOLERANCE_PA:.16g}</pressure_tolerance>
                <hamaker_constant>{hamaker_eff:.16g}</hamaker_constant>
                <specific_surface>{SPECIFIC_SURFACE_MASS:.16g}</specific_surface>
                <micro_solid_density_reference>{RHO_SOLID_REFERENCE:.16g}</micro_solid_density_reference>
                <micro_solid_volume_fraction_reference>{n_s:.16g}</micro_solid_volume_fraction_reference>
                <micro_liquid_density_reference>{MICRO_LIQUID_DENSITY_REFERENCE:.16g}</micro_liquid_density_reference>
                <micro_liquid_density_a>{MICRO_LIQUID_DENSITY_A:.16g}</micro_liquid_density_a>
                <micro_liquid_density_b>{MICRO_LIQUID_DENSITY_B:.16g}</micro_liquid_density_b>
                <initial_micro_water_content>{n_l0:.16g}</initial_micro_water_content>
                <local_nonlinear_solve_mode>scalar_notebook_mass_storage</local_nonlinear_solve_mode>
                <potential_role_mapping>notebook_roles</potential_role_mapping>
                <fd_jacobian_for_exchange>false</fd_jacobian_for_exchange>
                <micro_potential_convention>negative_attractive</micro_potential_convention>
                <micro_water_content_swelling_slope>{MICRO_SWELLING_SLOPE:.16g}</micro_water_content_swelling_slope>
            </vk_potential_exchange>
            <constitutive_relation>
                <type>LinearElasticIsotropic</type>
                <youngs_modulus>YoungModulus</youngs_modulus>
                <poissons_ratio>PoissonRatio</poissons_ratio>
            </constitutive_relation>
            <process_variables>
                <pressure>pressure</pressure>
                <displacement>displacement</displacement>
            </process_variables>
            <secondary_variables>
                <secondary_variable name="sigma"/>
                <secondary_variable name="swelling_stress"/>
                <secondary_variable name="saturation"/>
                <secondary_variable name="porosity"/>
                <secondary_variable name="dry_density_solid"/>
                <secondary_variable name="micro_pressure"/>
                <secondary_variable name="micro_saturation"/>
                <secondary_variable name="vk_micro_water_content"/>
                <secondary_variable name="vk_micro_porosity"/>
                <secondary_variable name="vk_micro_exchange_source"/>
            </secondary_variables>
            <specific_body_force>0 0</specific_body_force>
            <initial_stress>sigma0</initial_stress>
            <mass_lumping>true</mass_lumping>
        </process>
    </processes>
    <media>
        <medium>
            <phases>
                <phase>
                    <type>AqueousLiquid</type>
                    <properties>
                        <property><name>viscosity</name><type>Constant</type><value>1e-3</value></property>
                        <property><name>density</name><type>Constant</type><value>1e3</value></property>
                    </properties>
                </phase>
                <phase>
                    <type>Solid</type>
                    <properties>
                        <property><name>density</name><type>Constant</type><value>{RHO_SOLID:.16g}</value></property>
                        <property>
                            <name>swelling_stress_rate</name>
                            <type>SaturationDependentSwelling</type>
                            <swelling_pressures>0 0 0</swelling_pressures>
                            <exponents>1 1 1</exponents>
                            <lower_saturation_limit>0</lower_saturation_limit>
                            <upper_saturation_limit>1</upper_saturation_limit>
                        </property>
                    </properties>
                </phase>
            </phases>
            <properties>
                <property>
                    <name>biot_coefficient</name>
                    <type>Constant</type>
                    <value>1.0</value>
                </property>
                <property>
                    <name>permeability</name>
                    <type>KozenyCarman</type>
                    <initial_permeability>IntrinsicPermeability0</initial_permeability>
                    <initial_porosity>phi0</initial_porosity>
                </property>
                <property>
                    <name>porosity</name>
                    <type>PorosityFromMassBalance</type>
                    <initial_porosity>phi0</initial_porosity>
                    <minimal_porosity>0</minimal_porosity>
                    <maximal_porosity>1</maximal_porosity>
                </property>
                <property>
                    <name>reference_temperature</name>
                    <type>Constant</type>
                    <value>293.15</value>
                </property>
                <property>
                    <name>relative_permeability</name>
                    <type>Constant</type>
                    <value>1</value>
                </property>
                <property>
                    <name>saturation</name>
                    <type>SaturationTuller</type>
                    <area_factor_tuller>{AREA_FACTOR_TULLER:.16g}</area_factor_tuller>
                    <pore_area_shapefactor_tuller>{PORE_AREA_SHAPE_FACTOR_TULLER:.16g}</pore_area_shapefactor_tuller>
                    <characteristic_pore_size>{CHARACTERISTIC_PORE_SIZE:.16g}</characteristic_pore_size>
                    <surface_tension>{SURFACE_TENSION:.16g}</surface_tension>
                    <residual_liquid_saturation>0</residual_liquid_saturation>
                    <residual_gas_saturation>0</residual_gas_saturation>
                </property>
                <property>
                    <name>bishops_effective_stress</name>
                    <type>BishopsSaturationCutoff</type>
                    <cutoff_value>1</cutoff_value>
                </property>
            </properties>
        </medium>
    </media>
    <time_loop>
        <processes>
            <process ref="RM">
                <nonlinear_solver>basic_newton</nonlinear_solver>
                <convergence_criterion>
                    <type>PerComponentDeltaX</type>
                    <norm_type>NORM2</norm_type>
                    <abstols>5e-8 1e-13 1e-13</abstols>
                </convergence_criterion>
                <time_discretization><type>BackwardEuler</type></time_discretization>
                <time_stepping>
                    <type>FixedTimeStepping</type>
                    <t_initial>0</t_initial>
                    <t_end>{TIME_END_S}</t_end>
                    <timesteps><pair><repeat>120</repeat><delta_t>86400</delta_t></pair></timesteps>
                </time_stepping>
            </process>
        </processes>
        <output>
            <type>VTK</type>
            <prefix>{prefix}</prefix>
            <suffix>_ts_{{:timestep}}_t_{{:time}}</suffix>
            <fixed_output_times>{TIME_END_S}</fixed_output_times>
            <variables>
                <variable>pressure</variable>
                <variable>sigma</variable>
                <variable>swelling_stress</variable>
                <variable>saturation</variable>
                <variable>porosity</variable>
                <variable>dry_density_solid</variable>
                <variable>micro_pressure</variable>
                <variable>micro_saturation</variable>
                <variable>vk_micro_water_content</variable>
                <variable>vk_micro_porosity</variable>
                <variable>vk_micro_exchange_source</variable>
            </variables>
        </output>
    </time_loop>
    <parameters>
        <parameter><name>sigma0</name><type>Function</type><expression>0</expression><expression>0</expression><expression>0</expression><expression>0</expression></parameter>
        <parameter><name>YoungModulus</name><type>Constant</type><value>52e6</value></parameter>
        <parameter><name>PoissonRatio</name><type>Constant</type><value>0.3</value></parameter>
        <parameter><name>phi0</name><type>Constant</type><value>{case.phi0:.16g}</value></parameter>
        <parameter><name>IntrinsicPermeability0</name><type>Constant</type><value>{case.intrinsic_permeability:.16g}</value></parameter>
        <parameter><name>displacement0</name><type>Constant</type><values>0 0</values></parameter>
        <parameter><name>zero</name><type>Constant</type><value>0.0</value></parameter>
        <parameter><name>pressure_ic</name><type>Constant</type><value>{PRESSURE_IC_PA:.16g}</value></parameter>
        <parameter><name>pressure_bc_scale</name><type>Constant</type><value>1</value></parameter>
        <parameter><name>pressure_bc</name><type>CurveScaled</type><curve>pressure_release</curve><parameter>pressure_bc_scale</parameter></parameter>
    </parameters>
    <process_variables>
        <process_variable>
            <name>displacement</name>
            <components>2</components>
            <order>1</order>
            <initial_condition>displacement0</initial_condition>
            <boundary_conditions>
                <boundary_condition><mesh>square_1x1_quad_1e0_left</mesh><type>Dirichlet</type><component>0</component><parameter>zero</parameter></boundary_condition>
                <boundary_condition><mesh>square_1x1_quad_1e0_right</mesh><type>Dirichlet</type><component>0</component><parameter>zero</parameter></boundary_condition>
                <boundary_condition><mesh>square_1x1_quad_1e0_bottom</mesh><type>Dirichlet</type><component>1</component><parameter>zero</parameter></boundary_condition>
                <boundary_condition><mesh>square_1x1_quad_1e0_top</mesh><type>Dirichlet</type><component>1</component><parameter>zero</parameter></boundary_condition>
            </boundary_conditions>
        </process_variable>
        <process_variable>
            <name>pressure</name>
            <components>1</components>
            <order>1</order>
            <initial_condition>pressure_ic</initial_condition>
            <boundary_conditions>
                <boundary_condition><mesh>square_1x1_quad_1e0_left</mesh><type>Dirichlet</type><parameter>pressure_bc</parameter></boundary_condition>
                <boundary_condition><mesh>square_1x1_quad_1e0_right</mesh><type>Dirichlet</type><parameter>pressure_bc</parameter></boundary_condition>
                <boundary_condition><mesh>square_1x1_quad_1e0_bottom</mesh><type>Dirichlet</type><parameter>pressure_bc</parameter></boundary_condition>
                <boundary_condition><mesh>square_1x1_quad_1e0_top</mesh><type>Dirichlet</type><parameter>pressure_bc</parameter></boundary_condition>
            </boundary_conditions>
        </process_variable>
    </process_variables>
    <nonlinear_solvers>
        <nonlinear_solver><name>basic_newton</name><type>Newton</type><max_iter>60</max_iter><linear_solver>general_linear_solver</linear_solver></nonlinear_solver>
    </nonlinear_solvers>
    <linear_solvers>
        <linear_solver><name>general_linear_solver</name><eigen><solver_type>SparseLU</solver_type><scaling>true</scaling></eigen></linear_solver>
    </linear_solvers>
    <curves>
        <curve><name>pressure_release</name><coords>0 {TIME_END_S}</coords><values>{PRESSURE_IC_PA:.16g} 0</values></curve>
    </curves>
</OpenGeoSysProject>
"""
    project_path.write_text(xml)
    return {
        "hamaker_effective_J": hamaker_eff,
        "n_l0": n_l0,
        "micro_solid_volume_fraction_reference": n_s,
    }


def run_ogs(ogs_bin: Path, project_path: Path) -> None:
    """Execute OGS for a single temporary project."""
    subprocess.run(
        [str(ogs_bin), str(project_path)],
        cwd=ROOT,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )


def run_native_notebook_case(ogs_bin: Path, case: Case) -> dict:
    case_tag = f"dd{int(case.dry_density)}_native_notebook"
    project_path = ROOT / f"{case_tag}.prj"
    meta = write_native_notebook_project(case, project_path)
    try:
        run_ogs(ogs_bin, project_path)
        last_vtu = extract_last_vtu(case_tag)
        pressure_mpa = mean_total_stress_mpa(last_vtu)
    finally:
        cleanup_runtime(case_tag, project_path)
    return {"pressure_mpa": pressure_mpa, **meta}


def load_mfront_calibrated_curve(path: Path) -> tuple[np.ndarray, np.ndarray] | None:
    if not path.exists():
        return None
    x = []
    y = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            x.append(float(row["dry_density_kg_m3"]))
            y.append(float(row["mfront_calibrated_MPa"]))
    if not x:
        return None
    order = np.argsort(np.asarray(x))
    return np.asarray(x)[order], np.asarray(y)[order]


def main() -> None:
    """Run the dense sweep and write the CSV, JSON, and comparison plots."""
    parser = argparse.ArgumentParser(
        description=(
            "Villar-style dry-density sweep for the native notebook branch "
            "(micro-enabled native implementation)."
        )
    )
    parser.add_argument("--native-ogs", type=Path, default=DEFAULT_NATIVE_OGS)
    parser.add_argument(
        "--native-source",
        type=Path,
        default=NATIVE_NOTEBOOK_SOURCE,
        help="Native OGS source tree for commit-hash provenance in summary JSON.",
    )
    parser.add_argument("--dd-min", type=float, default=1400.0)
    parser.add_argument("--dd-max", type=float, default=1800.0)
    parser.add_argument("--dd-step", type=float, default=25.0)
    parser.add_argument(
        "--mfront-calibration-csv",
        type=Path,
        default=DEFAULT_MFRONT_CALIBRATION_CSV,
        help="Optional MFront calibration curve for overlay.",
    )
    args = parser.parse_args()

    dd_values = np.arange(args.dd_min, args.dd_max + 0.5 * args.dd_step, args.dd_step)
    cases = [Case(float(dd)) for dd in dd_values]

    rows = []
    for case in cases:
        result = run_native_notebook_case(args.native_ogs, case)
        target = case.villar_target_swelling_mpa
        native_p = result["pressure_mpa"]
        print(
            f"dd={case.dry_density:.0f} kg/m3: "
            f"Villar={target:.3f} MPa, "
            f"native-notebook={native_p:.3f} MPa"
        )
        rows.append(
            {
                "dry_density_kg_m3": case.dry_density,
                "dry_density_g_cm3": case.dry_density_g_cm3,
                "phi0": case.phi0,
                "phi_micro_assumed": case.phi_micro_assumed,
                "phi_macro_assumed": case.phi_macro_assumed,
                "k0_m2": case.intrinsic_permeability,
                "hamaker_literature_J": HAMAKER_LITERATURE,
                "hamaker_effective_J": result["hamaker_effective_J"],
                "specific_surface_mass_m2_g": SPECIFIC_SURFACE_MASS,
                "n_l0": result["n_l0"],
                "micro_solid_volume_fraction_reference": result[
                    "micro_solid_volume_fraction_reference"
                ],
                "target_villar_MPa": target,
                "native_notebook_branch_MPa": native_p,
                "native_minus_target_MPa": native_p - target,
            }
        )

    rows = sorted(rows, key=lambda r: r["dry_density_kg_m3"])

    out_csv = ROOT / "villar_dense_dd_native_notebook_branch.csv"
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "ogs_repo_hash": git_short_hash(RICHARDS_DATA_ROOT.parents[2]),
        "native_notebook_source_hash": git_short_hash(args.native_source),
        "native_notebook_source_path": str(args.native_source),
        "native_ogs_version": subprocess.check_output(
            [str(args.native_ogs), "--version"], text=True
        ),
        "total_initial_suction_MPa": TOTAL_SUCTION_MPA,
        "macro_initial_suction_MPa": MACRO_SUCTION_MPA,
        "micro_initial_suction_MPa": MICRO_SUCTION_MPA,
        "hamaker_literature_J": HAMAKER_LITERATURE,
        "specific_surface_mass_m2_g": SPECIFIC_SURFACE_MASS,
        "results": rows,
    }
    out_summary = ROOT / "villar_dense_dd_native_notebook_branch_summary.json"
    out_summary.write_text(json.dumps(summary, indent=2))

    x = np.array([r["dry_density_kg_m3"] for r in rows])
    y_target = np.array([r["target_villar_MPa"] for r in rows])
    y_native = np.array([r["native_notebook_branch_MPa"] for r in rows])

    plt.figure(figsize=(8.0, 5.2))
    plt.plot(x, y_target, "k-", linewidth=2.0, label="Villar Eq. (7) target")
    plt.plot(
        x,
        y_native,
        color="#d62728",
        marker="s",
        linestyle=":",
        linewidth=1.6,
        markersize=3.8,
        label="Native notebook branch",
    )

    mfront_curve = load_mfront_calibrated_curve(args.mfront_calibration_csv)
    if mfront_curve is not None:
        x_m, y_m = mfront_curve
        plt.plot(
            x_m,
            y_m,
            color="black",
            marker="o",
            linestyle="-",
            linewidth=1.4,
            markersize=3.4,
            label="Villar fit reference",
        )

    plt.xlabel("Dry density (kg/m$^3$)")
    plt.ylabel("Swelling pressure at full saturation (MPa)")
    plt.grid(True, alpha=0.35)
    plt.legend(frameon=False)
    plt.tight_layout()
    out_png = ROOT / "villar_dense_dd_native_notebook_branch_vs_villar.png"
    plt.savefig(out_png, dpi=220)
    plt.close()

    print(f"Wrote: {out_csv}")
    print(f"Wrote: {out_summary}")
    print(f"Wrote: {out_png}")


if __name__ == "__main__":
    main()
