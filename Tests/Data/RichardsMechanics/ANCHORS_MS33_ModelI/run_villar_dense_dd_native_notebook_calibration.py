#!/usr/bin/env python3
"""Calibrate the native notebook dense dry-density sweep against Villar data.

This workflow mirrors the MFront calibration logic so the two implementations
can be compared with the same dry-density sampling and output schema.
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


def write_native_notebook_project(
    case: Case, project_path: Path, vdw_multiplier: float, n_l0_fixed: float
) -> dict:
    """Create a temporary native notebook project for one multiplier trial."""
    prefix = project_path.stem
    hamaker_eff = HAMAKER_LITERATURE * vdw_multiplier
    n_l0 = n_l0_fixed
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


def run_native_notebook_case(
    ogs_bin: Path, case: Case, multiplier: float, n_l0_fixed: float, idx: int
) -> dict:
    case_tag = f"dd{int(case.dry_density)}_native_notebook_{idx:02d}"
    project_path = ROOT / f"{case_tag}.prj"
    meta = write_native_notebook_project(case, project_path, multiplier, n_l0_fixed)
    try:
        run_ogs(ogs_bin, project_path)
        last_vtu = extract_last_vtu(case_tag)
        pressure_mpa = mean_total_stress_mpa(last_vtu)
    finally:
        cleanup_runtime(case_tag, project_path)
    return {"pressure_mpa": pressure_mpa, "multiplier": multiplier, **meta}


def calibrate_multiplier_for_case(
    ogs_bin: Path,
    case: Case,
    n_l0_fixed: float,
    target_mpa: float,
    rel_tol: float = 0.02,
    max_iter: int = 18,
) -> dict:
    """Solve for the multiplier that matches the Villar target."""
    target_scale = max(target_mpa, 1e-12)
    max_multiplier = 1e18

    def rel_err(run: dict) -> float:
        pressure = run["pressure_mpa"]
        if not math.isfinite(pressure):
            return float("inf")
        return abs(pressure - target_mpa) / target_scale

    # Baseline.
    run_id = 0
    run_lo = run_native_notebook_case(ogs_bin, case, 1.0, n_l0_fixed, run_id)
    run_id += 1
    best = run_lo
    if rel_err(run_lo) < rel_tol:
        return run_lo

    # Bracket from above with geometric expansion.
    m_lo = 1.0
    m_hi = max(2.0, target_mpa / max(abs(run_lo["pressure_mpa"]), 1e-12))
    m_hi = float(np.clip(m_hi, 1e-6, max_multiplier))
    if m_hi <= m_lo:
        m_hi = 2.0

    run_hi = run_native_notebook_case(ogs_bin, case, m_hi, n_l0_fixed, run_id)
    run_id += 1
    if rel_err(run_hi) < rel_err(best):
        best = run_hi
    p_hi = run_hi["pressure_mpa"] if math.isfinite(run_hi["pressure_mpa"]) else float("inf")

    expand_limit = max_iter + 16
    while p_hi < target_mpa and m_hi < max_multiplier and run_id < expand_limit:
        m_lo, run_lo = m_hi, run_hi
        m_hi = min(max_multiplier, m_hi * 10.0)
        if m_hi <= m_lo:
            break
        run_hi = run_native_notebook_case(ogs_bin, case, m_hi, n_l0_fixed, run_id)
        run_id += 1
        if rel_err(run_hi) < rel_err(best):
            best = run_hi
        p_hi = (
            run_hi["pressure_mpa"]
            if math.isfinite(run_hi["pressure_mpa"])
            else float("inf")
        )

    if p_hi < target_mpa:
        return best

    # Log-space bisection inside bracket [m_lo, m_hi].
    for _ in range(max_iter):
        m_mid = math.sqrt(m_lo * m_hi)
        if m_mid <= m_lo * (1.0 + 1e-12) or m_mid >= m_hi * (1.0 - 1e-12):
            break
        run_mid = run_native_notebook_case(ogs_bin, case, m_mid, n_l0_fixed, run_id)
        run_id += 1
        err_mid = rel_err(run_mid)
        if err_mid < rel_err(best):
            best = run_mid
        if err_mid < rel_tol:
            return run_mid

        p_mid = (
            run_mid["pressure_mpa"]
            if math.isfinite(run_mid["pressure_mpa"])
            else float("inf")
        )
        if p_mid < target_mpa:
            m_lo, run_lo = m_mid, run_mid
        else:
            m_hi, run_hi = m_mid, run_mid

    return best


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
            "Dense dry-density calibration of effective micro vdW multiplier "
            "for the native notebook branch against Villar Eq.(7)."
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
    parser.add_argument("--rel-tol", type=float, default=0.02)
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
        target = case.villar_target_swelling_mpa
        # Keep n_l0 fixed per dry density from literature Hamaker baseline,
        # matching the MFront dense calibration workflow.
        n_l0_fixed = n_l0_from_micro_suction(case.phi0, HAMAKER_LITERATURE)
        baseline = run_native_notebook_case(args.native_ogs, case, 1.0, n_l0_fixed, 90)
        calibrated = calibrate_multiplier_for_case(
            args.native_ogs, case, n_l0_fixed, target, rel_tol=args.rel_tol
        )
        print(
            f"dd={case.dry_density:.0f} kg/m3: "
            f"target={target:.3f} MPa, "
            f"native_baseline={baseline['pressure_mpa']:.3f} MPa, "
            f"native_calibrated={calibrated['pressure_mpa']:.3f} MPa, "
            f"mult={calibrated['multiplier']:.3e}"
        )
        rows.append(
            {
                "dry_density_kg_m3": case.dry_density,
                "dry_density_g_cm3": case.dry_density_g_cm3,
                "phi0": case.phi0,
                "phi_micro_assumed": case.phi_micro_assumed,
                "phi_macro_assumed": case.phi_macro_assumed,
                "k0_m2": case.intrinsic_permeability,
                "target_villar_MPa": target,
                "native_baseline_MPa": baseline["pressure_mpa"],
                "native_calibrated_MPa": calibrated["pressure_mpa"],
                "vdw_multiplier": calibrated["multiplier"],
                "hamaker_literature_J": HAMAKER_LITERATURE,
                "hamaker_effective_J": calibrated["hamaker_effective_J"],
                "specific_surface_mass_m2_g": SPECIFIC_SURFACE_MASS,
                "n_l0": calibrated["n_l0"],
                "micro_solid_volume_fraction_reference": calibrated[
                    "micro_solid_volume_fraction_reference"
                ],
                "macro_suction_MPa": MACRO_SUCTION_MPA,
                "micro_suction_MPa": MICRO_SUCTION_MPA,
                "native_baseline_minus_target_MPa": baseline["pressure_mpa"] - target,
                "native_calibrated_minus_target_MPa": calibrated["pressure_mpa"] - target,
            }
        )

    rows = sorted(rows, key=lambda r: r["dry_density_kg_m3"])

    out_csv = ROOT / "villar_dense_dd_native_notebook_calibration.csv"
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    rel_err_baseline = [
        abs(r["native_baseline_minus_target_MPa"]) / max(r["target_villar_MPa"], 1e-12)
        for r in rows
    ]
    rel_err_calibrated = [
        abs(r["native_calibrated_minus_target_MPa"]) / max(r["target_villar_MPa"], 1e-12)
        for r in rows
    ]

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
        "mean_relative_error_baseline_percent": 100.0 * float(np.mean(rel_err_baseline)),
        "mean_relative_error_calibrated_percent": 100.0 * float(np.mean(rel_err_calibrated)),
        "max_relative_error_calibrated_percent": 100.0 * float(np.max(rel_err_calibrated)),
        "results": rows,
    }
    out_summary = ROOT / "villar_dense_dd_native_notebook_calibration_summary.json"
    out_summary.write_text(json.dumps(summary, indent=2))

    x = np.array([r["dry_density_kg_m3"] for r in rows])
    y_target = np.array([r["target_villar_MPa"] for r in rows])
    y_calibrated = np.array([r["native_calibrated_MPa"] for r in rows])
    y_mult = np.array([r["vdw_multiplier"] for r in rows])

    plt.figure(figsize=(8.2, 5.2))
    plt.plot(x, y_target, "k-", linewidth=2.0, label="Villar Eq. (7) target")
    plt.plot(
        x,
        y_calibrated,
        color="#2ca02c",
        marker="D",
        linestyle="--",
        linewidth=1.8,
        markersize=4.2,
        label="Native notebook calibrated",
    )
    mfront_curve = load_mfront_calibrated_curve(args.mfront_calibration_csv)
    if mfront_curve is not None:
        x_m, y_m = mfront_curve
        plt.plot(
            x_m,
            y_m,
            color="#d62728",
            marker="s",
            linestyle="-.",
            linewidth=1.5,
            markersize=3.6,
            label="MFront calibrated (existing run)",
        )
    plt.xlabel("Dry density (kg/m$^3$)")
    plt.ylabel("Swelling pressure at full saturation (MPa)")
    plt.grid(True, alpha=0.35)
    plt.legend(frameon=False)
    plt.tight_layout()
    out_cmp = ROOT / "villar_dense_dd_native_notebook_calibration_comparison.png"
    plt.savefig(out_cmp, dpi=220)
    plt.close()

    plt.figure(figsize=(8.2, 5.2))
    plt.semilogy(
        x,
        y_mult,
        color="#9467bd",
        marker="^",
        linestyle="-",
        linewidth=1.8,
        markersize=4.3,
    )
    plt.xlabel("Dry density (kg/m$^3$)")
    plt.ylabel("Native notebook effective vdW multiplier (-)")
    plt.grid(True, which="both", alpha=0.35)
    plt.tight_layout()
    out_mult = ROOT / "villar_dense_dd_native_notebook_vdw_multiplier.png"
    plt.savefig(out_mult, dpi=220)
    plt.close()

    print(f"Wrote: {out_csv}")
    print(f"Wrote: {out_summary}")
    print(f"Wrote: {out_cmp}")
    print(f"Wrote: {out_mult}")


if __name__ == "__main__":
    main()
