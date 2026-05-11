#!/usr/bin/env python3
"""Calibrate the DSM_native pure vdW model against the Villar dry-density curve.

Calibrates the effective Hamaker constant A_eff = m_vdW * A_lit at each dry
density so that the model swelling pressure matches the Villar relation
P_s = exp(6.77*rho_d - 9.07) [MPa, rho_d in g/cm³].

IMPORTANT — LinearElastic pressure cap:
  This script uses LinearElasticIsotropic (E=52 MPa, ν=0.3) which gives
  K_bulk = 43.33 MPa.  The maximum achievable swelling pressure is:
    P_sw_max = K_bulk × slope × (φ₀ - n_l0_lit) ≈ 43.33 × 0.1 × 0.40 ≈ 1.7 MPa

  The Villar target exceeds this cap for ρ_d > ~1420 kg/m³. To reproduce the
  full mfront Villar calibration curve (which uses ModCamClay_semiExpl_constE
  with pressure-dependent K_bulk ≈ 2000 MPa at 9 MPa), a native MCC model is
  required.  For the BEACON inflow simulation (target 604 kPa), use the
  mfront-calibrated parameters (A_eff = 2.065e-13 J) directly.

n_l0 convention:
  Initial micro water content n_l0 is computed from A_lit (the literature
  Hamaker constant), NOT from A_eff.  This matches the mfront calibration
  workflow which fixes n_l0 from HAMAKER_LITERATURE independent of the
  calibrated multiplier.
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

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
RHO_SOLID = 2780.0          # grain density [kg/m³]
RHO_LR_REF = 1000.0         # liquid reference density [kg/m³]

MS33_K0_REF = 5.6e-21       # MS33 reference intrinsic permeability [m²]
MS33_PHI_REF = 0.42         # MS33 reference porosity [-]

HAMAKER_LITERATURE = 5.1e-21   # J  (literature Hamaker for smectite)
SPECIFIC_SURFACE = 523.0        # m²/g  (specific surface area of MX-80)

TOTAL_SUCTION_MPA = 100.0
MACRO_SUCTION_MPA = 1.0
MICRO_SUCTION_MPA = TOTAL_SUCTION_MPA - MACRO_SUCTION_MPA
PRESSURE_IC_PA = -MACRO_SUCTION_MPA * 1e6

TIME_END_S = 120 * 86400       # 120 days
PRESSURE_TOLERANCE_PA = 1e-12
MASS_EXCHANGE_COEFFICIENT = 1e-13
MICRO_SWELLING_SLOPE = 0.1

# Micro liquid density parameters — kept at O(1e-6) so n_l is effectively a
# rescaled thermodynamic variable rather than a physical volume fraction.
# Matches the native-dsm-micromacro calibration convention in ogs-dsm_mfront.
MICRO_LIQUID_DENSITY_REFERENCE = 1e-6  # rho_l0
MICRO_LIQUID_DENSITY_A = 1e-16         # a_rho
MICRO_LIQUID_DENSITY_B = 1.0           # b_rho

# Tuller saturation model constants
AREA_FACTOR_TULLER = 1.0
PORE_AREA_SHAPE_FACTOR_TULLER = 0.8584073464102069
CHARACTERISTIC_PORE_SIZE = 1e-5   # m
SURFACE_TENSION = 0.0715          # N/m

# Mechanical parameters (LinearElastic; at high p_c0 MCC ≈ LinearElastic)
YOUNG_MODULUS = 52e6   # Pa
POISSON_RATIO = 0.3

ROOT = Path(__file__).resolve().parent
DATA_ROOT = ROOT.parent

DEFAULT_OGS = Path(
    "/Users/vinaykumar/git/build/ogs-worktrees/build/dsm_native-release/bin/ogs"
)
DEFAULT_LIB = Path(
    "/Users/vinaykumar/git/build/ogs-worktrees/build/dsm_native-release/lib"
)


@dataclass(frozen=True)
class Case:
    dry_density: float  # kg/m³

    @property
    def rho_d_g_cm3(self) -> float:
        return self.dry_density / 1000.0

    @property
    def phi0(self) -> float:
        return 1.0 - self.dry_density / RHO_SOLID

    @property
    def n_s(self) -> float:
        return 1.0 - self.phi0

    @property
    def permeability(self) -> float:
        p = self.phi0
        return (
            MS33_K0_REF
            * ((1.0 - MS33_PHI_REF) ** 2 / MS33_PHI_REF**3)
            * (p**3 / (1.0 - p) ** 2)
        )

    @property
    def villar_target_mpa(self) -> float:
        return math.exp(6.77 * self.rho_d_g_cm3 - 9.07)


def n_l0_from_suction(phi0: float) -> float:
    """Initial micro water content from suction split using literature Hamaker.

    Always uses HAMAKER_LITERATURE, matching the mfront calibration workflow
    (n_l0_fixed = n_l0_from_micro_suction(phi0, HAMAKER_LITERATURE)).
    A_eff only enters the dynamics (potential formula), not the initial state.
    """
    mu_abs = MICRO_SUCTION_MPA * 1e6 / RHO_LR_REF
    n_s = 1.0 - phi0
    prefactor = (
        HAMAKER_LITERATURE
        * (SPECIFIC_SURFACE * n_s * RHO_SOLID) ** 3
        / (6.0 * math.pi)
    )
    return max(1e-12, (prefactor / mu_abs) ** (1.0 / 3.0))


def read_mean_pressure_mpa(vtu_path: Path) -> float:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(vtu_path))
    reader.Update()
    grid = reader.GetOutput()
    sigma = vtk_to_numpy(grid.GetPointData().GetArray("sigma"))
    return float((-sigma[:, 0] - sigma[:, 1] - sigma[:, 2]).mean() / 3.0e6)


def last_vtu(prefix: str) -> Path:
    candidates = sorted(ROOT.glob(f"{prefix}_ts_*_t_*.vtu"))
    if not candidates:
        raise FileNotFoundError(f"No VTU outputs for prefix {prefix!r}")
    def _key(p: Path):
        m = re.search(r"_ts_(\d+)_t_([-0-9eE+.]+)$", p.stem)
        return (int(m.group(1)), float(m.group(2))) if m else (-1, -1.0)
    return sorted(candidates, key=_key)[-1]


def cleanup(prefix: str, prj: Path) -> None:
    for pat in (f"{prefix}.pvd", f"{prefix}_ts_*_t_*.vtu"):
        for p in ROOT.glob(pat):
            p.unlink(missing_ok=True)
    prj.unlink(missing_ok=True)


def write_project(
    case: Case,
    prj_path: Path,
    hamaker_eff: float,
    n_l0: float,
) -> None:
    prefix = prj_path.stem
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
            <micro_porosity>
                <mass_exchange_coefficient>{MASS_EXCHANGE_COEFFICIENT:.16g}</mass_exchange_coefficient>
                <nonlinear_solver>
                    <maximum_iterations>100</maximum_iterations>
                    <residuum_tolerance>1e-8</residuum_tolerance>
                    <increment_tolerance>1e-20</increment_tolerance>
                </nonlinear_solver>
            </micro_porosity>
            <potential_exchange>
                <enabled>true</enabled>
                <mode>full_potential</mode>
                <pressure_tolerance>{PRESSURE_TOLERANCE_PA:.16g}</pressure_tolerance>
                <hamaker_constant>{hamaker_eff:.16g}</hamaker_constant>
                <specific_surface>{SPECIFIC_SURFACE:.16g}</specific_surface>
                <micro_solid_density_reference>{RHO_SOLID:.16g}</micro_solid_density_reference>
                <micro_solid_volume_fraction_reference>{case.n_s:.16g}</micro_solid_volume_fraction_reference>
                <micro_liquid_density_reference>{MICRO_LIQUID_DENSITY_REFERENCE:.16g}</micro_liquid_density_reference>
                <micro_liquid_density_a>{MICRO_LIQUID_DENSITY_A:.16g}</micro_liquid_density_a>
                <micro_liquid_density_b>{MICRO_LIQUID_DENSITY_B:.16g}</micro_liquid_density_b>
                <initial_micro_water_content>{n_l0:.16g}</initial_micro_water_content>
                <local_nonlinear_solve_mode>scalar_micro_macro_mass_storage_mode</local_nonlinear_solve_mode>
                <fd_jacobian_for_exchange>false</fd_jacobian_for_exchange>
                <micro_potential_convention>negative_attractive</micro_potential_convention>
                <micro_water_content_swelling_slope>{MICRO_SWELLING_SLOPE:.16g}</micro_water_content_swelling_slope>
            </potential_exchange>
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
                <secondary_variable name="micro_water_content"/>
                <secondary_variable name="micro_exchange_source"/>
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
                <property><name>biot_coefficient</name><type>Constant</type><value>1.0</value></property>
                <property>
                    <name>permeability</name>
                    <type>KozenyCarman</type>
                    <initial_permeability>k0</initial_permeability>
                    <initial_porosity>phi0</initial_porosity>
                </property>
                <property>
                    <name>porosity</name>
                    <type>PorosityFromMassBalance</type>
                    <initial_porosity>phi0</initial_porosity>
                    <minimal_porosity>0</minimal_porosity>
                    <maximal_porosity>1</maximal_porosity>
                </property>
                <property><name>reference_temperature</name><type>Constant</type><value>293.15</value></property>
                <property><name>relative_permeability</name><type>Constant</type><value>1</value></property>
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
                <property><name>bishops_effective_stress</name><type>BishopsPowerLaw</type><exponent>1</exponent></property>
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
                <variable>micro_water_content</variable>
                <variable>micro_exchange_source</variable>
            </variables>
        </output>
    </time_loop>
    <parameters>
        <parameter><name>sigma0</name><type>Function</type><expression>0</expression><expression>0</expression><expression>0</expression><expression>0</expression></parameter>
        <parameter><name>YoungModulus</name><type>Constant</type><value>{YOUNG_MODULUS:.16g}</value></parameter>
        <parameter><name>PoissonRatio</name><type>Constant</type><value>{POISSON_RATIO:.16g}</value></parameter>
        <parameter><name>phi0</name><type>Constant</type><value>{case.phi0:.16g}</value></parameter>
        <parameter><name>k0</name><type>Constant</type><value>{case.permeability:.16g}</value></parameter>
        <parameter><name>displacement0</name><type>Constant</type><values>0 0</values></parameter>
        <parameter><name>zero</name><type>Constant</type><value>0.0</value></parameter>
        <parameter><name>pressure_ic</name><type>Constant</type><value>{PRESSURE_IC_PA:.16g}</value></parameter>
        <parameter><name>p_bc_scale</name><type>Constant</type><value>1</value></parameter>
        <parameter><name>pressure_bc</name><type>CurveScaled</type><curve>pressure_release</curve><parameter>p_bc_scale</parameter></parameter>
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
        <nonlinear_solver><name>basic_newton</name><type>Newton</type><max_iter>60</max_iter><linear_solver>ls</linear_solver></nonlinear_solver>
    </nonlinear_solvers>
    <linear_solvers>
        <linear_solver><name>ls</name><eigen><solver_type>SparseLU</solver_type><scaling>true</scaling></eigen></linear_solver>
    </linear_solvers>
    <curves>
        <curve>
            <name>pressure_release</name>
            <coords>0 {TIME_END_S}</coords>
            <values>{PRESSURE_IC_PA:.16g} 0</values>
        </curve>
    </curves>
</OpenGeoSysProject>
"""
    prj_path.write_text(xml)


def run_ogs(ogs_bin: Path, lib_path: Path, prj: Path) -> None:
    env = os.environ.copy()
    if lib_path.exists():
        env["DYLD_LIBRARY_PATH"] = str(lib_path)
    subprocess.run(
        [str(ogs_bin), str(prj)],
        cwd=ROOT,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
        env=env,
    )


def run_trial(
    ogs_bin: Path, lib_path: Path, case: Case, multiplier: float, idx: int
) -> dict:
    tag = f"_anchors_purevdw_dd{int(case.dry_density)}_{idx:03d}"
    prj = ROOT / f"{tag}.prj"
    hamaker_eff = HAMAKER_LITERATURE * multiplier
    n_l0 = n_l0_from_suction(case.phi0)
    write_project(case, prj, hamaker_eff, n_l0)
    try:
        run_ogs(ogs_bin, lib_path, prj)
        p_mpa = read_mean_pressure_mpa(last_vtu(tag))
    finally:
        cleanup(tag, prj)
    return {"multiplier": multiplier, "hamaker_eff": hamaker_eff, "n_l0": n_l0, "pressure_mpa": p_mpa}


def calibrate(
    ogs_bin: Path,
    lib_path: Path,
    case: Case,
    target: float,
    rel_tol: float = 0.02,
    max_iter: int = 18,
) -> dict:
    def rel_err(p: float) -> float:
        return abs(p - target) / max(target, 1e-12)

    run_id = 0
    r_lo = run_trial(ogs_bin, lib_path, case, 1.0, run_id)
    run_id += 1
    best = r_lo
    if rel_err(r_lo["pressure_mpa"]) < rel_tol:
        return r_lo

    m_lo = 1.0
    m_hi = max(2.0, target / max(abs(r_lo["pressure_mpa"]), 1e-12))
    m_hi = float(np.clip(m_hi, 1e-6, 1e18))
    if m_hi <= m_lo:
        m_hi = 2.0

    r_hi = run_trial(ogs_bin, lib_path, case, m_hi, run_id)
    run_id += 1
    if rel_err(r_hi["pressure_mpa"]) < rel_err(best["pressure_mpa"]):
        best = r_hi

    expand_limit = max_iter + 16
    while r_hi["pressure_mpa"] < target and m_hi < 1e18 and run_id < expand_limit:
        m_lo, r_lo = m_hi, r_hi
        m_hi = min(1e18, m_hi * 10.0)
        r_hi = run_trial(ogs_bin, lib_path, case, m_hi, run_id)
        run_id += 1
        if rel_err(r_hi["pressure_mpa"]) < rel_err(best["pressure_mpa"]):
            best = r_hi

    if r_hi["pressure_mpa"] < target:
        return best

    for _ in range(max_iter):
        m_mid = math.sqrt(m_lo * m_hi)
        if m_mid <= m_lo * (1 + 1e-12) or m_mid >= m_hi * (1 - 1e-12):
            break
        r_mid = run_trial(ogs_bin, lib_path, case, m_mid, run_id)
        run_id += 1
        if rel_err(r_mid["pressure_mpa"]) < rel_err(best["pressure_mpa"]):
            best = r_mid
        if rel_err(r_mid["pressure_mpa"]) < rel_tol:
            return r_mid
        if r_mid["pressure_mpa"] < target:
            m_lo, r_lo = m_mid, r_mid
        else:
            m_hi, r_hi = m_mid, r_mid

    return best


def main() -> None:
    ap = argparse.ArgumentParser(description="Pure vdW Villar dense dd calibration")
    ap.add_argument("--ogs-bin", type=Path, default=DEFAULT_OGS)
    ap.add_argument("--lib-path", type=Path, default=DEFAULT_LIB)
    ap.add_argument("--dd-min", type=float, default=1400.0)
    ap.add_argument("--dd-max", type=float, default=1800.0)
    ap.add_argument("--dd-step", type=float, default=25.0)
    ap.add_argument("--rel-tol", type=float, default=0.02)
    args = ap.parse_args()

    lib_path = args.lib_path
    dd_vals = np.arange(args.dd_min, args.dd_max + 0.5 * args.dd_step, args.dd_step)
    cases = [Case(float(d)) for d in dd_vals]

    rows = []
    for case in cases:
        target = case.villar_target_mpa
        result = calibrate(args.ogs_bin, lib_path, case, target, args.rel_tol)
        print(
            f"dd={case.dry_density:.0f} kg/m³: "
            f"target={target:.4f} MPa  "
            f"calibrated={result['pressure_mpa']:.4f} MPa  "
            f"m_vdW={result['multiplier']:.4e}  "
            f"n_l0={result['n_l0']:.4e}"
        )
        rows.append({
            "dry_density_kg_m3": case.dry_density,
            "dry_density_g_cm3": case.rho_d_g_cm3,
            "phi0": case.phi0,
            "n_s_ref": case.n_s,
            "k0_m2": case.permeability,
            "target_villar_MPa": target,
            "calibrated_MPa": result["pressure_mpa"],
            "vdw_multiplier": result["multiplier"],
            "hamaker_literature_J": HAMAKER_LITERATURE,
            "hamaker_effective_J": result["hamaker_eff"],
            "specific_surface": SPECIFIC_SURFACE,
            "n_l0": result["n_l0"],
            "delta_MPa": result["pressure_mpa"] - target,
        })

    out_csv = ROOT / "villar_dense_dd_native_purevdw_calibration.csv"
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "hamaker_literature_J": HAMAKER_LITERATURE,
        "specific_surface": SPECIFIC_SURFACE,
        "n_density_points": len(rows),
        "mean_rel_error_percent": float(
            100 * np.mean([abs(r["delta_MPa"]) / max(r["target_villar_MPa"], 1e-12) for r in rows])
        ),
        "max_rel_error_percent": float(
            100 * np.max([abs(r["delta_MPa"]) / max(r["target_villar_MPa"], 1e-12) for r in rows])
        ),
        "results": rows,
    }
    (ROOT / "villar_dense_dd_native_purevdw_calibration_summary.json").write_text(
        json.dumps(summary, indent=2)
    )

    x = np.array([r["dry_density_kg_m3"] for r in rows])
    y_tgt = np.array([r["target_villar_MPa"] for r in rows])
    y_cal = np.array([r["calibrated_MPa"] for r in rows])
    y_mult = np.array([r["vdw_multiplier"] for r in rows])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x, y_tgt, "k-", lw=2, label="Villar target")
    ax.plot(x, y_cal, "o--", color="#1f77b4", ms=4, lw=1.6, label="Native pure vdW calibrated")
    ax.set_xlabel("Dry density (kg/m³)")
    ax.set_ylabel("Swelling pressure (MPa)")
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(ROOT / "villar_dense_dd_native_purevdw_swelling_pressure.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(x, y_mult, "^-", color="#9467bd", ms=4, lw=1.8)
    ax.set_xlabel("Dry density (kg/m³)")
    ax.set_ylabel("vdW multiplier m_vdW (-)")
    ax.grid(which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(ROOT / "villar_dense_dd_native_purevdw_multiplier.png", dpi=220)
    plt.close(fig)

    print(f"\nWrote: {out_csv}")
    print(f"Wrote: {ROOT / 'villar_dense_dd_native_purevdw_calibration_summary.json'}")


if __name__ == "__main__":
    main()
