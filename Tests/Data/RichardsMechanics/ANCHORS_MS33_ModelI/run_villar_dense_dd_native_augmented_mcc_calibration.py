#!/usr/bin/env python3
"""Calibrate the DSM_native augmented vdW + ModCamClay model against the Villar curve.

Same strategy as run_villar_dense_dd_native_augmented_calibration.py but uses
MFront/ModCamClay_semiExpl_constE for the skeleton constitutive law.

Motivation
----------
LinearElasticIsotropic caps the achievable swelling pressure at
  P_sw_max = K_bulk × slope × φ₀ ≈ 43 MPa × 0.1 × 0.5 ≈ 2 MPa
regardless of the augmentation prefactor K, because the eigenstrain
Δε_sw = slope × Δn_l ≤ slope × φ₀ and K_bulk is constant.

With ModCamClay_semiExpl_constE and low p_c0 (here 0.2 MPa), the confined
swelling test engages MCC plastic hardening: once the effective mean stress
exceeds p_c0, plastic compression (ε_p_vol > 0) develops under the yield cap.
Because displacement is fixed (ε_vol = 0), this plastic compression adds to
the swelling eigenstrain in building compressive effective stress:
  p' = K_bulk × (ε_sw_vol + ε_p_vol)
The preconsolidation pressure p_c hardens and the stress state climbs along
the normal consolidation line, reaching pressures well above the elastic cap.

MCC parameters (fixed — from mfront_restart_part1.prj / MX-80 bentonite):
  E   = 52 MPa   ν = 0.3   M = 1.2   κ = 6.6e-3   λ = 7.7e-2   p_c0 = 0.2 MPa

The initial Bishop effective stress is ~0 at start (large suction → S_l → 0
→ χ_eff = S_l → 0), so p_c0 = 0.2 MPa puts the initial state safely inside
the yield surface.  VolumeRatio v0 = ρ_solid / ρ_d varies with dry density.

OMP parallelism: OMP_NUM_THREADS=12, OGS_ASM_THREADS=12.
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
# Physical constants — same as augmented calibration
# ---------------------------------------------------------------------------
RHO_SOLID = 2780.0
RHO_LR_REF = 1000.0
MS33_K0_REF = 5.6e-21
MS33_PHI_REF = 0.42
HAMAKER_LITERATURE = 5.1e-21
SPECIFIC_SURFACE = 523.0
TOTAL_SUCTION_MPA = 100.0
MACRO_SUCTION_MPA = 1.0
MICRO_SUCTION_MPA = TOTAL_SUCTION_MPA - MACRO_SUCTION_MPA
PRESSURE_IC_PA = -MACRO_SUCTION_MPA * 1e6
TIME_END_S = 120 * 86400
PRESSURE_TOLERANCE_PA = 1e-12
MASS_EXCHANGE_COEFFICIENT = 1e-13
MICRO_SWELLING_SLOPE = 0.1

MICRO_LIQUID_DENSITY_REFERENCE = 1e-6
MICRO_LIQUID_DENSITY_A = 1e-16
MICRO_LIQUID_DENSITY_B = 1.0
AREA_FACTOR_TULLER = 1.0
PORE_AREA_SHAPE_FACTOR_TULLER = 0.8584073464102069
CHARACTERISTIC_PORE_SIZE = 1e-5
SURFACE_TENSION = 0.0715

# ---------------------------------------------------------------------------
# ModCamClay parameters (fixed, from mfront_restart_part1.prj)
# ---------------------------------------------------------------------------
YOUNG_MODULUS = 52e6          # Pa — constant E for MCC_semiExpl_constE
POISSON_RATIO = 0.3
MCC_M = 1.2                   # CriticalStateLineSlope
MCC_KAPPA = 6.6e-3            # SwellingLineSlope
MCC_LAMBDA = 7.7e-2           # VirginConsolidationLineSlope
MCC_PC0 = 2e5                 # CharacteristicPreConsolidationPressure [Pa]

# OMP parallelism
OMP_THREADS = 12

ROOT = Path(__file__).resolve().parent
DEFAULT_OGS = Path(
    "/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs"
)
DEFAULT_LIB = Path(
    "/Users/vinaykumar/git/build/release-omp-mfront/lib"
)


@dataclass(frozen=True)
class Case:
    dry_density: float

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
    def volume_ratio_v0(self) -> float:
        """Specific volume v0 = 1 + e0 = rho_solid / rho_d."""
        return RHO_SOLID / self.dry_density

    @property
    def permeability(self) -> float:
        p = self.phi0
        return (
            MS33_K0_REF
            * ((1 - MS33_PHI_REF) ** 2 / MS33_PHI_REF**3)
            * (p**3 / (1 - p) ** 2)
        )

    @property
    def villar_target_mpa(self) -> float:
        return math.exp(6.77 * self.rho_d_g_cm3 - 9.07)


def n_l0_from_suction(phi0: float) -> float:
    mu_abs = MICRO_SUCTION_MPA * 1e6 / RHO_LR_REF
    n_s = 1.0 - phi0
    prefactor = HAMAKER_LITERATURE * (SPECIFIC_SURFACE * n_s * RHO_SOLID) ** 3 / (6 * math.pi)
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
        raise FileNotFoundError(f"No VTU outputs for {prefix!r}")

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
    K: float,
    lam: float,
    n_l0: float,
) -> None:
    prefix = prj_path.stem
    v0 = case.volume_ratio_v0
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
                <pressure_tolerance>{PRESSURE_TOLERANCE_PA:.16g}</pressure_tolerance>
                <hamaker_constant>{HAMAKER_LITERATURE:.16g}</hamaker_constant>
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
                <vdw_augmentation_prefactor>{K:.16g}</vdw_augmentation_prefactor>
                <vdw_augmentation_decay_length>{lam:.16g}</vdw_augmentation_decay_length>
            </potential_exchange>
            <constitutive_relation>
                <type>MFront</type>
                <behaviour>ModCamClay_semiExpl_constE</behaviour>
                <material_properties>
                    <material_property name="YoungModulus" parameter="YoungModulus"/>
                    <material_property name="PoissonRatio" parameter="PoissonRatio"/>
                    <material_property name="CriticalStateLineSlope" parameter="CriticalStateLineSlope"/>
                    <material_property name="SwellingLineSlope" parameter="SwellingLineSlope"/>
                    <material_property name="VirginConsolidationLineSlope" parameter="VirginConsolidationLineSlope"/>
                    <material_property name="CharacteristicPreConsolidationPressure" parameter="InitialPreConsolidationPressure"/>
                </material_properties>
                <initial_values>
                    <state_variable name="PreConsolidationPressure" parameter="InitialPreConsolidationPressure"/>
                    <state_variable name="VolumeRatio" parameter="InitialVolumeRatio"/>
                </initial_values>
            </constitutive_relation>
            <process_variables>
                <pressure>pressure</pressure>
                <displacement>displacement</displacement>
            </process_variables>
            <secondary_variables>
                <secondary_variable name="ElasticStrain"/>
                <secondary_variable name="EquivalentPlasticStrain"/>
                <secondary_variable name="PreConsolidationPressure"/>
                <secondary_variable name="PlasticVolumetricStrain"/>
                <secondary_variable name="VolumeRatio"/>
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
                <variable>ElasticStrain</variable>
                <variable>PreConsolidationPressure</variable>
                <variable>PlasticVolumetricStrain</variable>
                <variable>VolumeRatio</variable>
            </variables>
        </output>
    </time_loop>
    <parameters>
        <parameter><name>sigma0</name><type>Function</type><expression>0</expression><expression>0</expression><expression>0</expression><expression>0</expression></parameter>
        <parameter><name>YoungModulus</name><type>Constant</type><value>{YOUNG_MODULUS:.16g}</value></parameter>
        <parameter><name>PoissonRatio</name><type>Constant</type><value>{POISSON_RATIO:.16g}</value></parameter>
        <parameter><name>CriticalStateLineSlope</name><type>Constant</type><value>{MCC_M:.16g}</value></parameter>
        <parameter><name>SwellingLineSlope</name><type>Constant</type><value>{MCC_KAPPA:.16g}</value></parameter>
        <parameter><name>VirginConsolidationLineSlope</name><type>Constant</type><value>{MCC_LAMBDA:.16g}</value></parameter>
        <parameter><name>InitialPreConsolidationPressure</name><type>Constant</type><value>{MCC_PC0:.16g}</value></parameter>
        <parameter><name>InitialVolumeRatio</name><type>Constant</type><value>{v0:.16g}</value></parameter>
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


def run_ogs(ogs_bin: Path, lib_path: Path, prj: Path) -> bool:
    """Run OGS. Returns True on success, False if OGS exits non-zero."""
    env = os.environ.copy()
    if lib_path.exists():
        env["DYLD_LIBRARY_PATH"] = str(lib_path)
    env["OMP_NUM_THREADS"] = str(OMP_THREADS)
    env["OGS_ASM_THREADS"] = str(OMP_THREADS)
    result = subprocess.run(
        [str(ogs_bin), str(prj)],
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
        env=env,
    )
    return result.returncode == 0


def run_trial(
    ogs_bin: Path, lib_path: Path, case: Case, K: float, lam: float, n_l0: float, idx: int
) -> dict:
    """Returns dict with keys K, pressure_mpa, converged."""
    tag = f"_anchors_mcc_dd{int(case.dry_density)}_{idx:03d}"
    prj = ROOT / f"{tag}.prj"
    write_project(case, prj, K, lam, n_l0)
    try:
        ok = run_ogs(ogs_bin, lib_path, prj)
        if ok:
            p_mpa = read_mean_pressure_mpa(last_vtu(tag))
        else:
            p_mpa = 0.0
    finally:
        cleanup(tag, prj)
    return {"K": K, "pressure_mpa": p_mpa, "converged": ok}


def calibrate_K(
    ogs_bin: Path,
    lib_path: Path,
    case: Case,
    lam: float,
    n_l0: float,
    target: float,
    rel_tol: float = 0.02,
    max_iter: int = 20,
) -> dict:
    def rel_err(p: float) -> float:
        return abs(p - target) / max(target, 1e-12)

    run_id = 0
    r_lo = run_trial(ogs_bin, lib_path, case, 0.0, lam, n_l0, run_id)
    run_id += 1
    best = r_lo
    if rel_err(r_lo["pressure_mpa"]) < rel_tol:
        return r_lo

    deficit = max(target - r_lo["pressure_mpa"], 1e-12)
    K_hi = deficit * RHO_LR_REF * 10.0

    # Find a K_hi that both converges and brackets the target from above.
    # If OGS diverges (MCC integrator overloaded), halve K_hi until it converges.
    MAX_SHRINK = 30
    for _ in range(MAX_SHRINK):
        r_hi = run_trial(ogs_bin, lib_path, case, K_hi, lam, n_l0, run_id)
        run_id += 1
        if r_hi["converged"]:
            break
        K_hi *= 0.5
    else:
        return best  # no converging K_hi found

    if rel_err(r_hi["pressure_mpa"]) < rel_err(best["pressure_mpa"]):
        best = r_hi

    expand_limit = max_iter + 12
    K_lo = 0.0
    while r_hi["pressure_mpa"] < target and K_hi < 1e30 and run_id < expand_limit:
        K_lo, r_lo = K_hi, r_hi
        K_hi = K_hi * 5.0
        r_hi = run_trial(ogs_bin, lib_path, case, K_hi, lam, n_l0, run_id)
        run_id += 1
        # If expansion diverges, treat it as "K too large" — stop expanding
        if not r_hi["converged"]:
            break
        if rel_err(r_hi["pressure_mpa"]) < rel_err(best["pressure_mpa"]):
            best = r_hi

    if r_hi["pressure_mpa"] < target:
        return best

    for _ in range(max_iter):
        K_mid = 0.5 * (K_lo + K_hi)
        if K_mid <= K_lo * (1 + 1e-12) or abs(K_hi - K_lo) < 1e-30:
            break
        r_mid = run_trial(ogs_bin, lib_path, case, K_mid, lam, n_l0, run_id)
        run_id += 1
        # Diverged mid-point: shrink upper bracket
        if not r_mid["converged"]:
            K_hi = K_mid
            continue
        if rel_err(r_mid["pressure_mpa"]) < rel_err(best["pressure_mpa"]):
            best = r_mid
        if rel_err(r_mid["pressure_mpa"]) < rel_tol:
            return r_mid
        if r_mid["pressure_mpa"] < target:
            K_lo, r_lo = K_mid, r_mid
        else:
            K_hi, r_hi = K_mid, r_mid

    return best


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Augmented vdW + ModCamClay Villar dense dd calibration"
    )
    ap.add_argument("--ogs-bin", type=Path, default=DEFAULT_OGS)
    ap.add_argument("--lib-path", type=Path, default=DEFAULT_LIB)
    ap.add_argument("--dd-min", type=float, default=1400.0)
    ap.add_argument("--dd-max", type=float, default=1800.0)
    ap.add_argument("--dd-step", type=float, default=25.0)
    ap.add_argument(
        "--lam",
        type=float,
        default=1e-6,
        help="Decay length lambda [code units]. Default 1e-6 "
             "(= physical 1 nm; Sa uses m²/g convention so "
             "lambda_code = lambda_physical × 1000).",
    )
    ap.add_argument("--rel-tol", type=float, default=0.02)
    args = ap.parse_args()

    lam = args.lam
    lib_path = args.lib_path
    dd_vals = np.arange(args.dd_min, args.dd_max + 0.5 * args.dd_step, args.dd_step)
    cases = [Case(float(d)) for d in dd_vals]

    print(
        f"MCC parameters: E={YOUNG_MODULUS/1e6:.0f} MPa  ν={POISSON_RATIO}  "
        f"M={MCC_M}  κ={MCC_KAPPA}  λ={MCC_LAMBDA}  p_c0={MCC_PC0/1e6:.1f} MPa"
    )
    print(f"OMP_NUM_THREADS={OMP_THREADS}  OGS_ASM_THREADS={OMP_THREADS}")
    print()

    rows = []
    for case in cases:
        n_l0 = n_l0_from_suction(case.phi0)
        target = case.villar_target_mpa
        xi0 = n_l0 / (lam * case.n_s * RHO_SOLID * SPECIFIC_SURFACE)
        result = calibrate_K(args.ogs_bin, lib_path, case, lam, n_l0, target, args.rel_tol)
        converged = result.get("converged", True)
        print(
            f"dd={case.dry_density:.0f} kg/m³: "
            f"target={target:.4f} MPa  "
            f"calibrated={result['pressure_mpa']:.4f} MPa  "
            f"K={result['K']:.4e} J/kg  "
            f"v0={case.volume_ratio_v0:.4f}"
            + ("" if converged else "  [diverged]")
        )
        rows.append(
            {
                "dry_density_kg_m3": case.dry_density,
                "dry_density_g_cm3": case.rho_d_g_cm3,
                "phi0": case.phi0,
                "n_s_ref": case.n_s,
                "volume_ratio_v0": case.volume_ratio_v0,
                "target_villar_MPa": target,
                "calibrated_MPa": result["pressure_mpa"],
                "K_calibrated_J_kg": result["K"],
                "converged": converged,
                "lambda_m": lam,
                "xi0": xi0,
                "n_l0": n_l0,
                "hamaker_literature_J": HAMAKER_LITERATURE,
                "specific_surface": SPECIFIC_SURFACE,
                "mcc_M": MCC_M,
                "mcc_kappa": MCC_KAPPA,
                "mcc_lambda": MCC_LAMBDA,
                "mcc_pc0_Pa": MCC_PC0,
                "delta_MPa": result["pressure_mpa"] - target,
            }
        )

    tag = f"mcc_lam{lam:.0e}".replace("-", "n").replace("+", "")
    out_csv = ROOT / f"villar_dense_dd_native_augmented_{tag}_calibration.csv"
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "constitutive_model": "ModCamClay_semiExpl_constE",
        "mcc_M": MCC_M,
        "mcc_kappa": MCC_KAPPA,
        "mcc_lambda": MCC_LAMBDA,
        "mcc_pc0_Pa": MCC_PC0,
        "lambda_m": lam,
        "hamaker_literature_J": HAMAKER_LITERATURE,
        "specific_surface": SPECIFIC_SURFACE,
        "n_density_points": len(rows),
        "mean_rel_error_percent": float(
            100
            * np.mean(
                [abs(r["delta_MPa"]) / max(r["target_villar_MPa"], 1e-12) for r in rows]
            )
        ),
        "max_rel_error_percent": float(
            100
            * np.max(
                [abs(r["delta_MPa"]) / max(r["target_villar_MPa"], 1e-12) for r in rows]
            )
        ),
        "results": rows,
    }
    (ROOT / f"villar_dense_dd_native_augmented_{tag}_summary.json").write_text(
        json.dumps(summary, indent=2)
    )

    x = np.array([r["dry_density_kg_m3"] for r in rows])
    y_tgt = np.array([r["target_villar_MPa"] for r in rows])
    y_cal = np.array([r["calibrated_MPa"] for r in rows])
    y_K = np.array([r["K_calibrated_J_kg"] for r in rows])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x, y_tgt, "k-", lw=2, label="Villar target")
    ax.plot(
        x,
        y_cal,
        "s--",
        color="#d62728",
        ms=4,
        lw=1.6,
        label=f"MCC calibrated (λ={lam:.0e})",
    )
    ax.set_xlabel("Dry density (kg/m³)")
    ax.set_ylabel("Swelling pressure (MPa)")
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(
        ROOT / f"villar_dense_dd_native_augmented_{tag}_swelling_pressure.png", dpi=220
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(x, y_K + 1.0, "v-", color="#2ca02c", ms=4, lw=1.8)
    ax.set_xlabel("Dry density (kg/m³)")
    ax.set_ylabel("Augmentation prefactor K (J/kg)")
    ax.set_title(f"MCC  λ = {lam:.0e}")
    ax.grid(which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(
        ROOT / f"villar_dense_dd_native_augmented_{tag}_K_curve.png", dpi=220
    )
    plt.close(fig)

    print(f"\nWrote: {out_csv}")
    print(f"mean_rel_error = {summary['mean_rel_error_percent']:.2f} %")
    print(f"max_rel_error  = {summary['max_rel_error_percent']:.2f} %")


if __name__ == "__main__":
    main()
