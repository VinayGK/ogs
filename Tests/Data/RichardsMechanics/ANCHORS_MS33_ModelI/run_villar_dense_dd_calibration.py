#!/usr/bin/env python3
"""Calibrate the MFront dense dry-density sweep against the Villar curve.

The workflow runs one OGS case per dry density, solves for an effective vdW
multiplier, and writes the calibration table, JSON summary, and plots.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent

RHO_SOLID = 2780.0  # kg/m^3
RHO_SOLID_MFRONT = 2780.0  # kg/m^3
RHO_LR_REF = 1000.0  # kg/m^3

MS33_K0_REF = 5.6e-21
MS33_PHI_REF = 0.42
MS33_PB = 27e6
MS33_M = 0.45

HAMAKER_LITERATURE = 5.1e-21  # J, montmorillonite-water-montmorillonite
SPECIFIC_SURFACE_MASS = 523.0  # m^2/g

TOTAL_SUCTION_MPA = 100.0
MACRO_SUCTION_MPA = 1.0
MICRO_SUCTION_MPA = TOTAL_SUCTION_MPA - MACRO_SUCTION_MPA

NOTEBOOK_SWELLING_SLOPE = 0.1
MASS_EXCHANGE_COEFFICIENT = 1e-13
LIQUID_VISCOSITY = 1e-3
NOTEBOOK_SATURATION_MODE = 1
# 0 => coupled notebook mass-storage branch in RichardsMechanicsNotebookBridge_MCC
# (aligns with native scalar_notebook_mass_storage solve intent).
NOTEBOOK_LOCAL_SOLVE_MODE = 0
MICRO_POTENTIAL_CONVENTION = 1

PORE_AREA_SHAPE_FACTOR_TULLER = 0.8584073464102069
AREA_FACTOR_TULLER = 1.0
CHARACTERISTIC_PORE_SIZE = 1e-5
SURFACE_TENSION = 0.0715

DENSITY_A = 0.0
DENSITY_B = 1.0
RHO_L0 = 1e-6
RHO_LR0 = 1000.000001

PRESSURE_IC_PA = -MACRO_SUCTION_MPA * 1e6
TIME_END_S = 120 * 86400


@dataclass(frozen=True)
class Case:
    dry_density: float  # kg/m^3

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
    def initial_volume_ratio(self) -> float:
        return 1.0 / (1.0 - self.phi0)

    @property
    def intrinsic_permeability(self) -> float:
        return (
            MS33_K0_REF
            * ((1.0 - MS33_PHI_REF) ** 2 / MS33_PHI_REF**3)
            * (self.phi0**3 / (1.0 - self.phi0) ** 2)
        )

    @property
    def specific_surface_volumetric(self) -> float:
        # Convert m^2/g -> m^2/m^3 using dry density (kg/m^3 = 1000 g/L).
        return SPECIFIC_SURFACE_MASS * self.dry_density * 1000.0

    @property
    def villar_target_swelling_mpa(self) -> float:
        # Lloret/Villar (2007), Eq. (7): Ps [MPa], rho_d [g/cm^3].
        return math.exp(6.77 * self.dry_density_g_cm3 - 9.07)


def read_grid(path: Path) -> vtk.vtkUnstructuredGrid:
    """Read a VTU snapshot into a VTK unstructured grid."""
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    return reader.GetOutput()


def get_array(grid: vtk.vtkUnstructuredGrid, name: str):
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
    # Output name pattern: <prefix>_ts_<i>_t_<time>.vtu
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


def cleanup_runtime(prefix: str) -> None:
    """Remove transient outputs from a calibration trial."""
    for pattern in (f"{prefix}.pvd", f"{prefix}_ts_*_t_*.vtu"):
        for p in ROOT.glob(pattern):
            p.unlink(missing_ok=True)


def n_l0_from_micro_suction(phi0: float, hamaker_eff: float) -> float:
    # micro suction -> chemical potential target: mu = p / rho
    mu_abs = MICRO_SUCTION_MPA * 1e6 / RHO_LR_REF
    n_s = 1.0 - phi0
    prefactor = (
        abs(hamaker_eff)
        * (SPECIFIC_SURFACE_MASS * n_s * RHO_SOLID_MFRONT) ** 3
        / (6.0 * math.pi)
    )
    n_l0 = (prefactor / mu_abs) ** (1.0 / 3.0)
    return max(1e-12, n_l0)


def write_mfront_project(
    case: Case, multiplier: float, n_l0_fixed: float, project_path: Path
) -> dict:
    """Create a temporary MFront project for one multiplier trial."""
    hamaker_eff = HAMAKER_LITERATURE * multiplier
    prefix = project_path.stem

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
            <constitutive_relation>
                <type>MFrontRichardsMechanics</type>
                <behaviour>RichardsMechanicsNotebookBridge_MCC</behaviour>
                <material_properties>
                    <material_property name="YoungModulus" parameter="YoungModulus"/>
                    <material_property name="PoissonRatio" parameter="PoissonRatio"/>
                    <material_property name="CriticalStateLineSlope" parameter="CriticalStateLineSlope"/>
                    <material_property name="SwellingLineSlope" parameter="SwellingLineSlope"/>
                    <material_property name="VirginConsolidationLineSlope" parameter="VirginConsolidationLineSlope"/>
                    <material_property name="CharacteristicPreConsolidationPressure" parameter="InitialPreConsolidationPressure"/>
                    <material_property name="ResidualLiquidSaturation" parameter="ResidualLiquidSaturation"/>
                    <material_property name="ResidualGasSaturation" parameter="ResidualGasSaturation"/>
                    <material_property name="BubblePressure" parameter="BubblePressure"/>
                    <material_property name="VanGenuchtenExponent_m" parameter="VanGenuchtenExponent_m"/>
                    <material_property name="NotebookSaturationMode" parameter="NotebookSaturationMode"/>
                    <material_property name="NotebookLocalSolveMode" parameter="NotebookLocalSolveMode"/>
                    <material_property name="MicroPotentialConvention" parameter="MicroPotentialConvention"/>
                    <material_property name="SwellingSlope" parameter="NotebookSwellingSlope"/>
                    <material_property name="MassExchangeCoefficient" parameter="MassExchangeCoefficient"/>
                    <material_property name="MacroViscosity" parameter="MacroViscosity"/>
                    <material_property name="ReferenceLiquidDensityMacro" parameter="ReferenceLiquidDensityMacro"/>
                    <material_property name="ReferenceLiquidDensityMicro" parameter="ReferenceLiquidDensityMicro"/>
                    <material_property name="ReferenceDensitySolid" parameter="ReferenceDensitySolid"/>
                    <material_property name="MicroLiquidDensityA" parameter="MicroLiquidDensityA"/>
                    <material_property name="MicroLiquidDensityB" parameter="MicroLiquidDensityB"/>
                    <material_property name="HamakerConstant" parameter="HamakerConstant"/>
                    <material_property name="SpecificSurface" parameter="SpecificSurface"/>
                    <material_property name="AreaFactorTuller" parameter="AreaFactorTuller"/>
                    <material_property name="PoreAreaShapeFactorTuller" parameter="PoreAreaShapeFactorTuller"/>
                    <material_property name="CharacteristicPoreSize" parameter="CharacteristicPoreSize"/>
                    <material_property name="SurfaceTension" parameter="SurfaceTension"/>
                    <material_property name="InitialPorosity" parameter="phi0"/>
                </material_properties>
                <initial_values>
                    <state_variable name="PreConsolidationPressure" parameter="InitialPreConsolidationPressure"/>
                    <state_variable name="VolumeRatio" parameter="InitialVolumeRatio"/>
                    <state_variable name="n_l" parameter="n_l0"/>
                    <state_variable name="rho_lR" parameter="rho_lR0"/>
                    <state_variable name="epsilon_sw" parameter="epsilon_sw0"/>
                </initial_values>
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
                <secondary_variable name="n_l"/>
                <secondary_variable name="phi_m"/>
                <secondary_variable name="phi_M"/>
                <secondary_variable name="mu_lR"/>
            </secondary_variables>
            <specific_body_force>0 0</specific_body_force>
            <initial_stress>sigma0</initial_stress>
            <explicit_hm_coupling_in_unsaturated_zone>false</explicit_hm_coupling_in_unsaturated_zone>
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
                    <name>saturation</name>
                    <type>SaturationTuller</type>
                    <area_factor_tuller>{AREA_FACTOR_TULLER}</area_factor_tuller>
                    <pore_area_shape_factor_tuller>{PORE_AREA_SHAPE_FACTOR_TULLER}</pore_area_shape_factor_tuller>
                    <characteristic_pore_size>{CHARACTERISTIC_PORE_SIZE}</characteristic_pore_size>
                    <surface_tension>{SURFACE_TENSION}</surface_tension>
                </property>
                <property>
                    <name>relative_permeability</name>
                    <type>Constant</type>
                    <value>1</value>
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
                    <reltols>1e-10 1e-10 1e-9</reltols>
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
                <variable>n_l</variable>
                <variable>phi_m</variable>
                <variable>phi_M</variable>
                <variable>mu_lR</variable>
            </variables>
        </output>
    </time_loop>
    <parameters>
        <parameter><name>sigma0</name><type>Function</type><expression>0</expression><expression>0</expression><expression>0</expression><expression>0</expression></parameter>
        <parameter><name>YoungModulus</name><type>Constant</type><value>52e6</value></parameter>
        <parameter><name>PoissonRatio</name><type>Constant</type><value>0.3</value></parameter>
        <parameter><name>CriticalStateLineSlope</name><type>Constant</type><value>1.2</value></parameter>
        <parameter><name>SwellingLineSlope</name><type>Constant</type><value>6.6e-3</value></parameter>
        <parameter><name>VirginConsolidationLineSlope</name><type>Constant</type><value>7.7e-2</value></parameter>
        <parameter><name>InitialPreConsolidationPressure</name><type>Constant</type><value>1e10</value></parameter>
        <parameter><name>InitialVolumeRatio</name><type>Constant</type><value>{case.initial_volume_ratio:.16g}</value></parameter>
        <parameter><name>ResidualLiquidSaturation</name><type>Constant</type><value>0.0</value></parameter>
        <parameter><name>ResidualGasSaturation</name><type>Constant</type><value>0.0</value></parameter>
        <parameter><name>BubblePressure</name><type>Constant</type><value>{MS33_PB:.16g}</value></parameter>
        <parameter><name>VanGenuchtenExponent_m</name><type>Constant</type><value>{MS33_M:.16g}</value></parameter>
        <parameter><name>NotebookSaturationMode</name><type>Constant</type><value>{NOTEBOOK_SATURATION_MODE}</value></parameter>
        <parameter><name>NotebookLocalSolveMode</name><type>Constant</type><value>{NOTEBOOK_LOCAL_SOLVE_MODE}</value></parameter>
        <parameter><name>MicroPotentialConvention</name><type>Constant</type><value>{MICRO_POTENTIAL_CONVENTION}</value></parameter>
        <parameter><name>NotebookSwellingSlope</name><type>Constant</type><value>{NOTEBOOK_SWELLING_SLOPE:.16g}</value></parameter>
        <parameter><name>MassExchangeCoefficient</name><type>Constant</type><value>{MASS_EXCHANGE_COEFFICIENT:.16g}</value></parameter>
        <parameter><name>MacroViscosity</name><type>Constant</type><value>{LIQUID_VISCOSITY:.16g}</value></parameter>
        <parameter><name>ReferenceLiquidDensityMacro</name><type>Constant</type><value>{RHO_LR_REF:.16g}</value></parameter>
        <parameter><name>ReferenceLiquidDensityMicro</name><type>Constant</type><value>{RHO_L0:.16g}</value></parameter>
        <parameter><name>ReferenceDensitySolid</name><type>Constant</type><value>{RHO_SOLID_MFRONT:.16g}</value></parameter>
        <parameter><name>MicroLiquidDensityA</name><type>Constant</type><value>{DENSITY_A:.16g}</value></parameter>
        <parameter><name>MicroLiquidDensityB</name><type>Constant</type><value>{DENSITY_B:.16g}</value></parameter>
        <parameter><name>HamakerConstant</name><type>Constant</type><value>{hamaker_eff:.16g}</value></parameter>
        <parameter><name>SpecificSurface</name><type>Constant</type><value>{SPECIFIC_SURFACE_MASS:.16g}</value></parameter>
        <parameter><name>AreaFactorTuller</name><type>Constant</type><value>{AREA_FACTOR_TULLER:.16g}</value></parameter>
        <parameter><name>PoreAreaShapeFactorTuller</name><type>Constant</type><value>{PORE_AREA_SHAPE_FACTOR_TULLER:.16g}</value></parameter>
        <parameter><name>CharacteristicPoreSize</name><type>Constant</type><value>{CHARACTERISTIC_PORE_SIZE:.16g}</value></parameter>
        <parameter><name>SurfaceTension</name><type>Constant</type><value>{SURFACE_TENSION:.16g}</value></parameter>
        <parameter><name>phi0</name><type>Constant</type><value>{case.phi0:.16g}</value></parameter>
        <parameter><name>IntrinsicPermeability0</name><type>Constant</type><value>{case.intrinsic_permeability:.16g}</value></parameter>
        <parameter><name>n_l0</name><type>Constant</type><value>{n_l0_fixed:.16g}</value></parameter>
        <parameter><name>rho_lR0</name><type>Constant</type><value>{RHO_LR0:.16g}</value></parameter>
        <parameter><name>epsilon_sw0</name><type>Constant</type><value>0.0</value></parameter>
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
        "hamaker_effective": hamaker_eff,
        "n_l0": n_l0_fixed,
    }


def write_native_project(case: Case, project_path: Path) -> None:
    """Create the equivalent native project for comparison runs."""
    prefix = project_path.stem
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
            </secondary_variables>
            <specific_body_force>0 0</specific_body_force>
            <initial_stress>sigma0</initial_stress>
            <explicit_hm_coupling_in_unsaturated_zone>false</explicit_hm_coupling_in_unsaturated_zone>
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
                            <swelling_pressures>2e6 2e6 2e6</swelling_pressures>
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
                    <name>saturation</name>
                    <type>SaturationTuller</type>
                    <area_factor_tuller>{AREA_FACTOR_TULLER}</area_factor_tuller>
                    <pore_area_shape_factor_tuller>{PORE_AREA_SHAPE_FACTOR_TULLER}</pore_area_shape_factor_tuller>
                    <characteristic_pore_size>{CHARACTERISTIC_PORE_SIZE}</characteristic_pore_size>
                    <surface_tension>{SURFACE_TENSION}</surface_tension>
                </property>
                <property>
                    <name>relative_permeability</name>
                    <type>Constant</type>
                    <value>1</value>
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
                <convergence_criterion><type>PerComponentDeltaX</type><norm_type>NORM2</norm_type><reltols>1e-10 1e-10 1e-9</reltols></convergence_criterion>
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


def run_ogs(ogs_bin: Path, project_path: Path) -> None:
    """Execute OGS for a single temporary project."""
    subprocess.run(
        [str(ogs_bin), str(project_path)],
        cwd=ROOT,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )


def git_short_hash(repo: Path) -> str:
    """Return a short Git hash for provenance tracking."""
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=repo, text=True
        ).strip()
    except Exception:
        return ""
    return out


def run_mfront_case(
    ogs_bin: Path, case: Case, multiplier: float, n_l0_fixed: float, idx: int
) -> dict:
    """Run one MFront calibration trial and return the extracted metric."""
    case_tag = f"dd{int(case.dry_density)}_mfront_{idx:02d}"
    project_path = ROOT / f"{case_tag}.prj"
    meta = write_mfront_project(case, multiplier, n_l0_fixed, project_path)

    try:
        run_ogs(ogs_bin, project_path)
        last_vtu = extract_last_vtu(case_tag)
        pressure_mpa = mean_total_stress_mpa(last_vtu)
    finally:
        cleanup_runtime(case_tag)
        project_path.unlink(missing_ok=True)

    return {
        "case_tag": case_tag,
        "multiplier": multiplier,
        "pressure_mpa": pressure_mpa,
        **meta,
    }


def run_native_case(ogs_bin: Path, case: Case) -> float:
    """Run the native reference project and return the final swelling pressure."""
    case_tag = f"dd{int(case.dry_density)}_native"
    project_path = ROOT / f"{case_tag}.prj"
    write_native_project(case, project_path)

    try:
        run_ogs(ogs_bin, project_path)
        last_vtu = extract_last_vtu(case_tag)
        pressure_mpa = mean_total_stress_mpa(last_vtu)
    finally:
        cleanup_runtime(case_tag)
        project_path.unlink(missing_ok=True)

    return pressure_mpa


def calibrate_multiplier_for_case(
    ogs_bin: Path,
    case: Case,
    target_mpa: float,
    rel_tol: float = 0.02,
    max_iter: int = 12,
) -> dict:
    """Solve for the multiplier that matches the Villar target."""
    n_l0_fixed = n_l0_from_micro_suction(case.phi0, HAMAKER_LITERATURE)

    # 1) Baseline at multiplier=1.
    run0 = run_mfront_case(ogs_bin, case, 1.0, n_l0_fixed, 0)
    p0 = run0["pressure_mpa"]

    # 2) First scaled guess from near-linear response.
    m1 = target_mpa / max(abs(p0), 1e-9)
    m1 = float(np.clip(m1, 1e-6, 1e18))
    run1 = run_mfront_case(ogs_bin, case, m1, n_l0_fixed, 1)
    p1 = run1["pressure_mpa"]

    best = min(
        [run0, run1],
        key=lambda r: abs(r["pressure_mpa"] - target_mpa) / max(target_mpa, 1e-12),
    )

    if abs(p1 - target_mpa) / max(target_mpa, 1e-12) < rel_tol:
        return best

    m_prev, p_prev = 1.0, p0
    m_curr, p_curr = m1, p1
    i = 2

    while i <= max_iter:
        if abs(p_curr - p_prev) > 1e-12:
            m_new = m_curr + (target_mpa - p_curr) * (m_curr - m_prev) / (
                p_curr - p_prev
            )
        else:
            m_new = m_curr * (target_mpa / max(abs(p_curr), 1e-9))
        m_new = float(np.clip(m_new, 1e-6, 1e18))

        run_new = run_mfront_case(ogs_bin, case, m_new, n_l0_fixed, i)
        p_new = run_new["pressure_mpa"]
        err_new = abs(p_new - target_mpa) / max(target_mpa, 1e-12)
        if err_new < abs(best["pressure_mpa"] - target_mpa) / max(target_mpa, 1e-12):
            best = run_new
        if err_new < rel_tol:
            return run_new

        m_prev, p_prev = m_curr, p_curr
        m_curr, p_curr = m_new, p_new
        i += 1

    return best


def main() -> None:
    """Run the full dense calibration sweep and write all summaries."""
    parser = argparse.ArgumentParser(
        description=(
            "Dense dry-density calibration of effective micro vdW multiplier "
            "against Lloret/Villar swelling-pressure curve."
        )
    )
    parser.add_argument(
        "--ogs-bin",
        type=Path,
        default=Path("/Users/vinaykumar/git/build/release-mfront-tpm/bin/ogs"),
        help="Path to OGS binary (must support MFront).",
    )
    parser.add_argument("--dd-min", type=float, default=1400.0)
    parser.add_argument("--dd-max", type=float, default=1800.0)
    parser.add_argument("--dd-step", type=float, default=25.0)
    parser.add_argument("--rel-tol", type=float, default=0.02)
    args = parser.parse_args()

    dd_values = np.arange(args.dd_min, args.dd_max + 0.5 * args.dd_step, args.dd_step)
    cases = [Case(float(dd)) for dd in dd_values]

    rows = []
    for case in cases:
        target = case.villar_target_swelling_mpa
        calibrated = calibrate_multiplier_for_case(
            args.ogs_bin, case, target, rel_tol=args.rel_tol
        )
        native_pressure = run_native_case(args.ogs_bin, case)
        print(
            f"dd={case.dry_density:.0f} kg/m3: "
            f"target={target:.3f} MPa, "
            f"mfront={calibrated['pressure_mpa']:.3f} MPa, "
            f"native={native_pressure:.3f} MPa, "
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
                "specific_surface_mass_m2_g": SPECIFIC_SURFACE_MASS,
                "specific_surface_volumetric_m2_m3": case.specific_surface_volumetric,
                "target_villar_MPa": target,
                "mfront_calibrated_MPa": calibrated["pressure_mpa"],
                "native_MPa": native_pressure,
                "vdw_multiplier": calibrated["multiplier"],
                "hamaker_literature_J": HAMAKER_LITERATURE,
                "hamaker_effective_J": calibrated["hamaker_effective"],
                "n_l0": calibrated["n_l0"],
                "macro_suction_MPa": MACRO_SUCTION_MPA,
                "micro_suction_MPa": MICRO_SUCTION_MPA,
                "mfront_minus_target_MPa": calibrated["pressure_mpa"] - target,
                "native_minus_target_MPa": native_pressure - target,
            }
        )

    rows = sorted(rows, key=lambda r: r["dry_density_kg_m3"])

    csv_path = ROOT / "villar_dense_dd_calibration.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "ogs_repo_hash": git_short_hash(ROOT.parents[3]),
        "materialmodels_repo_hash": git_short_hash(
            Path("/Users/vinaykumar/Documents/GitHub/materialmodels")
        ),
        "hamaker_literature_J": HAMAKER_LITERATURE,
        "specific_surface_mass_m2_g": SPECIFIC_SURFACE_MASS,
        "total_initial_suction_MPa": TOTAL_SUCTION_MPA,
        "macro_initial_suction_MPa": MACRO_SUCTION_MPA,
        "micro_initial_suction_MPa": MICRO_SUCTION_MPA,
        "notes": (
            "Effective multiplier is applied as A_eff = multiplier * A_literature. "
            "n_l0 is recomputed for each (density, multiplier) to preserve the "
            "micro suction split in the initial state."
        ),
        "results": rows,
    }
    (ROOT / "villar_dense_dd_calibration_summary.json").write_text(
        json.dumps(summary, indent=2)
    )

    x = np.array([r["dry_density_kg_m3"] for r in rows])
    y_target = np.array([r["target_villar_MPa"] for r in rows])
    y_mfront = np.array([r["mfront_calibrated_MPa"] for r in rows])
    y_mult = np.array([r["vdw_multiplier"] for r in rows])

    plt.figure(figsize=(8.0, 5.2))
    plt.plot(x, y_target, "k-", linewidth=2.0, label="Villar Eq. (7) target")
    plt.plot(
        x,
        y_mfront,
        color="#d62728",
        marker="s",
        linestyle="--",
        linewidth=1.7,
        markersize=4.5,
        label="MFront calibrated (effective vdW)",
    )
    plt.xlabel("Dry density (kg/m$^3$)")
    plt.ylabel("Swelling pressure at full saturation (MPa)")
    plt.grid(True, alpha=0.35)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(ROOT / "villar_dense_dd_swelling_pressure_comparison.png", dpi=200)
    plt.close()

    plt.figure(figsize=(8.0, 5.2))
    plt.semilogy(
        x,
        y_mult,
        color="#2ca02c",
        marker="D",
        linestyle="-",
        linewidth=1.8,
        markersize=4.5,
    )
    plt.xlabel("Dry density (kg/m$^3$)")
    plt.ylabel("Effective micro vdW multiplier (-)")
    plt.grid(True, which="both", alpha=0.35)
    plt.tight_layout()
    plt.savefig(ROOT / "villar_dense_dd_vdw_multiplier.png", dpi=200)
    plt.close()

    print(f"Wrote: {csv_path}")
    print(f"Wrote: {ROOT / 'villar_dense_dd_calibration_summary.json'}")
    print(f"Wrote: {ROOT / 'villar_dense_dd_swelling_pressure_comparison.png'}")
    print(f"Wrote: {ROOT / 'villar_dense_dd_vdw_multiplier.png'}")


if __name__ == "__main__":
    main()
