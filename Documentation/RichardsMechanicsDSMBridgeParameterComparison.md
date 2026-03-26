# RichardsMechanics DSM Native-vs-Bridge Parameter Comparison

## Scope

This note compares the current reduced one-element parity pair:

- native shell:
  [mfront_parity_1element_native.prj](../Tests/Data/RichardsMechanics/mfront_parity_1element_native.prj);
- bridge shell:
  [mfront_parity_1element_bridge.prj](../Tests/Data/RichardsMechanics/mfront_parity_1element_bridge.prj).

Current verdict:

- the shared OGS problem setup is aligned well enough for the existing reduced parity CTest;
- the two parameter decks are not literally equivalent;
- the pair is not an exact same-parameter MCC plasticity parity test.

## Shared Setup And Active Overlap Parameters

These items are the same, or intentionally mapped with the same numerical value, on both sides of the reduced one-element parity pair.

| Item | Native shell | Bridge shell | Status | Notes |
| --- | --- | --- | --- | --- |
| Meshes | `square_1x1_quad_1e0` + side meshes | same | same and active | Same geometry and boundary subsets. |
| Axial symmetry | `true` | `true` | same and active | Same kinematic setting. |
| Time stepping | `t = 0..4`, `4 x dt=1` | same | same and active | Same output cadence too. |
| Pressure history | `pressure_ramp_parity` scaled by `1000` | same | same and active | Same ramp and same `pressure_bc_scale`. |
| Pressure initial condition | `0` | `0` | same and active | Same initial macro pressure. |
| Initial stress `sigma0` | `(-5e3, -5e3, -5e3, 0)` | same | same and active | Same initial total stress field. |
| `YoungModulus` | `52e6` | `52e6` | same and active | Consumed by both constitutive laws. |
| `PoissonRatio` | `0.3` | `0.3` | same and active | Consumed by both constitutive laws. |
| Mass exchange | `<micro_porosity><mass_exchange_coefficient>5e-15` | `MassExchangeCoefficient = 5e-15` | same value, different hookup | Native uses the RM process micro-porosity block; bridge uses an MFront material property. |
| `phi0` | `0.432` | `0.432` | same and active | Native uses it as total porosity IC; bridge also feeds `InitialPorosity`. |
| `phi_tr0` | `0.332` | `0.332` | same and active | Same transport-porosity IC. |
| `biot_coefficient` | `0.6` | `0.6` | same and active | Same medium property. |
| `permeability` | `2e-21` | `2e-21` | same and active | Same medium property. |
| `relative_permeability` | `1` | `1` | same and active | Same medium property. |
| Macro `saturation` law | van Genuchten, exponent `0.4`, `p_b = 15e6` | same | same and active | Same macro saturation law. |
| `bishops_effective_stress` | `BishopsSaturationCutoff`, cutoff `1` | same | same and active | Same effective-stress reduction law. |
| Liquid viscosity | `1e-3` | `1e-3` | same and active | Same medium property. |
| Liquid density | `1e3` | `1e3` | same and active | Same medium property. |
| Solid density in medium | `2780` | `2780` | same and active | Same solid phase density property. |

## Native MCC Parameters Versus Bridge Status

The native reduced shell uses an MCC constitutive law:
[mfront_parity_1element_native.prj](../Tests/Data/RichardsMechanics/mfront_parity_1element_native.prj)
binds `ModCamClay_semiExpl_constE`.

The bridge shell does not bind an MCC law:
[mfront_parity_1element_bridge.prj](../Tests/Data/RichardsMechanics/mfront_parity_1element_bridge.prj)
binds `RichardsMechanicsNotebookBridge`.

| MCC-related item | Native shell | Bridge shell | Status | Notes |
| --- | --- | --- | --- | --- |
| Constitutive behaviour | `ModCamClay_semiExpl_constE` | `RichardsMechanicsNotebookBridge` | not the same | Different constitutive models. |
| `CriticalStateLineSlope` | `1.2`, active | `1.2`, defined only | same value, inactive on bridge | Not consumed by `RichardsMechanicsNotebookBridge`. |
| `SwellingLineSlope` | `6.6e-3`, active | `6.6e-3`, defined only | same value, inactive on bridge | Native MCC parameter only. |
| `VirginConsolidationLineSlope` | `7.7e-2`, active | `7.7e-2`, defined only | same value, inactive on bridge | Native MCC parameter only. |
| `InitialPreConsolidationPressure` | `2e5`, active | `2e5`, defined only | same value, inactive on bridge | Native hardening state only. |
| `InitialVolumeRatio` | `1.78571428571428571429`, active | same value, defined only | same value, inactive on bridge | Native hardening state only. |
| Plastic state outputs | `EquivalentPlasticStrain`, `PlasticVolumetricStrain`, `PreConsolidationPressure`, `VolumeRatio` requested | not requested and not available | native only | Bridge behaviour does not expose matching MCC state. |

## Bridge-Only Active Reduced-Microstate Parameters

These parameters are active only in the bridge behaviour described in
[RichardsMechanicsNotebookBridge.mfront](../MaterialLib/SolidModels/MFront/RichardsMechanicsNotebookBridge.mfront).

| Bridge-only item | Value | Status | Notes |
| --- | --- | --- | --- |
| `SwellingSlope` | `0.1` | active on bridge only | Bridge uses isotropic swelling from micro-porosity change. |
| `ReferenceLiquidDensityMacro` | `1000.0` | active on bridge only | Consumed by bridge EOS / chemical potential update. |
| `ReferenceLiquidDensityMicro` | `2072.8234319102588` | active on bridge only | Consumed by bridge EOS in the current reduced parity shell. |
| `ReferenceDensitySolid` | `2470.0` | active on bridge only | Internal bridge microstate parameter; distinct from medium solid density `2780`. |
| `MicroLiquidDensityA` | `1.3` | active on bridge only | Bridge EOS parameter. |
| `MicroLiquidDensityB` | `1.0` | active on bridge only | Bridge EOS parameter. |
| `HamakerConstant` | `-6e-20` | active on bridge only | Bridge micro potential parameter. |
| `SpecificSurface` | `100.0` | active on bridge only | Bridge micro potential parameter. |
| `AreaFactorTuller` | `1.0` | active on bridge only | Bridge saturation update parameter. |
| `PoreAreaShapeFactorTuller` | `0.8584073464102069` | active on bridge only | Bridge saturation update parameter. |
| `CharacteristicPoreSize` | `1e-5` | active on bridge only | Bridge saturation update parameter. |
| `SurfaceTension` | `0.0715` | active on bridge only | Bridge saturation update parameter. |
| `n_l0` | `0.1` | active on bridge only | Initial bridge micro-liquid content. |
| `rho_lR0` | `2072.8234319102588` | active on bridge only | EOS-consistent initial bridge micro-liquid density. |
| `epsilon_sw0` | `0.0` | active on bridge only | Initial bridge swelling strain state. |

For the native part-1 benchmark pressure level `pressure_ic = -5e3`, these
one-element bridge initial-state values are not pressure-universal. Solving the
current bridge microstate equations at the benchmark initial pressure gives a
different equilibrium anchor:

| Benchmark-only bridge IC item | Value | Status | Notes |
| --- | --- | --- | --- |
| `n_l0` | `0.012069019712402708` | benchmark-equilibrium only | Pressure-consistent bridge micro-liquid content for the native part-1 initial pressure. |
| `rho_lR0` | `2267.4495975433856` | benchmark-equilibrium only | Pressure-consistent bridge micro-liquid density used by the tracked benchmark bridge shell. |

Those values are now part of the tracked run-level benchmark-shell solution.
Dedicated material-point tests confirm that the bridge accepts them at
`pressure_ic = -5e3` for `dt = 0`, for `dt = 1000` at the exact RM-aligned
first-step anchor, and at the exact former process-failure state. The tracked
benchmark bridge deck also now runs on the native part-1 load/time scale.
Benchmark-shell parity is still not closed, however, because the benchmark
`ts_1` fields still differ materially from the native deck and the bridge
benchmark deck is still the reduced `RichardsMechanicsNotebookBridge` law, not
the native MCC benchmark law.

## Bridge-Defined But Currently Inactive Parameters

These entries are present in the bridge project file but are not consumed by
[RichardsMechanicsNotebookBridge.mfront](../MaterialLib/SolidModels/MFront/RichardsMechanicsNotebookBridge.mfront).

| Bridge parameter | Value | Status | Notes |
| --- | --- | --- | --- |
| `CriticalStateLineSlope` | `1.2` | inactive | Present only to mirror the native deck numerically. |
| `SwellingLineSlope` | `6.6e-3` | inactive | Present only to mirror the native deck numerically. |
| `VirginConsolidationLineSlope` | `7.7e-2` | inactive | Present only to mirror the native deck numerically. |
| `InitialPreConsolidationPressure` | `2e5` | inactive | Present only to mirror the native deck numerically. |
| `InitialVolumeRatio` | `1.78571428571428571429` | inactive | Present only to mirror the native deck numerically. |
| `SaturationPressureScale` | `1e3` | inactive | Leftover parameter; not wired into the bridge behaviour. |

## Plasticity Coverage Status

Native RM already has MCC-capable test coverage in this repository:

- [double_porosity_swelling_RM.prj](../Tests/Data/RichardsMechanics/DoubleStructureBenchmark/double_porosity_swelling_RM.prj)
  is registered in
  [Tests.cmake](../ProcessLib/RichardsMechanics/Tests.cmake)
  and uses `ModCamClay_semiExpl`.

The current reduced one-element native-vs-bridge parity pair is still not an MCC plasticity parity test:

- the native one-element shell binds MCC and requests plastic outputs;
- the bridge behaviour stores only `n_l`, `rho_lR`, and `epsilon_sw`, and updates stress by an elastic law plus swelling contribution;
- the current one-element pressure path is small (`pressure_bc_scale = 1000`) relative to the native `InitialPreConsolidationPressure = 2e5`, so this shell should be treated as a reduced non-yielding overlap test unless dedicated plastic-state outputs prove otherwise.

## Consequence For CTests

The current green compare CTest is valid as a reduced DSM overlap test for:

- `displacement`
- `pressure`
- `sigma`
- `epsilon`
- `saturation`
- `swelling_stress`

It is not valid as proof of exact same-parameter MCC plasticity parity.

No CTest was changed here to force plastic yielding, because that would not create an exact same-model comparison and would predictably fail for constitutive-model reasons rather than infrastructure reasons.

## Plan For Exact Same-Parameter MCC Plasticity Parity

1. Replace or extend
   [RichardsMechanicsNotebookBridge.mfront](../MaterialLib/SolidModels/MFront/RichardsMechanicsNotebookBridge.mfront)
   so that it consumes the active MCC parameter set:
   `CriticalStateLineSlope`, `SwellingLineSlope`,
   `VirginConsolidationLineSlope`, `InitialPreConsolidationPressure`,
   and `InitialVolumeRatio`.
2. Expose matching plastic state variables and secondary outputs on the bridge side:
   `EquivalentPlasticStrain`, `PlasticVolumetricStrain`,
   `PreConsolidationPressure`, and `VolumeRatio`.
3. Build a yield-driving one-element shell with the same mesh, BCs, time stepping,
   and active parameter list on both sides.
4. Only after steps 1 to 3 are green, extend the compare CTest to diff the plastic
   state outputs in addition to the existing overlap fields.
