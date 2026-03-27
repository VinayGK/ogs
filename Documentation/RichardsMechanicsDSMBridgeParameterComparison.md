# RichardsMechanics DSM Native-vs-Bridge Parameter Comparison

## What this note answers

This note answers a simple question:

Do the current native and bridge parity tests use the same parameters?

Short answer:

- the shared OGS problem setup is aligned well
- the active constitutive parameter sets are not literally the same
- the current green compare test is therefore a reduced overlap test, not exact
  same-parameter MCC plasticity parity

There is now also a separate aligned unsaturated elastic pair where both sides
use the same Tuller saturation law. That pair is exact, but it is only an
elastic saturation-law parity check, not MCC parity.

There is also now a second bridge surface, the hybrid
`RichardsMechanicsNotebookBridge_MCC`, which keeps the verified MCC carrier
surface but also stores notebook auxiliary state. That new hybrid surface is a
staging step toward notebook-driven parity; it is not yet the final
notebook-derived constitutive closure.

In the current verified hybrid step:

- notebook microstate now feeds back into the returned effective stress through
  the swelling correction
- the returned saturation law is still intentionally kept equal to the
  verified MCC carrier surface
- the hybrid now also exposes the notebook support-state outputs `phi`, `n_S`,
  `n_L`, `rho_LR`, `omega_l`, `delta_epsilon_sw`, and `sigma_S`

## Files compared

- native shell:
  [mfront_parity_1element_native.prj](../Tests/Data/RichardsMechanics/mfront_parity_1element_native.prj)
- bridge shell:
  [mfront_parity_1element_bridge.prj](../Tests/Data/RichardsMechanics/mfront_parity_1element_bridge.prj)

For the aligned unsaturated elastic check:

- native shell:
  [mfront_parity_1element_unsat_native.prj](../Tests/Data/RichardsMechanics/mfront_parity_1element_unsat_native.prj)
- bridge shell:
  [mfront_parity_1element_unsat_bridge.prj](../Tests/Data/RichardsMechanics/mfront_parity_1element_unsat_bridge.prj)

## What is the same

These items are the same, or intentionally mapped to the same numerical value,
in the reduced one-element parity pair.

| Item | Native shell | Bridge shell | Meaning |
| --- | --- | --- | --- |
| Mesh and boundary meshes | same | same | Same geometry and same boundary subsets |
| Axial symmetry | `true` | `true` | Same kinematics |
| Time stepping | `t = 0..4`, `4 x dt=1` | same | Same time grid |
| Pressure history | same ramp | same ramp | Same loading path |
| Initial macro pressure | `0` | `0` | Same macro pressure start |
| Initial stress `sigma0` | `(-5e3, -5e3, -5e3, 0)` | same | Same total initial stress |
| `YoungModulus` | `52e6` | `52e6` | Same elastic modulus |
| `PoissonRatio` | `0.3` | `0.3` | Same elastic Poisson ratio |
| Mass exchange coefficient | `5e-15` | `5e-15` | Same numerical value, different hookup |
| `phi0` | `0.432` | `0.432` | Same initial total porosity |
| `phi_tr0` | `0.332` | `0.332` | Same initial transport porosity |
| `biot_coefficient` | `0.6` | `0.6` | Same Biot coupling |
| `permeability` | `2e-21` | `2e-21` | Same hydraulic conductivity input |
| `relative_permeability` | `1` | `1` | Same reduced permeability |
| Macro saturation law | same van Genuchten law | same | Same macro saturation curve |
| Bishop law | same cutoff law | same | Same effective-stress reduction law |
| Liquid viscosity | `1e-3` | `1e-3` | Same viscosity |
| Liquid density | `1e3` | `1e3` | Same macro liquid density |
| Medium solid density | `2780` | `2780` | Same medium property |

## What is not the same

The main difference is the constitutive law itself.

| Item | Native shell | Bridge shell | Meaning |
| --- | --- | --- | --- |
| Constitutive law | `ModCamClay_semiExpl_constE` | `RichardsMechanicsNotebookBridge` | Not the same model |
| Active plastic state | MCC state variables are active | no matching MCC plastic state | Native has plasticity path; bridge does not expose the same one |
| Active microstate inputs | process-side DSM + native law inputs | reduced bridge microstate inputs | Different parameter surfaces |

## Native MCC parameters vs bridge status

These parameters are active on the native side.
Some of them are written into the bridge project file too, but they are not
used by the current bridge law.

| MCC-related item | Native shell | Bridge shell | Status |
| --- | --- | --- | --- |
| `CriticalStateLineSlope` | active | defined only | inactive on bridge |
| `SwellingLineSlope` | active | defined only | inactive on bridge |
| `VirginConsolidationLineSlope` | active | defined only | inactive on bridge |
| `InitialPreConsolidationPressure` | active | defined only | inactive on bridge |
| `InitialVolumeRatio` | active | defined only | inactive on bridge |
| Plastic outputs | requested and available | not available | native only |

So matching numbers in the project files do not automatically mean matching
constitutive behavior.

## Bridge-only active parameters

These parameters are active only in the current bridge law
[RichardsMechanicsNotebookBridge.mfront](../MaterialLib/SolidModels/MFront/RichardsMechanicsNotebookBridge.mfront).

| Bridge-only item | Value | Meaning |
| --- | --- | --- |
| `SwellingSlope` | `0.1` | swelling from micro-porosity change |
| `ReferenceLiquidDensityMacro` | `1000.0` | macro liquid density reference |
| `ReferenceLiquidDensityMicro` | `2072.8234319102588` | micro liquid density reference |
| `ReferenceDensitySolid` | `2470.0` | bridge microstate solid density |
| `MicroLiquidDensityA` | `1.3` | bridge EOS parameter |
| `MicroLiquidDensityB` | `1.0` | bridge EOS parameter |
| `HamakerConstant` | `-6e-20` | bridge micro-potential parameter |
| `SpecificSurface` | `100.0` | bridge micro-potential parameter |
| `AreaFactorTuller` | `1.0` | bridge saturation parameter |
| `PoreAreaShapeFactorTuller` | `0.8584073464102069` | bridge saturation parameter |
| `CharacteristicPoreSize` | `1e-5` | bridge saturation parameter |
| `SurfaceTension` | `0.0715` | bridge saturation parameter |
| `n_l0` | `0.1` | initial bridge micro-liquid content |
| `rho_lR0` | `2072.8234319102588` | initial bridge micro-liquid density |
| `epsilon_sw0` | `0.0` | initial bridge swelling strain |

## Benchmark-only bridge initial state

The benchmark shell uses a different pressure level than the reduced
one-element test, so the bridge needs a different pressure-consistent initial
microstate there.

| Benchmark bridge IC item | Value | Meaning |
| --- | --- | --- |
| `n_l0` | `0.012069019712402708` | pressure-consistent micro-liquid content at benchmark pressure |
| `rho_lR0` | `2267.4495975433856` | pressure-consistent micro-liquid density at benchmark pressure |

These values are part of the run-level benchmark-shell fix.
They do not make the benchmark pair an exact same-model parity pair.

## What exact same-parameter MCC parity would mean

For MCC, the yield function is usually written as

\[
f(q,p,p_c) = q^2 + M^2 p (p - p_c),
\]

with mean effective stress and deviatoric stress measure

\[
p = -\frac{1}{3}\mathrm{tr}(\boldsymbol{\sigma}_{eff}),
\qquad
q = \sqrt{\frac{3}{2}\,\boldsymbol{s}:\boldsymbol{s}}.
\]

An exact same-parameter MCC parity test would therefore need:

- the same active MCC parameters on both sides
- the same active hardening state on both sides
- the same stress meaning on both sides
- the same plastic-state outputs on both sides

In practice that means the bridge would need to consume and expose at least:

- `CriticalStateLineSlope`
- `SwellingLineSlope`
- `VirginConsolidationLineSlope`
- `InitialPreConsolidationPressure`
- `InitialVolumeRatio`
- `EquivalentPlasticStrain`
- `PlasticVolumetricStrain`
- `PreConsolidationPressure`
- `VolumeRatio`

The current bridge does not do that yet.

## What the current green CTest is still good for

The current green compare CTest is still valuable.

It is a valid reduced DSM overlap test for:

- `displacement`
- `pressure`
- `sigma`
- `epsilon`
- `saturation`
- `swelling_stress`

It is just important to describe it honestly:

- good description:
  “reduced native-vs-bridge overlap parity”
- bad description:
  “exact same-parameter MCC plasticity parity”

## Next step for exact MCC parity

1. Extend or replace the current bridge law so it uses the active MCC parameter
   set.
2. Expose matching plastic state variables and outputs on the bridge side.
3. Build a yield-driving shell where both sides use the same active parameter
   list.
4. Only then extend the compare test to include plastic-state outputs.

## What the new hybrid bridge changes

The repository now has a new pressure-coupled bridge:

- [RichardsMechanicsNotebookBridge_MCC.mfront](../MaterialLib/SolidModels/MFront/RichardsMechanicsNotebookBridge_MCC.mfront)

This hybrid bridge uses the active MCC carrier parameters

- `YoungModulus`
- `PoissonRatio`
- `CriticalStateLineSlope`
- `SwellingLineSlope`
- `VirginConsolidationLineSlope`
- `CharacteristicPreConsolidationPressure`
- `ResidualLiquidSaturation`
- `ResidualGasSaturation`
- `BubblePressure`
- `VanGenuchtenExponent_m`
- `NotebookSaturationMode`

and also stores notebook auxiliary-state inputs

- `SwellingSlope`
- `MassExchangeCoefficient`
- `ReferenceLiquidDensityMacro`
- `ReferenceLiquidDensityMicro`
- `ReferenceDensitySolid`
- `MicroLiquidDensityA`
- `MicroLiquidDensityB`
- `HamakerConstant`
- `SpecificSurface`
- `AreaFactorTuller`
- `PoreAreaShapeFactorTuller`
- `CharacteristicPoreSize`
- `SurfaceTension`
- `InitialPorosity`

The widened hybrid bridge now exposes these notebook support-state outputs too:

- `phi`
- `n_S`
- `n_L`
- `rho_LR`
- `omega_l`
- `delta_epsilon_sw`
- `sigma_S`

Those support-state outputs are locked against the committed notebook baseline
histories in
[RichardsMechanicsNotebookBridgeMCC.cpp](../Tests/MaterialLib/MFront/RichardsMechanicsNotebookBridgeMCC.cpp).
- `n_l0`
- `rho_lR0`
- `epsilon_sw0`

So the hybrid bridge now has

\[
\eta_{hybrid} = (\eta_{MCC}, \eta_{NB}),
\]

but in the current verified step it still returns the same constitutive
surface as the pressure-coupled MCC bridge when notebook coupling is neutral.

That means the remaining open work is no longer “add MCC state to the bridge”.
That step is already done. The open work is to decide how notebook state should
feed back into the returned stress, tangents, and possibly saturation without
breaking the verified MCC carrier surface.

That decision is now partially implemented:

- stress-side notebook feedback is active in the hybrid bridge
- saturation-side notebook feedback is now also available in an optional
  reduced-shell Tuller mode
- the same notebook Tuller saturation branch is now also benchmark-compared on
  the part-1 shell
- the neutral benchmark shell for the hybrid bridge is now tested in
  [mfront_restart_part1_notebook_mcc_bridge.prj](../Tests/Data/RichardsMechanics/mfront_restart_part1_notebook_mcc_bridge.prj)

The reduced same-law MCC + notebook saturation pair is now also tested in

- [mfront_parity_1element_notebook_mcc_tuller_native.prj](../Tests/Data/RichardsMechanics/mfront_parity_1element_notebook_mcc_tuller_native.prj)
- [mfront_parity_1element_notebook_mcc_tuller_bridge.prj](../Tests/Data/RichardsMechanics/mfront_parity_1element_notebook_mcc_tuller_bridge.prj)

with the compare CTest

- `ogs-RichardsMechanics_mfront_parity_1element_notebook_mcc_tuller_compare`

The benchmark same-law MCC + notebook saturation pair is now also tested in

- [mfront_restart_part1_notebook_mcc_tuller_native.prj](../Tests/Data/RichardsMechanics/mfront_restart_part1_notebook_mcc_tuller_native.prj)
- [mfront_restart_part1_notebook_mcc_tuller_bridge.prj](../Tests/Data/RichardsMechanics/mfront_restart_part1_notebook_mcc_tuller_bridge.prj)

with the compare CTest

- `ogs-RichardsMechanics_mfront_restart_part1_notebook_mcc_tuller_compare`

So the narrower remaining open question is now:

- which notebook-owned non-neutral constitutive branches should be widened
  beyond the now-tested Tuller saturation branch

## What the new project files test

| File | Purpose | Exact test meaning |
| --- | --- | --- |
| [mfront_parity_1element_notebook_mcc_tuller_native.prj](../Tests/Data/RichardsMechanics/mfront_parity_1element_notebook_mcc_tuller_native.prj) | Reduced native MCC shell with native MPL Tuller saturation | Native reduced reference for notebook-Tuller saturation on the one-element shell |
| [mfront_parity_1element_notebook_mcc_tuller_bridge.prj](../Tests/Data/RichardsMechanics/mfront_parity_1element_notebook_mcc_tuller_bridge.prj) | Reduced hybrid bridge shell with `NotebookSaturationMode = 1` | Bridge reduced reference for the same notebook-Tuller saturation branch |
| [mfront_restart_part1_notebook_mcc_tuller_native.prj](../Tests/Data/RichardsMechanics/mfront_restart_part1_notebook_mcc_tuller_native.prj) | Native benchmark part-1 MCC shell with native MPL Tuller saturation | Native benchmark reference for the notebook-Tuller saturation branch on the full part-1 load scale |
| [mfront_restart_part1_notebook_mcc_tuller_bridge.prj](../Tests/Data/RichardsMechanics/mfront_restart_part1_notebook_mcc_tuller_bridge.prj) | Hybrid bridge benchmark part-1 shell with `NotebookSaturationMode = 1` and neutral notebook swelling/mass exchange | Bridge benchmark reference for the same notebook-Tuller saturation branch |

For the benchmark Tuller pair, the same native part-1 geometry and load range
are kept, but the branch is substepped with fixed `\Delta t = 1` and only the
start and end states are compared. That keeps the benchmark-scale load path
while making the non-neutral saturation branch loadable on both sides.
