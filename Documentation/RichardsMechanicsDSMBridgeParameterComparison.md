# RichardsMechanics DSM Native-vs-Bridge Parameter Comparison

## What this note answers

This note answers a simple question:

Do the current native and bridge parity tests use the same parameters?

Short answer:

- the shared OGS problem setup is aligned well
- the active constitutive parameter sets are not literally the same
- the current green compare test is therefore a reduced overlap test, not exact
  same-parameter MCC plasticity parity

## Files compared

- native shell:
  [mfront_parity_1element_native.prj](../Tests/Data/RichardsMechanics/mfront_parity_1element_native.prj)
- bridge shell:
  [mfront_parity_1element_bridge.prj](../Tests/Data/RichardsMechanics/mfront_parity_1element_bridge.prj)

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
