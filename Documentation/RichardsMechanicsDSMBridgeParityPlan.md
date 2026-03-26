# RichardsMechanics DSM Native-vs-MFront Parity Plan

## Goal

Establish a native-vs-MFront comparison for the RichardsMechanics double-structure model (DSM) using the same CTest slice. The finished state is:

- one native DSM run;
- one MFront bridge DSM run;
- one binary-vs-binary comparison CTest over the same output fields and timesteps;
- agreement at zero tolerance where numerically possible, with any residual tolerance justified at machine-precision scale only.

## Current Verified State

As of 2026-03-25 on branch `dsm-nb-mfront-transition`:

- the repository configures and builds with `OGS_USE_MFRONT=ON`;
- the following targeted test slice passes:
  - `ogs-RichardsMechanics/mfront_restart_part1`
  - `ogs-RichardsMechanics/mfront_restart_part2`
  - `ogs-RichardsMechanics_mfront_restart_part1_rm_bridge`
  - `ogs-RichardsMechanics_mfront_parity_1element_native`
  - `ogs-RichardsMechanics_mfront_parity_1element_bridge`
  - `ogs-RichardsMechanics_mfront_parity_1element_compare`
- the bridge restart smoke input is stable enough to be part of the default CTest slice;
- the one-element parity inputs both run successfully and produce five VTU outputs (`ts_0` through `ts_4`).
- the parity comparison CTest runs both one-element shells and compares the generated VTUs directly in the build output directory.
- the bridge restart shell itself is still only a reduced smoke case at
  `t_end = 1`, `pressure_ic = -100`, `top_pressure = -100`;
- an exploratory benchmark-shaped native-vs-bridge pair was rerun on the
  native part-1 load and time scale
  (`t_end = 1000`, `pressure_ic = -5e3`, `top_pressure = -1e5`);
- with a pressure-consistent bridge initial microstate
  `n_l0 = 0.012069019712402708`, `rho_lR0 = 3004.336830222012`,
  the benchmark-shaped pair matches at `ts_0` up to machine precision;
- the same benchmark-shaped pair still diverges at `ts_1`, so the benchmark
  parity task is not closed.
- the RM process now skips duplicate secondary-variable registrations, which
  removes the `swelling_stress` name clash between process-owned output and
  pressure-coupled MFront internal variables without regressing the reduced
  one-element parity slice.

## Verified Parity Findings

The current one-element parity pair now compares as follows:

- exact match for `displacement`;
- exact match for `pressure`;
- exact match for `epsilon`;
- exact match for `saturation`;
- exact match for `swelling_stress`;
- machine-precision mismatch only for `sigma`.

Observed `sigma` mismatch after aligning the bridge shell:

- present already at `t = 0`;
- persists through `t = 4`;
- absolute maximum norm is about `1.8e-12` per normal component;
- relative maximum norm is about `3.6e-16`.

This is no longer a model-parity blocker. It is only floating-point residue from two numerically equivalent paths not rounding bit-for-bit identically.

## Current Parity Test Shape

The parity slice now contains:

- one native one-element DSM shell:
  [mfront_parity_1element_native.prj](../Tests/Data/RichardsMechanics/mfront_parity_1element_native.prj);
- one bridge one-element DSM shell:
  [mfront_parity_1element_bridge.prj](../Tests/Data/RichardsMechanics/mfront_parity_1element_bridge.prj);
- one dedicated binary-vs-binary compare test implemented via:
  [CompareRichardsMechanicsMFrontParity.cmake](../scripts/cmake/test/CompareRichardsMechanicsMFrontParity.cmake).

Important implementation notes:

- the compare CTest is a plain `add_test()` entry in
  [Tests.cmake](../ProcessLib/RichardsMechanics/Tests.cmake);
- it runs the native shell, runs the bridge shell, and then executes direct
  binary-vs-binary `vtkdiff` checks;
- `sigma` uses `5e-12` absolute tolerance and zero relative tolerance;
- all other compared fields use `1e-15` absolute tolerance and zero relative tolerance.

The small `sigma` tolerance is justified by the measured machine-precision residue only. It is not masking a model discrepancy.

## Parameter Comparison And Plasticity Status

For a parameter-by-parameter comparison of the current reduced parity pair, see:

- [RichardsMechanicsDSMBridgeParameterComparison.md](./RichardsMechanicsDSMBridgeParameterComparison.md).

Important consequence:

- the current one-element parity slice is a reduced shared-overlap DSM parity test;
- it is not yet an exact same-parameter MCC plasticity parity test.

## Benchmark-Shell Status

The benchmark-shell blocker is now narrower and better characterized than it
was at the start of this effort.

What is now verified:

- the bridge can be stabilized on the native part-1 load and time scale if its
  initial bridge microstate is made pressure-consistent for
  `pressure_ic = -5e3`;
- with that equilibrium initial state, the benchmark-shaped native and bridge
  shells match at `ts_0` for `displacement`, `pressure`, `epsilon`,
  `saturation`, and `swelling_stress`, with `sigma` differing only at
  machine-precision scale;
- the earlier large `ts_0` stress jump was caused by starting the benchmark
  bridge shell from a one-element parity microstate that is not in equilibrium
  at the benchmark initial pressure.
- replaying the native process-owned micro-porosity path
  (`micro_porosity`, `saturation_micro`, `swelling_stress_rate`,
  `transport_porosity`) inside the bridge benchmark shell still leaves the
  `ts_1` gap, so that process-owned DSM update is no longer the primary
  suspect by itself.
- a pressure-coupled MCC bridge attempt using `ModCamClay_semiExpl_coupled`
  reproduces `saturation`, `micro_saturation`, `EquivalentPlasticStrain`,
  `PreConsolidationPressure`, `PlasticVolumetricStrain`, and
  `swelling_stress_ip` exactly at `ts_1`, while `pressure`, `micro_pressure`,
  `displacement`, `sigma`, `ElasticStrain`, `porosity`,
  `transport_porosity`, and `VolumeRatio` still diverge.

What still fails at `ts_1`:

- `displacement` absolute maximum norm about `8.31e-4`;
- `pressure` absolute maximum norm about `3.22e5`;
- `sigma` absolute maximum norm about `5.81e4`;
- `epsilon` absolute maximum norm about `8.31e-4`;
- `swelling_stress` absolute maximum norm about `5.65e-2`;
- `saturation` still matches exactly.

Current interpretation:

- this is no longer an initial-condition bug;
- the process-owned micro-porosity path is no longer the only plausible
  blocker, because replaying that path in the bridge shell keeps `ts_0`
  exact but still leaves the `ts_1` gap;
- the exact agreement in the plastic state variables, saturation variables,
  and `swelling_stress_ip` narrows the remaining mismatch to the
  stress/pressure coupling boundary rather than to MCC plasticity itself;
- the most likely remaining ownership gap is now the RichardsMechanics
  pressure-coupled solid contract, especially how
  [IntegrationPointData.h](../ProcessLib/RichardsMechanics/IntegrationPointData.h),
  [PressureCoupledSolidData.h](../ProcessLib/RichardsMechanics/ConstitutiveRelations/PressureCoupledSolidData.h),
  [RichardsMechanicsFEM-impl.h](../ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h),
  and
  [MFrontRichardsMechanics.h](../MaterialLib/SolidModels/MFront/MFrontRichardsMechanics.h)
  interpret the returned stress and pressure-coupling blocks.

## Why The Remaining Gap Looks Interface-Owned

At continuum level the quasi-static balance of momentum is

\[
\nabla\cdot\boldsymbol{\sigma}_\mathrm{tot} + \rho\mathbf{b} = \mathbf{0}.
\]

The native RichardsMechanics benchmark is assembled around the split

\[
p_\mathrm{cap} = -p_L,
\qquad
p_F^R = -\chi(S_L)\,p_\mathrm{cap} = \chi(S_L)\,p_L,
\]

with total and effective stress related by

\[
\boldsymbol{\sigma}_\mathrm{tot}
=
\boldsymbol{\sigma}_\mathrm{eff}
- \alpha\,p_F^R\,\mathbf{I}.
\]

Here:

- \(p_L\) is the liquid pressure;
- \(p_\mathrm{cap}\) is the capillary pressure;
- \(\chi(S_L)\) is Bishop's effective-stress factor;
- \(\alpha\) is the Biot coefficient;
- \(\mathbf{I}\) is the identity tensor.

So the continuum stress closure expected by RM is

\[
\boldsymbol{\sigma}_\mathrm{tot}
=
\boldsymbol{\sigma}_\mathrm{eff}
- \alpha\,\chi(S_L)\,p_L\,\mathbf{I}.
\]

The native DSM microstate update is also process-owned. In reduced form, the
benchmark advances

\[
(\Delta \phi_m,\ \Delta \varepsilon_{sw},\ \Delta p_L^m,\ \Delta
\boldsymbol{\sigma}_{sw})
=
\mathcal{G}(S_L^m,\ \dot{\boldsymbol{\sigma}}_{sw},\ k_{ex},\dots),
\]

where \(\phi_m\) is the micro-porosity, \(\varepsilon_{sw}\) is the swelling
strain, \(p_L^m\) is the micro-liquid pressure, and \(k_{ex}\) is the mass
exchange coefficient. In OGS this update is carried by
`computeMicroPorosity(...)`.

The pressure-coupled MFront bridge contract is different in shape. It returns a
constitutive package

\[
\left(
\boldsymbol{\sigma}^{\star},
S_L^{\star},
\frac{\partial \boldsymbol{\sigma}^{\star}}{\partial p_L},
\frac{\partial S_L^{\star}}{\partial p_L}
\right),
\]

which RichardsMechanics then consumes through the pressure-coupled solid path.
The benchmark only closes if both sides agree that
\(\boldsymbol{\sigma}^{\star}\) means

- effective stress,
- total stress,
- or a mixed stress already containing some pore-pressure contribution.

The formal continuum contract required by RM is that

\[
\boldsymbol{\sigma}^{\star} = \boldsymbol{\sigma}_\mathrm{eff},
\qquad
\boldsymbol{\sigma}_\mathrm{tot}
=
\boldsymbol{\sigma}^{\star}
- \alpha\,\chi(S_L^\star)\,p_L\,\mathbf{I},
\]

with pressure derivative

\[
\frac{\partial \boldsymbol{\sigma}_\mathrm{tot}}{\partial p_L}
=
\frac{\partial \boldsymbol{\sigma}^{\star}}{\partial p_L}
- \alpha
\left[
\chi(S_L^\star)
+ p_L\,\chi'(S_L^\star)\,
\frac{\partial S_L^\star}{\partial p_L}
\right]
\mathbf{I}.
\]

If the bridge instead returns total stress while RM still interprets it as
effective stress, the pore-pressure contribution is subtracted twice.

The current diagnostics say the plastic law itself is probably not the dominant
problem. For the coupled MCC benchmark attempt, the following quantities match
exactly at `ts_1`:

\[
\Lambda_p,\qquad
p_c,\qquad
\varepsilon_v^p,\qquad
S_L,\qquad
S_L^m,\qquad
\boldsymbol{\sigma}_{sw}.
\]

Those are the key state variables for the MCC yield and hardening part:

\[
f(q,p,p_c) = q^2 + M^2\,p\,(p-p_c),
\qquad
p = -\frac{1}{3}\operatorname{tr}(\boldsymbol{\sigma}_\mathrm{eff}),
\qquad
q = \sqrt{\frac{3}{2}\,\boldsymbol{s}:\boldsymbol{s}}.
\]

If the dominant bug were in MCC plasticity itself, one would expect
\(\Lambda_p\), \(p_c\), or \(\varepsilon_v^p\) to drift as well. They do not.

What still mismatches at `ts_1` is the hydro-mechanical side:

\[
p_L,\qquad
p_L^m,\qquad
\mathbf{u},\qquad
\boldsymbol{\sigma},\qquad
\boldsymbol{\varepsilon}^{el},\qquad
\phi,\qquad
\phi_{tr},\qquad
v^r.
\]

That pattern is much more consistent with a contract mismatch in how the RM
pressure-coupled path interprets
\(\boldsymbol{\sigma}^{\star}\) and
\(\partial\boldsymbol{\sigma}^{\star}/\partial p_L\), or with a mismatch in
how the returned pressure-coupled stress feeds the momentum equation and the
pressure Jacobian.

## Work Packages

### WP1: Keep The Reduced One-Element Parity Slice Green

This work package is complete for the reduced one-element shell.

Current verified scope:

- same process path and same governing equations;
- same geometry and same mesh;
- same axial/plane-strain setting;
- same time stepping;
- same pressure history;
- same initial stress;
- same porosity and transport porosity setup;
- same saturation laws;
- same Biot and effective stress settings;
- same swelling law semantics;
- same densities and other material parameters;
- same output variables.

Acceptance criteria:

- both project files represent the same DSM problem;
- both runs remain green as standalone CTests;
- field mismatch analysis is meaningful because the inputs are genuinely comparable.

### WP2: Extend The Bridge To The Benchmark Shell

The remaining real work is no longer the reduced one-element `sigma`
residue. It is the benchmark-shell evolution gap after the first accepted
step.

Primary inspection points for that gap:

- [IntegrationPointData.h](../ProcessLib/RichardsMechanics/IntegrationPointData.h);
- [PressureCoupledSolidData.h](../ProcessLib/RichardsMechanics/ConstitutiveRelations/PressureCoupledSolidData.h);
- [RichardsMechanicsFEM-impl.h](../ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h).
- [MFrontRichardsMechanics.h](../MaterialLib/SolidModels/MFront/MFrontRichardsMechanics.h);

Most likely causes to check:

- pressure-consistent benchmark initialization is now known and should be kept
  fixed while diagnosing the remaining step-1 mismatch;
- if the benchmark path stays process-owned, the bridge/model interface must be
  widened instead of compensating with deck-only tuning;
- verify the stress convention expected by the RM pressure-coupled path
  against the stress convention returned by the candidate MFront behaviour;
- verify the ownership split for `dStress_dLiquidPressure` and
  `dSaturation_dLiquidPressure` against the native RM assembly path.

Recommended debugging sequence:

- keep the benchmark bridge initial state at
  `n_l0 = 0.012069019712402708`, `rho_lR0 = 3004.336830222012`;
- compare the first accepted benchmark step at `t = 1000`;
- keep the process-owned micro-porosity replay available as a control case,
  because it already proves that matching that path alone does not close
  `ts_1`;
- trace the RM pressure-coupled solid path from constitutive integration to
  assembly, especially the returned stress, `dStress_dLiquidPressure`, and
  `dSaturation_dLiquidPressure`;
- decide whether the benchmark path belongs inside a dedicated
  RichardsMechanics-specific MFront bridge behaviour or in a widened
  process-to-bridge interface;
- only after that, add or extend unit coverage around the chosen ownership
  boundary.

Acceptance criteria:

- the benchmark-shaped native and bridge shells both run on the native part-1
  load and time scale;
- `displacement`, `pressure`, `sigma`, `epsilon`, `saturation`, and
  `swelling_stress` match at the benchmark output timesteps with only
  machine-precision residue where justified.

### WP3: Add The Benchmark Comparison CTest

After the benchmark-shell evolution gap is closed, add the actual
benchmark-shell native-vs-bridge compare test.

Recommended implementation:

- keep the two run tests as separate smoke/run prerequisites;
- add a dedicated comparison test that depends on both run tests;
- invoke `vtkdiff` directly on binary outputs for:
  - `displacement`
  - `pressure`
  - `sigma`
  - `epsilon`
  - `saturation`
  - `swelling_stress`
- compare timesteps `0` through `4`.

Implementation note:

- do not use `OgsTest` for this comparison;
- do not use `AddTest(TESTER vtkdiff)` directly for this comparison;
- use a small custom `add_test()` wrapper or a CMake script that performs binary-vs-binary `vtkdiff` calls.

Status:

- completed for the one-element parity shell.
- not yet implemented for the native part-1 benchmark shell.

Acceptance criteria:

- the compare CTest passes in the MFront-enabled build;
- the compare CTest is part of the same RichardsMechanics CTest slice as the two run tests;
- the smoke tests remain, and the compare test is the parity gate for this reduced benchmark.

## Suggested Execution Order

1. Keep the new compare CTest green while changing the bridge or the native DSM shell.
2. Keep the benchmark-shell bridge initialization pressure-consistent while debugging the native part-1 load case.
3. Treat the reduced one-element `sigma` tolerance question as a floating-point reproducibility issue, not as the remaining benchmark blocker.
4. Reuse the existing binary-vs-binary compare pattern for the benchmark shell only after the step-1 evolution gap is closed.

## Useful Commands

Configure:

```bash
cmake -S /Users/vinaykumar/git/ogs-TPM_Swelling_MCC_Coupled \
  -B /Users/vinaykumar/git/build/release-mfront-tpm \
  -DCMAKE_BUILD_TYPE=Release \
  -DOGS_USE_MFRONT=ON \
  -DCPM_SOURCE_CACHE=/Users/vinaykumar/git/.cpm-cache \
  -DTFELHOME=/Users/vinaykumar/.cache/CPM/_ext/TFEL/895fd7874cfa49079f7711831b3fa10069d8eee57388ab29ff17c265f8ffe777 \
  -DMFRONT=/Users/vinaykumar/.cache/CPM/_ext/TFEL/895fd7874cfa49079f7711831b3fa10069d8eee57388ab29ff17c265f8ffe777/bin/mfront \
  -DMFRONT_QUERY=/Users/vinaykumar/.cache/CPM/_ext/TFEL/895fd7874cfa49079f7711831b3fa10069d8eee57388ab29ff17c265f8ffe777/bin/mfront-query
```

Targeted test slice:

```bash
ctest --test-dir /Users/vinaykumar/git/build/release-mfront-tpm \
  --output-on-failure \
  -R 'ogs-RichardsMechanics/mfront_restart_part1$|ogs-RichardsMechanics/mfront_restart_part2$|ogs-RichardsMechanics_mfront_restart_part1_rm_bridge$|ogs-RichardsMechanics_mfront_parity_1element_native$|ogs-RichardsMechanics_mfront_parity_1element_bridge$|ogs-RichardsMechanics_mfront_parity_1element_compare$'
```

Manual parity check pattern:

```bash
/Users/vinaykumar/git/build/release-mfront-tpm/bin/vtkdiff \
  native.vtu bridge.vtu -a sigma -b sigma --abs 0 --rel 0
```

## Definition Of Done

The reduced one-element DSM native-vs-MFront parity task is done when all of the following are true:

- the parity pair is a true native DSM vs MFront bridge DSM comparison for the reduced shell;
- `displacement`, `pressure`, `epsilon`, `saturation`, and `swelling_stress` match exactly;
- `sigma` matches within machine-precision absolute tolerance only;
- the comparison is enforced by a dedicated CTest in the same RichardsMechanics slice;
- the compare test passes in the MFront-enabled build without manual intervention.

The benchmark-shell parity task is done only when all of the following are true:

- the bridge runs on the native part-1 load and time scale without reducing the
  shell back to the smoke-only boundary;
- the bridge uses a pressure-consistent initial microstate for the benchmark
  pressure level;
- `displacement`, `pressure`, `sigma`, `epsilon`, `saturation`, and
  `swelling_stress` match for the benchmark outputs with only justified
  machine-precision residue;
- the benchmark comparison is enforced by a dedicated CTest, not just by manual
  `vtkdiff` checks.
