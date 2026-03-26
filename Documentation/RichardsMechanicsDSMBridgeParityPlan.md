# RichardsMechanics DSM Native-vs-MFront Parity Plan

## Goal

Establish a native-vs-MFront comparison for the RichardsMechanics double-structure model (DSM) using the same CTest slice. The finished state is:

- one native DSM run;
- one MFront bridge DSM run;
- one binary-vs-binary comparison CTest over the same output fields and timesteps;
- agreement at zero tolerance where numerically possible, with any residual tolerance justified at machine-precision scale only.

## Current Verified State

As of 2026-03-26 on branch `dsm-nb-mfront-transition`:

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
  the full benchmark-shaped bridge run still fails before any accepted output
  with `Notebook bridge microstate Newton line search failed`;
- direct material-point bridge tests at the same benchmark pressure do pass for
  `dt = 0`, both with the raw benchmark stress state and with the
  RM-equivalent effective stress state, so the negative-pressure initial state
  itself is not the remaining blocker;
- the current benchmark blocker is therefore the first nonzero coupled
  microstate solve under the native benchmark load/time scale, not the
  zero-step material-point response;
- the RM process now skips duplicate secondary-variable registrations, which
  removes the `swelling_stress` name clash between process-owned output and
  pressure-coupled MFront internal variables without regressing the reduced
  one-element parity slice;
- the RM pressure-coupled assembly now keeps the process-owned Biot/Bishop
  pressure contribution active even when a pressure-coupled solid returns
  `dStress_dLiquidPressure`, and that change keeps the focused CTest slice
  green.

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

- the bridge material itself accepts the benchmark pressure level
  `pressure_ic = -5e3` at `dt = 0` when its initial bridge microstate is made
  pressure-consistent for that pressure;
- the same direct material-point benchmark response remains finite when the
  previous stress is given either as the raw benchmark initial stress
  `(-5e3, -5e3, -5e3, 0)` or as the RM-equivalent effective stress
  `(-1e4, -1e4, -1e4, 0)`;
- the full benchmark-shaped process run still fails before producing `ts_0`,
  and the failure now surfaces the underlying bridge exception text:
  `Notebook bridge microstate Newton line search failed`;
- the reduced one-element parity slice still stays green after keeping the
  process-owned Biot/Bishop contribution active in the pressure-coupled RM
  assembly path.

What is therefore still blocked:

- the first nonzero benchmark constitutive solve inside the full
  RichardsMechanics process;
- specifically, the bridge microstate Newton solve under the coupled benchmark
  load/time step, not the zero-step material-point evaluation.

Current interpretation:

- this is no longer a blind `status -1` failure; the verified constitutive
  failure mode is the bridge microstate Newton line search under the benchmark
  process step;
- this is not explained by the negative-pressure initial state alone, because
  dedicated material-point tests now pass at that pressure for the two stress
  conventions that matter here;
- the remaining ownership gap is still centered on the RichardsMechanics
  pressure-coupled path and how it drives the first nonzero bridge solve,
  especially through
  [IntegrationPointData.h](../ProcessLib/RichardsMechanics/IntegrationPointData.h),
  [PressureCoupledSolidData.h](../ProcessLib/RichardsMechanics/ConstitutiveRelations/PressureCoupledSolidData.h),
  [RichardsMechanicsFEM-impl.h](../ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h),
  and
  [MFrontRichardsMechanics.h](../MaterialLib/SolidModels/MFront/MFrontRichardsMechanics.h);
- if the benchmark shell is revisited next, the first check should be the
  bridge microstate nonlinear solve under the benchmark `dt = 1000` step
  rather than further deck-only tuning of the zero-step initial state.

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
residue. It is the benchmark-shell failure at the first nonzero coupled
microstate solve.

Primary inspection points for that gap:

- [IntegrationPointData.h](../ProcessLib/RichardsMechanics/IntegrationPointData.h);
- [PressureCoupledSolidData.h](../ProcessLib/RichardsMechanics/ConstitutiveRelations/PressureCoupledSolidData.h);
- [RichardsMechanicsFEM-impl.h](../ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h).
- [MFrontRichardsMechanics.h](../MaterialLib/SolidModels/MFront/MFrontRichardsMechanics.h);

Most likely causes to check:

- pressure-consistent benchmark initialization is now known and should be kept
  fixed while diagnosing the remaining nonzero-step failure;
- if the benchmark path stays process-owned, the bridge/model interface must be
  widened instead of compensating with deck-only tuning;
- verify the stress convention expected by the RM pressure-coupled path
  against the stress convention returned by the candidate MFront behaviour;
- verify the ownership split for `dStress_dLiquidPressure` and
  `dSaturation_dLiquidPressure` against the native RM assembly path;
- inspect the bridge microstate Newton solve itself under the benchmark
  `dt = 1000` step, because the current verified failure mode is a line-search
  breakdown there rather than an opaque MGIS status.

Recommended debugging sequence:

- keep the benchmark bridge initial state at
  `n_l0 = 0.012069019712402708`, `rho_lR0 = 3004.336830222012`;
- keep the new direct material-point tests as the zero-step negative-pressure
  guardrail while changing the full benchmark path;
- reproduce and inspect the first nonzero bridge solve at the benchmark
  `dt = 1000` step;
- trace the RM pressure-coupled solid path from constitutive integration to
  assembly, especially the returned stress, `dStress_dLiquidPressure`, and
  `dSaturation_dLiquidPressure`;
- if needed, instrument or regularize the bridge microstate Newton solve rather
  than re-tuning the pressure-consistent initial state that already passes at
  `dt = 0`;
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
4. Use the new direct material-point negative-pressure tests to separate zero-step constitutive behavior from the full benchmark process failure.
5. Reuse the existing binary-vs-binary compare pattern for the benchmark shell only after the first nonzero benchmark bridge solve is stable.

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

Focused bridge material-point guardrail:

```bash
/Users/vinaykumar/git/build/release-mfront-tpm/bin/testrunner \
  --gtest_filter='MaterialLib_RichardsMechanicsNotebookBridgeMFront.PlaneStrainFactoryPathZeroDtBenchmarkPressureConsistentNegativePressureResponse:MaterialLib_RichardsMechanicsNotebookBridgeMFront.PlaneStrainFactoryPathZeroDtBenchmarkPressureConsistentRMeffectiveStressResponse'
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
