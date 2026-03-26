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
- the bridge restart run input is stable enough to be part of the default CTest slice;
- the one-element parity inputs both run successfully and produce five VTU outputs (`ts_0` through `ts_4`).
- the parity comparison CTest runs both one-element shells and compares the generated VTUs directly in the build output directory.
- the tracked bridge restart shell now runs on the native part-1 load and time
  scale
  (`t_end = 1000`, `pressure_ic = -5e3`, `top_pressure = -1e5`);
- that run-level benchmark-shell fix has three verified pieces:
  - a pressure-consistent bridge initial microstate
    `n_l0 = 0.012069019712402708`, `rho_lR0 = 2267.4495975433856`;
  - a bracketed fallback for the bridge microstate solve when the original
    Newton line search stalls;
  - damped global Newton in the benchmark bridge deck
    (`damping = 0.1`, `damping_reduction = 20`);
- direct material-point bridge tests now pass for:
  - `dt = 0` at the benchmark pressure with the raw benchmark stress state;
  - `dt = 0` at the benchmark pressure with the RM-equivalent effective stress
    state;
  - `dt = 1000` at the benchmark pressure/state anchor;
  - the exact former process-failure state extracted from the coupled run;
- a clean native-vs-bridge benchmark-shell rerun shows exact agreement at
  `ts_0` for `displacement`, `pressure`, `epsilon`, `saturation`, and
  `swelling_stress`, with only machine-precision `sigma` residue
  (`1.8e-12` absolute, `3.6e-16` relative);
- that same clean rerun still shows a material benchmark-shell mismatch at
  `ts_1`, with maximum absolute differences
  `|Δu_y|_∞ ≈ 1.51e-4`, `|Δp|_∞ ≈ 6.80e3`, `|Δσ|_∞ ≈ 9.21e3`,
  `|Δε_yy|_∞ ≈ 1.51e-4`, and `|ΔS_L|_∞ ≈ 1.87e-1`;
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

The benchmark-shell loadability problem is solved for the current tracked
bridge benchmark deck.

What is now verified:

- the bridge material itself accepts the benchmark pressure level
  `pressure_ic = -5e3` at `dt = 0` when its initial bridge microstate is made
  pressure-consistent for that pressure;
- the same direct material-point benchmark response remains finite when the
  previous stress is given either as the raw benchmark initial stress
  `(-5e3, -5e3, -5e3, 0)` or as the RM-equivalent effective stress
  `(-1e4, -1e4, -1e4, 0)`;
- the exact former process-failure state now also passes as a direct bridge
  regression;
- the tracked benchmark bridge deck now reaches `ts_0` and `ts_1` on the true
  native part-1 load/time scale;
- the reduced one-element parity slice still stays green after the benchmark
  loadability fixes.

What is therefore still blocked:

- quantitative native-vs-bridge benchmark-shell parity at `ts_1`;
- a benchmark compare CTest on the true part-1 shell.

Current interpretation:

- the solved part is benchmark-shell run stability, not benchmark-shell
  quantitative parity;
- the verified run-level solution is a combination of the pressure-consistent
  bridge initial microstate, the bracketing fallback in the bridge microstate
  solve, and damped global Newton in the benchmark bridge deck;
- the clean `ts_0` match shows that shell alignment and initialization are now
  good enough to compare the first real benchmark step;
- the remaining `ts_1` mismatch cannot yet be interpreted as a pure RM
  interface bug, because the tracked bridge benchmark deck still binds the
  reduced `RichardsMechanicsNotebookBridge` behaviour while the native
  benchmark deck binds `ModCamClay_semiExpl_constE`;
- the continuum pressure-coupled contract below therefore remains the correct
  audit baseline for future interface work, but the current tracked
  benchmark-shell mismatch is still a combined reduced-law-plus-process gap,
  not an isolated same-model parity result.

## Why Benchmark-Shell Parity Is Still Open

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

The continuum contract above is still the right interpretation surface for
future RM pressure-coupled audit work, but the current tracked benchmark-shell
comparison does not yet isolate an interface bug by itself.

What is now solved is the run-level benchmark-shell failure. The bridge deck
reaches the native part-1 load/time scale and produces both `ts_0` and `ts_1`
after:

- using the pressure-consistent initial bridge microstate
  `n_l0 = 0.012069019712402708`,
  `rho_lR0 = 2267.4495975433856`;
- falling back from the original microstate Newton line search to a bracketed
  reduced solve when necessary;
- damping the global Newton process solve.

The clean benchmark-shell comparison is now easy to summarize:

- at `ts_0`, `displacement`, `pressure`, `epsilon`, `saturation`, and
  `swelling_stress` match exactly, while `sigma` differs only by the same
  machine-scale residue already seen in the reduced one-element shell;
- at `ts_1`, the benchmark-shell mismatch is still material:
  `|Δu_y|_∞ ≈ 1.51e-4`, `|Δp|_∞ ≈ 6.80e3`,
  `|Δσ|_∞ ≈ 9.21e3`, `|Δε_yy|_∞ ≈ 1.51e-4`,
  `|ΔS_L|_∞ ≈ 1.87e-1`.

That mismatch is not yet a pure same-model parity signal, because the tracked
native and bridge benchmark decks still do not bind the same constitutive law:

- the native deck uses `ModCamClay_semiExpl_constE`;
- the bridge deck uses `RichardsMechanicsNotebookBridge`.

So the current benchmark-shell result is:

- loadability on the native part-1 shell is solved;
- quantitative parity on that shell is still open;
- a benchmark compare CTest would still be misleading until the constitutive
  surface is aligned more strictly or the remaining mismatch is intentionally
  scoped as reduced-law mismatch rather than native-vs-native-equivalent
  parity.

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

### WP2: Keep The Benchmark Shell Loadable

This work package is complete for benchmark-shell run stability on the tracked
reduced bridge deck.

Primary inspection points for that gap:

- [IntegrationPointData.h](../ProcessLib/RichardsMechanics/IntegrationPointData.h);
- [PressureCoupledSolidData.h](../ProcessLib/RichardsMechanics/ConstitutiveRelations/PressureCoupledSolidData.h);
- [RichardsMechanicsFEM-impl.h](../ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h).
- [MFrontRichardsMechanics.h](../MaterialLib/SolidModels/MFront/MFrontRichardsMechanics.h);

What solved it:

- pressure-consistent benchmark initialization
  `n_l0 = 0.012069019712402708`,
  `rho_lR0 = 2267.4495975433856`;
- a bracketed fallback inside the bridge microstate solve;
- damped global Newton in the benchmark bridge deck.

Regression guardrails that should stay:

- keep the benchmark bridge initial state at
  `n_l0 = 0.012069019712402708`,
  `rho_lR0 = 2267.4495975433856`;
- keep the direct material-point tests for the benchmark pressure anchor,
  `dt = 1000` first step, and the exact former process-failure state;
- keep the tracked bridge deck damping unless a later parity fix proves it is
  no longer needed;
- do not regress the reduced one-element parity slice while changing the
  benchmark shell.

Acceptance criteria:

- the benchmark-shaped native and bridge shells both run on the native part-1
  load and time scale;
- the direct bridge unit tests cover the former process-failure state and stay
  green;
- the tracked bridge restart shell stays green in the focused CTest slice.

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
5. Treat benchmark-shell loadability as closed on the tracked reduced bridge deck, and only reuse the existing binary-vs-binary compare pattern after the constitutive surface is aligned tightly enough for the result to be interpretable.

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
  --gtest_filter='MaterialLib_RichardsMechanicsNotebookBridgeMFront.PlaneStrainFactoryPathZeroDtBenchmarkPressureConsistentExactRMStateResponse:MaterialLib_RichardsMechanicsNotebookBridgeMFront.PlaneStrainFactoryPathBenchmarkPressureConsistentExactRMStateFirstStepResponse:MaterialLib_RichardsMechanicsNotebookBridgeMFront.PlaneStrainFactoryPathBenchmarkPressureConsistentProcessFailureStateResponse'
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
- the bridge keeps the microstate fallback and benchmark-shell damping needed to
  maintain that run-level boundary, unless later parity work proves they can be
  removed;
- `displacement`, `pressure`, `sigma`, `epsilon`, `saturation`, and
  `swelling_stress` match for the benchmark outputs with only justified
  machine-precision residue;
- the benchmark comparison is enforced by a dedicated CTest, not just by manual
  `vtkdiff` checks.
