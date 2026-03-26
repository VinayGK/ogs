# RichardsMechanics DSM Native-vs-MFront Parity Plan

## Goal

Build an honest native-vs-MFront comparison for the RichardsMechanics
double-structure model (DSM) inside the same CTest slice.

The finished state should contain:

- one native DSM run
- one MFront bridge DSM run
- one native-vs-bridge comparison test on the same output fields and time steps
- zero tolerance where that is realistic, and only machine-scale tolerance
  where exact bitwise equality is not realistic

## Short version

Two things are already true:

- the reduced one-element native-vs-bridge compare works and is green
- the benchmark bridge shell now runs on the native part-1 load/time scale

One thing is still open:

- the benchmark native-vs-bridge fields still differ at `ts_1`

So the remaining problem is no longer “can the benchmark bridge run?”.
The remaining problem is “why do the native and bridge benchmark fields still
separate after the first real step?”.

## What is already verified

As of 2026-03-26 on branch `dsm-nb-mfront-transition`:

- the repository configures and builds with `OGS_USE_MFRONT=ON`
- this focused test slice passes:
  - `ogs-RichardsMechanics/mfront_restart_part1`
  - `ogs-RichardsMechanics/mfront_restart_part2`
  - `ogs-RichardsMechanics_mfront_restart_part1_rm_bridge`
  - `ogs-RichardsMechanics_mfront_parity_1element_native`
  - `ogs-RichardsMechanics_mfront_parity_1element_bridge`
  - `ogs-RichardsMechanics_mfront_parity_1element_compare`
- the RM pressure-coupled carrier now keeps the bridge saturation derivative
  with respect to strain `dS/d\varepsilon`, and the pressure-equation
  displacement block now consumes that term instead of dropping it

## Reduced one-element parity status

The reduced one-element pair is the current parity gate.

Compared fields:

- `displacement`
- `pressure`
- `sigma`
- `epsilon`
- `saturation`
- `swelling_stress`

Result:

- `displacement`: exact
- `pressure`: exact
- `epsilon`: exact
- `saturation`: exact
- `swelling_stress`: exact
- `sigma`: machine-scale residue only

Observed `sigma` residue:

- absolute max norm about `1.8e-12`
- relative max norm about `3.6e-16`

So the reduced one-element compare is closed for practical purposes.
The remaining `sigma` difference is a floating-point issue, not a model-parity
issue.

## Benchmark-shell status

### What is solved

The benchmark bridge shell now runs on the native part-1 load/time scale.

That run-level result currently depends on three things:

- pressure-consistent bridge initial microstate
  - `n_l0 = 0.012069019712402708`
  - `rho_lR0 = 2267.4495975433856`
- a bracketed fallback in the bridge microstate solve
- damped global Newton in the benchmark bridge project

Direct bridge tests now also cover:

- `dt = 0` at benchmark pressure with the raw benchmark stress state
- `dt = 0` at benchmark pressure with the RM-equivalent effective stress state
- `dt = 1000` at the benchmark first-step anchor
- the exact former process-failure state

### What is still open

A clean native-vs-bridge benchmark rerun now shows:

- at `ts_0`
  - `displacement`, `pressure`, `epsilon`, `saturation`, and
    `swelling_stress` match exactly
  - `sigma` differs only at machine scale
- at `ts_1`
  - the mismatch is still material

Representative `ts_1` differences:

- `|Δu_y|_∞ ≈ 1.51e-4`
- `|Δp|_∞ ≈ 6.80e3`
- `|Δσ|_∞ ≈ 9.21e3`
- `|Δε_yy|_∞ ≈ 1.51e-4`
- `|ΔS_L|_∞ ≈ 1.87e-1`

So benchmark-shell loadability is closed, but benchmark-shell parity is not.

One useful negative result is now also known:

- the missing bridge saturation-strain tangent was a real carrier omission, and
  it is now fixed
- but the benchmark `ts_1` field differences stay at essentially the same scale
  after that fix

So the remaining benchmark gap is not explained by that dropped tangent alone.

## Why the benchmark gap is still hard to interpret

The benchmark decks still do not use the same constitutive law:

- native benchmark law:
  `ModCamClay_semiExpl_constE`
- bridge benchmark law:
  `RichardsMechanicsNotebookBridge`

That means the current benchmark mismatch is not yet a pure same-model
interface signal.

It may still contain two effects mixed together:

- a real RM pressure-coupled interface issue
- a real constitutive difference between the native and bridge laws

## The pressure-coupled contract in plain language

The key question is simple:

When the bridge returns stress and pressure derivatives, does RM interpret them
the same way the bridge means them?

### Continuum form

The quasi-static momentum balance is

\[
\nabla \cdot \boldsymbol{\sigma}_{tot} + \rho \mathbf{b} = \mathbf{0}.
\]

The liquid mass balance is

\[
\dot{m}_L + \nabla \cdot \mathbf{w}_L = q_L.
\]

In RichardsMechanics, total stress is built from effective stress by

\[
\boldsymbol{\sigma}_{tot}
=
\boldsymbol{\sigma}_{eff}
- \alpha \chi(S_L) p_L \mathbf{I}.
\]

So if the bridge returns effective stress, RM should add the pore-pressure
contribution once.

If the bridge already returns total stress, but RM still applies the same
pressure correction, the pore-pressure term is counted twice.

### Weak form

The weak form of momentum is

\[
\int_\Omega \boldsymbol{\varepsilon}(\delta \mathbf{u}) :
\boldsymbol{\sigma}_{tot}\, d\Omega
- \int_\Omega \delta \mathbf{u}\cdot \rho \mathbf{b}\, d\Omega
- \int_{\Gamma_t} \delta \mathbf{u}\cdot \bar{\mathbf{t}}\, d\Gamma = 0.
\]

The weak form of liquid mass balance is

\[
\int_\Omega \delta p \, \dot{m}_L \, d\Omega
+ \int_\Omega \nabla \delta p \cdot \mathbf{w}_L \, d\Omega
- \int_\Omega \delta p \, q_L \, d\Omega = 0.
\]

The important point is that the momentum equation depends on pressure through
`\sigma_tot`.
So a wrong stress meaning or a wrong pressure derivative changes the coupled
`u-p` Jacobian.

### Semi-discrete form

With standard FE interpolation,

\[
\mathbf{u}_h = N_u \mathbf{d},
\qquad
p_h = N_p \mathbf{p},
\qquad
\boldsymbol{\varepsilon}(\mathbf{u}_h) = B \mathbf{d}.
\]

The momentum residual is

\[
R_u(\mathbf{d},\mathbf{p})
=
\int_\Omega B^T \boldsymbol{\sigma}_{tot}\, d\Omega - f_{ext}.
\]

The pressure-to-momentum Jacobian block is

\[
K_{up}
=
\frac{\partial R_u}{\partial \mathbf{p}}
=
\int_\Omega B^T
\left(
\frac{\partial \boldsymbol{\sigma}_{tot}}{\partial p_L}
\right)
N_p\, d\Omega.
\]

The pressure equation also has a displacement-coupling block. If the liquid
storage contains a term like
\[
m_L \supset \phi \rho_{LR} S_L,
\]
then linearizing with respect to displacement gives the saturation-strain part
\[
K_{pu}^{(S)}
=
\frac{\partial R_p}{\partial \mathbf{d}}
\supset
\int_\Omega N_p^T \phi \rho_{LR}
\left(
\frac{\partial S_L}{\partial \boldsymbol{\varepsilon}}
\right)
B\, d\Omega.
\]

That first `dS/d\varepsilon` term was a real RM-side carrier omission. It is
now assembled.

This gives one useful negative result:

- the benchmark `ts_1` gap stays at essentially the same scale after restoring
  that term

So the remaining benchmark gap can still come from a wrong meaning of:

- returned stress
- `dS/d\varepsilon`
- `dSigma/dp_L`
- `dS_L/dp_L`

even if the local bridge solve itself is finite and stable.

## What the current parity gate does and does not prove

### What it proves

The current one-element compare proves that the native and bridge paths agree
on the reduced shared shell that is actually being tested.

### What it does not prove

It does not prove:

- benchmark-shell parity
- exact same-parameter MCC plasticity parity
- full notebook-to-native constitutive equivalence

## Parameter status

For a parameter-by-parameter comparison of the reduced one-element pair, see:

- [RichardsMechanicsDSMBridgeParameterComparison.md](./RichardsMechanicsDSMBridgeParameterComparison.md)

Short version:

- the shared OGS setup is aligned well enough for the reduced parity CTest
- the active constitutive parameter sets are not literally the same
- the current green compare test is not an exact same-parameter MCC plasticity
  proof

## Work packages

### WP1: Keep the reduced one-element parity slice green

Status:

- complete

Acceptance criteria:

- native one-element run stays green
- bridge one-element run stays green
- compare CTest stays green
- only machine-scale tolerance is used where exact equality is not realistic

### WP2: Keep the benchmark shell loadable

Status:

- complete for run-level loadability

Guardrails:

- keep the pressure-consistent benchmark bridge IC
- keep the bridge microstate fallback
- keep the benchmark damping unless tests show it is safe to remove
- keep the direct benchmark-pressure bridge unit tests green

### WP3: Close the benchmark quantitative gap

Status:

- open

Practical task:

- compare clean native and bridge benchmark outputs
- separate constitutive differences from RM interface differences
- only then tighten the interpretation of the remaining `ts_1` mismatch

### WP4: Add the benchmark comparison CTest

Status:

- open

Do this only after WP3 is honest enough to interpret.

Recommended shape:

- keep the two run tests
- add a separate compare test
- compare:
  - `displacement`
  - `pressure`
  - `sigma`
  - `epsilon`
  - `saturation`
  - `swelling_stress`

## Definition of done

### Reduced one-element task

Done means:

- the native and bridge one-element shells both run
- the compare CTest is green
- all compared fields match exactly except for justified machine-scale residue

### Benchmark task

Done means:

- the native and bridge benchmark shells both run on the native part-1 scale
- the benchmark fields match up to only justified machine-scale residue
- the benchmark comparison is enforced by a dedicated CTest

## Useful commands

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

Focused CTest slice:

```bash
ctest --test-dir /Users/vinaykumar/git/build/release-mfront-tpm \
  --output-on-failure \
  -R 'ogs-RichardsMechanics/mfront_restart_part1$|ogs-RichardsMechanics/mfront_restart_part2$|ogs-RichardsMechanics_mfront_restart_part1_rm_bridge$|ogs-RichardsMechanics_mfront_parity_1element_native$|ogs-RichardsMechanics_mfront_parity_1element_bridge$|ogs-RichardsMechanics_mfront_parity_1element_compare$'
```

Direct benchmark-pressure bridge guardrail:

```bash
/Users/vinaykumar/git/build/release-mfront-tpm/bin/testrunner \
  --gtest_filter='MaterialLib_RichardsMechanicsNotebookBridgeMFront.PlaneStrainFactoryPathZeroDtBenchmarkPressureConsistentExactRMStateResponse:MaterialLib_RichardsMechanicsNotebookBridgeMFront.PlaneStrainFactoryPathBenchmarkPressureConsistentExactRMStateFirstStepResponse:MaterialLib_RichardsMechanicsNotebookBridgeMFront.PlaneStrainFactoryPathBenchmarkPressureConsistentProcessFailureStateResponse'
```
