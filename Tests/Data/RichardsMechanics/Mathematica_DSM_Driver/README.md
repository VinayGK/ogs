# Mathematica DSM Driver Project

This is a clean, reversible Mathematica (Wolfram Language) driver project for
the DSM micro-macro constitutive replay used in OGS transition work.

It is intended for:
- single-integration-point constitutive checks,
- microstate evolution replay (`n_l`, `rho_lR`, `phi_m`, `phi_M`),
- swelling-strain/stress updates,
- dry-density sweep studies with a vdW multiplier curve.

It is not a full FE/process replacement for OGS BEACON/ANCHORS runs.

## What is transferred here

The project includes the recent fine-tuning logic:
- finite microstate ceiling (`n_l` never unbounded),
- `n_l <= phi_total` at every step,
- `phi_m + phi_M == phi_total`,
- total porosity update by kinematics
  `phi_trial = (phi_prev + delta_eps_v)/(1 + delta_eps_v)` when valid,
- robust fallback when local Newton does not converge.

## Files

- `DsmMicromacroDriver.wl`
  - reusable package with model equations and driver functions.
- `run_dsm_driver_demo.wl`
  - runnable script that generates:
    - pressure-path history,
    - strain-coupled history,
    - dry-density sweep.
- `_outputs/`
  - generated CSV/JSON files (created by the demo script).

## Run

From this directory:

```bash
wolframscript -file run_dsm_driver_demo.wl
```

Expected outputs:

- `_outputs/overlap_history.csv`
- `_outputs/overlap_history.json`
- `_outputs/strain_coupled_history.csv`
- `_outputs/strain_coupled_history.json`
- `_outputs/dry_density_sweep.csv`
- `_outputs/dry_density_sweep.json`

## Scope and comparability

Comparable to OGS at constitutive-driver level:
- same local microstate closure form,
- same porosity bounds logic,
- same exchange-driven storage update structure,
- same swelling increment form from microstate change.

For full benchmark-level parity (`1a01`, `1b`, `1c`, EPFL/BGR, ANCHORS FE
geometry and BCs), OGS process equations and FE discretization remain required.
