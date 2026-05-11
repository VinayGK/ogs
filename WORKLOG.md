# DSM_native BEACON + ANCHORS work log

## Goal

Port the BEACON/ANCHORS calibrated benchmark from ogs-dsm_mfront to ogs-dsm_native and
extend it with the exponential augmented vdW model (K·exp(-h/λ)) to test whether the
augmentation can replace the implausibly large vdW multiplier (m_vdW ∼ 10⁷–10⁸) required
by the pure vdW model.

## Background (from benchmark note)

See `materialmodels/src/TPM/THMDSMRichardsRM_BEACON_ANCHORS_benchmark_note.tex`.

Experimental targets:
| Case          | Observable                      | Target                          |
|---------------|---------------------------------|---------------------------------|
| BEACON 1a01   | Axial swelling pressure         | 604 kPa                         |
| BEACON 1a01   | Radial swelling pressure        | 994 kPa                         |
| BEACON 1a01   | Stage-1 dry density             | 1655 kg/m³                      |
| BEACON 1b     | Dry density                     | 1520 kg/m³                      |
| ANCHORS MS33  | SP vs dry density (Villar)      | P_s = exp(6.77ρ_d − 9.07) MPa  |
| P2-1 A-B′     | Vertical confined SP            | 2.71 MPa (BGR) / 3.1–3.5 (exp) |

## Key finding from prior calibration (ogs-dsm_mfront)

The pure vdW model requires an effective Hamaker multiplier m_vdW that grows from
∼1×10⁵ at ρ_d = 1400 kg/m³ to ∼9×10⁸ at ρ_d = 1800 kg/m³ — several orders of
magnitude above the literature value A_lit = 5.1×10⁻²¹ J.

The augmented exponential model (K·exp(-h/λ)) may explain part of this range more
physically by capturing thin-film surface forces that dominate at high compaction.

## What is in this repo (DSM_native) vs what was ported

### Already present (smoke tests only):
- `Tests/Data/RichardsMechanics/beacon_1a01_dsm_micromacro_smoke.prj`
- `Tests/Data/RichardsMechanics/beacon_1b_dsm_micromacro_smoke.prj`
- `Tests/Data/RichardsMechanics/beacon_1c_dsm_micromacro_smoke.prj`
- `Tests/Data/RichardsMechanics/beacon_1a01_dsm_micromacro_inflow.prj` (exploratory)
- `Tests/Data/RichardsMechanics/beacon_1a01_dsm_micromacro_stressprobe.prj`
- `Tests/Data/RichardsMechanics/beacon_1a01_dsm_micromacro_inflowprobe.prj`

### Ported / created in this session:
- `Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/`
  - `run_villar_dense_dd_native_purevdw_calibration.py` – pure vdW, LinearElastic
  - `run_villar_dense_dd_native_augmented_calibration.py` – augmented vdW, LinearElastic
  - `compare_purevdw_vs_augmented.py` – comparison plots
- `Tests/Data/RichardsMechanics/beacon_1a01_dsm_micromacro_calibrated_inflow.prj` – calibrated pure vdW
- `Tests/Data/RichardsMechanics/beacon_1a01_dsm_micromacro_augmented_inflow.prj` – augmented vdW

## OGS binary

Binary: `/Users/vinaykumar/git/build/ogs-worktrees/build/dsm_native-release/bin/ogs`
Library path: `/Users/vinaykumar/git/build/ogs-worktrees/build/dsm_native-release/lib`
Version: `vdw-baseline-2026-05-08`
MFront: OFF (not available in this build)

## Important adaptation: LinearElastic vs MFront MCC

The mfront calibration uses ModCamClay_semiExpl_constE with InitialPreConsolidationPressure=1e10 Pa.
Even at this very high p_c0 (no plastic yielding during the confined hydration test),
the MCC *elastic* tangent bulk modulus is K_e = ((1+e)/κ)·p' which depends on stress.
At p' ≈ 0.6 MPa (BEACON 1a01 stage 1) K_e ≈ 186 MPa; at p' ≈ 9 MPa (Villar target)
K_e ≈ 2 GPa. The LinearElasticIsotropic stand-in uses K_bulk = E/(3(1-2ν)) = 43.33 MPa
(E=52 MPa, ν=0.3), a *fixed* and much smaller value.

Consequences:
- LinearElastic pressure cap: p_sw,max ≈ K_bulk · α_sw · φ₀ ≈ 1.73 MPa.
- Villar calibration valid only for ρ_d ≲ 1418 kg/m³ (where Villar target < cap).
- BEACON 1a01 native (LinearElastic) gives ≈ 890 kPa, MFront MCC gives 603.7 kPa, target 604 kPa.
  The 890/604 gap reflects the bulk-modulus difference and is *not* a port error.

To get exactly 604 kPa with LinearElastic, A_eff would need a separate native-specific calibration.

## Workflow to run calibrations

```bash
OGS=/Users/vinaykumar/git/build/ogs-worktrees/build/dsm_native-release/bin/ogs
DYLD_LIBRARY_PATH=/Users/vinaykumar/git/build/ogs-worktrees/build/dsm_native-release/lib
export DYLD_LIBRARY_PATH

cd Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI

# Pure vdW calibration — restrict to ρ_d ≤ 1418 kg/m³ (LinearElastic cap)
python3 run_villar_dense_dd_native_purevdw_calibration.py --ogs-bin $OGS --dd-min 1400 --dd-max 1418

# Augmented calibration at physical 1 nm decay (default; see lambda unit note below)
python3 run_villar_dense_dd_native_augmented_calibration.py --ogs-bin $OGS --lam 1e-6

# Quick demo (3 points, low-density window only)
python3 run_villar_dense_dd_native_purevdw_calibration.py --ogs-bin $OGS --dd-min 1380 --dd-max 1418 --dd-step 20
python3 run_villar_dense_dd_native_augmented_calibration.py --ogs-bin $OGS --lam 1e-6 --dd-min 1380 --dd-max 1418 --dd-step 20

# Comparison plot (after both calibrations complete)
python3 compare_purevdw_vs_augmented.py
```

## Status

- [x] ANCHORS_MS33_ModelI directory created
- [x] Pure vdW calibration script written
- [x] Augmented calibration script written
- [x] Comparison script written
- [x] Calibrated BEACON 1a01 project file written
- [x] Augmented BEACON 1a01 project file written (placeholder K; update after calibration)
- [ ] Pure vdW calibration run completed (17 points)
- [ ] Augmented calibration run completed (17 points)
- [ ] BEACON 1a01 calibrated inflow run and compared to 604 kPa target

## Physical parameters used (ANCHORS calibration cell)

- Grain density ρ_s = 2780 kg/m³
- Specific surface Sa = 523 m²/g (from Villar/ANCHORS MS33)
- Literature Hamaker A_lit = 5.1×10⁻²¹ J
- Suction split: ψ_tot = −100 MPa → ψ_macro = −1 MPa + ψ_micro = −99 MPa
- Initial macro pressure p_ic = −1 MPa, released to 0 over 120 days
- MS33 permeability: Kozeny-Carman with k₀=5.6×10⁻²¹ m², φ_ref=0.42
- Saturation: Tuller model (char. pore = 10 μm, γ = 0.0715 N/m)
- Confinement: all boundaries fixed (constant volume)
- Time horizon: 120 days

## Parameters used for BEACON 1a01 calibrated project

The Villar calibration (9.22 MPa target) and the BEACON-specific calibration
(604 kPa short-time target) give *different* A_eff. Do not interchange them.

BEACON-specific calibration (from mfront bridge report,
`beacon_1a01_axial_fit_native_bridge_report.json`):
- hamaker_constant = 1.9606596248030517×10⁻¹⁵ J  (m_vdW ≈ 384 — much smaller than Villar)
- n_l0 = 1.53385355×10⁻³ (from BEACON mfront fit)
- t_end = 1×10⁵ s = 27.8 h (matches mfront bridge fit duration)
- specific_surface = 523
- micro_solid_volume_fraction_reference = 0.6
- local_nonlinear_solve_mode = scalar_micro_macro_mass_storage_mode
- micro_liquid_density_reference = 1e-6, _a = 1e-16, _b = 1.0

## Lambda unit convention (augmented model)

`specific_surface` is in m²/g (not converted to SI m²/kg in the code).
Therefore the code-internal film thickness `h_code = n_l/(nS·ρ_SR·Sa)`
is 1000× the physical SI value, and:

    λ_code = 1000 · λ_physical

Use λ_code = 1e-6 to represent a physical 1 nm decay length (default).
λ_code = 1e-9 (previous default) represented physical 1 pm — sub-atomic.

## n_l0 convention

Always use HAMAKER_LITERATURE (A_lit) to compute n_l0 from the micro-suction
split, *never* A_eff. This matches the mfront calibration reference
(`n_l0_from_micro_suction(case.phi0, HAMAKER_LITERATURE)`). Using A_eff
gives a grossly wrong initial state (n_l0 ≈ 0.087 at ρ_d=1500 with Villar
A_eff) and leaves no swelling headroom.
