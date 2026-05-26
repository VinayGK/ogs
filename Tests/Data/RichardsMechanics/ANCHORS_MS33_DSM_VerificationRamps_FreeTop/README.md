# DSM Verification Ramps — FREE-TOP single-element saturation/desaturation

Companion to the confined `ANCHORS_MS33_DSM_VerificationRamps/` suite.
Same triangle pressure ramp at the three Model-I-calibrated dry densities,
same native DSM hierarchical binary, **but with the top displacement
boundary released** so the column free-swells in y while remaining
laterally confined and bottom-supported. The initial stress is set to
the Biot-balancing isotropic effective stress
sigma0_eff = alpha * chi(p_L_init) * p_L_init (~ -33 MPa at p=-100 MPa,
chi ~ 0.33) so the column starts in mechanical equilibrium under the
high initial suction.

Used to visualise the kinematic-partition behaviour of the rigid
hierarchical split that is invisible in the confined runs: total
porosity phi(t) now evolves nontrivially, the algebraic ceiling
n_l <= phi(t) rises with mechanical dilation, and phi_M does not
remain constant under free swelling (the rigid algebraic split forces
dot(phi) to distribute into both phi_M and phi_m). The cancellation
identity in paper_DSM.tex S2.2 guarantees that mass conservation is
preserved.

## Setup

- **Template PRJ**: `../ANCHORS_MS33_StrictParity/ms33_dsm_parity_native.prj`
  (DSM hierarchical path active; `<potential_exchange>`, vdW Hamaker, mass
  exchange wired; `micro_*` secondary variables in the output list).
- **Per-density override**: `phi0` and `IntrinsicPermeability0` substituted
  from `../ANCHORS_MS33_ModelI/ms33_model_i_dd{1400,1600,1800}.prj`.
- **Ramp**: pressure BC triangle `-100 MPa → +0.1 MPa → -100 MPa` over
  120 days (dt = 1 day, 120 steps). Shape matches the notebook (saturation
  half then desaturation half); range covers the suction regime over which
  the Model-I calibration is active.
- **Mesh**: reuses `../square_1x1_quad_1e0.vtu` (single axisymmetric quad).
- **Binary**: native DSM hierarchical
  `~/git/build/native-release-omp-sharedcache/bin/ogs`
  (worktree `dsm_native_hierarchical`).

## Run

```bash
python3 run_verification_ramps.py             # generate, run, post-process, plot
python3 run_verification_ramps.py --only-gen  # regenerate PRJ files only
python3 run_verification_ramps.py --skip-run  # skip OGS, re-post-process from existing VTU
```

## Outputs

For each density:

| File                                       | Content                              |
|--------------------------------------------|--------------------------------------|
| `dsm_ramp_dd{dd}.prj`                      | generated project file                |
| `dsm_ramp_dd{dd}.pvd`, `_ts_*.vtu`         | OGS simulation results (not tracked) |
| `dsm_ramp_dd{dd}.csv`                      | cell-averaged scalar time series      |
| `dsm_ramp_dd{dd}_content_vs_time.png`      | water-content components vs time      |
| `dsm_ramp_dd{dd}_content_vs_potential.png` | water content vs suction              |
| `dsm_ramp_dd{dd}_hysteresis.png`           | macro saturation hysteresis loop      |
| `dsm_ramp_dd{dd}_stress_vs_time.png`       | mean stress vs time (swelling trace)  |

## Verified outcome (2026-05-25)

All three sims complete 120 steps without rejected steps. dd1400 example:

| t (d) | p (MPa) | n_l   | φ_M   | S_L  | σ_mean (MPa) |
|-------|---------|-------|-------|------|--------------|
|  0    | −100    | 0.210 | 0.362 | 0.33 | −0.91        |
| 25    | −58.3   | 0.274 | 0.307 | 0.48 | −1.18        |
| 55    | −8.2    | 0.496 | 0.000 | 0.95 | −2.15 (peak) |
| 85    | −41.6   | 0.334 | 0.244 | 0.59 | −1.44        |
| 115   | −91.7   | 0.219 | 0.355 | 0.35 | −0.94        |

The `φ_M → 0` collapse at peak saturation when `n_l = φ_0` is the rigid
hierarchical-split kinematic limit documented in
`tex/dsm-bgr-paper/draft/paper_DSM.tex` (kinematic-scope section).
