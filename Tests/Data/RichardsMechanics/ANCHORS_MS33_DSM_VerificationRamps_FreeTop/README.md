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

## Verified outcome (2026-05-26, after active_nS fix)

All three sims complete 120 steps without rejected steps. dd1400 example:

| t (d) | p (MPa) | n_l   | φ_M   | S_L  | σ_mean (MPa) |
|-------|---------|-------|-------|------|--------------|
|  0    | −100    | 0.248 | 0.330 | 0.33 | −33.97       |
| 25    | −58.3   | 0.283 | 0.528 | 0.48 | −13.00       |
| 60    | +0.1    | 0.670 | 0.000 | 1.00 | −13.58 (peak)|
| 85    | −41.6   | 0.306 | 0.513 | 0.59 | −13.04       |
| 115   | −91.7   | 0.253 | 0.546 | 0.35 | −12.95       |

The free-top dd1400 column starts under the Biot-balancing initial
effective stress `sigma0_eff = alpha * chi(p_L_init) * p_L_init` 
(≈ −34 MPa with chi ≈ 0.33 at p = −100 MPa). The hydration arc raises
n_l from 0.25 to 0.67 (= phi(t) at peak) and lifts the mean effective
stress to ≈ −13.6 MPa. The `φ_M → 0` collapse at peak saturation when
`n_l = φ(t)` is the rigid hierarchical-split kinematic limit documented
in `tex/dsm-bgr-paper/draft/paper_DSM.tex` (kinematic-scope section).

**Note (2026-05-26):** the pre-2026-05-26 numbers in earlier revisions
of this table (n_l=0.210 at t=0, σ=−2.15 MPa at peak, etc.) were
produced before the `active_nS` physics fix (native commit
`8192021299`, mfront commit `340b928856`) which corrected the
denominator of `omega_l` from `(1 − phi_M)` to `(1 − n_l)`. The new
values in the table above reflect the post-fix physics.
