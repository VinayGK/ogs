# DSM element Newton tangent against a central finite difference of the assembled residual (2026-10-13)

Status words: **computed** (a run of the test or script named here), **derived** (follows from computed
numbers by the step shown), **open**. No statement is made about convergence rates or about the step
counts of the production runs.

Branch `dsm_z_weight_2026-10_jac5`, cut from `dsm_z_weight_2026-10` (tip cc5e71e12a). No existing source
file, deck or reference frame is modified. Material:

| file | content |
|---|---|
| `Tests/ProcessLib/RichardsMechanics/DSMElementJacobianFD.cpp` | Google test `RichardsMechanics.DSMElementJacobianFD` (the harness, the FD, a regression guard) |
| `Tests/Python/dsm_jacobian_fd_run.sh` | driver: runs the test with the dumps on and writes the tables |
| `Tests/Python/dsm_jacobian_fd_analysis.py` | tables and the attribution fits from the dumps |
| `docs/dsm_jacobian_fd_verification_2026-10_data/` | the dumps of the run reported here: `results.csv` (every block x state x switch set x step), `results.csv.state.csv` (integration-point state and its central differences), `mat/` (analytic and FD element matrices, 12 x 12, row-major CSV), `tables.md` (all tables, also those for h = 1e-2 and 1e-7), `testrunner.log` |

Build: Release, `-DOGS_BUILD_TESTING=ON -DOGS_USE_MFRONT=OFF -DOGS_USE_PIP=OFF`, gcc 13.3, build directory
outside the source tree. Three dependencies could not be fetched through the network policy of the build
machine and were taken from mirrors: tclap (github.com/mirror/tclap, branch 1.4), xmlpatch
(github.com/ufz/xmlpatch, v0.4.2), xdmf (Ubuntu source `xdmf_3.0+git20190531`, library targets renamed
to `OgsXdmf`/`OgsXdmfCore`). None of them is on the path of this test (no XDMF output, no CLI parsing).
`build/bin/testrunner --gtest_filter='RichardsMechanics.DSM*'`: the 39 tests present on the branch
PASS (computed). The new test runs in about 1.2 s.

## 1. Method (computed)

**Harness.** The deck `Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj` is read
by the regular project reader (`Simulation::initializeDataStructures`, the code path of `bin/ogs`), so the
element (the 1 cm axisymmetric Quad4 `square_1e-2_quad_1e0.vtu`, 2 x 2 Gauss points), the material set
and every `potential_exchange` switch of the deck are used as shipped: `micro_ceiling_treatment = kkt`,
`micro_ceiling_pu_tangent = kkt_active`, `micro_ceiling_sw_tangent = kkt_active`,
`micro_ceiling_saturation_gate = bishop_relperm`, `macro_balance_drops_micro_biot_term = true`,
`darcy_relative_permeability_mobility = kirchhoff_element_mean`, `macro_storage_uses_macro_porosity =
true`, `macro_storage_exact_time_levels = true`, `micro_mass_strain_term_eulerian = true`,
`swelling_stress_k_level = true`, linear elastic skeleton (E = 52 MPa, nu = 0.3), Biot b = 1. The
switch-set variants are patched copies of the deck written to a temporary directory:

| switch set | change against the deck |
|---|---|
| `deck_step_klevel` | none (step form of the swelling stress, K-level) |
| `level` | `+ swelling_stress_form = level` |
| `level_z` | `+ swelling_stress_form = level`, `+ micro_potential_partner_z_weight = true` |
| `diag_deck_pu_all_exchange`, `diag_level_pu_all_exchange`, `diag_level_z_pu_all_exchange` | the three above with `micro_ceiling_pu_tangent = all_exchange` instead of `kkt_active`; attribution only, not a shipped setting |

**States.** The deck's own time loop drives the element (all four pressure nodes are Dirichlet nodes of
the suction ramp -100 MPa -> 0 over 20 d; all displacements are prescribed). At each state the local
assembler of the element is called directly with x_prev = the last accepted solution and x = the first
Newton iterate of the next step (x_prev with the pressure Dirichlet value at t + dt), dt = the next time
step of the deck's stepper:

| state | t_n [s] | dt [s] | p at x [MPa] | n_l | phi_M | S_L (macro) | KKT-active IPs | eps_v |
|---|---|---|---|---|---|---|---|---|
| `dry` | 0 | 10 | -100.0 | 0.00118 | 0.424 | 1.8e-8 | 0/4 | 0 |
| `mid` (micro level filling) | 864000 | 76114 | -45.6 | 0.414 | 0.0178 | 8.4e-8 | 0/4 | 0 |
| `mid_eps_step` | as `mid`, x carries in addition u_z = -0.01 (z - z_min): a 1 % compressive volumetric strain increment in this step | | | 0.413 | 0.0088 | | 0/4 | -0.01 |
| `mid_eps_hist` | 864000, deck run with top u_z = -1e-4 m from t = 0 (eps_v = -0.01 in x and x_prev) | 76114 | -45.6 | 0.411 | 0.0127 | | 0/4 | -0.01 |
| `wet_ceiling` (on the micro ceiling) | 8640000 | 86400 | 0 | 0.4245 | 0 | 1 | 4/4 | 0 |

(values of `deck_step_klevel`; the strained states of the level and Z sets differ in n_l and phi_M,
see `results.csv.state.csv`.)

**FD.** J_fd(:, j) = [R(x + h_j e_j) - R(x - h_j e_j)] / (2 h_j), R = -local_rhs, for all 12 DOFs (4 p,
4 u_r, 4 u_z). Before every evaluation the integration-point states (`current_states_`, `output_data_`)
are restored to the snapshot S0 taken before the analytic evaluation, i.e. the FD is that of the residual
R(x; S0) at the fixed lagged state the Newton iteration starts from (the convention of the in-assembler
`micro_ceiling_fd_check`). Steps: `abs` h_p = 1 Pa, h_u = 1e-9 m; `rel1e-2` and `rel1e-7` h = 1e-2 |x_j|
resp. 1e-7 |x_j|, with |x_j| replaced by 1e8 Pa (p) or 1e-4 m (u) where x_j = 0; `abs+relag5` = abs with
the perturbed residual evaluated six times in sequence from S0 (lagged values re-converged at the
perturbed point).

**Error measures** per block B (pp 4x4, pu 4x8, up 8x4, uu 8x8): max abs = max |J_an - J_fd|; max rel =
max abs / max(max |J_an|, max |J_fd|) (block-normalised); ‖·‖ rel = ‖J_an - J_fd‖_F / ‖J_fd‖_F with the p
rows and the u rows each scaled by their largest |J_fd| entry. Units: pp [kg/(s Pa)], pu [kg/(s m)],
up [m^2], uu [N/m] (residual of the axisymmetric element, 2 pi r included).

## 2. Error table (computed; h = abs, fixed lag)

| switch set | state | kkt IPs | pp abs | pp rel | pu abs | pu rel | up abs | up rel | uu abs | uu rel | ‖·‖ rel |
|---|---|---|---|---|---|---|---|---|---|---|---|
| deck_step_klevel | dry | 0 | 5.4e-22 | 1.0e-05 | 1.4e-06 | 1.0e+00 | 1.4e-08 | 1.0e+00 | 4.3e+02 | 1.8e-04 | 6.6e-01 |
| deck_step_klevel | mid | 0 | 6.9e-25 | 1.7e-08 | 1.9e-07 | 1.0e+00 | 2.9e-06 | 1.0e+00 | 2.5e+03 | 1.0e-03 | 6.6e-01 |
| deck_step_klevel | mid_eps_step | 0 | 8.3e-25 | 2.1e-08 | 1.9e-07 | 1.0e+00 | 3.0e-06 | 1.0e+00 | 2.4e+03 | 1.0e-03 | 6.6e-01 |
| deck_step_klevel | wet_ceiling | 4 | 0.0e+00 | 0.0e+00 | 1.4e-16 | 8.0e-11 | 4.1e-13 | 2.9e-09 | 1.3e-03 | 5.5e-10 | 3.6e-10 |
| deck_step_klevel | mid_eps_hist | 0 | 5.4e-25 | 1.3e-08 | 1.9e-07 | 1.0e+00 | 2.9e-06 | 1.0e+00 | 2.3e+03 | 9.4e-04 | 6.6e-01 |
| level | dry | 0 | 5.4e-22 | 1.0e-05 | 1.4e-06 | 1.0e+00 | 1.4e-08 | 1.0e+00 | 4.5e+02 | 1.9e-04 | 6.6e-01 |
| level | mid | 0 | 6.9e-25 | 1.7e-08 | 1.9e-07 | 1.0e+00 | 9.8e-06 | 1.0e+00 | 5.4e+04 | 2.3e-02 | 6.6e-01 |
| level | mid_eps_step | 0 | 8.3e-25 | 2.1e-08 | 1.9e-07 | 1.0e+00 | 9.7e-06 | 1.0e+00 | 4.6e+05 | 2.0e-01 | 6.8e-01 |
| level | wet_ceiling | 4 | 0.0e+00 | 0.0e+00 | 1.4e-16 | 8.0e-11 | 4.1e-13 | 2.9e-09 | 1.3e-03 | 5.9e-10 | 4.8e-10 |
| level | mid_eps_hist | 0 | 5.4e-25 | 1.3e-08 | 1.9e-07 | 1.0e+00 | 9.8e-06 | 1.0e+00 | 5.4e+04 | 2.2e-02 | 6.6e-01 |
| level_z | dry | 0 | 5.4e-22 | 1.0e-05 | 1.4e-06 | 1.0e+00 | 1.4e-08 | 1.0e+00 | 4.5e+02 | 1.9e-04 | 6.6e-01 |
| level_z | mid | 0 | 6.9e-25 | 1.7e-08 | 1.9e-08 | 1.0e+00 | 9.8e-06 | 1.0e+00 | 7.5e+03 | 3.1e-03 | 6.6e-01 |
| level_z | mid_eps_step | 0 | 2.2e-20 | 5.5e-04 | 1.4e-08 | 1.0e+00 | 9.7e-06 | 1.0e+00 | 5.1e+05 | 2.1e-01 | 6.8e-01 |
| level_z | wet_ceiling | 4 | 0.0e+00 | 0.0e+00 | 1.4e-16 | 8.0e-11 | 4.1e-13 | 2.9e-09 | 1.3e-03 | 5.9e-10 | 4.8e-10 |
| level_z | mid_eps_hist | 0 | 2.2e-20 | 5.5e-04 | 1.6e-08 | 1.0e+00 | 9.7e-06 | 1.0e+00 | 6.9e+03 | 2.9e-03 | 6.6e-01 |

Block magnitudes max |J_fd| (all sets): pp 4e-17 to 1.5e-16, pu 1.4e-8 to 1.8e-6, up 1.4e-8 (dry),
3e-6 to 1e-5 (mid), 1.4e-4 (wet), uu 2.1e6 to 2.4e6 (`tables.md`).

Attribution runs (`micro_ceiling_pu_tangent = all_exchange`, computed):

| switch set | state | kkt IPs | pp abs | pp rel | pu abs | pu rel | up abs | up rel | uu abs | uu rel | ‖·‖ rel |
|---|---|---|---|---|---|---|---|---|---|---|---|
| diag_deck_pu_all_exchange | dry | 0 | 5.4e-22 | 1.0e-05 | 2.6e-09 | 1.8e-03 | 1.4e-08 | 1.0e+00 | 4.3e+02 | 1.8e-04 | 1.2e-03 |
| diag_deck_pu_all_exchange | mid | 0 | 6.9e-25 | 1.7e-08 | 6.8e-08 | 2.6e-01 | 2.9e-06 | 1.0e+00 | 2.5e+03 | 1.0e-03 | 2.3e-01 |
| diag_deck_pu_all_exchange | mid_eps_step | 0 | 8.3e-25 | 2.1e-08 | 6.6e-08 | 2.6e-01 | 3.0e-06 | 1.0e+00 | 2.4e+03 | 1.0e-03 | 2.3e-01 |
| diag_deck_pu_all_exchange | wet_ceiling | 4 | 0.0e+00 | 0.0e+00 | 1.4e-16 | 8.0e-11 | 4.1e-13 | 2.9e-09 | 1.3e-03 | 5.5e-10 | 3.6e-10 |
| diag_deck_pu_all_exchange | mid_eps_hist | 0 | 5.4e-25 | 1.3e-08 | 6.7e-08 | 2.6e-01 | 2.9e-06 | 1.0e+00 | 2.3e+03 | 9.4e-04 | 2.3e-01 |
| diag_level_pu_all_exchange | dry | 0 | 5.4e-22 | 1.0e-05 | 2.6e-09 | 1.8e-03 | 1.4e-08 | 1.0e+00 | 4.5e+02 | 1.9e-04 | 1.2e-03 |
| diag_level_pu_all_exchange | mid | 0 | 6.9e-25 | 1.7e-08 | 6.8e-08 | 2.6e-01 | 9.8e-06 | 1.0e+00 | 5.4e+04 | 2.3e-02 | 2.4e-01 |
| diag_level_pu_all_exchange | mid_eps_step | 0 | 8.3e-25 | 2.1e-08 | 6.6e-08 | 2.6e-01 | 9.7e-06 | 1.0e+00 | 4.6e+05 | 2.0e-01 | 2.8e-01 |
| diag_level_pu_all_exchange | wet_ceiling | 4 | 0.0e+00 | 0.0e+00 | 1.4e-16 | 8.0e-11 | 4.1e-13 | 2.9e-09 | 1.3e-03 | 5.9e-10 | 4.8e-10 |
| diag_level_pu_all_exchange | mid_eps_hist | 0 | 5.4e-25 | 1.3e-08 | 6.7e-08 | 2.6e-01 | 9.8e-06 | 1.0e+00 | 5.4e+04 | 2.2e-02 | 2.3e-01 |
| diag_level_z_pu_all_exchange | dry | 0 | 5.4e-22 | 1.0e-05 | 2.6e-09 | 1.8e-03 | 1.4e-08 | 1.0e+00 | 4.5e+02 | 1.9e-04 | 1.2e-03 |
| diag_level_z_pu_all_exchange | mid | 0 | 6.9e-25 | 1.7e-08 | 9.3e-09 | 3.3e-01 | 9.8e-06 | 1.0e+00 | 7.5e+03 | 3.1e-03 | 3.3e-01 |
| diag_level_z_pu_all_exchange | mid_eps_step | 0 | 2.2e-20 | 5.5e-04 | 9.7e-09 | 4.0e-01 | 9.7e-06 | 1.0e+00 | 5.1e+05 | 2.1e-01 | 4.7e-01 |
| diag_level_z_pu_all_exchange | wet_ceiling | 4 | 0.0e+00 | 0.0e+00 | 1.4e-16 | 8.0e-11 | 4.1e-13 | 2.9e-09 | 1.3e-03 | 5.9e-10 | 4.8e-10 |
| diag_level_z_pu_all_exchange | mid_eps_hist | 0 | 2.2e-20 | 5.5e-04 | 9.9e-09 | 3.8e-01 | 9.7e-06 | 1.0e+00 | 6.9e+03 | 2.9e-03 | 4.0e-01 |

## 3. Blocks exact to round-off (computed)

Criterion: max rel <= 1e-6 for at least one of the three fixed-lag steps, and insensitive to the step
beyond round-off/truncation (section 6).

- **pp** at every state of `deck_step_klevel` and `level` (1.2e-8 to 4.9e-7 depending on h; the 1e-5 at
  `dry` with h = 1 Pa is round-off: the same block is 1.2e-8 with h = 1 MPa), and of `level_z` at
  eps_v = 0 (`dry`, `mid`).
- **pp, pu, up, uu** at `wet_ceiling` (all four IPs KKT-active) for every switch set: pp 0, pu 8e-11,
  up 2.9e-9, uu 5.5e-10 to 5.9e-10 (h = abs).
- The regression guard of the test asserts exactly this list (EXPECT_LT 1e-6 on the smallest of the three
  fixed-lag errors).

## 4. Blocks with error above 1e-6 relative, with the code site (computed, derived)

Every site below is in `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h` (line numbers of
cc5e71e12a). The attribution uses two independent measurements: (i) the error matrix is fitted
(least squares) to the volumetric patterns P = int B^T m N_p w (up), P^T (pu) and Q = int B^T m m^T B w
(uu) of this element, computed independently in `dsm_jacobian_fd_analysis.py` (P is verified against the
analytic up block at `wet_ceiling`, where it is exactly -P: ‖J_up + P‖/‖P‖ = 4.4e-16); (ii) the mean
swelling stress and n_l of the integration-point state are differenced under a uniform pressure shift
(+-100 Pa) and a uniform volumetric strain (+-1e-7), one evaluation from S0 each. Fit table (h = abs;
c is the coefficient of the pattern, "resid" the relative residual of the fit):

| switch set | state | c_pu = fit err_pu / P^T | resid | c_up = fit err_up / P [Pa/Pa] | resid | -dsigma_sw/dp (IP FD) | c_uu = fit err_uu / Q [Pa] | resid | dsigma_sw/deps_v (IP FD) | dn_l/deps_v | (dsigma_sw/dp)/(dn_l/dp) dn_l/deps_v |
|---|---|---|---|---|---|---|---|---|---|---|---|
| deck_step_klevel | dry | -1.0336e-02 | 5.9e-10 | 1.0324e-04 | 2.0e-07 | 1.0324e-04 | 1.0185e+04 | 1.1e-08 | -1.0185e+04 | 1.6288e-04 | -1.0690e+04 |
| deck_step_klevel | mid | -1.3584e-03 | 2.6e-09 | 2.1070e-02 | 2.1e-07 | 2.1070e-02 | 5.8891e+04 | 3.8e-07 | -5.8891e+04 | 5.6073e-02 | -3.8810e+05 |
| deck_step_klevel | mid_eps_step | -1.3378e-03 | 2.5e-09 | 2.1557e-02 | 2.4e-07 | 2.1557e-02 | 5.8029e+04 | 4.9e-07 | -5.8029e+04 | 5.5330e-02 | -3.8998e+05 |
| deck_step_klevel | wet_ceiling | -4.7099e-13 | 9.1e-01 | 6.1241e-10 | 9.8e-01 | -0.0000e+00 | -8.2100e-03 | 9.5e-01 | -2.7590e+06 | 5.7554e-01 | 0.0000e+00 |
| deck_step_klevel | mid_eps_hist | -1.3571e-03 | 3.7e-09 | 2.1083e-02 | 1.7e-07 | 2.1083e-02 | 5.4142e+04 | 5.5e-07 | -5.4142e+04 | 5.6004e-02 | -3.8704e+05 |
| level | dry | -1.0336e-02 | 5.9e-10 | 1.0329e-04 | 1.5e-07 | 1.0329e-04 | 1.0695e+04 | 8.5e-09 | 3.8153e+04 | 1.6288e-04 | -1.0695e+04 |
| level | mid | -1.3584e-03 | 2.6e-09 | 7.0255e-02 | 7.4e-08 | 7.0255e-02 | 1.2941e+06 | 2.8e-08 | -1.5685e+06 | 5.6073e-02 | -1.2941e+06 |
| level | mid_eps_step | -1.3378e-03 | 2.5e-09 | 6.9624e-02 | 5.9e-08 | 6.9624e-02 | -1.1079e+07 | 2.8e-09 | -1.5181e+06 | 5.5330e-02 | -1.2595e+06 |
| level | wet_ceiling | -4.7099e-13 | 9.1e-01 | 6.1241e-10 | 9.8e-01 | -0.0000e+00 | -8.2780e-03 | 9.6e-01 | -1.3543e+07 | 5.7554e-01 | 0.0000e+00 |
| level | mid_eps_hist | -1.3571e-03 | 3.7e-09 | 6.9913e-02 | 5.4e-08 | 6.9913e-02 | 1.2835e+06 | 2.0e-08 | -1.5844e+06 | 5.6004e-02 | -1.2835e+06 |
| level_z | dry | -1.0329e-02 | 6.0e-10 | 1.0329e-04 | 1.5e-07 | 1.0329e-04 | 1.0688e+04 | 8.4e-09 | 3.8160e+04 | 1.6277e-04 | -1.0688e+04 |
| level_z | mid | -1.3405e-04 | 2.2e-08 | 7.0255e-02 | 7.4e-08 | 7.0255e-02 | 1.7806e+05 | 1.7e-07 | -4.5244e+05 | 7.7155e-03 | -1.7806e+05 |
| level_z | mid_eps_step | -1.0370e-04 | 4.4e-08 | 6.9394e-02 | 6.9e-08 | 6.9394e-02 | -1.2177e+07 | 1.8e-09 | -4.3771e+05 | 7.1273e-03 | -1.6202e+05 |
| level_z | wet_ceiling | -4.7099e-13 | 9.1e-01 | 6.1241e-10 | 9.8e-01 | -0.0000e+00 | -8.2780e-03 | 9.6e-01 | -1.3543e+07 | 5.7554e-01 | 0.0000e+00 |
| level_z | mid_eps_hist | -1.1598e-04 | 3.2e-08 | 6.9357e-02 | 8.8e-08 | 6.9357e-02 | 1.6590e+05 | 1.7e-07 | -5.1656e+05 | 7.2691e-03 | -1.6590e+05 |
| diag_deck_pu_all_exchange | dry | 1.8510e-05 | 3.3e-07 | 1.0324e-04 | 2.0e-07 | 1.0324e-04 | 1.0185e+04 | 1.1e-08 | -1.0185e+04 | 1.6288e-04 | -1.0690e+04 |
| diag_deck_pu_all_exchange | mid | 4.8362e-04 | 7.4e-09 | 2.1070e-02 | 2.1e-07 | 2.1070e-02 | 5.8891e+04 | 3.8e-07 | -5.8891e+04 | 5.6073e-02 | -3.8810e+05 |
| diag_deck_pu_all_exchange | mid_eps_step | 4.7123e-04 | 7.2e-09 | 2.1557e-02 | 2.4e-07 | 2.1557e-02 | 5.8029e+04 | 4.9e-07 | -5.8029e+04 | 5.5330e-02 | -3.8998e+05 |
| diag_deck_pu_all_exchange | wet_ceiling | -4.7099e-13 | 9.1e-01 | 6.1241e-10 | 9.8e-01 | -0.0000e+00 | -8.2100e-03 | 9.5e-01 | -2.7590e+06 | 5.7554e-01 | 0.0000e+00 |
| diag_deck_pu_all_exchange | mid_eps_hist | 4.7870e-04 | 1.0e-08 | 2.1083e-02 | 1.7e-07 | 2.1083e-02 | 5.4142e+04 | 5.5e-07 | -5.4142e+04 | 5.6004e-02 | -3.8704e+05 |
| diag_level_pu_all_exchange | dry | 1.8510e-05 | 3.3e-07 | 1.0329e-04 | 1.5e-07 | 1.0329e-04 | 1.0695e+04 | 8.5e-09 | 3.8153e+04 | 1.6288e-04 | -1.0695e+04 |
| diag_level_pu_all_exchange | mid | 4.8362e-04 | 7.4e-09 | 7.0255e-02 | 7.4e-08 | 7.0255e-02 | 1.2941e+06 | 2.8e-08 | -1.5685e+06 | 5.6073e-02 | -1.2941e+06 |
| diag_level_pu_all_exchange | mid_eps_step | 4.7123e-04 | 7.2e-09 | 6.9624e-02 | 5.9e-08 | 6.9624e-02 | -1.1079e+07 | 2.8e-09 | -1.5181e+06 | 5.5330e-02 | -1.2595e+06 |
| diag_level_pu_all_exchange | wet_ceiling | -4.7099e-13 | 9.1e-01 | 6.1241e-10 | 9.8e-01 | -0.0000e+00 | -8.2780e-03 | 9.6e-01 | -1.3543e+07 | 5.7554e-01 | 0.0000e+00 |
| diag_level_pu_all_exchange | mid_eps_hist | 4.7870e-04 | 1.0e-08 | 6.9913e-02 | 5.4e-08 | 6.9913e-02 | 1.2835e+06 | 2.0e-08 | -1.5844e+06 | 5.6004e-02 | -1.2835e+06 |
| diag_level_z_pu_all_exchange | dry | 1.8497e-05 | 3.3e-07 | 1.0329e-04 | 1.5e-07 | 1.0329e-04 | 1.0688e+04 | 8.4e-09 | 3.8160e+04 | 1.6277e-04 | -1.0688e+04 |
| diag_level_z_pu_all_exchange | mid | 6.6544e-05 | 4.5e-08 | 7.0255e-02 | 7.4e-08 | 7.0255e-02 | 1.7806e+05 | 1.7e-07 | -4.5244e+05 | 7.7155e-03 | -1.7806e+05 |
| diag_level_z_pu_all_exchange | mid_eps_step | 6.9293e-05 | 6.7e-08 | 6.9394e-02 | 6.9e-08 | 6.9394e-02 | -1.2177e+07 | 1.8e-09 | -4.3771e+05 | 7.1273e-03 | -1.6202e+05 |
| diag_level_z_pu_all_exchange | wet_ceiling | -4.7099e-13 | 9.1e-01 | 6.1241e-10 | 9.8e-01 | -0.0000e+00 | -8.2780e-03 | 9.6e-01 | -1.3543e+07 | 5.7554e-01 | 0.0000e+00 |
| diag_level_z_pu_all_exchange | mid_eps_hist | 7.0597e-05 | 5.2e-08 | 6.9357e-02 | 8.8e-08 | 6.9357e-02 | 1.6590e+05 | 1.7e-07 | -5.1656e+05 | 7.2691e-03 | -1.6590e+05 |

### 4.1 up block, every switch set, every state off the ceiling: error 100 % (omission)

- computed: the analytic up block is identically zero at `dry` and `mid*` (the Bishop term
  alpha chi(S_L) with S_L <= 1e-7 is the only up entry assembled); the FD up block is c_up P with fit
  residual <= 2.4e-7, and c_up equals -d sigma_sw,mean/dp of the integration-point state to 4 to 5
  digits (1.03e-4 Pa/Pa dry; 0.0211 step form, 0.0703 level form at `mid`).
- derived: the missing up entry is the whole pressure derivative of the swelling eigenstress (its n_l(p)
  channel), i.e. an omission, not an approximation.
- sites: `constexpr bool enable_dsm_swelling_up_jacobian = false;` (line 8087), declared at 8081:
  "u-side swelling Jacobian blocks (K[u,p]/K[u,u]) kept in code but OFF by default — enabling them
  singularizes the assembled tangent on stiff/dense MS33 cases (dd1800, ModelIII gap2mm: SparseLU
  compute() failure) per the at-scale comparison 2026-06-09; flip to true to opt in."; for the level form
  also line 7813: "NOT included (also omitted by the shipped tangent): the n_l channel of sigma_sw (u-p,
  enable_dsm_swelling_up_jacobian = false), the dependence of L on n_S = 1 - phi_M through eps_v
  (dL/dn_S dn_S/deps_v), and a suction-dependent skeleton."
- On the ceiling (all IPs KKT-active, n_l pinned, beta_LR = 0) d sigma_sw/dp = 0 and the block is exact.

### 4.2 pu block, every switch set, every state off the ceiling: error 100 % with the shipped kkt_active (omission), 0.2 % to 40 % without it (open)

- computed: with `kkt_active` the analytic pu block is 3 to 4 orders of magnitude smaller than the FD
  one (dry: max 2.5e-10 against 1.4e-6). With `all_exchange` the error drops to 1.8e-3 (dry) and 0.26
  (`mid*`, step and level form) resp. 0.33 to 0.40 (`level_z`).
- derived: the dominant part is the erasure of the exchange p-u entries at KKT-inactive points by the
  Q9 assignment, an omission by design. Site lines 8493-8520: "Q9 (DESIGN.md D1):
  micro_ceiling_pu_tangent. overwritten (default, the shipped line, bitwise): every exchange p-u entry
  accumulated above is erased. all_exchange: Kpu/dt is ADDED to the accumulated exchange entries (repo
  rule 4.1, accumulators use +=), which also revives the Maxwell, film and live-K p-u entries at inactive
  points. kkt_active / all_exchange: the KKT-active entries are added after the assignment." and line
  6982: "effective only at micro_ceiling_pu_tangent = all_exchange (otherwise erased by the Kpu/dt
  assignment below)".
- computed: the remainder at `all_exchange` is purely volumetric (fit to P^T with residual <= 3e-7),
  c_pu = 1.85e-5 (dry), 4.8e-4 (`mid*`, step and level), 6.7e-5 to 7.1e-5 (`level_z`) [kg/(m^3 s), coefficient of P^T]. open: which strain channel of the pressure residual it is
  (candidates by the code comments: the n_l(eps_v) and phi_M(eps_v) dependences of the macro storage and
  of the exchange; line 8564 lists "phi_s tangents in J_pp, J_pu" among the known inexact terms). Not
  attributed by this measurement.

### 4.3 uu block

Step form (`deck_step_klevel`), error 1.8e-4 (dry) to 1.0e-3 (`mid*`): omission.
- computed: err_uu = c_uu Q with residual <= 5.5e-7 and c_uu = -d sigma_sw,mean/d eps_v of the state to
  5 digits (1.0185e4 Pa dry, 5.889e4 Pa `mid`).
- derived: the whole strain derivative of the step-form swelling stress at fixed lag is missing from the
  tangent. Site: the same compiled-out block, `enable_dsm_swelling_up_jacobian = false` (8087; it holds
  the film-ON K[u,u] term "d sigma_sw/d eps_v = +(1-phi_M)*n_l*b*K_drained"), and line 7700: "the step
  rule's partial omits it (shipped, kept bitwise when the level form is off)".
- computed: with the lag re-converged (`abs+relag5`) the uu error grows to 2.8e-2: the lagged drain
  p_conf adds a strain response the fixed-lag tangent does not carry (line 7699: "p_conf lagged, a
  constant here").

Level form (`level`, `level_z`), states with no strain increment in the step (`dry`, `mid`,
`mid_eps_hist`), error 1.9e-4 to 2.3e-2: omission of the n_l channel.
- computed: c_uu equals, to 4 to 5 digits, -(d sigma_sw/dp)/(dn_l/dp) * dn_l/deps_v from the state
  differences (`level` mid: 1.2941e6 Pa against 1.2941e6 Pa; `level_z` mid: 1.7806e5 against 1.7806e5).
- derived: the level-form tangent carries everything of d sigma_sw/d eps_v except the implicit n_l
  channel, as declared at line 7813 ("NOT included ...: the n_l channel of sigma_sw").
- The Z weight lowers this error (2.3e-2 -> 3.1e-3 at `mid`) because the Z partner lowers dn_l/deps_v
  (0.056 -> 0.0077; computed), not because the channel is added.

Level form at `mid_eps_step` (strain increment -0.01 in the step), error 0.20 (`level`) and 0.21
(`level_z`): an additional term (approximation of the gate).
- computed: c_uu = -1.1079e7 Pa (`level`), -1.2177e7 Pa (`level_z`); after removal of the n_l channel
  (+1.2595e6, +1.6202e5) an extra -1.234e7 Pa remains in both; the FD response of the state
  (d sigma_sw/d eps_v) is the same as at `mid`.
- derived: the extra equals, to 0.5 % (L approximated by sigma_sw, L_0 neglected), the n_S-channel term
  ds/deps_v = (L/n_S)/(1 + c) (-1/(1 - n_l)) (alpha - phi)/(1 + w) = 1.228e7 Pa, i.e. at this state the
  block "n_S channel (2026-10-04, review should-fix 6, JACOBIAN-ONLY)" (lines 7854-7915) does not
  contribute while the residual still carries the dependence. The block zeroes its dphi/deps_v when its
  clamp detector `clamp_active_lv` (line 7897) fires. open: why the detector fires at a strain increment
  (phi is not clamped there; the recomputed phi_unclamped and the porosity of the constitutive setting
  differ by more than 1e-12 relative); not instrumented here.

### 4.4 pp block with the Z-weighted partner at eps_v != 0: error 5.5e-4 (open)

- computed: `level_z` pp max rel 5.5e-4 at `mid_eps_step` and `mid_eps_hist` for all three steps (4.5e-4
  with the lag re-converged); 1.7e-8 at eps_v = 0; `level` (drained-line partner) pp 1.3e-8 to 2.1e-8 at
  the same strained states; unchanged by `all_exchange`.
- derived: the error belongs to the Z-weighted partner; it vanishes where the partner vanishes
  (eps_v = 0). open: which of its inputs or tangents (d/dn_l total with n_S held, d/drho_lR, n_S and
  s_0 - L_0 read from the states, design note section 6) is responsible.

## 5. The Z-weighted partner: do its Jacobian sites make pu and pp exact? (computed, derived)

The partner is dispatched at the residual fold (line 1008) and at three Jacobian sites
(`computeMechanicalMicroPotentialPartner` at 5372 in `assemble`, 6928 macro exchange in
`assembleWithJacobian`, 7047 the p-u tangent).
- **pp: not exact.** Exact at eps_v = 0 (partner zero), 5.5e-4 relative at eps_v = -0.01 (section 4.4),
  where the drained-line partner is exact to 2e-8.
- **pu: not exact.** With the shipped `kkt_active` the p-u entries written by the sites 6928/7047 are
  erased at every KKT-inactive point by the Q9 assignment (section 4.2), so the block is 100 % off at
  every state off the ceiling for all switch sets alike. With `all_exchange` (entries kept) the block is
  still 33 % (`mid`), 40 % (`mid_eps_step`), 38 % (`mid_eps_hist`) off, against 26 % for `level`;
  the absolute pu error is smaller with Z (9.3e-9 against 6.8e-8 at `mid`) because the pu block itself
  is smaller (max |J_fd| 1.9e-8 against 1.9e-7).
- On the ceiling (all IPs KKT-active) the partner does not enter the exchange and every block is exact.

## 6. Sensitivity to the step h (computed)

Max rel per block for the three fixed-lag steps and the lag-re-converged FD (shipped switch sets):

| switch set | state | block | h = 1 Pa / 1e-9 m | h = 1e-2 |x| | h = 1e-7 |x| | 1 Pa / 1e-9 m, lag re-converged |
|---|---|---|---|---|---|---|
| all | dry | pp | 1.0e-5 | 1.2e-8 | 1.7e-6 | 1.0e-5 |
| all | mid* | pp (eps_v = 0 or no Z) | 1.3e-8 to 2.1e-8 | 4.8e-7 to 4.9e-7 | 6.8e-9 to 7.2e-9 | same as h = abs |
| all | wet_ceiling | pu / up / uu | 8e-11 / 2.9e-9 / 5.5e-10 | 1.7e-13 / 2.3e-15 / 1.2e-8 to 2.3e-8 | 2.6e-8 / 3.2e-10 / 3.2e-8 to 4.1e-8 | as h = abs |
| deck_step_klevel | mid* | uu | 9.4e-4 to 1.0e-3 | same | same | 2.7e-2 to 2.8e-2 |
| level | mid / mid_eps_hist / mid_eps_step | uu | 2.3e-2 / 2.2e-2 / 0.20 | same | same | same |
| level_z | mid / mid_eps_hist / mid_eps_step | uu | 3.1e-3 / 2.9e-3 / 0.21 | same | same | 2.5e-3 / 2.2e-3 / 0.21 |
| level_z | mid_eps_* | pp | 5.5e-4 | 5.5e-4 | 5.5e-4 | 4.5e-4 |
| all | off the ceiling | pu, up | 1.0 | 1.0 | 1.0 | 1.0 |

- computed: every error above 1e-6 is independent of h to two digits over five decades of h (1e-2 |x|
  to 1 Pa / 1e-7 |x|): these are structural (missing or extra terms), not FD artefacts.
- computed: the exact blocks show the usual FD behaviour: at `dry` the pp block (|J| ~ 5e-17) is
  round-off limited at h = 1 Pa (1e-5) and 10 Pa (1.7e-6) and exact at h = 1 MPa (1.2e-8); at `mid`
  h = 0.46 MPa gives a truncation error of 4.9e-7.
- computed: re-converging the lag changes only the step-form uu (x 28: lagged drain) and slightly the
  Z-partner pp/uu (x 0.8): every other block is lag-independent at these states.
- The full per-case table is `docs/dsm_jacobian_fd_verification_2026-10_data/tables.md` and
  `results.csv` (also max entrywise relative errors).

## 7. Overall error norm (computed)

‖J_an - J_fd‖_F / ‖J_fd‖_F with row scaling (h = abs): 0.66 at every state off the ceiling for the
three shipped switch sets (0.68 at `mid_eps_step` with the level form), dominated by the pu and up blocks
that are 100 % off; 3.6e-10 to 4.8e-10 at `wet_ceiling`. With `all_exchange`: 1.2e-3 (dry) and 0.23 to
0.47 (`mid*`).

## 8. Scope and what stays open (open)

- One element, uniform fields, the first Newton iterate of the next step at four times of one deck
  (dd1600, Model I); a state after a rejected step is not covered. The FD is at fixed lagged state;
  the lag-re-converged FD is reported as a second column only.
- The remaining pu error at `all_exchange` (section 4.2), the cause of the n_S-channel gate at a strain
  increment (4.3) and the source of the Z-partner pp error (4.4) are measured, not explained.
- To reproduce: `Tests/Python/dsm_jacobian_fd_run.sh <build> <outdir>` (needs numpy).
