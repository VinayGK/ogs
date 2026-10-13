# One-element closed loop of the level-form pair (sigma_sw, mu_lR): STATUS NOTE (2026-10-13)

Status words: computed (a run or a script), derived (follows by the step shown), open.
Branch `dsm_z_weight_2026-10_loop1`, cut from `dsm_z_weight_2026-10` (cc5e71e12a).

## 0. Status of this session: BLOCKED AT THE BUILD (open)

No OGS run was made, so this note has **no closure check, no budget table and no ctest result**.
Nothing in it is a computed loop number. What stopped it:

1. The cloud network policy refuses CONNECT to `git.code.sf.net` (tclap), `gitlab.opengeosys.org`
   (xmlpatch, xdmflib) and `gitlab.kitware.com` (the VTK patch). CPM cannot fetch these three
   dependencies, so `cmake` cannot configure (recorded in the proxy status: `connect_rejected`).
2. Swapping in sources from GitHub mirrors (`mirror/tclap` branch 1.4, `ufz/xmlpatch` v0.4.2, and
   the upstream xdmf from `live-clones/xdmf` patched locally to the `OgsXdmf`/`OgsXdmfCore` target
   names) got `cmake` to configure, but the session's safety policy refused to compile third-party
   code from unofficial mirrors. That was not worked around.
3. Also seen: with CMake 4.0 the prompt's `-B build` (build dir inside the source tree) fails at
   generate ("Base directories in file set cannot be subdirectories of each other"); an
   out-of-tree build dir (`../build/release`, as AGENTS.md says) avoids this. CMake 3.28 (the
   system one) is below the required 3.31.

What would unblock it: run on a machine with the CPM sources reachable (or a populated
`CPM_SOURCE_CACHE`), or allow `git.code.sf.net` and `gitlab.opengeosys.org` in the environment's
network policy.

## 1. Design of the loop (as implemented in the generator; not yet run)

Generator `Tests/Data/RichardsMechanics/DSM_closed_loop_2026-10/make_loop_decks.py` reads the parent
`ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj` (unchanged) and writes one deck per cycle and switch.
Changes against the parent; every material and numerical value is the parent's:

- `<swelling_stress_form>level</swelling_stress_form>` under `<potential_exchange>`; scalar K and
  `swelling_stress_k_level = true` kept; ON decks add `<micro_potential_partner_z_weight>true</...>`.
- Top face (y = 0.01 m) driven in y by curve `strain` x 1e-4 m (= +1 % of the height). Left, right
  (radial) and bottom (axial) stay fixed, so eps_v = eps_zz.
- Suction by curve `suction` (p_L in Pa) on top and bottom, i.e. on all four nodes.
- Schedule [d]: dry hold H0 = 400 at p_L = -1e8 (the loop starts at its END, so that an initial
  n_l out of equilibrium with 100 MPa suction, the suspected cause of the earlier 0.00115 vs 0.00035
  mismatch, is relaxed before the loop), then four legs, each a 200 d linear ramp + a hold of H = 400:
  cycle A wet, strain, dry, unstrain; cycle B strain, wet, unstrain, dry. t_end = 2800 d.
  H0 and H are placeholders until a pilot run shows every state variable stationary to 1e-8 at each
  hold end (checked by `loop_budget.py --hold-ends`; not assumed).
- Output every 0.2 d in ramps and 2 d in holds (5000 frames); `--ramp-dt 0.1 --hold-dt 1` gives the
  halved set for the integration-error estimate. `micro_liquid_density` and `epsilon` are written in
  addition to the parent's fields.
- `--test` adds the `<test_definition>` of the brief (displacement 1e-9/1e-7, porosity 1e-8/1e-7,
  micro_water_content 1e-8/1e-7, sigma and swelling_stress 1e3/1e-2, pressure 1e3/1e-2) on the final
  frame. Not used yet: no reference frame exists, and nothing is registered in Tests.cmake.

The four decks in that directory are the generator's default output (H0 = H = 400 d), untested.

## 2. Sign conventions (read from the code; to be confirmed by D >= 0 on a run)

- `micro_pressure` = p_L_m = -rho_pi mu_lR, rho_pi = micro liquid density when
  `use_micro_liquid_density_for_micro_pressure` (true in the parent): `RichardsMechanicsFEM-impl.h`,
  the `.p_L_m = -p_L_m_density * micro_potential.mu_lR` return (search `p_L_m_density`).
- Level form: `L = -n_S n_l [Pi - b p_conf]`, `d sigma_sw = L_curr - L_prev`, closed form
  `L (1 + c) = -n_S n_l Pi - c (m_hat + s_prev - L_prev)`, block `if (level_form_sw)` of the same file.
- `micro_liquid_density` is an output name of the process (`ConstitutiveRelations/MicroLiquidDensity.h`).

## 3. Evaluation tool (computed on a synthetic loop only)

`tools/loop_budget.py RUN.pvd --t-start H0 [--legs ...] [--hold-ends ...] [--json f]`: nodal mean
over the nodes, trapezoidal closed integrals of W_sigma = int sigma_sw d eps_v (eps_v from
d phi/(1-phi), and from the written epsilon as a check), W_mu = int mu_lR d m_mic (mu_lR =
-micro_pressure/rho_lR, m_mic = rho_lR (1-phi_M) n_l), W_mu_nl = int (1-phi_M) rho_lR mu_lR d n_l,
R_pair = W_sigma + W_mu, D = int (mu_LR - mu_lR) d m_mic (mu_LR = p_L/1000), W_chem = int mu_LR
d m_tot (m_tot = 1000 S_L phi_M + m_mic), budget = W_sigma + W_chem - D; start/end differences;
per-leg integrals; the sign of every D increment; stationarity at hold ends.
Check (computed): on a synthetic conservative pair Psi = a eps n + k eps^2/2 + c n^2 over a closed
(eps, n) rectangle, 401 frames, |W_sigma| = |W_mu| = 9.0e4 J/m3 and R_pair = 7.5e-5 J/m3. This checks
the tool, not the model.

## 4. Open

Everything the brief asks to compute: the closure check, the OFF/ON budget tables for cycles A and B,
the halved-interval error estimate, the reference frame, the Tests.cmake entry and its ctest.
