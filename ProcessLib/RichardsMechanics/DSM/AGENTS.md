# AGENTS.md — ANCHORS MS33 DSM (streamlined)

## Scope
- Models: I / III / IV / VII (plus V_LE/V_MCC where present).
- Core implementation: `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`.
- Branch: `dsm_native_hierarchical`.
- Binary: `/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs`.

## Roadmap (one-line commit refs)
- Step 1 REV-scale storage + split consistency: `0d7a9edd64`.
- Step 2 thermodynamic swelling stress + K recalibration: `88d42c98fd`.
- Step 3 Pi-path Gibbs–Duhem consistency + flag cleanup: `c4888b6db4`, `ce9178fa96`.
- Step 5 vdW dimensional fix (`/rho_lR`) + literature A lock: `0d579e8aeb`.
- Step 6 DSM hardening (viscosity guards, micro-pressure density default true): `66b782afa1`.
- Step 7 dead-code removal (compatibility overload/unused flag): `4d47efff55`, `ce9178fa96`.
- Step 8 DSM micro-macro test refactor (13/13 passing): `3ac6b7de1f`.

## Key physics/implementation invariants
- Porosity split: `phi = phi_M + phi_m`, with micro state carried by `n_l`.
- Storage is REV-scale: `phi_m * rho_lR`.
- Swelling stress uses thermodynamic Pi-path tied to `rho_LR` for Gibbs–Duhem consistency.
- vdW potential terms are additive; never replace additive update with assignment.
- `hamaker_constant` is literature/material-fixed; calibration target is K (not A).

## Execution instructions
- Keep committed runs benchmark-spec compliant.
- After physics changes, require:
  1. Model-I Villar target check within tolerance,
  2. canonical LE reruns with zero rejected steps,
  3. `./bin/testrunner --gtest_filter='*DSMMicroMacro*'` passing.

## Current summary
- Production path stable under latest DSM fixes.
- Canonical LE outcomes unchanged in accepted/rejected-step sense.
- Open benchmark-quality work is primarily calibration/interpretation side (not immediate solver-break state).

## DSM_NATIVE_HIERARCHICAL_PATCH_RECIPE.md maintenance rule

`ProcessLib/RichardsMechanics/DSM/DSM_NATIVE_HIERARCHICAL_PATCH_RECIPE.md` is the reconstruction
recipe for this branch from a fresh `master`. It must stay current.

**Update DSM_NATIVE_HIERARCHICAL_PATCH_RECIPE.md before committing whenever:**
- Any hunk in `RichardsMechanicsFEM-impl.h` or `PotentialExchange.h` changes.
- The DSMMicroMacro unit tests change (step 8 section + passing count).
- Any PRJ `hamaker_constant`, `potential_augmentation_prefactor` (K), or pre-consolidation
  pressure (`pc_char_mcc` / `InitialPreConsolidationPressure`, the MCC cap `pc`) value changes.
- A new benchmark model is added to the canonical LE set.
- Build flags or the verification `ctest` invocation changes.
- A new step beyond Step 8 is added (add a new numbered section).

Do not mark a step done in AGENTS.md unless DSM_NATIVE_HIERARCHICAL_PATCH_RECIPE.md already reflects it.

## Known limitations (logged 2026-05-27)

### Forgotten Maxwell pair — mean stress absent from μ_lR  (NEXT IMPLEMENTATION GOAL)

**Location:** `PotentialExchange.h` — `μ_lR = μ_lR(n_l, ρ_lR)` (aggregate fraction
`1−n_l`, no stress/strain); swelling enters σ one-way as the Π-eigenstress
`σ_sw = −φ_m·Π(n_l)` (`RichardsMechanicsFEM-impl.h`, commit `72f4f3a192`).

**Issue:** from one `Ψ(ε,n_l)`, `∂σ/∂n_l = ∂μ_lR/∂ε` (Maxwell). The Π-eigenstress
gives `∂σ/∂n_l ≠ 0` (mean-stress, isotropic); `μ_lR` has no stress dependence so
`∂μ_lR/∂ε = 0` → **broken pair** → not derivable from a single `Ψ` →
non-conservative (`∮≠0`) past the gate. (= the "unlicensed equipresence deletion"
of the T&N lecture's Maxwell-pairs unit.)

**Gate:** exact **iff `σ_n < Π`** — below the disjoining pressure the missing
mechanical-expulsion term is physically zero; swelling-driven monotonic loading is
in-domain. Breaks past `Π` (over-compaction): `Π(n_l)=σ_n` is never closed (the
exchange does only the chemical balance `ψ_M=ψ_m`) → `φ_M→0` crash. Current claim:
restricted-domain admissibility (A).

**Right fix = the next goal:** give `μ_lR` a mean-effective-stress dependence — the
Maxwell partner of the swelling eigenstress, derived from one `Ψ` (coefficient
fixed by the existing swelling closure, **no new constant**; re-verify `K`). Full
spec — diagnosis, code plan, tiers (B1 sharp / B2 smear `Π`), verification — in
**`MAXWELL_PAIR_RESTORATION.md`**.

### Hydraulic-side double-counting of suction in Darcy flux

**Location:** `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h:2627` (and
the parallel branches at 3368 and 3991) — assembly of the macro Darcy flux:

```cpp
laplace_p += dNdx_p^T · (k_intr·k_rel·ρ_LR/μ) · dNdx_p · w
// drives:  q_L = -(k_intr·k_rel/μ) · ∇p_L
```

**Issue:** `p_L` is the primary `pressure` process variable, which the boundary
conditions impose as *total* suction (lab-measured, capillary + molecular).
The full ∇p_L therefore drives the macro Darcy flux, including the molecular
(disjoining-pressure Π) component. Physically the molecular component should
drive the micro→macro mass exchange via the DSM source term ρ̇_micro→macro,
NOT bulk advection in the macro pore.

**Manifestation:** boundary suction ramps of order 100 MPa create huge ∇p_L
gradients that propagate the wetting front much faster than the material
physically would. Affects all transient III/IV/VII results at finite t. Does
NOT affect Model I (no spatial gradient) and does NOT affect t→∞ equilibria
of III/IV/VII (saturated swelling pressure, asymptotic gap closure, equilibrium
void ratio).

**Mechanical-side companion (already fixed):** the same total p_L was being
fed into Bishop's effective stress via the `BishopsPowerLaw(exponent=1)` path
(χ = S_L → χ·p_L = molecular component leaks into σ_eff during unsaturated
phase). Fixed by switching all 4 MS33 PRJs to `BishopsSaturationCutoff(cutoff=1)`
so χ=0 below S_L=1; DSM Π then carries all swelling-source work in the
unsaturated regime.

**Right fix (open):** split p_L into p_L_macro (capillary, bounded by
Young–Laplace ~3 MPa for compacted bentonite) and p_L_micro (disjoining, Π).
Use ∇p_L_macro only in Darcy; route the (p_L − p_L_macro) residual through
the existing DSM micro↔macro exchange. Requires a new constitutive law for
the split and a process-variable refactor in `RichardsMechanicsProcess`.
Estimate: 1–2 weeks careful work + tests.

### MCC tension-apex non-convergence (ModelVII AND ModelIV) — use LE

**Location:** `MaterialLib/SolidModels/MFront/ModCamClay_semiExpl_constE.mfront`
(yield `f = q² + M²·p·(p−pc)`, `p = −trace(σ)/3 + pamb`).

**Issue:** in the free-swelling ModelVII, differential swelling at the wetting
front drives integration points to the tensile APEX of the cam-clay ellipse —
the failing Gauss point sits at `p ≈ 0` with a residual deviatoric `q ≈ 0.16 MPa`
and `pc ≈ 12 MPa` (healthy). At p=0 the ellipse pinches to the single point
(0,0); a state with q≠0 has no admissible return → MFront `status -1`.

**ModelIV joins this class (2026-05-29).** Once the clay–pellet ModelIV uses the
*physical* soft pellet modulus (`E_mcc_pellet = 9.2549 MPa = C·ρ_d³` at ρ_d=900,
parity with the LE variant), the compliant pellet lets differential swelling drive
the pellet/clay interface to the apex (`p → 0.03 MPa`, `q → 4.8 MPa`, `eqpl=0`,
`status -1` at ~22 d). The earlier MCC ModelIV that "completed" (1977 steps) used an
unphysically stiff pellet (`E = 52 MPa = E_clay`) that suppressed the interface
deformation. Per the user decision (physical params preferred), **ModelIV is now
submitted with the LE variant** (`ANCHORS_MS33_LE_PER_DD/ModelIV`, soft pellet,
converges: clay 13.9 / pellet ~0 MPa). ModelIII does NOT fail because its soft gap
is a per-material LinearElastic zone (id=1) — no MCC integration there, no apex.

**Tested (2026-05-29) and rejected as the fix:** added a gated `TensionCutoff`
(`pt_cutoff`) @Parameter (default −1 = OFF; converging models I/III/IV stay
byte-identical, verified dd1400 elastic 4.908 MPa / eqpl=0). `pt_cutoff=0`
(strict no-tension) routes p<0 elastic predictions to a volumetric-only return
that caps the mean stress at 0. It works (nodal min mean −63 → +52 kPa) but VII
still fails at the same step: the cutoff caps p but leaves q, so `f = q² > 0`
persists at the apex. Completing the Jacobian did not change the outcome →
structural, not a numerical bug. The earlier "negative mean stress" reading was
a nodal-extrapolation artifact; the true obstruction is the p→0 / q≠0 apex
coincidence. `pc_min` and `pamb` were separately ruled out (see
`project_dsm_mcc_bishop_cutoff.md` memory).

**Verdict:** the gated cutoff is left disabled in the .mfront as the recorded
experiment, NOT a production lever. ModelVII free swelling and ModelIV clay–pellet
both use the LE variant (`ANCHORS_MS33_LE_PER_DD/ModelVII` and `.../ModelIV`).
A genuine cure would relax BOTH p→0 and q→0 (full apex/fissuring collapse), which
zeroes wetting-front stiffness and is not globally solvable in OGS RM.

**Pragmatic interim:** submit current results with explicit caveat in the
deliverable (deck frame "Known Limitation — Hydraulic Side of OGS RM-DSM").

**2026-06-01 -- full-Pi closure; beta_sw retired; EMDD=rho_d calibration.**
The disjoining-pressure eigenstress (Delta sigma_sw = n_S(n_l_prev Pi_prev - n_l Pi) I,
Pi = -rho_lR psi_Micro) is the implemented micro-swelling closure; the legacy
micro_water_content_swelling_slope (beta_sw) branch and its
<accumulate_swelling_contributions> PRJ tag are removed. K re-fit (dt-converged P*
basis) to the Dixon (2023) MX-80 anchor under the EMDD=rho_d ANCHORS-groups
agreement: targets 4.92/14.16/40.86 MPa, K=35625.4/85312.6/224610 J/kg at
dd1400/1600/1800. CLAUDE.md 12.1 updated; 12.2 provenance synced across the LE
suite; 12.2 blocks added to the ModelVII experimental BC variants. Verified:
ogs+testrunner build clean; 14/14 DSM unit tests pass (incl. corrected
active_nS=1-n_l, section-2 incident); full MS33 suite (10/11 to t_end; dd1600
documented corner crash; endpoints ~2% high, first-order dt error).

**2026-06-08 -- K(rho_d): augmentation prefactor as a function of dry density.**
Branch `dsm_native_pdisj_maxwell_kofdd` (off `dsm_native_pdisj_maxwell`); build
`/Users/vinaykumar/git/build/kofdd_20260608/bin/ogs`. Implemented Option A
(parse-time table): `potential_augmentation_prefactor` may be set by a
`<potential_augmentation_prefactor_vs_dry_density>` block (child lists
`<dry_densities>`/`<prefactors>`) evaluated at each material's `<dry_density>`
(rho_d). K is resolved to a scalar at parse time -> *initial/target* rho_d,
constant in time, so NO Jacobian/tangent change downstream (Vinay's call
2026-06-08). The shared table inherits into per-`<medium id>` overrides via the
existing defaults mechanism; the scalar key and the table are mutually exclusive
(OGS_FATAL if both); a table requires a `<dry_density>`; `getValue` clamps
outside the node range. Files: `PotentialExchangeParameters.h` (+2 fields,
forward-decl MathLib::PiecewiseLinearInterpolation), `CreateRichardsMechanics
Process.cpp` (parse + resolve). DONE.
- Test (pellet block, Model IV): curve-vs-scalar equivalence. New PRJs
  `ms33_modelIV_pellets_kref100x.prj` (scalar K) and `..._kofdd.prj` (table K),
  both at k0 x100 spec (test acceleration, Vinay 2026-06-08; rate-only, endpoint
  unchanged), differing ONLY in how K is specified. Table nodes {(900, 20600),
  (1600, 103879.0)} J/kg carried verbatim from `ms33_modelIV_pellets.prj` so the
  table reproduces each material's existing per-material K. VERIFIED 2026-06-08:
  both run to t=200 d (ts_689); vtkdiff abs-max diff = 0 on all 14 output fields
  -> bit-for-bit identical. Registered both run-only in `Tests.cmake`.
- Back-compat VERIFIED: unmodified `ms33_modelI_dd1600.prj` runs clean through the
  refactored parser (scalar path).
- OPEN (deferred, Vinay): (i) closed-form vs table interpolation shape for
  intermediate rho_d (linear/log) is a modelling choice, not yet decided;
  (ii) the *current/evolving* rho_d variant (K riding porosity, needs a tangent
  term, double-count risk vs the exp(-xi) porosity dependence) is NOT
  implemented -- this delivery is initial/target rho_d only.

## Strained-film disjoining law h(w_m, eps_v) (2026-06-09, branch dsm_native_h_of_eps)

Goal (Vinay): reversible de-swelling/expulsion under load — film thickness
varies with compressive stress so the potential reverses; the dissipative
residual stays a future flow rule. Design + decision record:
**STRAINED_FILM_IMPLEMENTATION.md** (this directory).

- DONE 2026-06-09: enums + params (film_strain_coupling off|kinematic|
  equilibrium, film_strain_kappa aggregate|unity), computeStrainedFilmState +
  invertDisjoiningPressure (PotentialExchange.h), fold-point rewiring
  (applyFilmPressureMicroPotential REPLACES the shipped integrable partner when
  ON — D3 provisional), eigenstress threading (eps_v sentinel args), PRJ
  parsing, unit tests (StrainedFilmPotential.cpp). Off = bit-for-bit.
- STRUCTURAL FINDING [D]: a pure geometric squeeze of any repulsive Pi(h) can
  never reverse the potential (Pi'(h)<0 ⇒ imbibition); the reversal lives in
  the Derjaguin load term +b*p_conf/rho_lR, made h-live here. Emergent gate
  b*K_drained > 3*kappa*Pi(h) (kinematic) / min()-branch at p_conf = Pi(w_m)
  (equilibrium) — no bolted-on Macaulay gate.
- HONESTY NOTE: implemented cut = operational Derjaguin form, NOT yet
  Maxwell-exact; exact one-Psi closed forms derived in the design doc §9a,
  AWAITING Vinay's review before coding. Do not cite the branch as
  "Maxwell-exact".
- TODO: build + unit tests + dd1400 off-mode regression (in progress
  2026-06-09); §9a exact forms; confined expulsion probe; K re-calibration
  [PRED: saturated swelling-pressure equilibrium shifts in both modes].

## 2026-06-11 — LIVE K(rho_d) variant (Vinay: "K(rho_d) try it")

- DONE (2026-06-11): live (evolving dry-density) K(rho_d) implemented per
  Vinay's order; see K_OF_RHO_D_LIVE.md. New PRJ bool
  `potential_augmentation_prefactor_live_dry_density` (default false =
  parse-time freeze, bit-for-bit); helper effectiveAugmentationPrefactor
  (PotentialExchangeParameters.h) evaluates K_table(rho_SR*(1-phi)) at all
  FEM sites with porosity in scope (context phi / new defaulted
  total_porosity arg on the swelling increment / assembly phi); scalar
  fallback where no phi exists. Endpoint-clamped (getValue endpoint hold).
- Verified: 31/31 RichardsMechanics unit tests (3 new
  RichardsMechanicsLiveKOfRhoD, structural knots); dd1400 off-mode
  regression sigma_zz = -4.9218 MPa = recorded baseline
  (runs/2026-06-10_0841_dsm_native_h_of_eps_successful).
- PROVISIONAL: linear knot interpolation (shape undecided); dK tangent
  OMITTED first cut [PRED: extra Newton iterations, not benchmarked]; no
  live-mode production run yet (behavior under live K = predicted only).

## Equipresent Pi(n_l, eps_v) + compressible-liquid carrier D2 (2026-06-11, branch dsm_native_Pi_fofnlev)

Goal (Vinay, 2026-06-11): make the load coupling energetically compliant —
"Pi does not carry p_conf; the p_conf is still a bolt-on to the micro, not an
energetically compliant bolt-on; equipresence says Pi(n_l, <mech. state>)."
Equipresent argument is eps_v (configuration), not p_conf (force). Two
deliverables, PRJ-selectable, default off (bit-for-bit):
(E) exact one-Psi film pair (closes STRAINED_FILM_IMPLEMENTATION.md §9a) and
(L) compressible-liquid carrier (D2 proper — over-pressure from Psi_liq with
confined K_liq instead of the bolted +b*p_conf/rho_lR).

Design + decision record: **PI_OF_NL_EV_IMPLEMENTATION.md** (this directory) —
complete implementation/test/docs/beamer plan, written for an implementing
agent; decision queue Q1–Q5 (Vinay) inside.

- DONE 2026-06-11: branch + worktree created off 7ff8861847; design doc
  written; memory file project_dsm_pi_fofnlev + MEMORY.md pointer; beamer
  maxwell_from_psi.tex Step 19/23 "in design" markers.
- DONE 2026-06-11 (was: blocked on Q1; resolved by Vinay's "implement that
  now"): (E) implemented — FilmEnergyRoute enum/param/predicate,
  computeStrainedFilmEnergyPair (x_over_kappa-stable closed forms), exact
  fold branch (bare at TRUE n_l + one-Psi partner, g-cutoff product rule,
  eigenstress site unchanged), film_energy_route parsing + mode-matrix
  OGS_FATAL, ExactFilmEnergyPair.cpp (8 tests: 6 active + 2 Q3/Q4 skips).
  VERIFIED: 36/36 unit tests; T-5 loop measured |∮|/scale 8.4e-9 (exact) vs
  0.93 (operational — §9a "small" prediction corrected by measurement);
  T-1 dd1400 off-mode bitwise-identical (12/12 VTUs, parent-head binary vs
  new, one input). Build ~/git/build/Pi_fofnlev_20260611. UNCOMMITTED.
- PARKED 2026-06-11 (Vinay: "park it, this is getting very intricate"):
  (L) computeMicroLiquidCompression + tests T-6/T-8. Self-contained decision
  brief (routes (a)/(b), K_liq candidate + magnitudes, Q3<->Q4 coupling,
  unpark conditions) is in PI_OF_NL_EV_IMPLEMENTATION.md §8 — read THAT
  before any (L) work; do not unpark without Vinay's Q3+Q4 answers.
- DONE 2026-06-11 (partial): STRAINED_FILM_IMPLEMENTATION.md §9a annotated
  (measured correction); beamer Step 21/23/7b updated to implemented+measured
  status. STILL TODO: "Step 24 first numbers" frame (gated on T-8 / MS33 VII
  runs); Doxygen tag doc for film_energy_route (joint TODO with the
  undocumented film_strain tags).

## Form (a) vs Form (b) paired comparison (2026-06-12) — READ BEFORE TOUCHING THE FILM ROUTES

Definitions (the two FORMS of the micro potential; beamer maxwell_from_psi.tex Step 21b):
- FORM (a): mu(Pi(n_l), eps_v) — Pi frozen at n_l, strain a separate argument
  (the ch.1 maxwell partner; default/off mode; maxwell_conjugate lineage).
- FORM (b): mu(Pi(n_l, eps_v)) — strain enters THROUGH the film state
  (film_strain_coupling=kinematic + film_energy_route=exact on THIS branch).
- [D] (a) = (b) Taylor-truncated at eps_v->0. Mutually exclusive at runtime
  (create-time guard); coexistent in one build as limit test + baseline.

Paired runs (common pre-recalibration K base, two binaries:
mc_20260608 d98f5f8324 for (a), pi_fofnlev_20260611 @4c7a5d03b2 for (b)),
record: ogs/formAB_2026-06-12/FORMAB_RESULTS.md + eurad-anchors snapshot. MEASURED:
1. dd1600 control (eps_v~0): BITWISE identical across forms AND binaries —
   the truncation identity holds exactly in running code.
2. VII discriminator: e_end 1.4995 (a) vs 1.3530 (b) — self-relaxing drive
   removes 0.146 of the over-swell; exact within 3e-4 of operational here.
3. OPPOSITE PULL: Task-13 DD 1.4592 (a) vs 1.4839 (b) vs expt 1.4139 —
   (b) helps VII but hurts Task-13. No reversible form satisfies both =>
   Task-13 residual is the IRREVERSIBLE channel (see TASK13_MCC_BLOCKAGE.md
   + MCC_INTERNAL_SWELLING_IMPLEMENTATION.md design).
4. (b) also better-conditioned on the Task-13 wetting front (1256 steps/138 s
   vs 2673/713 s) and ends force-balanced (Pi +0.55 MPa ~ load vs -1.85).
Variant PRJs committed next to the base files (suffix _formB_piexact_2026-06-12).

## Live-K(rho_d) analytic tangent completion (2026-06-12)

- DONE 2026-06-12: analytic dK/dphi = -rho_SR*(table segment slope) tangent
  (Vinay-approved completion of the live-K first cut; Jacobian-only,
  residual untouched) wired into the live p-u augmentation Jacobian block.
  New `AugmentationPrefactorTable::getSegmentSlope` (exact clamped
  piecewise-linear slope; zero outside/at edges, left slope at interior
  knots), `effectiveAugmentationPrefactorPhiDerivative`, mu-level exact
  K-partials `dmu_lR_dK`/`ddmu_lR_dnl_dK` (mu_aug linear in K), and the
  `PorosityFromMassBalance` dphi/deps_v = (alpha-phi)/(1+w) chain. Details +
  measured verification: K_OF_RHO_D_LIVE.md "Analytic tangent completion".
- MEASURED: 41 RM unit tests (39 pass + 2 designed skips; 2 new FD-vs-
  analytic tangent tests); dd1400 off-mode bitwise-identical vs the
  pre-tangent h_of_eps_20260609 binary (12/12 VTUs); truncated 1a_robin_A_Kl
  live-K sanity converges with iterations equal to before (17/2/2).
- CURE TEST (task42 1b *_Kl step-1 singularity): NOT CURED — both 1b_A_Kl
  and 1b_B_Kl still die in step #1, but the Newton trajectory measurably
  changed (contraction to |dx|_uz=1.09 over 7 its before the it.8 blow-up,
  vs pre-tangent monotonic divergence). Second mechanism suspected
  (hypothesis, not verified); evidence in task42_case1_2026-06-12/
  out_1b_{A,B}_Kl/run.log + _diagnostics_1bKl/README_DIAG.md addendum.
  OPEN: Vinay's call on the next probe; no further patching done.

## Ultracode-review fixes (2026-06-14) — branch dsm_native_Pi_fofnlev_review_fixes_2026-06-14

Implements the fixes in DSM/ULTRACODE_REVIEW_2026-06-14.md, one commit per
fix. Branch off dsm_native_Pi_fofnlev tip 9795f252e1. Build:
~/git/build/pi_fofnlev_fixes_20260614. Baseline: 39 RM unit tests pass + 2
designed skips. Every Jacobian-only fix re-verified the off-mode dd1400
final-VTU SHA256 = 91f404a5...577 (the STEP-0 reference) BYTE-FOR-BYTE.

- DONE 2026-06-14 (H2, JAC-only, commit 7e2eee0190): exact-route fold passed
  the parse-time scalar K into computeStrainedFilmEnergyPair while the bare
  out.mu_lR used the live effectiveAugmentationPrefactor(phi); g_cut corrupted
  under live K. Now uses effectiveAugmentationPrefactor(params, local_context
  .phi) at L730 (matches :768/:1364/:2084). Scalar-mode bit-for-bit; off-mode
  bitwise verified.
- DONE 2026-06-14 (M1, JAC-only, commit 6708ef1d98): the p-u Maxwell exchange
  tangent always used the integrable-partner dmu_lR_mech_deps_v; now dispatches
  dmu_lR/deps_v on the route (Off -> integrable; operational-strained ->
  d(bare(w_eff))/deps_v + b*(dp_conf/deps_v)/rho with dp_conf/deps_v=-K_drained;
  exact -> g_cut*pair.dmu_mech_deps_v). Off-mode bitwise verified.
- DONE 2026-06-14 (M2+L2, JAC-only, commit cbe3e11ed6): wired the displacement-
  side live-K swelling-eigenstress tangent d(delta_sigma_sw)/dK*dK/dphi*
  (dphi/deps_v, dphi/dp) into K[u,u]/K[u,p] (the 1b compliant-top cure
  candidate). SCOPE DECISION (announced): wired ONLY the live-K chain, NOT the
  pre-existing swelling u-p/u-u term Vinay set OFF 2026-06-01 (enable_dsm_
  swelling_up_jacobian left at its default; independent term). Gated on
  film_pressure_coupling && dK/dphi != 0 && PorosityFromMassBalance. Off-mode
  bitwise verified.
- DONE 2026-06-14 (L1, JAC-only, commit 78a71ae6df): exact-fold dg_dnl mixed
  live-nS out.dmu_lR_dnl with frozen-nS pair.dmu_bare_dnl_pre under
  CurrentPorositySplit; recompute the bare-pre derivative with the caller's
  dnS_dnl. Reference mode bit-for-bit. Off-mode bitwise verified.
- DONE 2026-06-14 (N1, JAC-only, commit 6a03a58213): zero the live-K
  dphi/deps_v=(alpha-phi)/(1+w) chain when PorosityFromMassBalance clamps phi
  (detect via unclamped-vs-stored phi; bounds are private). Applied at the
  live-K p-u block and the new M2 block. Off-mode bitwise verified.
- SUPERSEDED 2026-06-14 (L3, DOC-only, commit f23f69c5b4): documented the
  live-K dn_l/dK local-solve strain channel in ScalarReferenceMassStorage mode
  as a DELIBERATE PARTIAL TANGENT (not wired). Judgment call (prompt-authorized):
  wiring it would risk the converged forward solve for a LOW-severity gap; the
  dominant cure is M2. Off-mode bitwise verified. NOTE: the "risks the converged
  forward solve" concern was the basis for not wiring; it is resolved below by
  keeping the wiring strictly Jacobian-only (the local solve / residual / the
  computeImplicitNlDpL return are all untouched). See the WIRED entry next.
- DONE 2026-06-14 (L3, JAC-only, WIRES the above): the L3 implicit-n_l(K) strain
  channel is now wired into the global displacement Jacobian, RESIDUAL-SAFE.
  - New sibling helper computeImplicitNlDK (RichardsMechanicsFEM-impl.h, right
    after computeImplicitNlDpL): returns dn_l/dK = -(dr/dK)/(dr/dn_l) for
    ScalarReferenceMassStorage (0 in every other mode and at dt<=0). It rebuilds
    dr_dn_l by the SAME REV-mass reduction as computeImplicitNlDpL and uses
    dr/dK = -dt*exchange.drho_l_hat_dmu_lR*micro_potential.dmu_lR_dK (K enters r
    only through the exchange/mu_lR; the mass term carries no K). Reads ONLY the
    already-converged (n_l, rho_lR, micro_potential, exchange) the caller
    threads in -- it never re-solves and never mutates forward state.
  - Wired at the M2 swelling-eigenstress site as the implicit half of M2's
    explicit-K chain: d(delta_sigma_sw)/dn_l * dn_l/dK * dK/dphi * dphi/d(.),
    folded into the SAME dsig_sw_deps_v_scalar / dsig_sw_dp_scalar that ride the
    existing C*C_el^-1*identity2 map into K[u,u]/K[u,p]. d(delta_sigma_sw)/dn_l
    = -n_S*(Pi_curr - n_l*rho_curr*dmu_lR/dw_eval*dw_eval/dn_l) (dw_eval/dn_l=1
    on the OFF branch, film_state.dw_eff_dnl on the strained branch), matching
    the residual increment's argument chain and its dnS_dnl=0.
  - SCOPE (residual-consistency): the implicit chain fires only where the
    telescoped eigenstress IS the assembled residual (OFF + film-ON-operational);
    on the EXACT route (H1, residual = one-Psi pair) it is gated OFF -> the
    cured 1b_B tangent is left bit-for-bit untouched. Off / frozen K / clamped
    table edge / non-mass-storage -> all factors 0 -> Jacobian bit-for-bit.
  - JACOBIAN-ONLY: the local forward n_l solve, the residual, and the
    computeImplicitNlDpL return value are all unchanged. The prior agent's note
    inside computeImplicitNlDpL is replaced with a WIRED note.
  - PREDICTED (not yet verified, §5; build is a later phase): this completes the
    live-K mass-storage displacement tangent and is the candidate cure for the
    1b_A form-(a) step-1 divergence the 1b cure verdict left open; the converged
    root is unaffected by construction. Off-mode bitwise: NOT re-run here (build
    deferred) -> the off/frozen/clamped/non-mass-storage gates make the change a
    no-op there by construction, to be confirmed at build.
- DONE 2026-06-14 (H1, RESIDUAL-CHANGING, Vinay-authorized, commit 6391e357a2):
  under film_energy_route=Exact, source the eigenstress increment from the
  one-Psi pair.sigma_sw_m (drained-line, telescoped curr-prev) instead of the
  operational Pi(w_eff)-b*p_conf. VALIDATED by the new assembled loop-closure
  probe (AssembledExactPairClosesOperationalSigmaDoesNot): |W|/scale = 8.41e-09
  with exact sigma (post-H1, CLOSES) vs 0.0675 with operational sigma (pre-H1,
  does NOT close). H1 is CORRECT (not reverted). PHYSICS TRADEOFF (predicted):
  drained-line p_conf vs actual GP p_conf; they agree on the drained line,
  differ off it (e.g. fully confined). Off-mode unaffected (film off -> branch
  not entered): bitwise verified.
- DONE 2026-06-14 (NEW TEST §8, commit 9732b46498): added
  RichardsMechanicsLiveKOfRhoD.AssembledDisplacementTangentExactKinematicLiveK
  (StrainedFilmPotential.cpp) — helper-level FD-vs-analytic identity on the
  exact-route mu_lR eps_v tangent (H2/M1) and the live-K eigenstress eps_v
  tangent (M2), anchor (d), scale-derived tolerances, no Vinay expected value.
  Scope-noted: no run-level assembleWithJacobian harness exists in the unit-
  test dir, so this checks the tangent FORMULAE the assembly reconstructs.
- DONE 2026-06-14 (N2/N3, DOC-only, commit 9c8423766c): N2 labelled the loop-
  test 100x separation as the MEASURED conservative floor per §5.1 (defect/
  bound ~3.0e3 at N=400, grows like N^2; not a derived 100). N3 re-confirmed
  computeMaxwellConjugateMicroPotential fully dead (zero live callers) and
  strengthened the RETIRED banner; kept on disk per §6.3.
- L5 — NOT FIXED (deliberate, guardrail §1.1/§12.2). The two formB_piexact
  PRJs on disk (ms33_modelI_dd1600_formB_piexact_2026-06-12.prj,
  ms33_modelVII_freeswelling_formB_piexact_2026-06-12.prj) carry TODO(Vinay)
  §12.2 provenance locators (E/nu, micro-EOS, Tuller geom, specific_surface,
  lambda) inherited verbatim from the base PRJs. They are NOT in Tests.cmake
  and STAY OUT pending Vinay's cited source locators (cannot be invented).
  Their calibration K=103879 J/kg is §12.1-clean (Dixon 2023 Fig.1).

OPEN (for Vinay): (1) the 1b cure verdict (M2) — run task42 1b_*_Kl with the
new binary; (2) whether to also flip enable_dsm_swelling_up_jacobian (the
pre-existing OFF swelling tangent), kept untouched here; (3) §12.2 locators to
register the piexact PRJs.

## 1b cure verdict (2026-06-14, the M2 headline)

Ran task42 1b_A_Kl.prj and 1b_B_Kl.prj (live K(rho_d) table, film_pressure_
coupling OFF, ScalarReferenceMassStorage; 1b_A = form (a) Off coupling, 1b_B =
form (b) kinematic + exact route) with the review-fixed binary (commit
936488482c), scratch ~/git/build/pi_fofnlev_fixes_20260614/cure_1b_test.
Pre-fix baseline (out_1b_*_Kl/run.log in task42): BOTH died in time step #1
("Newton: the linear solver failed"). MEASURED with the new binary:
 - 1b_B_Kl (exact route, H1+H2+M1+M2 all active): CURED — passes step #1 and
   keeps stepping (reached Time step #200+ with ZERO step-1 failures and zero
   nonlinear failures; observed VTUs out to t=2.59e6 s in an earlier run).
 - 1b_A_Kl (form (a), Off film_strain_coupling): NOT cured — still fails in
   time step #1 (terminated with error). The M2 live-K eigenstress tangent now
   fires for it (gate = live-K flag), so this confirms the review §3 "second,
   still unidentified mechanism" is NOT (only) the M2 eigenstress tangent for
   the form-(a) ScalarReferenceMassStorage case. Strong candidate: the L3
   dn_l/dK local-solve strain channel left as a DOCUMENTED partial tangent
   (f23f69c5b4) — Vinay's call whether to wire it next.
VERDICT: M2 (+H1/H2/M1) CURES the exact-route 1b_B; the operational-route 1b_A
remains a step-1 failure (open, points at L3 / a distinct mechanism).

UPDATE 2026-06-14 (L3 now WIRED): the L3 dn_l/dK local-solve strain channel —
the strong candidate this verdict flagged for 1b_A — has since been wired into
the displacement Jacobian (Jacobian-only; see the L3 WIRED worklog entry above).
It is the implicit-n_l(K) partner of the M2 eigenstress tangent and fires on the
1b_A OFF/operational regime (the M2 gate fired the explicit-K half there but the
local-solve n_l(K) half was missing). PREDICTED candidate cure for the 1b_A
step-1 divergence; NOT yet re-run (build is the next phase) — re-run task42
1b_A_Kl.prj with the new binary to confirm or refute (§5 predicted, not
verified). If 1b_A still fails after L3, the residual points to a genuinely
distinct mechanism beyond the live-K tangent gaps.

---

## 2026-06-23 — MS33 gating ctests: vtkdiff regression baseline added (DONE) + '--' parse-fix

DONE (commit c5615fcfd2, dsm_native_maxwell_conjugate, pushed 4 remotes):
the 6 registered MS LE gating ctests (ANCHORS_MS33_Model{I_dd1400/1600/1800,
III_gap2mm, IV_pellets, VII_freeswelling}) now carry an in-PRJ
`<test_definition>` vtkdiff regression baseline (CLAUDE.md §3, anchor (f);
Vinay-approved 2026-06-23). Mechanism (Applications/ApplicationsLib/
TestDefinition.cpp): `<regex>` is matched against the REFERENCE dir, then
vtkdiff compares output<->reference of IDENTICAL filename. Reference = the
single committed final-time output VTU per model (ts at t_end), regex pins the
end time. Two-tier tolerances (Vinay's choice), documented in a comment block
in each PRJ:
  TIGHT (geometry/equilibrium-pinned): displacement, saturation, porosity,
    transport_porosity, micro_porosity, micro_water_content, dry_density_solid.
  LOOSE (calibration/stress-path-sensitive; rel=1e-2, abs=1e3 Pa — INTENTIONAL,
    tolerates the documented ~1% Model-I non-reproducibility): sigma,
    swelling_stress, pressure, micro_pressure.
abs floors derived from extracted final-state problem scale. VERIFIED 6/6 PASS
on build maxwell_floor_20260619 (ogs @ 71366ac0d3; vtkdiff target built there);
adversarial teeth-check passed (perturbed reference -> vtkdiff FAIL; ref-vs-ref
PASS). CAVEAT: adaptive IterationNumberBasedTimeStepping pins ts_NNN to the
canonical toolchain -> a different CI toolchain may need reference regeneration.
Tests are LARGE-labelled (run via `ctest -L Large` or by name, NOT default ctest).

INCIDENT folded into the same commit: the prior provenance commit 30a1d0570f
wrote '--' (illegal double-hyphen) inside XML comments, breaking parse of 7
committed PRJs (the 3 gating Model I + confined_expulsion x2, dd1600 formB,
dd900). Fixed ' -- ' -> ' - '. The §6.7 gate did not catch it (it checks
header==live, not parse); default CI did not catch it (LARGE-excluded). Lesson +
proposed §6.7 parse-check: see memory incident_xml_double_hyphen_comment_parse.

---

## Full consistent tangent for the Maxwell local Jacobian (2026-06-09, branch dsm_maxwell_jac_parallel)

Tangent-only gap-closing (Vinay's AceGen derivation, `THM_DSM_Richards_maxwell_web.wl`); residual UNCHANGED.

- (a) Analytic micro 2×2 Jacobian in `solveReferenceMassStorageCoupledState`
  (RichardsMechanicsFEM-impl.h): replaced the 4 FD `evaluate()` calls with a
  closed-form `J = d(mass_res, dens_res)/d(n_l, ρ_lR)` (`evaluate_analytic_jacobian`
  lambda). Reuses the helper-chain μ_lR derivatives; recovers the live-nS chain by
  re-running the vdW helper with the right `dnS_dnl`. Branch on
  `fd_jacobian_for_exchange` (default analytic, FD path kept as fallback). J22=1 exact.
- (b) Enabled the u-side swelling Jacobian: `enable_dsm_swelling_up_jacobian`
  false→true (~L4394). The film-ON K[u,p]/K[u,u] block already matched the Maxwell
  identity (`d σ_sw/d ε_v = +(1-φ_M)·n_l·b·K_drained`; `d σ_sw/d n_l = -(1-φ_M)·(p_film+n_l·Π')`),
  no formula fix needed.
- VERIFIED 2026-06-09 (Model I dd1400, Maxwell film path; maxjac vs pre-edit floor binary):
  (i) solution-unchanged max rel diff = 6.748e-15 (PASS); (ii) local 2×2 cascade
  quadratic, `|R_{k+1}|/|R_k|²` ~0.24 const (4.84e-2 → 5.57e-4 → 7.23e-8 → 1.3e-15),
  analytic = FD to round-off; (iii) global iters/step identical 1.338 (max 2) pre vs
  post (this benchmark's global problem is near-linear per step, `a=1e-16` EOS-bypass —
  global count cannot move; gain is in the local cascade). SPLICE B not needed.
  Full numbers + table in `AUDIT_maxwell_local_jacobian_2026-06-09.md`.
- OPEN: a two-way-coupled benchmark (EOS active, `a` not bypassed) to measure the
  predicted global iters/step reduction — not yet exercised.
