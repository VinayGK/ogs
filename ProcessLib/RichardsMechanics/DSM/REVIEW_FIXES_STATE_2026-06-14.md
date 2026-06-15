# STATE — dsm_native_Pi_fofnlev review-fixes effort (2026-06-14)

Single-page status for the ultracode-review fix campaign. Source of truth for
the per-fix detail is `DSM/AGENTS.md` (worklog) and `DSM/ULTRACODE_REVIEW_2026-06-14.md`
(the review that scoped these fixes); regression evidence is
`task42_review_regression_2026-06-14/REGRESSION_RESULTS.md`.

## Branch / worktree

- Worktree : `/Users/vinaykumar/git/ogs-worktrees/Pi_fofnlev_review_fixes_wt`
- Branch   : `dsm_native_Pi_fofnlev_review_fixes_2026-06-14`
- Base     : `dsm_native_Pi_fofnlev` tip `9795f252e1`
- State    : **12 commits over `9795f252e1`, NOT pushed.** Tip `02d31b1cfd`.
- Build    : `~/git/build/pi_fofnlev_fixes_20260614/bin/ogs`
- Baseline test count : 39 RM unit tests pass + 2 designed skips.

## (1) The 12 fixes — one line each, with commit

| # | Fix | Class | Commit |
|---|---|---|---|
| H2 | Exact-route fold fed parse-time scalar K into `computeStrainedFilmEnergyPair` while the bare `out.mu_lR` used live `effectiveAugmentationPrefactor(phi)`; now feeds live K at L730 (matches :768/:1364/:2084). | Jacobian-only | `7e2eee0190` |
| M1 | p-u Maxwell exchange tangent always used the integrable-partner `dmu_lR_mech_deps_v`; now dispatches `dmu_lR/deps_v` on the film route (Off→integrable; operational-strained→d(bare)/deps_v − b·K_drained/rho; exact→g_cut·`pair.dmu_mech_deps_v`). | Jacobian-only | `6708ef1d98` |
| M2 + L2 | Wired the displacement-side live-K swelling-eigenstress tangent `d(delta_sigma_sw)/dK · dK/dphi · (dphi/deps_v, dphi/dp)` into K[u,u]/K[u,p] (the 1b compliant-top cure candidate). SCOPE: live-K chain ONLY; the pre-existing `enable_dsm_swelling_up_jacobian` term Vinay set OFF 2026-06-01 left untouched. Gated on `film_pressure_coupling && dK/dphi != 0 && PorosityFromMassBalance`. | Jacobian-only | `cbe3e11ed6` |
| L1 | Exact-fold `dg_dnl` mixed live-nS `out.dmu_lR_dnl` with frozen-nS `pair.dmu_bare_dnl_pre` under `CurrentPorositySplit`; recompute the bare-pre derivative with the caller's `dnS_dnl`. | Jacobian-only | `78a71ae6df` |
| N1 | Zero the live-K `dphi/deps_v=(alpha-phi)/(1+w)` chain when `PorosityFromMassBalance` clamps phi (detect via unclamped-vs-stored phi; bounds private). Applied at the live-K p-u block and the new M2 block. | Jacobian-only | `6a03a58213` |
| L3 | DOC-only: documented the live-K `dn_l/dK` local-solve strain channel in `ScalarReferenceMassStorage` mode as a DELIBERATE PARTIAL TANGENT (not wired). Prompt-authorized judgment: wiring risks the converged forward solve for a LOW-severity gap; dominant cure is M2. | Doc-only | `f23f69c5b4` |
| H1 | RESIDUAL-CHANGING (Vinay-authorized): under `film_energy_route=Exact`, source the eigenstress increment from the one-Psi `pair.sigma_sw_m` (drained-line, telescoped curr−prev) instead of the operational `Pi(w_eff)−b·p_conf`. | Residual-changing | `6391e357a2` |
| NEW TEST §8 | Added `RichardsMechanicsLiveKOfRhoD.AssembledDisplacementTangentExactKinematicLiveK` (StrainedFilmPotential.cpp) — helper-level FD-vs-analytic identity on the exact-route mu_lR eps_v tangent (H2/M1) and the live-K eigenstress eps_v tangent (M2). Anchor (d) symmetry/derived identity, scale-derived tolerances, no Vinay expected value. | Test | `9732b46498` |
| N2 + N3 | DOC-only: N2 labelled the loop-test 100× separation as MEASURED conservative floor per §5.1 (defect/bound ~3.0e3 at N=400, grows ~N²; not a derived 100). N3 re-confirmed `computeMaxwellConjugateMicroPotential` fully dead, strengthened RETIRED banner; kept on disk per §6.3. | Doc-only | `9c8423766c` |

Plus three worklog/follow-up commits on the branch (not fixes themselves):
`4ab7ce80da` (AGENTS.md worklog for the campaign), `936488482c`
(M2 gate: fire the live-K eigenstress tangent whenever live K is on),
`02d31b1cfd` (AGENTS.md 1b cure verdict).

Note on H1/H2 numbering as requested (H1, H2, M1, M2/L2, L1, L3-doc, N1, N2,
N3, new test): all ten are present above; H1 lands later in the log than H2
because it was the one residual-changing fix and was sequenced last.

## (2) Regression verdict (from REGRESSION_RESULTS.md)

Set (Vinay): MS33, EBS Task13, MGR23, MGR27, EPFL. PRE-fix binary
`~/git/build/Pi_fofnlev_20260611/bin/ogs` (HEAD `9795f252e1`) vs POST-fix
`~/git/build/pi_fofnlev_fixes_20260614/bin/ogs`.

- **MS33 LE standard suite (§12.3): residual-invariant CONFIRMED — zero regression.**
  - ModelI dd1400/1600/1800: rc=0, sigma drift **1e-14** scale (3.6e-14% / 1.3e-14% / 0%),
    **identical step sequence** (dd1400 #308 to t=1.728e7).
  - ModelIII gap2mm: 14 fields **BIT-IDENTICAL (diff 0.0)**, end t=1.728e7.
  - ModelIV pellets: 14 fields **BIT-IDENTICAL (diff 0.0)**, end t=1.728e7.
  - ModelVII freeswell: 15 fields **BIT-IDENTICAL (diff 0.0)**, end t=2.074e7.
  - **off-mode dd1400 final-VTU SHA256 identical** (`91f404a5...577`, the STEP-0
    reference) byte-for-byte across all Jacobian-only fixes.
- **MCC (Task13) + MGR23 orthogonal/unchanged:** Task13 12a/12b rc=1, MFront
  status −1 @ step #1 = pre-fix MCC integrator blockage; MGR23 robust rc=1,
  Eigen solver fail @ step #131, unchanged by the fixes.
- **MGR27 le_calk0 (LE):** POST rc=0 / PRE rc=0, head-to-head clean.

The Jacobian-only fixes change no converged result; this is the standard-suite
preservation gate that release-blocks (§12.3) and it holds.

## (3) M2 partial cure (1b compliant-top, the headline)

The M2 displacement-side live-K eigenstress tangent was the candidate cure for
the task42 1b live-K step-#1 singularity. MEASURED outcome:

- **1b_B_Kl (exact route, H1+H2+M1+M2 all active): CURED.** Was step-#1 death
  pre-fix; post-fix passes step #1 and keeps stepping (reached time step #200+,
  zero step-1/nonlinear failures; VTUs observed out to t=2.59e6 s).
- **1b_A_Kl (form (a), Off film_strain_coupling): NOT cured.** Still fails in
  time step #1 with the new binary. The M2 live-K eigenstress tangent fires for
  it (gate = live-K flag), so this confirms the review §3 "second, still
  unidentified mechanism" is NOT (only) the M2 eigenstress tangent for the
  form-(a) `ScalarReferenceMassStorage` case. Strong candidate: the **L3
  `dn_l/dK` local-solve strain channel** left as a DOCUMENTED partial (`f23f69c5b4`).
- **MGR27 kinematic_livek:** POST fail #125 / PRE fail #120 — M2 +5 steps,
  **marginal, not cured** — consistent with the same picture.

VERDICT: M2 (+H1/H2/M1) **cures the exact-route 1b_B**; the operational-route
1b_A and the live-K MGR27 family are **helped but not cured**, pointing at L3 /
a distinct mechanism as the next lever.

## (4) Open items

1. **L3 wiring (this round's deferred lever).** The live-K `dn_l/dK` local-solve
   strain channel in `ScalarReferenceMassStorage` mode is documented as a
   deliberate partial tangent, not wired (`f23f69c5b4`). It is the strongest
   candidate for the still-failing 1b_A_Kl / form-(a) live-K case. Vinay's call
   whether to wire it next; would risk the converged forward solve for a
   LOW-severity gap, hence deferred.
2. **L5 piexact provenance — NOT FIXED (deliberate, §1.1/§12.2 guardrail).** The
   two on-disk formB_piexact PRJs
   (`ms33_modelI_dd1600_formB_piexact_2026-06-12.prj`,
   `ms33_modelVII_freeswelling_formB_piexact_2026-06-12.prj`) carry `TODO(Vinay)`
   §12.2 provenance locators (E/nu, micro-EOS, Tuller geom, specific_surface,
   lambda) **inherited verbatim** from the base PRJs. They are NOT in
   Tests.cmake and STAY OUT pending Vinay's cited source locators — these
   **cannot be invented (needs Vinay's locators)**. Their calibration K=103879
   J/kg is §12.1-clean (Dixon 2023 Fig.1).
3. **Merge decision (Vinay).** The branch is 12 commits over `9795f252e1`,
   **not pushed**. Decide whether to merge into `dsm_native_Pi_fofnlev`. The
   Jacobian-only fixes are residual-invariant and bit-verified on the standard
   suite; the one residual-changing fix (H1) only affects the exact-route branch
   (off-mode unaffected, bit-verified). Pre-merge: §6.7 provenance gate, and the
   two follow-up calls below.
   - Whether to also flip `enable_dsm_swelling_up_jacobian` (the pre-existing OFF
     swelling tangent Vinay set off 2026-06-01) — kept untouched here.
   - §12.2 locators (item 2) to register the piexact PRJs.

## Cross-reference / provenance notes

- **H1 loop closure (validated):** the new assembled loop-closure probe
  `AssembledExactPairClosesOperationalSigmaDoesNot` reports |W|/scale = **8.41e-09**
  with exact sigma (post-H1, CLOSES) vs 0.0675 with operational sigma (pre-H1,
  does NOT close). Consistent with the parked feature's helper-level T-5 number
  (|∮|/scale 8.4e-9 exact vs 0.93 operational). H1 is CORRECT, not reverted.
- **§5 (predicted ≠ verified):** the residual-invariance and bit-identity of the
  standard suite are MEASURED (REGRESSION_RESULTS.md, both binaries head-to-head).
  The H1 physics tradeoff (drained-line p_conf vs actual GP p_conf — agree on the
  drained line, differ off it, e.g. fully confined) is **predicted**, not yet
  verified by a confined campaign. The 1b_A / MGR27 root cause attributed to L3
  is a **strong candidate (predicted)**, not confirmed by a wiring-and-rerun.
- **§6.7 PROVENANCE CAVEAT (from REGRESSION_RESULTS.md):** the POST-fix binary
  (linked 18:54:05) predates the final M2 gate commit `936488482c` (18:54:55, 50s
  later). The tested binary therefore contains the fixes through M2/L2 **wiring**
  (`cbe3e11ed6`) but NOT the always-on gate; the RESULTS section at the foot of
  REGRESSION_RESULTS.md states all runs settled with a post-fix binary at tip
  incl. the M2 gate — read the live-K/MGR signal against that caveat.
- **§7:** no co-author trailer on any of the 12 commits (verified: commit
  messages use "review fix" / "as specified" phrasing, no `Co-Authored-By`).
