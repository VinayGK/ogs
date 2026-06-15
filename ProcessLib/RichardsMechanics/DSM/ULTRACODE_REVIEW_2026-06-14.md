I'll synthesize the report directly from the confirmed findings. No file inspection needed — the verifiers already did the source-level work and the locators are in the payload. Let me organize, merge duplicates, and write.

The confirmed list has heavy duplication along the "exact route + K" axis (findings tagged bug1b/cpp/physics/jacobian all describe the same line-730 defect from different angles). I'll merge those and preserve the distinct mechanisms.

# Ultracode Review — `dsm_native_Pi_fofnlev` (d98f5f8324..9795f252e1)
DSM strained-film + exact one-Psi route + live-K(ρ_d) + Jacobian tangent

## 1. Verdict

**Sound to keep building on. No MS LE / §12.3 release-blocker today.** Every confirmed defect is either (a) latent — it fires only under feature combinations no *committed, registered* PRJ enables (`exact` route × live-K, `CurrentPorositySplit` × exact, `PorosityFromMassBalance` clamp × live-K), or (b) a Jacobian/test-coverage gap that degrades convergence without moving the converged root. The standard MS LE suite (Model I/III/IV/VII, LE, Operational route, scalar K) is untouched by all of them; the strained-film and live-K code paths are bit-for-bit identical to baseline in their off/frozen/clamped-edge modes.

Two things gate *building further on the headline combination* (strained-film exact + live-K), not the current suite:
- the exact route is wired with the **wrong K** (scalar, not live) — corrupts the residual the moment live-K is switched on;
- the **eigenstress half of the exact route was never wired at all** — so the one-Ψ conservation the branch exists to deliver is not present in the assembled residual even at scalar K.

Both are real now, invisible now (no PRJ exercises them), and would silently produce wrong physics the day someone combines the features the design explicitly supports.

---

## 2. Findings by severity

### HIGH

**H1 — Exact route never wired the eigenstress half; the one-Ψ pair is not in the assembled residual.**
`RichardsMechanicsFEM-impl.h:2056-2164` (eigenstress) vs `:703-746` (μ fold).
`applyFilmPressureMicroPotential` folds the exact closed-form `pair.mu_mech` into μ, but `computeReferenceMicroPorositySwellingStressIncrement` **never branches on `film_energy_route`** — it always builds σ_sw from the *operational* telescoped form `Π(w_eff) − b·p_conf`. The μ-half and σ-half are therefore gradients of two different functionals, so `dσ_sw/dn_l == (1−φ_M)·ρ_lR·dμ_mech/dε_v` (the Maxwell pair the branch exists for) does **not** hold in the assembled residual. The shipped `ms33_modelI_dd1600_formB_piexact` is fully confined (4-side Dirichlet u=0), so `p_conf` grows to the full swelling pressure while ε_v≈0 — off the drained line `p_conf=−K_d·ε_v` the mismatch is essentially the entire `nS·n_l·b·p_conf` swelling-stress drain (worked example: 87% eigenstress mismatch). The §9a operational-loop defect (|W|/scale = 0.93) the exact route was meant to cure is *not cured at assembly level*.
**Fix:** add the route switch at the eigenstress site (PI_OF_NL_EV §4 step 5): under `Exact`, source `delta_sigma_sw` from `computeStrainedFilmEnergyPair`'s `sigma_sw_m`/`dsigma_sw_dnl` (one helper call per GP, both halves from one Ψ). The doc's "eigenstress site unchanged by design" is incorrect — the plan itself required wiring `sigma_sw_m`.
*Interaction:* H1 is the σ-channel sibling of H2/M1 (μ-channel). Fix both at the same site or the pair stays inconsistent.

**H2 — Exact route feeds scalar K (not live K) into `computeStrainedFilmEnergyPair`.**
`RichardsMechanicsFEM-impl.h:730` (consolidates bug1b/cpp/physics duplicates).
Line 730 passes `potential_exchange_params.potential_augmentation_prefactor` (parse-time scalar) while the bare `out.mu_lR` it folds against was built with `effectiveAugmentationPrefactor(...,phi)` (live K) at `:1364`, and the operational route at `:768` and eigenstress at `:2084` both use live K. This is the **only remaining raw-scalar prefactor feeding a potential evaluation** in the file. Consequence under live-K: `g_cut = out.mu_lR / pair.mu_bare_pre` (`:738`) divides a live-K bare by a scalar-K bare, so it stops being the macro-floor cutoff factor (numeric check: g_cut → 3.46 instead of 1.0, 246% corruption at K_live=103879 vs K_scalar=29999), and `out.mu_lR += g_cut·pair.mu_mech` (`:743`) writes a corrupted residual; `dg_dnl` (`:739`) inherits it. The two tracked piexact PRJs use scalar K, so `g_cut==1` exactly today → latent in CI, but the branch was rebased onto live-K and K_OF_RHO_D_LIVE.md treats live-K as branch-wide.
**Fix:** replace line 730 with `effectiveAugmentationPrefactor(potential_exchange_params, local_context.phi)`; **or** add a create-time `OGS_FATAL` forbidding `film_energy_route=exact` + `potential_augmentation_prefactor_live_dry_density=true` (the create-time validator currently checks only `exact ⇒ kinematic`, not this coupling). Add a test exercising exact × live-K table (none exists).

**H3 — Exact energy route has no run-level/ctest coverage; 6 new piexact PRJs unregistered.**
`Tests/.../ExactFilmEnergyPair.cpp:1-331`.
All 8 exact tests call `computeStrainedFilmEnergyPair`/`computeIntegrableMechanicalMicroPotential` directly; none drives a mesh/solve. The assembled fold (`:703-748`), the eigenstress site, and the global Jacobian wiring are covered by **no** automated test. The 6 `ms33_*_formB_piexact_2026-06-12.prj` are added with zero cmake change. T-1 "OffModeBitForBit" is documented as run-level but realised only as `EXPECT_EQ(params.film_energy_route, Operational)` — a default-value check, not a field-data regression.
**Fix:** register one exact-route PRJ (dd1600 formB_piexact) as a run-only ctest; add the design's T-1 off-mode bitwise regression to CI; or `GTEST_SKIP`-document that the assembled exact path is verified only by manual campaign.

### MEDIUM

**M1 — p-u Maxwell exchange tangent always uses the integrable-partner `dmu_lR_mech_deps_v`, even when the residual uses the strained/exact μ.**
`RichardsMechanicsFEM-impl.h:4427-4496`.
When `film_strain_coupling != Off`, the residual μ depends on ε_v through `w_eff/mu_load` (operational) or `pair.mu_mech` (exact), but the p-u block unconditionally builds `mech_pu = computeIntegrableMechanicalMicroPotential(...)` at unstrained n_l and never dispatches on the route. The exact partner's `dmu_mech/dε_v = −2μ_v/(1+x)³ + μ_a·E·(1−ξ0−ξ0·x) − b·K_d·ε/ρ` differs from the integrable partner's `−(Π+n_l·Π'+b·K_d·ε)/ρ` by O(κ·ε_v); they agree only in the T-4 κ·ε_v→0 limit. The piexact PRJs run kinematic+exact and swell to finite ε_v, so the analytic Jacobian (no `jacobian_assembler` tag → analytic is used) is finitely inconsistent with the residual. **Residual stays correct (converged answer right *if* it converges); convergence degrades — consistent with the III/IV stagnation in the run-status board.** Design doc §4 step 6 explicitly required `dmu_mech_deps_v` "wired into the same slots the operational mode fills" — it is not.
**Fix:** in the p-u block, when `film_strain_coupling != Off`, source `dmu_lR/dε_v` from the route: exact → `pair.dmu_mech_deps_v` (× cutoff g); operational-strained → `d/dε_v[bare(w_eff)] + b·dp_conf/dε_v/ρ`.

**M2 — Live-K swelling-eigenstress K(φ(ε_v)) dependence is in the residual but absent from the displacement Jacobian — the "second mechanism" behind 1b compliant-top divergence.** *(See §3 for the root-cause writeup.)*
`:2357-2365` (residual) vs `:4498-4564` / `:4665-4727` (Jacobian).

### LOW

**L1 — Exact-fold `dg_dnl` mixes live-nS `out.dmu_lR_dnl` with frozen-nS `pair.dmu_bare_dnl_pre` under `CurrentPorositySplit`.** `:739-741`. `out` carries the `dnS_dnl=−1` chain (`:1353`); `pair.dmu_bare_dnl_pre` is computed with `dnS_dnl=0` (`:732`). Value `g_cut` is unaffected (dnS-independent); only `dg_dnl` is wrong. Committed piexact PRJs default to `Reference` mode → dormant. **Fix:** recompute one leg with the caller's `dnS_dnl`, or assert exact-route ⇒ Reference mode.

**L2 — Swelling-stress u-u/u-p Jacobian is dead code (`enable_dsm_swelling_up_jacobian=false`) and its live-K chain is admittedly omitted.** `:4665-4800`. Residual eigenstress fully uses live-K (`:2127-2163`); the only strain tangent for it is constexpr-false dead code that explicitly carries no dK/dφ chain → residual/Jacobian inconsistent under live-K with no test. **Fix:** flip+wire (owner call) or add an FD-vs-analytic consistency test that currently fails, and flag the dead block in AGENTS.md as untested-not-silently-dead.

**L3 — In `ScalarReferenceMassStorage` mode (the 1b PRJ) n_l is solved locally with K(φ); the K-channel holds n_l fixed, omitting the implicit `dn_l/dK·dK/dφ·dφ/dε_v` chain.** `:907-916, 1085-1094, 4546-4562`. Consistent with the established fixed-n_l strain-channel approximation, but compounds the compliant-top inconsistency. **Fix:** add `dn_l/dK` via the existing local-solve sensitivity machinery, or document as deliberate partial tangent.

**L4 — `StrainedFilmEnergyPairData` strain-tangent/eigenstress fields (`dmu_mech_deps_v`, `sigma_sw_m`, `dsigma_sw_dnl`, `Psi_film`) are computed but consumed only by tests.** `PotentialExchange.h:795-808`. Proximate cause of M1/H1 — the exact route's Maxwell-consistent tangent is built then discarded. **Fix:** wire them (per M1/H1) or annotate explicitly test-only so the gap is documented, not silent.

**L5 — formB_piexact PRJ §12.2 provenance fields read `TODO(Vinay)`** (E/ν, micro-EOS, Tuller geom, specific_surface=523 m²/g, λ=7.5e-7). `ms33_modelI_dd1600_formB_piexact_2026-06-12.prj` header. *Inherited verbatim from the registered base PRJs* (literal-tag diff empty), and these files are **not** in `Tests.cmake`, so §12.3 CI-gating does not fire. Calibration K=103879 J/kg **is** §12.1-compliant (Dixon 2023 Fig.1, EMDD≡ρ_d, target 14.161 MPa) and the header K equals the live tag (§6.7 holds). *(See §4.)*

### NIT

**N1 — Live-K p-u `dφ/dε_v = (α−φ)/(1+w)` ignores the `PorosityFromMassBalance` min/max clamp** (`:4530-4543`); spurious nonzero tangent at clamped φ. Fires only when live-K active **and** φ interior to the table **and** φ clamped — narrow, converged result unaffected. **Fix:** zero the chain at `phi_min/phi_max`.
**N2 — Loop-test operational-defect floor hard-codes a ×100 separation over the quadrature bound** (`ExactFilmEnergyPair.cpp:298-313`); the 100 is fitted to this Pi-dominant loop, not derived. **Fix:** derive the O(Π·κ·ε_v) leading term, or label "measured" per §5.1.
**N3 — `computeMaxwellConjugateMicroPotential` is fully dead code** (`PotentialExchange.h:333-397`), already banner-marked RETIRED 2026-06-08; flagged for completeness only.

---

## 3. The 1b live-K compliant-top divergence — root cause (predicted, not verified)

**Mechanism.** In live-K mode the swelling eigenstress σ_sw depends on K(φ(ε_v)) through `K_aug_sw = effectiveAugmentationPrefactor(params, total_porosity)` in the *residual* (`:1977, 2108-2163`, entering R_u via `volumetric_mechanical_strain`). The analytic live-K tangent was wired into **only the pressure-displacement block** (`:4554`); the displacement-side swelling tangent is the constexpr-false dead block (L2) and explicitly carries no dK/dφ chain. So **`dσ_sw/dε_v|_K` is missing from K[u,u] and K[u,p] everywhere**, while the residual fully contains it.

**Why only compliant-top × live-K:** on a rigid lid ε_v≈0 → the term is inert (lidTop completes, 24 its); under scalar K σ_sw has no K(φ) dependence → residual==Jacobian (out_1b_A_Ks completes, 1005 steps); only live-K × compliant-top makes the missing mechanical tangent large. Magnitude estimate: `dσ_sw/dK ≈ (1−φ_M)·n_l·ρ·sign·exp(−ξ)`, `dK/dφ = −ρ_SR·slope ≈ −8.2e5 J/kg` (1b knots 45217/103879/266767 at ρ_d 1400/1600/1800, ρ_SR=2780), `dφ/dε_v ≈ 0.55` → an **O(1e8 Pa) eigenstress-strain coupling absent from K[u,u]**. Corroboration: out_1a_robin_A_Kl (live-K, mild) completes but 32517 vs 949 steps (same mechanism, mild form); the failing out_1b_A_Kl log shows the pressure component diverging monotonically (1.66e7→…→3.66e13) — an inconsistent global tangent, not overshoot (damping 4/10 did not rescue).

**Proposed minimal fix (predicted-not-verified):** add the symmetric displacement-side live-K tangent — `dR_u/dε_v` and `dR_u/dp` via `dσ_sw/dK · dK/dφ · (dφ/dε_v, dφ/dp)`, assembled as `B^T C C_el^{-1}(…)·identity2` into K[u,u]/K[u,p], hosting it in the existing OFF block (`:4672`) once the K-channel is added there and the block enabled. **Robust fallback:** freeze K per time step (parse-time lag / Picard on K), as K_OF_RHO_D_LIVE.md already proposes. *Predicted to restore Newton consistency on the compliant top; K_OF_RHO_D_LIVE.md records the cure test as NOT CURED with the current wiring, so this is a hypothesis, not a verified fix.* Note the residual is correct either way — this is purely a convergence repair.

**Honest residual uncertainty:** K_OF_RHO_D_LIVE.md attributes the 1b failure to "a second, still unidentified mechanism." M2 is the strongest-evidenced candidate (worked magnitude, completes-on-rigid-lid control, scalar-K control, mild-form corroboration), but the doc's own language stops short of confirming it, and no FD-vs-analytic Jacobian probe has been run to close it. Treat M2's identification as predicted.

---

## 4. Guardrail ledger

- **§12.2 (FLAG, not STOP) — L5.** The two/six formB_piexact PRJs carry `TODO(Vinay)` provenance locators (E/ν TH-G table, micro-EOS a/b primary source, Tuller geom, specific_surface, λ co-calibration note). §12.2 lists `TODO`/`unknown`/`standard` as non-compliant *for Tests.cmake registration*. **They are not registered** (`grep piexact *.cmake` empty; `Tests.cmake` diff empty), so §12.3 CI-gating does **not** fire. The TODOs are **inherited verbatim** from the registered base originals (literal-tag diff empty) — a propagated pre-existing gap, not a new violation. Calibration K=103879 J/kg is §12.1-clean (Dixon 2023 Fig.1, EMDD≡ρ_d). **Action:** resolve locators before any registration/submission; keep unregistered until then (current state). Committing/registering as-is would trip §6.7/§12.2.
- **§1.1 — clean.** No uncited material/expected/BC literal introduced; the one calibration literal (K) traces to Dixon 2023 with the §12.2 anchor block intact and header==live-value.
- **§2 (no fit-and-verify) — clean.** No test both fits K and asserts on a K-derived quantity; the new tests assert against analytical limits, FD of the residual, and loop-closure identities.
- **§4 (`+=` accumulators, units) — clean.** The fold uses `out.mu_lR += g_cut·pair.mu_mech` (`:743`, additive, correct); no `=` overwrite on an accumulator found. (The bug in H2/H1 is *which K* is used, not the accumulation operator.)
- **§5 (predicted≠verified) — observed in docs, not violated by code.** K_OF_RHO_D_LIVE.md correctly labels the live-K run as "no measured gain" and the 1b cure as "NOT CURED" — predicted/verified separation is honoured. The one place to watch: PI_OF_NL_EV_IMPLEMENTATION.md's "eigenstress site unchanged by design" is *factually wrong* (H1) rather than a mislabeled prediction — correct the statement, not the label.

---

## 5. Test-coverage gaps + proposed missing tests (structure only, no expected values)

Confirmed gaps, in priority order:

1. **Exact-route fold at assembly level (H3, M1, H2).** No test drives `applyFilmPressureMicroPotential` exact branch; g_cut recovery, the `+=`/dg_dnl folding, and the exact × live-K path are all unverified above the helper.
2. **Live-K assembled p-u tangent.** Only `effectiveAugmentationPrefactor`/`PhiDerivative`/`getSegmentSlope` and mu-level `dmu_lR_dK` are FD-tested; the assembled `dmu_lR_dK·dK/dφ·dφ/dε_v` insertion into `local_Jac` (sign, the (α−φ)/(1+w) chain, the n_l·ddmu telescoping) is not.
3. **Live-K under a free/compliant boundary with evolving ρ_d** — exactly the 1b case the feature targets, recorded NOT CURED, invisible to the suite.
4. **Floor / positivity-guard branches** of `computeStrainedFilmEnergyPair` and the vdW law — every sample sets `floor=0` and ε_v ≥ −0.03, so the `floored` and `clamped_f` derivative-freeze branches (PotentialExchange.h:842-921) get zero coverage despite the floor being a shipped, campaign-used (litfloor) feature. Design T-3 required "ON and OFF the guard" states — the ON states are absent.
5. **`ReplacementIsExclusiveAtZeroStrain`** uses only `EXPECT_NE` (`StrainedFilmPotential.cpp:284`) — cannot distinguish "partner replaced" from "partner added on top"; the no-double-count property is unverified. Fix: compare Kinematic-operational μ against the closed-form `bare(w_eff)+b·p_conf/ρ` it already builds at `:256-259`.
6. **Equilibrium loaded-branch derivatives** (`dw_eff_dnl`, `dw_eff_deps_v`, hardcoded 0) never FD-checked; inverter asserts only to a tolerance looser than its own. Partly inherent — at minimum document the intentional-zero tangent.
7. **Maxwell loop closure** is verified on a helper-reconstructed integrand (g==1, frozen K), not on the assembled residual that folds μ through g and live-K — the central conservation claim is verified one level below where it ships.

**Proposed test (per CLAUDE.md §3, structure only):**
```
Name:            LiveK_CompliantTop_EvolvingRhoD_JacobianConsistency
Physics anchor:  (d) symmetry / derived identity — assembled tangent vs central FD of R w.r.t. nodal u
Input config:    single axisymmetric element; free/Robin top; live-K table spanning the
                 swelling rho_d range; PorosityFromMassBalance; one Newton step
Probe:           rel. error of analytic local_Jac displacement-column vs central FD of
                 rho_L_hat (and, separately, run under film_energy_route=exact to also
                 catch M1/H1)
Expected:        TODO(Vinay) — owner sets threshold
Catches:         missing dK/dphi*dphi/deps_v in the mechanical (eigenstress) block (M2);
                 frozen-vs-live K in the exact fold (H2); integrable-vs-exact partner (M1)
Overlap:         none
```
Plus a fold-level companion: build params with `exact + kinematic + live-K table`, call `applyFilmPressureMicroPotential`, assert (a) g_cut recovers the macro-floor cutoff factor and (b) the augmentation channel uses live K — **fails today**, exposing H2/H1.

---

## 6. Checked and found clean (informative absences)

- **MS LE standard suite (Model I/III/IV/VII, LE, Operational route, scalar K):** no defect touches it. All flagged code is bit-for-bit identical to baseline in off/frozen/clamped-edge modes (off-mode bitwise regression intent confirmed, dd1400 bitwise per the design doc).
- **Live-K analytic chain, component-level:** independently re-derived clean — `dφ/dε_v = (α−φ)/(1+w)` is the exact derivative of `PorosityFromMassBalance` (verified by substitution; `α`=BiotData matches the porosity law, `β_SR` matches grain compressibility); `dK/dφ = ρ_SR·getSegmentSlope` is FD-verified to 1e-9 against the residual value. The defects are in *what's wired where* (M2, L2, L3), not the math.
- **`getSegmentSlope` C0 discontinuity at knots/edges** — *examined and dismissed.* It is the truthful, FD-verified derivative of the clamped piecewise-linear residual (correct-by-design, not an inconsistency); causality is inverted (it's the documented *fix* for the first-cut missing-tangent, not the cause); the measured cure test shows Newton contracts 7 iterations with it active. C1 smoothing is a flagged-undecided future option, not a correctness fix.
- **dφ/dε_v clamp/(1+w)-singularity escalated to "tangent wrong exactly where 1b drives φ to the clamp"** — *dismissed.* The faulted line is unreachable in every committed PRJ (live flag defaults false, no PRJ enables it), self-skips at the clamp via the coupled `ρ_d=ρ_SR(1−φ)` table-edge-zero, and self-mitigates at α=1 in the cited PRJ. The surviving narrow version is N1 (nit).
- **"Small-strain K-channel destabilizes at ε_v O(1–10)"** — *dismissed as a fabricated premise.* The MS33 column is 0.070 m, not 7.8 m; actual displacement norms are O(1e-5 m) → ε_v O(1e-4), four-to-five orders below the claimed magnitude; the partner form is the exact ½bK_dε² free-energy term, not a truncated Taylor series.

**Net residual uncertainty:** the two HIGH findings (H1/H2) are code-evident and unambiguous. M2 (the 1b root cause) is strongly evidenced but the upstream doc still calls the mechanism "unidentified" and no FD Jacobian probe has closed it — carry it as predicted. Everything below MEDIUM is latent or coverage, not active in any committed run.