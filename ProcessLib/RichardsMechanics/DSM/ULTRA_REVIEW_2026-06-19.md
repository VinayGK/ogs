<!-- Ultra review (multi-agent, 38 agents/2.9M tok) of maxwell tip 71366ac0d3, run wf_164b2bde-577, 2026-06-19. 31 findings, 28 verified, 3 refuted, 0 blockers. -->

# Adversarial Review — OGS Maxwell DSM merge (tip `71366ac0d3`)

Worktree: `/Users/vinaykumar/git/ogs-worktrees/dsm_native_maxwell_conjugate_wt`
Scope: new-physics code (exact one-Psi pair, live-K, strained-film), unit tests, MS33 benchmarks/ctests.

## 1. Verdict

The tip is **not confirmed CI-green** and **not confirmed releasable** — but no defect found here is a CI break or wrong-physics blocker. The new-physics math is correct (MATH-1/3/4, F4 mode-guard all pass as positive controls), the registered MS33 suite parses (floor tags present on all 9), and every flagged item is either a doc/comment imprecision, a provenance gap, or a test-coverage gap. **The single most important thing to fix: the tip is unbuilt.** The only `ogs` on the system (`maxwell-conjugate-20260602`, Jun-17) predates the Jun-19 macro-mandatory commit, so the mandatory-floor parse path and the floored MS33 run-to-completion have never been compiled or executed. Build the tip and run the 9 MS33 ctests before any "verified" claim — everything else (provenance headers, run-only test coverage) is owner-decision debt that does not gate the build.

## 2. BLOCKERS (CI-breaking or wrong-physics)

**None.** No finding in the verified set breaks the build, fails parse on a registered ctest, or encodes wrong physics in a residual-bearing term. The mandatory-floor guard (the headline risk of this merge) does **not** break the registered suite — see F7 below.

## 3. MAJOR findings

### M1 — Tip is unbuilt; mandatory-floor parse path + floored MS33 runs never compiled/executed
`ProcessLib/RichardsMechanics/CreateRichardsMechanicsProcess.cpp:706-745` (was JSON id F2/floor-dim)
HEAD `71366ac0d3` (2026-06-19 02:59:34) added the no-default `getConfigParameter<double>("macro_porosity_floor")` arm. The system `ogs` and `libRichardsMechanics.dylib` are dated Jun-17 17:16 — the dylib *postdates* the micro-mandatory commit (so the micro-floor path is compiled) but *predates* the macro commit (`: 0.0` default still in place at that state). `find … -newermt '2026-06-19 02:59:34'` returns nothing.
**What breaks:** the "missing-key now fatals" and "0.08 floor runs clean" claims are predicted, not verified; the §12.3 gating suite has never been run to completion under the mandatory floors.
**Fix:** build the tip, run the 9 registered MS33 OgsTests, confirm (a) parse under the mandatory floors and (b) completion without solver failure at `7.06e-5`/`0.08`. Until then label both claims predicted (§5).

### M2 — `micro_water_content_floor=7.06e-5` is an orphan literal on the registered suite
`Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1400.prj:112` (JSON MS33-PROV-01; overlaps floor-dim F3)
Line 112 `<micro_water_content_floor>7.06e-5</micro_water_content_floor>` carries **no inline source** and is **absent from the §12.2 header** (header block lines 37-61 lists 13 groups, neither floor). Its sibling `macro_porosity_floor=0.08` (line 113) *does* carry `(Vinay 2026-06-17: MS33=0.08)`. Same pattern: dd1600:112, dd1800:111, IIIgap2mm:118, IVpellets:141, VIIfreeswelling:129. It is a dimensioned clamp that shapes the dry-endpoint micro solve.
**Why major:** §1.1 STOP / §6.7 item-3 orphan-literal — the gate would fire at commit. (Not a parse/physics break; runs work.)
**Fix:** get a §1.1 source (most likely a Vinay approval like the 0.08 macro floor, or an in-file derivation from `n_l0`/`w0`); add an inline comment on the line in all 6 form-(a) PRJs; add a "disjoining floor" group to the §12.2 header. Also add `0.08`/`7.06e-5` to the §12.2 block mirroring the E=52 WARNING pattern (the macro floor is uncited-in-literature).

### M3 — Standard-suite ctests are run-only; §12.3 "regression baseline" asserts nothing on values
`ProcessLib/RichardsMechanics/Tests.cmake:31-47` (JSON F1/benchmark-dim; same as MS33-PROV-05)
All 9 ANCHORS_MS33 registrations use `OgsTest(PROJECTFILE … RUNTIME N)` — no `TESTER vtkdiff`, no `DIFF_DATA`; `OgsTest.cmake:70-72` runs `-r` only and `AddTestWrapper.cmake` has no diff stage. `find ANCHORS_MS33_Model{I,III,IV,VII} -name '*.vtu' | grep -iE 'ref|expect'` → 0 reference VTUs. Contrast the `AddTest` at Tests.cmake:56 (`double_porosity_swelling_dsm_micromacro_constbc`) which *does* carry `vtkdiff DIFF_DATA`.
**Why major:** §12.3 calls this the standard regression baseline and a CI break a release blocker, but a run-only test catches only crashes/`OGS_FATAL` — a converging sign-flip or unit error passes green. This compounds with the documented dd1400 K non-reproducibility (header lines 8-9: clean code → Ps 4.9218 at K=45217; dirty re-fit → 5.0516), exactly the drift a run-only test cannot catch.
**Fix (owner call):** either (a) document run-only as intentional in §12.3, or (b) freeze reference VTUs for the 6 form-(a) PRJs and convert I/III/IV/VII to `AddTest` with `vtkdiff`. **Do not auto-generate references** — that is fit-and-verify in one artifact (§2); the reference endpoint must be the Dixon/Villar target asserted independently of K. dd1400 needs special care given its documented divergence — freeze only from a run you certify, cite the commit.

## 4. MINOR / nits (terse)

- **MATH-2** `PotentialExchange.h:796-798` — helper comment says the FEM eigenstress site is "unchanged … with the actual p_conf"; under `Exact` it is **replaced** with the drained-line `sigma_sw_m` (`RichardsMechanicsFEM-impl.h:2162-2236`). Two comments contradict; reword the helper one. *(docs)*
- **F3 (impl-dim)** `RichardsMechanicsFEM-impl.h:2162` — exact-route eigenstress on drained line `p_conf=-K_d*eps_v` vs operational actual-GP `p_conf`; deliberate one-Psi tradeoff, correctly §5-labelled *predicted*. Keep predicted-vs-verified status in the run report until an MS33 exact-route confined-path re-run confirms it. *(maxwell)*
- **MS33-PROV-02** `ms33_modelIV_pellets.prj:139` — `K_pellet=20600` fit to user-set 0.400 MPa (Dixon extrapolation gives 0.35; ρ_d=900 below measured range). Header is transparent ("§12.1 extrapolation edge"); §6.7 header==live holds. Confirm the user target stands as §1.1 source; flag at commit per §12.5 that K is anchored to a user target, not a Dixon/Villar row. *(provenance)*
- **MS33-PROV-03 / F3(bench)** dd1400:56,58,59 (and III:66, IV:70/372, VII:66) — `TODO(Vinay)` placeholders for micro-EOS, Tuller geom/`char_pore_size=10um`, `specific_surface=523 m2/g`. §12.2-non-compliant strings, but params are deliverable-inert (BishopsSaturationCutoff=1; `micro_liquid_density_a=1e-16` bypasses EOS). Supply locators (SSA → Seiphoori 2014/EPFL per memory). *(provenance)*
- **F2 (MCC-dim)** `ANCHORS_MS33_MCC_ABSP_SUITE/.../ms33_modelI_dd1400.prj:63` — dead `accumulate_swelling_contributions` tag in 7 MCC-suite PRJs (none registered). Would `OGS_FATAL` if any are registered/copied. Annotate as historical/non-runnable per §6.3 before any reuse. *(parse, latent)*
- **F1 (impl-dim)** `RichardsMechanicsFEM-impl.h:842` — strained-film branch uses `out.mu_lR = …` (replacement) vs `+=` in sibling branches (L775, L879). Intentional supersede (correct, §4.1 carve-out), but add a "REPLACE not += (by design)" marker for auditability. *(accumulator)*
- **F2 (impl-dim)** `RichardsMechanicsFEM-impl.h:744` — `g_cut` denominator comment attributes nonzero-ness to the FATALs; the operative reason is sign-uniformity (core + augmentation share `potential_sign_factor`). Tighten the comment or add a defensive `OGS_FATAL` guard. *(guardrail/comment)*
- **F5 (merge-dim)** `RichardsMechanicsFEM-impl.h:2203` — eigenstress recomputes `active_nS` via default-constructed context; safe only because `computeActiveMicroSolidVolumeFraction` ignores `local_context`. Add an invariant comment or thread the live context. *(fragility)*
- **TQ-1** `Tests/ProcessLib/RichardsMechanics/ExactFilmEnergyPair.cpp:325` (and :406) — operational-defect separation factor `100.0` is a hand-picked literal (measured ~3.0e3 at N=400), not scale-derived like the exact-side `50/N^2`. Discriminates today; only a ~30× regression slips through. Derive from the leading `O(Pi*kappa*eps_v)` coeff or assert `|W_op|/Wabs_op ≥ c`. *(test tolerance)*
- **TQ-2** `Tests/ProcessLib/RichardsMechanics/StrainedFilmPotential.cpp:596` — the (A) H2/M1 leg's header claims floor-active discrimination, but the sample state keeps `g_cut==1` (cutoff inactive, floors default 0.0). The floor-active `g_cut!=1` path the comment advertises is untested anywhere. Add an (A)-leg instance with a nonzero floor, or correct the comment. *(coverage)*
- **F4 (bench-dim) / F7** positive controls — `isValidFilmEnergyRouteCombination` mode-guard correct (`PotentialExchangeParameters.h:124-129`); §12.3 standard set fully registered with both floor tags on all 9 PRJs, no dead tag (`Tests.cmake:31-47`). Optionally add the (Exact,off)→fatal / (Exact,kinematic)→ok unit test (no built binary exercises it — see M1).
- **F4/F5(bench), TQ-3/4/5, MATH-3/4** — positive controls, no action: E=52 WARNING present and well-formed (dd1400:285, header 39-47); experimental scratch PRJs (dd900, _calib, formB_piexact, v0..v7, …) orphaned but unguarded inside registered dirs (annotate/relocate per §6.8); test literals are derived structural probes or scale-relative parity floors; GTEST_SKIPs (T-6/T-8) correctly parked on Q3/Q4; live-K partials, `getSegmentSlope`, `invertDisjoiningPressure` all exact.

## 5. Coverage gaps — what this static review could NOT establish

This is a **source + on-disk static review on an unbuilt tip**. It did not, and could not, establish:

1. **That the tip compiles** (no build reflects `71366ac0d3`).
2. **That the mandatory-floor `OGS_FATAL` actually fires** on a missing key (parse arm unbuilt — M1).
3. **That the 9 registered MS33 ctests complete** under `7.06e-5`/`0.08` (suite never run on a mandatory-macro binary).
4. **Bit-for-bit invariance** of III/IV/VII vs the prior deliverable under the floor-pairing change (claimed in memory, not re-confirmed here).
5. **MATH-1 `kappa→0` cross-block reduction (UNCERTAIN — needs runtime/symbolic check):** ~5/6 of the exact-pair claim reproduces symbolically, but `dmu_mech_drho_lR` does **not** reduce to the integrable partner's `drho` block (the partner's `drho` is the inexact one; the pair's is exact). Impact on the residual is nil (FD Jacobian per the in-file F2 note), but the "reduces EXACTLY across all blocks" assertion is refuted for the `drho` tangent. Confirm with a direct FD `drho` cross-test before relying on the reduction claim in any doc.

Recommended commands to close gaps 1-4 (run from the build dir against this worktree):

```
# configure + build the tip
cmake --build /Users/vinaykumar/git/build/<tip-build-dir> --target ogs -j$(sysctl -n hw.physicalcpu)

# parse + run the 9 registered MS33 ctests (run-only; confirms parse + completion)
ctest --test-dir /Users/vinaykumar/git/build/<tip-build-dir> -R 'ANCHORS_MS33' --output-on-failure -j$(sysctl -n hw.physicalcpu)

# negative control: confirm the mandatory-floor fatal actually fires
#   temporarily strip <macro_porosity_floor> from a scratch copy of a DSM PRJ and run ogs on it;
#   expect OGS_FATAL "Key <macro_porosity_floor> has not been found"
```

## 6. Prioritized action list

1. **Build the tip + run `ctest -R ANCHORS_MS33`** — converts M1 and gaps 1-4 from predicted to verified. Add the missing-key negative control to prove the new fatal fires. *(blocker on "verified" status)*
2. **Source `7.06e-5`** (M2) — a Vinay approval analogous to the 0.08 floor, plus inline comment in all 6 form-(a) PRJs and a §12.2 "disjoining floor" group. Clears the §6.7 orphan-literal gate so the suite is commit-clean.
3. **Decide run-only vs vtkdiff for §12.3** (M3) — either document run-only intent in §12.3 or baseline Dixon/Villar-target reference VTUs (independent of K, not auto-generated). dd1400 last, with a certified run.
4. **Confirm MATH-1 `drho` reduction** with a direct FD cross-test (gap 5); correct the "reduces EXACTLY across all blocks" wording.
5. **Provenance sweep** — resolve `TODO(Vinay)` locators (MS33-PROV-03), confirm IV pellet user-target stance (MS33-PROV-02), add `0.08`/`7.06e-5` WARNING-style §12.2 entries.
6. **Comment/hygiene cleanups** (batchable, no functional change): MATH-2 helper comment, F1/F2/F5 impl markers, TQ-1/TQ-2 test fixes, dead-tag annotation (F2 MCC), scratch-PRJ status suffixes (F6).