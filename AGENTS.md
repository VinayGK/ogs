# AGENTS.md — repository root

This file is read by Claude / Codex agents at the repository root.

The repository contains several long-lived feature lineages. Branch-specific
goals are recorded here at the top so any agent that lands on this tree knows
*what this branch is for* before editing.

For DSM physics and the hierarchical micro–macro framework see
[ProcessLib/RichardsMechanics/DSM/AGENTS.md](ProcessLib/RichardsMechanics/DSM/AGENTS.md).

**Canonical design doc for this branch:**
`/Users/vinaykumar/tex/cc2024/VK_B35_Pinion_May_2026/tuller_macro_wrc.tex`
(the "Two-Branch (Tuller) Macro WRC" beamer, Kumar & Nagel, May 2026). That
deck is authoritative for the physics, the closure, the Truesdell–Noll
standing, and the four decision points reproduced below. Where this AGENTS.md
and the beamer disagree on a *branch name*, see the explicit inconsistency
note in "Implementation scope".

---

## Branch: `dsm_native_tuller_macro_film` — Option B (Tuller corner + film)

**Scope.** This branch is the home of **Option B**: macro-scale adsorption via
a Tuller corner-plus-film coexistence closure. It is a **deliberately parked
alternative**, not blocked work. The decision the beamer records is already
made:

- **Option A (cavitation cutoff) is ACTIVE** on `dsm_native_hierarchical` —
  the pragmatic, paper-consistent, MFront-free fix that simply enforces what
  paper §2.4 already states (macro pores drain into the micro structure at the
  cavitation pressure).
- **Option B (this branch) is preserved for later**, if/when macro-scale
  adsorption is ever actually required. It is physically richer and scores
  better under Truesdell–Noll, but it carries architectural cost (below).

**Native C++ only on this branch.** The work here is the OGS MPL property +
its factory + its unit test. The MFront port and the native↔bridge parity
are a **follow-on**, on the `dsm_mfront` lineage, and must land *before any
merge* — see "MFront follow-on". The MFront DSM bridge does **not** exist on
this native branch (no `RichardsMechanicsDSMMicroMacroBridge.mfront`, no
tuller parity ctest, no tuller PRJ here).

**Implementation discipline (Vinay, 2026-05-29): in-situ, no new files.**
Extend the existing class / factory / unit test in place. Do **not** add a
`SaturationTullerWithFilm` sibling class, a new gtest file, or new MFront
sources. Existing PRJ test data is re-used or modified in place; new PRJ test
data requires explicit clearance from Vinay first.

**Created.** 2026-05-29 (Vinay). No implementation commits yet — this
document is the entry point. Intent is to merge Option B back into the
hierarchical mainline once it has converged *and* MFront parity holds.

---

## The problem Option B fixes

(beamer §"The problem")

A relative humidity is imposed; the Kelvin equation converts it to a *total*
suction; the equilibrated probe sits at one total potential everywhere:

```
ψ_Macro = ψ_Micro = ψ_total ~ O(10^2) MPa     (driest bentonite: w ≈ 14%)
```

- **Micro is well-posed.** The adsorptive (vdW / electrochemical) potential
  `ψ_Micro(ω_l)` is defined down to a monolayer; the inversion `h(ψ_Micro)`
  never hits a wall.
- **Macro is *not* well-posed under a capillary-only interpretation.**
  Young–Laplace `ψ_cap = −2σ/(ρ^L r)` needs a meniscus. At ψ=100 MPa the
  implied radius is `r = 2σ/(ρ^L ψ) ≈ 1.4 nm`, far below any macro-pore
  feature. Capillarity *cannot* generate this suction; the macro pores are
  drained and the capillary potential at `S_Macro → 0` is undefined.

The current model evaluates a pure van Genuchten capillary curve (`p_b = 27`
MPa) at ~100 MPa — ~4× its scaling pressure, deep in its drying *tail*, where
the water it reports is physically adsorbed film, and where VG has no film
term to account for it. Option B replaces that artefact with an explicit
film branch.

---

## Option B closure (fully specified by the beamer)

(beamer §"Option B (alternative): the Tuller two-branch closure")

Macro retention = capillary **corner** water + adsorbed **film** water,
**coexisting at every potential**, both smooth functions of the single macro
potential ψ — **no air-entry threshold** (Tuller 1999, Or 1999):

```
S_Macro(ψ) = S_cap(ψ) + S_film(ψ),     ψ < 0
```

**Corner (capillary) branch.**
```
r(ψ) = −σ/(ρ^L ψ)        (defined for all ψ < 0)
```
Corner liquid area ∝ r², upscaled over the pore-size distribution f(L) gives
`S_cap(ψ)`. It shrinks *continuously* to 0 as ψ→−∞ — no snap-off, no cutoff.
**The existing single-branch Tuller `1 − exp(−C/ψ²)` *is* this upscaled corner
term** — i.e. the corner branch is already in `SaturationTuller.cpp`; Option B
*adds* the film branch to it.

**Film (adsorptive) branch.**
```
ψ_film = −A_H / (6π ρ^L h³)   ⇒   h(ψ),    S_film = s^a_Macro · h(ψ) / φ^Macro
```
Same vdW physics as the micro reservoir, applied on the lower-specific-surface
macro walls. Persists as ψ→−∞; it sets the film residual `S_r`. The macro
specific surface `s^a_Macro` is *distinct* from the micro specific surface —
the film term is a **partition** of adsorption onto the macro walls, not a
duplicate of the micro vdW term (this is the answer to the double-counting
worry; the real cost is architectural, see "Tradeoffs").

**Smoothness / monotonicity / well-posedness** (beamer §"Smoothness…"):
- Smoothness is automatic — *no switch point*. Both branches are smooth in ψ
  over the whole range; their sum has no kink (per-pore drainage events are
  integrated out by the pore-size distribution). In the pure coexistence
  reading there is **no C¹ blend / air-entry junction to engineer**.
- Monotone: `dS_Macro/dψ < 0` on both terms ⇒ summed WRC monotone (needed for
  a stable global solve and for `k_r(S_Macro)`).
- Analytic derivatives: the MPL property needs
  `∂S_Macro/∂p_cap = ∂(S_cap + S_film)/∂p_cap` in closed form. **Mirror the
  existing micro vdW `domega_l_d*` pattern** — the film term reuses it
  directly.
- Residual: `S_r` is now *physically* the film-held water, not a numerical
  floor — removes the "ride the VG asymptote" artefact.
- Well-posedness: the macro balance's primary variable *is* ψ, so we only ever
  evaluate *forward* ψ → (r, h) → S_Macro. Tuller parameterizes S_Macro(ψ),
  never ψ(S_Macro) — the inversion that would be undefined at S_Macro→0 never
  occurs. The IC can have ψ = total suction (~100 MPa) with S_Macro
  small-but-defined, carried by corner + film water.

---

## Truesdell–Noll standing (why Option B is the principled closure)

(beamer §"Through Truesdell–Noll")

| Principle | Behaviour | Verdict |
|---|---|---|
| 1. Determinism | `S_Macro(ψ)` algebraic in current ψ; no history (hysteresis not modelled) | Satisfied |
| 2. Local action | Pointwise; no gradients of ψ or S in the closure | Satisfied |
| 3. Frame-indifference | ψ, S scalars; ψ_ae, A_H, s^a material constants | Satisfied (trivially) |
| 4. Material symmetry | Isotropic pore-size / surface statistics; anisotropy out of scope | Satisfied under isotropy |
| **5. Equipresence** | **Improved.** Capillary and film branches share the *same* driving variable ψ. Pure VG appends no film branch and uses the capillary form past its validity. | Satisfied *non-trivially*; strictly better than VG |
| 6. Dimensional invariance | ψ [J/kg], h [m], s^a [m²/kg], A_H [J] consistent across both branches | Satisfied |
| **7. Thermodynamic admissibility** | **Improved.** Both branches are stationarity conditions of one pore-water free energy `F = σ A^LV(r) + ∫₀^h Π dh' A^SL − μ m^L`; the macro `p^disj` *descends* from F rather than being postulated | Satisfied *by construction* |

**Bottom line.** The Tuller split keeps 1–4, 6 and *upgrades* 5 and 7 from
"postulated VG fit" to "equipresent, free-energy-derived" — the same upgrade
proposed for the micro closure. VG/Brooks–Corey pass 1–4, 6 but are postulated
fits: no shared-variable film branch (fail 5 non-trivially) and no free-energy
derivation (fail 7 non-trivially). Tuller is the only common WRC framework
that passes 5 and 7 *by construction*, so the whole double-structure model
rests on one free energy.

**Caveat to state in the paper:** Tuller assumes instantaneous
film↔capillary equilibrium at each ψ. Fast transients need a kinetic exchange
term — which already exists as `ρ̂ = −ζ(ψ_Macro − ψ_Micro)`.

---

## Tradeoffs — the honest cost of Option B

(beamer §"Two ways…", Option B block + §"Is it big surgery?")

Option B is physically richer but **contradicts the paper's strict scale
separation** (macro = capillarity only, micro = vdW only). That premise
breaking ripples into:

1. **`k_r`** — relative permeability derived from the macro WRC changes when
   the film branch is added.
2. **A second `p^disj`** in the effective stress — the macro film introduces a
   disjoining-pressure contribution on the macro side, in addition to the
   micro one already carried in Δσ^sw. (On the *mechanical* side the existing
   model already neutralises the micro term via `BishopsSaturationCutoff(χ=0)`;
   on the *hydraulic* side a `p_L` macro/micro split is still open — see the
   MS33 run report. The macro film must be wired so the two `p^disj`
   contributions partition cleanly, not overlap.)
3. **The two-scale premise itself** — once adsorption lives on both scales the
   "macro = capillarity only" story in the paper must be rewritten, including
   the "no cross-contribution" claim and the figure caption.

The hard *engineering* cost is **physics parity across two hand-coded tangent
implementations** (MPL and MFront), not the algebra — but that cost is
deferred to the MFront follow-on, not incurred on this native branch.

---

## Implementation scope — native C++, in-situ, on this branch

(beamer §"Implementation scope"; scope here trimmed to the native MPL work
this branch owns)

**Already present on base** (scaffold to extend, *not* build from zero):

- [SaturationTuller.{h,cpp}](MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationTuller.h)
  — single-branch (`1 − exp(−C/p²)`, the upscaled *corner* term), uncapped
  (no `cavitation_pressure`). Constructor currently takes
  `(name, residual_liquid_saturation, maximum_liquid_saturation,
  area_factor_tuller, pore_area_shapefactor_tuller, characteristic_pore_size,
  surface_tension, pressure_tolerance)` with
  `coefficient_ = 4·pore_area_shapefactor·σ² / (area_factor·char_pore²)`.
  `value()`/`dValue()`/`d2Value()` are implemented for the capillary pressure
  only.
- [PotentialExchange.h](ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h)
  — `computeVanDerWaalsMicroPotential` (micro, the mirror target for the film
  term) **and** `computeYoungLaplaceMacroPotential` (capillary-branch
  scaffold).
- Registered in
  [CreateProperty.cpp:189](MaterialLib/MPL/CreateProperty.cpp) (also accepts
  the `TullerRetention` alias); CMake auto-globs (no CMake edit).

**Edits this branch lands** — in dependency order, each commit leaving the
tree compiling and previously-passing tests still passing (Strict Rule 5):

1. **Extend `SaturationTuller` in place with the film term.**
   - Take the film-closure parameters (`A_H`, `s^a_Macro`, `φ^Macro`, and
     whatever the chosen blend needs) as **additional constructor arguments**.
   - Add the film contribution `S_film = s^a_Macro·h(ψ)/φ^Macro` (with
     `h(ψ)` from `ψ_film = −A_H/(6π ρ^L h³)`) into `value`, `dValue`, and
     `d2Value`, mirroring the micro vdW `domega_l_d*` derivative pattern.
   - **Pure-Tuller numerical behaviour MUST be recovered exactly** when the
     film parameters are at their "off" sentinel. This is what makes the
     change Strict-Rule-5 safe — existing PRJs and the unit test continue to
     pass unchanged.
2. **Extend `CreateSaturationTuller.{h,cpp}` in place.**
   - Parse the new *optional* XML keys (default each to the "off" sentinel)
     and forward them to the constructor.
3. **Extend [TestMPLSaturationTuller.cpp](Tests/MaterialLib/TestMPLSaturationTuller.cpp)
   in place.**
   - Cases covering the film-augmented closure: monotonicity, ψ→−∞ limit
     (film residual), smoothness of `dValue`, and **exact recovery of the
     pure-Tuller numbers when the film parameters are off**.

**Beamer branch-name inconsistency — read this.** The beamer's
implementation-scope frames ("Where the surgery lands", "Is it big surgery?",
"Work estimate", "Decision points") name the surgery branch as
`dsm_native_tuller_macro_wrc` and list MFront work inline. That is
**inconsistent** with the beamer's own decision frame, which parks Option B on
`dsm_native_tuller_macro_film`. Per Vinay's confirmed assignment
(2026-05-29): **Option B lives on `dsm_native_tuller_macro_film`;
`dsm_native_tuller_macro_wrc` stays clean** (Option A is a separate effort).
So when reading those beamer frames, mentally substitute `_macro_film` for the
native MPL items and treat all MFront items as the follow-on below.

---

## MFront follow-on (NOT on this native branch; gates the merge)

The corner+film split must be ported to MFront and kept at bit-parity with the
native MPL property **before Option B can merge**. This work is on the
`dsm_mfront` lineage, where the bridge sources actually live:

- `MaterialLib/SolidModels/MFront/RichardsMechanicsDSMMicroMacroBridge.mfront`
  (saturation lambda; analytic `∂S_L/∂Δp`) **plus** the `_MCC` variant —
  3 code paths total (MPL, MFront bridge, MFront `_MCC`).
- Kept at parity by the
  `mfront_parity_1element_..._tuller_{native,bridge}` ctest pair (see
  [reference_dsm_parity_script.md](../../.claude/projects/-Users-vinaykumar-git-ogs/memory/reference_dsm_parity_script.md)
  — `python3 scripts/run_dsm_parity.py`). Reuse the existing
  `mfront_parity_1element_dsm_micromacro_mcc_tuller_{native,bridge}.prj` pair
  on `dsm_mfront` (modify in place to exercise the film parameters; do not add
  a new pair).
- This is the **highest-risk** part (hand-coded analytic Jacobian); the
  beamer flags it as the main cost driver.

---

## Decision points (from you, before the agent starts — beamer §"Decision points")

1. **Air-entry ψ_ae** — from the largest drainable macro pore `r_max`: a
   measured macro-pore size, or derived from `p_b = 27` MPa?
   *Note:* in the pure coexistence reading (corner shrinks continuously to 0,
   "no air-entry threshold"), ψ_ae and the blend are **moot** — they only
   matter if a corner cutoff / junction is reintroduced. Resolve which reading
   is intended.
2. **Film `h` as state** — track macro film thickness as an output/state, or
   compute pointwise and discard? (Decides plumbing cost; the beamer's
   "new process state variable?" is *likely no* but flags "verify if film h
   must be tracked".)
3. **Blend** — hard switch at ψ_ae (simpler, kink risk) vs. C¹ smoothing
   (Newton-safe, more algebra). Only relevant if a junction exists (see 1).
4. **Calibration target** — keep the Villar swelling-pressure target, or also
   fit the high-suction retention tail if data exist?

Record the choices in the Status log when made.

---

## Calibration & verification

- Re-anchor any swelling-pressure calibration so the macro film term does not
  double-count what the micro `vdw_augmentation_prefactor K` already supplies.
  Because `s^a_Macro ≠ s^a_Micro`, the film is a partition, not a duplicate —
  but the *calibrated* K must be re-fit so the total Δσ^sw still matches
  targets. Targets: Villar (Dixon EMDD = ρ_d), per the active group agreement
  (memory `feedback_emdd_dry_density_convention`).
- Run the four-pack MS33 LE/MCC suites on the film-augmented build and compare
  to the hierarchical baseline numbers archived under
  `tex/cc2024/VK_SB_EURAD_DSM/`.

## Documentation (Strict Rules 4 & 6)

On every code-level step above: update this AGENTS.md (Status log) and
`ProcessLib/RichardsMechanics/DSM/AGENTS.md`. Mark each item DONE with a date
when it lands; never erase items — annotate.

---

## Status log

- **2026-05-29** — Branch goal recorded. No implementation commits yet.
  Reconstruction note: the design name "Tuller Option A" first surfaces in
  `stash@{0}` ("auto-stashed 2026-05-28 before Tuller Option A"); the detailed
  Option A vs. Option B distinction was supplied directly by Vinay on
  2026-05-29 from his beamer / paper §2.4 reconciliation notes.
- **2026-05-29** — Implementation discipline pinned (Vinay): **in-situ edits
  only, no new files.** `SaturationTuller`, `CreateSaturationTuller`,
  `TestMPLSaturationTuller.cpp` are the only native edit targets. New PRJ test
  data requires explicit clearance before adding. This supersedes the earlier
  "new class vs. extend" open question — the decision is: extend in place.
- **2026-05-29 — SUPERSEDED (see next entry).** An earlier reading of *only*
  the EURAD-DSM deck (not the beamer) concluded that Option B double-counts
  the micro adsorption and that the architecture demanded Option A, and
  recorded this as a "BLOCKING — architecture decision pending Vinay" open
  question with implementation on hold. **This finding was wrong on two
  counts** and is retained here only as the historical record (Strict Rules
  3 & 7): (a) the decision was *already made* — Option A is active on
  `dsm_native_hierarchical`, Option B is a deliberately parked alternative on
  this branch, so nothing is "blocking"; (b) the film term is a *partition* of
  adsorption onto the distinct macro specific surface `s^a_Macro`, not a
  duplicate of the micro vdW term — the real cost of Option B is architectural
  (scale-separation premise, `k_r`, a second `p^disj`), not double-counting.
- **2026-05-29** — Read the authoritative beamer
  `tex/cc2024/VK_B35_Pinion_May_2026/tuller_macro_wrc.tex` in full (per
  Vinay's pointer) and **fully rewrote this AGENTS.md to match it.** Captured:
  the IC problem Option B fixes, the corner+film closure (corner = existing
  `1 − exp(−C/ψ²)`; film = `ψ_film = −A_H/(6π ρ^L h³) ⇒ h(ψ)`,
  `S_film = s^a_Macro·h(ψ)/φ^Macro`, mirroring the micro vdW `domega_l_d*`),
  the Truesdell–Noll standing (keeps 1–4,6; upgrades 5 & 7), the honest
  tradeoffs (scale-separation premise, `k_r`, second `p^disj`), the
  native-C++-only in-situ scope, and the four decision points. Scoped MFront
  to a **follow-on** on the `dsm_mfront` lineage that gates the merge (the
  bridge sources do not exist on this native branch). Flagged the beamer's
  internal branch-name inconsistency (its implementation-scope frames say
  `dsm_native_tuller_macro_wrc`; Option B belongs on `_macro_film`).
  Verified against the tree: scaffold exists
  ([PotentialExchange.h](ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h),
  [CreateProperty.cpp:189](MaterialLib/MPL/CreateProperty.cpp)); no MFront DSM
  bridge / tuller parity ctest / tuller PRJ on this branch. **No code changes
  made** — implementation still not started, now un-blocked and awaiting the
  decision-point answers (ψ_ae reading, film-h-as-state, blend, calibration
  target).

---
