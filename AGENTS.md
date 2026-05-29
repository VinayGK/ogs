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

> **ACTIVE DESIGN (2026-05-29, latest) — read this first.** The closure on this
> branch is now a **sharp cavitation handover** at the cavitation pressure
> `p_cav` (Frydman & Baker 2009), **not** the additive corner+film
> *coexistence sum* described in the older sections below. The two macro
> branches are **sequential in `p_cap`**, joined at `p_cav` — *not* summed.
> The fully agent-executable spec is the new section
> **"ACTIVE closure — sharp cavitation handover"** immediately below; the
> beamer's matching slides are its `\section{The active plan: sharp cavitation
> handover}` (the older additive-sum slides there now carry a "THEORETICAL"
> watermark). The additive-sum sections further down (and the shipped additive
> `value()` placeholder, see the final Status-log entry) are **superseded** and
> retained only as the historical/parked reading per Strict Rules 3 & 7.

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

## ACTIVE closure — sharp cavitation handover (2026-05-29 latest)

(beamer §"The active plan: sharp cavitation handover"; Frydman & Baker 2009,
after Tuller 1999 / Or 1999. **This SUPERSEDES the additive coexistence sum**
in the next section — keep that next section only as the parked theoretical
reading, Strict Rules 3 & 7.)

Macro liquid saturation is **piecewise** in the capillary pressure `p_cap`, the
two branches **sequential** (not summed), joined **sharply** at the cavitation
pressure `p_cav`:

```
                ┌ S_Lmax · (1 − exp(−C_T / p_cap²)),   0 ≤ p_cap ≤ p_cav   (capillary corner)
S_Macro(p_cap) =│
                └ C_film · p_cap^(−1/3),               p_cap > p_cav        (adsorbed film = residual)

C_T   = 4 F_γ σ² / (A_n L²)              [Pa²]      (the existing single-branch Tuller coefficient_)
C_film = (s^a_Macro / φ^Macro) · (A_H/6π)^(1/3)     [Pa^(1/3)]   (the existing computeFilmCoefficient)
```

- **Capillary regime** (`p_cap < p_cav`): the Tuller *corner* term carries the
  macro water — exactly the term already in `SaturationTuller`.
- **Cavitation** (`p_cap = p_cav`): macropore menisci cavitate, corner water
  drains.
- **Adsorptive regime** (`p_cap > p_cav`): only the vdW *film*
  `C_film·p_cap^(−1/3)` survives — it **is** the macro residual (so the constant
  floor `S_L_res` is dropped when the film is on).

**Junction (the downward jump).** At `p_cav`,
`ΔS_Macro = S_Lmax(1 − exp(−C_T/p_cav²)) − C_film·p_cav^(−1/3) ≥ 0`. It must be
*downward* (corner-just-wet ≥ film-just-dry) — a **parameter-sanity check, not a
clamp**. Vinay's physical reading: "smooth in saturation" (a vertical plateau in
`p_cap(S_Macro)`), "jumpy in `p_cap`" (a downward step in `S_Macro(p_cap)`).

**Derivatives (each branch smooth; reuse the existing code).**
```
p_cap < p_cav:  dS/dp_cap = −S_Lmax · (2 C_T / p_cap³) · exp(−C_T/p_cap²)   < 0
p_cap > p_cav:  dS/dp_cap = −(1/3) · C_film · p_cap^(−4/3)                  < 0
```
Within a Newton evaluation `p_cav` is fixed, so each branch has a well-defined
tangent; only the crossing is non-smooth.

### Agent-executable rewrite spec (native C++, in-situ, no new files)

Target files (extend in place — do **not** add a sibling class or new gtest):

- [SaturationTuller.{h,cpp}](MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationTuller.h)
- [CreateSaturationTuller.cpp](MaterialLib/MPL/Properties/CapillaryPressureSaturation/CreateSaturationTuller.cpp)
- [TestMPLSaturationTuller.cpp](Tests/MaterialLib/TestMPLSaturationTuller.cpp)

1. **Replace the shipped additive placeholder.** The current `value()` does
   `S = S_L_res_ + (S_Lmax−S_L_res_)(1−e^{−C/p²}); if (film_active_) S += C_film·cbrt(1/p); return std::clamp(...)`.
   That additive-plus-clamp form is the **known-wrong placeholder** and is
   removed. (`computeFilmCoefficient` and `coefficient_` themselves are correct
   and stay.)
2. **New `value()` — gated piecewise:**
   - `!film_active_` → `S_L_res_ + (S_Lmax − S_L_res_)(1 − e^{−C_T/p²})`
     — the unchanged single-branch Tuller, **byte-for-byte off-recovery**.
   - `film_active_` → piecewise on `p_cav`:
     `p ≤ p_cav → S_Lmax(1 − e^{−C_T/p²})`; `p > p_cav → C_film·p^(−1/3)`.
     **No `std::clamp`; no `S_L_res_` floor.**
3. **`dValue()` / `d2Value()`** — same gate, same split; reuse the per-branch
   corner and film derivative formulas above (the film `dS = −(C_film/3)p^(−4/3)`,
   `d²S = (4/9)C_film·p^(−7/3)` are already implemented — just move them under
   the `p > p_cav` branch). Carry the unit comment on each changed line (§4.2).
4. **New constructor parameter `cavitation_pressure` (`p_cav`)** on
   `SaturationTuller` + `CreateSaturationTuller`. **Default `+∞`** (inactive) ⇒
   the cut never triggers ⇒ uncut curve recovered. Parse it as an *optional* XML
   key. Use the **same parameter** as Option A's `cavitation_pressure` on
   `dsm_native_hierarchical` — one shared concept, not a second knob.

**One parameter, three behaviours (the unification — verify by test, do not
assert numerically here):**

| film (`s^a_Macro`) | `p_cav` | behaviour |
|---|---|---|
| off | `+∞` | pure single-branch Tuller (byte-for-byte) |
| off | finite | **Option A**: capped capillary; macro drains to micro via the exchange |
| on  | finite | **Option B** (this branch): capped capillary + macro film residual |

### Test plan (structure only; expected values `TODO(user)`, §3)

In [TestMPLSaturationTuller.cpp](Tests/MaterialLib/TestMPLSaturationTuller.cpp),
property-based, no expected-value literals:

1. **Off-recovery** — `s^a_Macro = 0` ⇒ value/dValue/d2Value equal the
   single-branch Tuller across a `p_cap` sweep. *Anchor: prior-approved
   regression baseline (the in-file pure-Tuller test).*
2. **Regime selection** — `p_cap < p_cav` uses the corner; `p_cap > p_cav` uses
   the film. *Anchor: analytical limit.*
3. **Monotonicity** — `dValue < 0` within each branch. *Anchor: analytical.*
4. **Junction** — `S_Macro(p_cav⁻) ≥ S_Macro(p_cav⁺)` (downward jump).
   *Anchor: analytical/physical.*
5. **FD vs analytic** derivative within each branch (same FD tolerances as the
   existing in-file pure-Tuller test). *Anchor: analytical.*

### Open §9 sub-decision (Vinay, formulation — NOT decided by the agent)

The junction treatment is the one remaining formulation choice. Options:
**(a)** genuine jump — sharp, the default, accept the `S_Macro` discontinuity at
`p_cav` (Vinay's "jumpy in `p_cap`" favours this); **(b)** continuity by picking
`p_cav` at the crossover `S_corner = S_film` (a slope kink only, no jump);
**(c)** a narrow regularisation band around `p_cav`. **Confirm with Vinay before
choosing (b) or (c).** Implement (a) unless told otherwise.

### Parameter sourcing & merge gates (this design)

- `p_cav` ≈ `1.4×10⁸` Pa — **Or & Tuller (2002)**, already the DSM cavitation
  pressure (memory `project_dsm_mcc_bishop_cutoff`). Confirm the same value
  governs the macro film cut before writing a literal in any PRJ.
- `A_H` ≈ `2.2×10⁻²⁰` J — **literature anchor** (Israelachvili & Adams 1978 /
  Mitchell & Soga 2005), **NOT a knob** (CLAUDE.md §1.1 incident).
- `s^a_Macro` — still needs a **FEBEX**-sourced value (CLAUDE.md §12.1 family #6)
  **before any `.prj` enables the film** → STOP and ASK at that point (§1.1).
- `φ^Macro` — from the PRJ.
- **MFront port + native↔bridge parity** on the `dsm_mfront` lineage gates the
  merge (see "MFront follow-on" below).
- **Recalibrate K** (`vdw_augmentation_prefactor`) against Villar targets
  (Dixon EMDD ≡ ρ_d) after the closure changes — *predicted* to shift, not yet
  verified (§5).

---

## Option B closure (SUPERSEDED — additive coexistence sum, parked theoretical)

> **SUPERSEDED 2026-05-29** by the sharp cavitation handover above. This section
> describes the *additive coexistence sum* `S_Macro = S_cap + S_film` (both
> branches present at every potential, summed). It is retained as the historical
> / parked theoretical reading (Strict Rules 3 & 7); the corresponding beamer
> slides now carry a "THEORETICAL" watermark. **Do not implement the additive
> sum** — implement the piecewise cavitation handover above.

(beamer §"Implementation discussion --- volume referencing of the sum",
now watermarked THEORETICAL)

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

## Decision points (beamer §"Decision points") — RESOLVED 2026-05-29 (Vinay)

> **Refined by the active design.** Decisions (1) "no air-entry, no junction"
> and (3) "no blend" below were made under the *additive coexistence* reading
> (both branches summed at every potential). The **active** sharp cavitation
> handover **reintroduces a single junction at `p_cav`** (the branches are now
> sequential, not summed). Its treatment is the one open **§9 formulation
> sub-decision** — see "ACTIVE closure → Open §9 sub-decision" above:
> (a) genuine jump [default], (b) continuity at crossover, (c) regularisation.
> Decisions (2) film-`h`-not-a-state and (4) Villar-only calibration are
> unchanged.

Original questions retained for the record (Strict Rule 6); each is annotated
with Vinay's ruling on 2026-05-29.

1. **Air-entry ψ_ae** — *from the largest drainable macro pore `r_max`: a
   measured macro-pore size, or derived from `p_b = 27` MPa?*
   **RESOLVED: no air-entry, no junction.** Pure coexistence reading — the
   corner shrinks continuously to 0 and the film persists. There is no ψ_ae to
   pick. (This is the beamer's "Resolution" frame, not its decision-points
   frame.)
2. **Film `h` as state** — *track macro film thickness as an output/state, or
   compute pointwise and discard?*
   **RESOLVED: `h` is NOT a state and NOT the target.** The closure is
   parameterized through **water content** — the measurable quantity the paper
   uses (`ψ ↔ ω_l`), exactly as the micro vdW potential already is. Film
   thickness `h` is hard to measure and is, at most, an internal pointwise
   intermediate; it is never tracked, never plumbed as a primary/secondary
   variable, and never a calibration target. **Be careful here:** the
   implementation must express `S_film` via the adsorption/water-content
   expression, not by carrying `h`.
3. **Blend** — *hard switch at ψ_ae vs. C¹ smoothing?*
   **RESOLVED: no blend.** Follows directly from decision 1 (no junction); both
   branches are smooth in ψ and their sum has no kink.
4. **Calibration target** — *Villar swelling pressure only, or also fit a
   high-suction retention tail?*
   **RESOLVED: Villar swelling pressure only.** No high-suction retention-tail
   data exist, and none are needed — the water-content/adsorption expression
   *produces* the high-suction behaviour automatically. Do not fit a separate
   tail.

**Parameter sourcing (resolved 2026-05-29).** `φ^Macro` is taken from the PRJ
files (point the implementer at the specific PRJ + key). Remaining FEBEX-bentonite
material parameters come from the newly-added **FEBEX** source family (CLAUDE.md
§12.1, family #6: Villar 2002 ENRESA PT; Lloret & Villar 2007, Phys. Chem.
Earth 32, 701–715) — "most parameters are there" (Vinay). `A_H` stays the
literature anchor (2.2e-20 J, Israelachvili & Adams 1978), not a knob; confirm
the same `A_H` is intended for the macro walls before writing it. `s^a_Macro`
still needs a cited value (FEBEX or another §12.1 source) before any literal is
written — STOP and ASK at that point (§1.1).

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
- **2026-05-29** — All four decision points RESOLVED by Vinay (see the
  Decision points section): (1) no air-entry, no junction — pure coexistence;
  (2) film `h` is not a state and not the target — closure parameterized via
  **water content** (`ψ ↔ ω_l`, the measurable, as in the micro vdW closure),
  `h` only an internal intermediate; (3) no blend (follows from no junction);
  (4) calibrate to Villar swelling pressure only, high-suction behaviour is
  emergent from the adsorption expression (no tail fit). Parameter sourcing:
  `φ^Macro` from the PRJ; FEBEX added as CLAUDE.md §12.1 source family #6
  (Villar 2002 ENRESA PT; Lloret & Villar 2007, PCE 32, 701–715) at Vinay's
  instruction; `A_H` stays the literature anchor; `s^a_Macro` still needs a
  cited value before any literal (§1.1 STOP/ASK). A standalone evaluation
  Q&A was written alongside the beamer at
  `tex/cc2024/VK_B35_Pinion_May_2026/tuller_macro_wrc_QA.md`. Still no code
  changes — the native MPL surgery can now begin once `s^a_Macro` is sourced.
- **2026-05-29 — DONE: native MPL film surgery implemented and unit-verified.**
  **(NOW SUPERSEDED — see the next entry: the shipped additive `value()` is the
  known-wrong placeholder; the active design is the sharp cavitation handover,
  not this additive sum + clamp.)**
  Extended in place (no new files):
  [SaturationTuller.h](MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationTuller.h),
  [SaturationTuller.cpp](MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationTuller.cpp),
  [CreateSaturationTuller.cpp](MaterialLib/MPL/Properties/CapillaryPressureSaturation/CreateSaturationTuller.cpp),
  [TestMPLSaturationTuller.cpp](Tests/MaterialLib/TestMPLSaturationTuller.cpp).
  Closure: from the vdW disjoining-pressure balance `p_cap = A_H/(6π h³)`
  (ρ^L cancels in pressure form; same `6π` as the micro vdW in
  PotentialExchange.h) → `h(p_cap) = (A_H/(6π p_cap))^(1/3)`,
  `S_film = (a_v/φ_M)·h = C_film·p_cap^(−1/3)` [dimensionless],
  `C_film = (a_v/φ_M)·(A_H/(6π))^(1/3)` [Pa^(1/3)]. Added into
  `value`/`dValue`/`d2Value` (`dS_film = −(C_film/3)p^(−4/3)`,
  `d²S_film = (4/9)C_film·p^(−7/3)`), with unit comments (§4.2).
  Three new constructor args `macro_specific_surface` (a_v),
  `hamaker_constant` (A_H), `macro_porosity` (φ_M), all defaulting to 0;
  `film_active_ = (a_v > 0)` gates the branch. **Off recovers pure Tuller
  byte-for-byte** (the `if (film_active_)` blocks are skipped), verified by a
  dedicated point-for-point regression test. `h` is an internal pointwise
  intermediate only — never stored, never a state/target; the closure is
  parameterized through saturation/water content as Vinay required. New XML
  keys are optional in CreateSaturationTuller (default off). Tests added
  (property-based, no expected-value literals): off-recovery, film-on
  monotonicity + analytic-vs-FD derivative consistency (same FD tolerances as
  the in-file pure-Tuller test), and a create-path parse check. Build (warm
  CPM cache, MFront OFF) + `testrunner` green: 5/5 SaturationTuller tests and
  all 83 MaterialPropertyLib tests pass. **Guardrails noted (§0.1):** §1.1 —
  film-on tests use `a_v`/`φ_M` as *synthetic fixtures* (no physical magnitude
  asserted; A_H cited to Israelachvili & Adams 1978); a FEBEX-sourced
  `s^a_Macro` is still required before any *PRJ* enables the film branch. §9
  (formulation, left to Vinay, NOT decided here): (a) whether the film should
  replace `S_L_res` (QA Q7: film "sets the residual"); (b) the near-saturation
  region where `S_cap + S_film` can hit the existing `S_L_max` clamp — the
  literal additive reading + existing clamp is what shipped; the film-on test
  samples the high-suction window where it is unclamped. **Remaining gates
  before merge:** MFront port + native↔bridge parity on the `dsm_mfront`
  lineage; a §12-compliant Tuller-film PRJ once `s^a_Macro` is sourced.
- **2026-05-29 — CLOSURE SUPERSEDED: additive sum → sharp cavitation handover.**
  Vinay's decision (recorded in memory `project_tuller_film_volume_referencing`
  and `project_tuller_macro_options`): the additive coexistence sum
  `S_Macro = S_cap + S_film` is **parked as theoretical**; the active closure is
  a **sharp cavitation handover** at `p_cav` (Frydman & Baker 2009) — the two
  macro branches are **sequential in `p_cap`**, joined at `p_cav`, **not summed**
  (capillary corner for `p_cap ≤ p_cav`; vdW film = residual for `p_cap > p_cav`).
  This unifies Options A and B through one shared `cavitation_pressure` param:
  film off + `p_cav = +∞` → pure Tuller; film off + finite `p_cav` → Option A;
  film on + finite `p_cav` → Option B. Documented as the **active design** in
  this AGENTS.md (new section "ACTIVE closure — sharp cavitation handover",
  with the agent-executable rewrite spec for `SaturationTuller` /
  `CreateSaturationTuller` / `TestMPLSaturationTuller`) and in the beamer
  `tex/cc2024/VK_B35_Pinion_May_2026/tuller_macro_wrc.tex` (new
  `\section{The active plan: sharp cavitation handover}`, 6 slides; the older
  additive-sum slides watermarked "THEORETICAL", kept per §6). The shipped
  additive `value()` (previous entry) is therefore the **known-wrong
  placeholder** to be replaced by the gated piecewise form. **No code changed in
  this step — documentation only.** Open **§9 formulation sub-decision (Vinay,
  not the agent):** junction treatment — (a) genuine jump [default], (b)
  continuity at the crossover (slope kink), or (c) a narrow regularisation band;
  implement (a) unless told otherwise. Sourcing unchanged: `p_cav` from Or &
  Tuller (2002), `A_H` literature anchor, `s^a_Macro` still needs a FEBEX value
  (§1.1 STOP/ASK) before any PRJ enables the film.

---
