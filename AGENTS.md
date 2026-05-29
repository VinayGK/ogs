# AGENTS.md — repository root

This file is read by Claude / Codex agents at the repository root.

The repository contains several long-lived feature lineages. Branch-specific
goals are recorded here at the top so any agent that lands on this tree knows
*what this branch is for* before editing.

For DSM physics and the hierarchical micro–macro framework see
[ProcessLib/RichardsMechanics/DSM/AGENTS.md](ProcessLib/RichardsMechanics/DSM/AGENTS.md).

---

## Branch: `dsm_native_tuller_macro_film` — Tuller coexistence (Option B)

**Scope.** Only this branch carries adsorption on the macro side. The
sibling branch `dsm_native_tuller_macro_wrc` does **not** carry this work;
it stays on the hierarchical baseline pending a separate Option A surgery
(cavitation cutoff on the capillary-only curve, see *Contrast* below).
The two are intentionally kept apart so the surgeries can be reviewed
independently. Intent is to merge Option B back into the hierarchical
mainline once it has converged.

**Created.** 2026-05-29 (Vinay). No implementation commits on this branch
yet — this document is the entry point.

### Option B (this branch)

Macro WRC = capillary corner **plus** adsorptive film, both active at all
suctions ψ. No drainage threshold. Above the cavitation pressure `p_cav`
the macro film carries the potential rather than draining to the micro pool.

This is the "Tuller coexistence" reading: pure Tuller has no air-entry
threshold, so the macro retention has support at all ψ via the adsorbed
film term, even when the capillary corner contribution is exhausted.

| What | Value |
|---|---|
| Macro WRC | capillary corner + adsorbed film, both at all ψ |
| Above `p_cav` | macro film carries the potential; no drainage threshold |
| Threshold | none (Tuller has no air entry) |
| Paper change required | rewrite "no cross-contribution" claim + Fig. caption |
| Surgery size | larger: new macro film branch in MPL + 2× MFront |

### Contrast — Option A (NOT this branch)

For the record (so future agents do not confuse the two):

| What | Value |
|---|---|
| Macro WRC | capillary only, capped at `p_cav` |
| Above `p_cav` | macro drains completely → water lives in micro, reached via ρ̂ = −ζ(ψ_M − ψ_m) |
| Threshold | yes — `p_cav`, the cavitation pressure (already in §2.4 of the paper) |
| Paper change required | none — implements what §2.4 already says |
| Surgery size | smaller: cap the existing curve + ensure exchange handles the drained regime |

Option A is paper-consistent (Frydman–Baker basis); Option B is the beamer
("Tuller coexistence") reading. The user reaffirmed on 2026-05-29 that this
branch (`_macro_film`) is Option B; Option A is a separate effort.

---

## Surgery scope (Option B)

**Hard constraint (Vinay, 2026-05-29): in-situ, no new files.** Extend
existing classes / existing test file / existing MFront sources in place.
Do not introduce a `SaturationTullerWithFilm` sibling class, a new
gtest file, or new MFront `.mfront` sources. Existing PRJ test data
should be re-used or modified in place; do not add new ones without
clearing it with Vinay first.

Concrete edits this branch needs to land — listed in dependency order so
each commit leaves the tree compiling and previously-passing tests still
passing (per repo Strict Rule 5).

1. **MPL: extend `SaturationTuller` in place with the film term**
   - Edit
     [SaturationTuller.h](MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationTuller.h)
     and
     [SaturationTuller.cpp](MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationTuller.cpp)
     in place to take the new film-closure parameters as additional
     constructor arguments and to add the film contribution into
     `value`, `dValue`, and `d2Value`. The pure-Tuller numerical
     behaviour MUST be recovered exactly when the film parameters are
     set to their "off" sentinel (this is what makes the change strict-
     rule-5 safe — existing PRJs and the unit test continue to pass).
   - Edit
     [CreateSaturationTuller.cpp](MaterialLib/MPL/Properties/CapillaryPressureSaturation/CreateSaturationTuller.cpp)
     and its header in place to parse the new optional XML keys (default
     to the "off" sentinel) and forward them to the constructor.
   - Closed form: `S_L(ψ) = S_cap(ψ) + S_film(ψ)`, with `S_cap` the
     existing Tuller capillary corner and `S_film(ψ)` a Hamaker-
     disjoining-pressure-derived film term, both supported at all ψ.
     Film closure form (Or-Tuller analytical vs. exponential matched
     to the micro vdW augmentation already on `n_l`) is an **open
     decision** — record the choice in the Status log when made.

2. **2× MFront bridge: extend in place**
   - Native side: thread the film closure parameters through the
     existing macro constitutive bridge files. No new bridge file.
   - MFront side: corresponding additions to the existing macro MFront
     source(s) so `dsm_native ↔ dsm_mfront` strict parity holds (see
     [reference_dsm_parity_script.md](../../.claude/projects/-Users-vinaykumar-git-ogs/memory/reference_dsm_parity_script.md)
     — `python3 scripts/run_dsm_parity.py`).
   - Reuse the existing
     `mfront_parity_1element_dsm_micromacro_mcc_tuller_{native,bridge}.prj`
     pair on `dsm_mfront` as the parity vehicle (modify in place to
     exercise the film parameters; do not add a new pair).

3. **Tests (extend in place — no new test files)**
   - Extend [TestMPLSaturationTuller.cpp](Tests/MaterialLib/TestMPLSaturationTuller.cpp)
     with additional cases covering the film-augmented closure:
     monotonicity, limits, smoothness of `dValue`, and exact recovery
     of pure-Tuller numbers when the film parameters are off.

4. **Calibration & verification**
   - Re-anchor any swelling-pressure calibration so the macro film term
     does not double-count what the micro `vdw_augmentation_prefactor K`
     is already supplying. Targets: Villar (Dixon EMDD = ρ_d), per the
     active group agreement (memory `feedback_emdd_dry_density_convention`).
   - Run the four-pack MS33 LE/MCC suites on the film-augmented build
     and compare to the hierarchical baseline numbers archived under
     `tex/cc2024/VK_SB_EURAD_DSM/`.

5. **Documentation**
   - On every code-level step above: update this AGENTS.md (Status log
     below) and `ProcessLib/RichardsMechanics/DSM/AGENTS.md` per repo
     Strict Rule 4. Mark each surgery item as DONE with a date when it
     lands (Strict Rule 6). Do not erase items — annotate.

## Open questions (decide before implementation)

- Film closure form: Or-Tuller analytical vs. exponential matched to the
  micro vdW augmentation? Both are defensible; the second has the
  advantage of using parameters already calibrated in the
  [ANCHORS_MS33_ModelI](Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/)
  workflow.
- **(BLOCKING — architecture decision pending Vinay.)** Does Option B
  even apply to *this* model, or does the paper's architecture already
  demand Option A? Reading of the EURAD-DSM sources on 2026-05-29:
  - The paper assumes a clean two-scale split — **macro = capillarity
    only** (Tuller 1999), **micro = van der Waals / adsorption only**
    (Frydman 2009, Baker 2009; `μ^Micro = A_Hamaker/(6πh³)`)
    (`tex/cc2024/VK_SB_EURAD_DSM/VK_SB_EURAD_DSM.tex:290–342`).
  - The cavitation cutoff `p_cav ≈ 140 MPa` (Or & Tuller 2002) plays a
    **physical drainage role**: "beyond `p_cav` the macro meniscus
    cavitates, the macro pore water drains into the micro-structure,
    and S_L freezes" (`VK_SB_EURAD_DSM_output.tex:786–795`). It is the
    macro→micro handoff, not a numerical guard.
  - Consequence: the "adsorption potential active at zero water content"
    motivating Option B is, in this model, **the micro continuum**.
    Adding a macro adsorbed film re-introduces on the macro side the
    adsorption already carried on the micro side — the same
    double-counting already corrected on the mechanical side
    (`BishopsSaturationCutoff(χ=0)`, report §"Bishop's effective
    stress") and still open on the hydraulic side (the `p_L`
    macro/micro split, report §"Hydraulic-side double-counting").
  - **Therefore:** as a bolt-on to the existing micro=vdW model, Option B
    double-counts and **Option A is the architecture-consistent choice**.
    Option B is coherent *only* as a deliberate re-architecture that
    moves/partitions adsorption macro↔micro so the two do not overlap.
    This is a model decision, not a coding detail — do not start the
    surgery until Vinay rules: (a) re-architect for Option B with an
    explicit macro/micro adsorption partition, or (b) convert this
    branch's goal to Option A (capillary + `p_cav` drain-to-micro).
- Film closure form (only if Option B / re-architecture is chosen):
  Or-Tuller analytical vs. exponential matched to the micro vdW
  augmentation? The second reuses parameters already calibrated in the
  [ANCHORS_MS33_ModelI](Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/)
  workflow.

Note: the cavitation cutoff is **not in the code** on these branches —
[SaturationTuller.cpp](MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationTuller.cpp)
is the 2-branch uncapped form (no `cavitation_pressure`), despite the
deck describing a 3-branch capped form as "implemented but not yet
wired." Whichever option is chosen, the cap/freeze branch would itself
be new in-situ code.

## Status log

- **2026-05-29** — Branch goal recorded. No implementation commits yet.
  Reconstruction note: the design name "Tuller Option A" first surfaces
  in `stash@{0}` ("auto-stashed 2026-05-28 before Tuller Option A"); the
  detailed Option A vs. Option B distinction above was supplied directly
  by Vinay on 2026-05-29 from his beamer / paper §2.4 reconciliation
  notes and is canonical for this branch.
- **2026-05-29** — Implementation discipline pinned (Vinay): **in-situ
  edits only, no new files.** `SaturationTuller`,
  `CreateSaturationTuller`, `TestMPLSaturationTuller.cpp`, and the
  existing macro MFront bridge sources are the only edit targets. New
  PRJ test data requires explicit clearance before adding. This
  supersedes the earlier "new class vs. extend" open question — the
  decision is: extend in place.
- **2026-05-29** — Read the EURAD-DSM paper/deck sources to settle the
  `p_cav` role (per Vinay's request). Finding: the paper's macro scale
  is capillarity-only and the micro scale carries all adsorption (vdW),
  and `p_cav` is the physical macro→micro drainage threshold — so a
  macro adsorbed film (Option B) double-counts the micro adsorption
  unless adsorption is re-partitioned across scales. Recorded as the
  BLOCKING open question above. **Implementation is on hold pending
  Vinay's architecture ruling (re-architect for Option B vs. convert
  this branch to Option A).** No code changes made.

---
