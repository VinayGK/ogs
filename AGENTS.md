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

Concrete edits this branch needs to land — listed in dependency order so
each commit leaves the tree compiling and previously-passing tests still
passing (per repo Strict Rule 5).

1. **MPL: new macro saturation property with film term**
   - Either a new class `SaturationTullerWithFilm` next to
     [SaturationTuller.h](MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationTuller.h),
     or extend `SaturationTuller` with an optional film contribution
     gated by a new constructor argument. **Open question** — keep the
     existing single-responsibility class untouched (new class) vs.
     in-place extension. Default lean: new class, to avoid touching the
     unit test invariants in [TestMPLSaturationTuller.cpp](Tests/MaterialLib/TestMPLSaturationTuller.cpp).
   - Closed form: `S_L(ψ) = S_cap(ψ) + S_film(ψ)`, with `S_cap` the existing
     Tuller capillary corner and `S_film(ψ)` a Hamaker-disjoining-pressure-
     derived film term, both supported at all ψ. Exact film closure (Or &
     Tuller form vs. a parametric exponential matched to the micro vdW
     augmentation already on `n_l`) is an **open decision** — record the
     choice in this file when made.
   - Implement `value`, `dValue` (w.r.t. capillary pressure), `d2Value`,
     and the constructor argument validation, matching the existing
     property's style.

2. **2× MFront bridge: macro film parameters on both sides**
   - Native side: thread the film closure's parameters through the macro
     constitutive bridge so the native RichardsMechanics loop sees them.
   - MFront side: corresponding additions to the macro MFront bridge so
     `dsm_native ↔ dsm_mfront` strict parity holds (see
     [reference_dsm_parity_script.md](../../.claude/projects/-Users-vinaykumar-git-ogs/memory/reference_dsm_parity_script.md)
     — `python3 scripts/run_dsm_parity.py`).
   - Add a parity suite entry mirroring the existing
     `mfront_parity_1element_dsm_micromacro_mcc_tuller_{native,bridge}.prj`
     pair on `dsm_mfront`, with the new film-augmented property.

3. **Tests**
   - Extend or fork [TestMPLSaturationTuller.cpp](Tests/MaterialLib/TestMPLSaturationTuller.cpp)
     to cover the film-augmented closure: monotonicity, limits, smoothness
     of `dValue`, and that turning the film term off recovers the pure-
     Tuller numbers exactly.
   - Add an integration smoke test in `Tests/Data/RichardsMechanics/`
     exercising the new property end-to-end on a single element.

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

- New `SaturationTullerWithFilm` class vs. in-place extension of
  `SaturationTuller`? (See Surgery step 1.)
- Film closure form: Or-Tuller analytical vs. exponential matched to the
  micro vdW augmentation? Both are defensible; the second has the
  advantage of using parameters already calibrated in the
  [ANCHORS_MS33_ModelI](Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/)
  workflow.
- Does the macro film term replace the cavitation cap entirely on this
  branch, or coexist with a soft `p_cav` cap as a numerical safeguard?
  Default reading from the table above: replace entirely (no threshold).

## Status log

- **2026-05-29** — Branch goal recorded. No implementation commits yet.
  Reconstruction note: the design name "Tuller Option A" first surfaces
  in `stash@{0}` ("auto-stashed 2026-05-28 before Tuller Option A"); the
  detailed Option A vs. Option B distinction above was supplied directly
  by Vinay on 2026-05-29 from his beamer / paper §2.4 reconciliation
  notes and is canonical for this branch.

---
