# Agent instructions — ANCHORS MS33 Model VII

Part of the **ANCHORS MS33** benchmark suite (Models I, III, IV, VII), DSM
native hierarchical model in `RichardsMechanics`.

**Authoritative MS33 instructions:** `../ANCHORS_MS33_ModelIV/AGENTS.md` — read
it for the suite-wide context, the negative-macro-porosity defect, and the
output-path provenance rule.

## This model — OPEN DEFECT (void-ratio mismatch)

Model VII (`ms33_modelVII_freeswelling.prj`): free-swelling void-ratio vs
axial-stress path.

- **Fixed:** stress-path interpretation. The benchmark compares **total** axial
  stress `sigma_total = sigma_effective + chi * p_L` with `chi = S`. Postprocess
  was corrected (macro-porosity priority `porosity -> dry_density_solid ->
  transport_porosity`) and the loading path now matches targets within
  ~0.01 MPa at all required steps (200, 205, ..., 240 d).
- **Open:** simulated void ratio during loading/unloading is `e ~ 2.81..2.96`,
  far above the benchmark reference band `0.4..1.2`. This is genuine model
  behaviour, **not** a postprocessing artefact. It must be resolved before
  Model VII can move from PARTIAL to PASS.
- Note `transport_porosity` is identically zero here after hydration (macropores
  fully closed but clamped at 0) — the same DSM mechanism that goes negative in
  Model IV. After the Model IV macroporosity fix lands, recheck Model VII.

## PRJ change (2026-05-20) — slope = 0 propagated

`micro_water_content_swelling_slope = 0` was set in the PRJ (Pi-path default).
K = 23423.8 J/kg (reference calibration at ρ_d = 1600 kg/m³, same as Models
III and IV). This change does not alter the on-hold status below.

**With current C++** (early-return bug in `computeReferenceMicroPorositySwellingStressIncrement`),
slope = 0 means sigma_sw = 0. Do not rerun until Action 3 is fixed and Model I
dd1600 K is verified. After Action 3: if K(dd1600) is recalibrated in Model I,
propagate the new K here too.

## Do not touch until Model IV fix is merged

Model VII is on hold. Do not modify the PRJ, postprocess script, or reference
VTUs until the Model IV macroporosity fix is committed and verified, and the
void-ratio root cause is understood.

Once the Model IV fix is confirmed and Model I dd1600 K is verified (Action 3 +
Action 5 from `../ANCHORS_MS33_ModelIV/AGENTS.md`): rerun Model VII with the
updated K, recheck the void-ratio trajectory, and reassess the PARTIAL status.

## Binding rule — record output and result paths

Every figure / number produced here must carry its full provenance chain
(result artefact, postprocess script, `.prj` / `.pvd` / `.vtu` / `.log`) in the
run summary and the postprocess docstring. The status deck must show a source
line on every slide that plots or pulls this data. See
`../ANCHORS_MS33_ModelIV/AGENTS.md` and
`materialmodels/src/TPM/VK_SB_EURAD_DSM/AGENTS.md`.

Model VII artefacts: `ms33_modelVII_freeswelling.prj` ->
`ms33_modelVII_freeswelling.pvd`; `ms33_postprocess_modelVII.py` ->
`ms33_modelVII_void_ratio_stress.pdf`; logs `anchors_modelVII_full.log`,
`ms33_modelVII_freeswelling.log`.
