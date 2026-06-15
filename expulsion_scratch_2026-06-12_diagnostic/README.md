# Expulsion scratch — saturated column, mechanical squeeze (2026-06-12) — DIAGNOSTIC

Scratch probe (NOT the §12.3 standard suite; suffix per §6.8 diagnostic
precedent of 2026-06-11). Template:
`Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/ms33_confined_expulsion_dd1400.prj`
(all BC/IC/material literals from that approved file, unchanged); variants
differ ONLY by appended mode tags. Branch `dsm_native_Pi_fofnlev` @ 4c7a5d03b2,
binary `~/git/build/Pi_fofnlev_20260611/bin/ogs`, 4 runs concurrent.

Sequence: suction 100→0 MPa over 20 d (saturate), then compressive top
displacement ramp to ~5% strain by 200 d (squeeze).

## Measured

- Wetting phase discriminates: operational/equilibrium (ungated +b·p_conf/ρ
  term active under confined swelling) imbibe less than off/exact
  (n_l at 20 d: 0.418 vs 0.463).
- Squeeze phase: water IS pressed out in all variants — n_l 0.485 → 0.4706
  (20→200 d), exiting through the suction-0 boundary.
- BUT n_l == phi exactly in ALL variants from ~20 d on (no
  macro_porosity_floor in this PRJ): the micro fills the entire pore space.
  The expulsion is VOLUME-CONSTRAINT-driven (cap shrinks with phi), identical
  across all four mu-routes to 6 decimals. The mu-physics stays alive only in
  p_L_m (off/exact ≈ 16.6 MPa vs operational/equilibrium ≈ 10.5–10.7 MPa at
  200 d) and σ_sw (±0.2 MPa) — it cannot move n_l at the ceiling.

## Conclusion (sign-only)

"Press water out of a saturated confined column" WORKS in every mode — but
this configuration is mode-blind for the expulsion question (same
ceiling-degeneracy as MS33 VII). A discriminating T-8 needs n_l OFF the cap
at squeeze start: partial-saturation hold, macro_porosity_floor > 0, or an
open/drained macro — parameter/BC choices that are Vinay's call (§1.1/§9).

## Run 2 (same day): macro_porosity_floor = 0.08 (value: Vinay, this scratch)

Same 4 variants + `<macro_porosity_floor>0.08</macro_porosity_floor>`; same
binary (Pi_fofnlev_20260611 — deliberately NOT rebuilt: the worktree carries
uncommitted live-K-Jacobian source edits from a parallel session).

n_l relative to the cap (phi-0.08)/0.92 at 200 d (5% squeeze, sigma_zz ~ -7.5 MPa):

| mode            | n_l @200d | n_l - cap | reading |
| :--             | --:       | --:       | :-- |
| off             | 0.44189   | +0.0174   | TRAPPED above the shrinking cap — frozen-Pi cannot expel |
| kin_exact       | 0.42795   | +0.0034   | expels, tracks the cap from just above (conservative) |
| kin_operational | 0.42058   | -0.0040   | expels BELOW the cap (mu-driven, ungated bolt-on) |
| equilibrium     | 0.42059   | -0.0040   | == operational to 5 decimals (same Derjaguin load term) |

Sign-only findings:
1. The floor un-pins n_l: the configuration now discriminates the mu-routes.
2. off (frozen-Pi): squeeze shrinks the cap faster than the model expels —
   water stays trapped in the micro (the no-expulsion pathology, visible).
3. All strain/load-coupled modes expel; exact is measurably gentler than the
   operational cut (0.428 vs 0.421) — consistent with the operational form's
   known over-aggressive non-integrable load term (|loop dW|/scale 0.93).
4. p_L_m -> ~0 MPa in all modes at 200 d (film equilibrated with boundary).
Magnitudes remain TODO(Vinay) (T-8 discipline); this is direction-of-effect only.
