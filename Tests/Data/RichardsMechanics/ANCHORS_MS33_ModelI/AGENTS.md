# Agent instructions — ANCHORS MS33 Model I

Part of the **ANCHORS MS33** benchmark suite (Models I, III, IV, VII), DSM
native hierarchical model in `RichardsMechanics`.

**Authoritative MS33 instructions:** `../ANCHORS_MS33_ModelIV/AGENTS.md` — read
it for the suite-wide context, the negative-macro-porosity defect, and the
output-path provenance rule.

## This model

Model I (`ms33_modelI_dd1400.prj`, `dd1600`, `dd1800`): 1-element constant-volume
dry-density calibration against the Villar swelling-pressure curve.

**Villar target:** `p_sw = exp(6.77 * rho_d[g/cm³] − 9.07)` MPa
→ 1.504 MPa (dd1400), 5.824 MPa (dd1600), 22.556 MPa (dd1800).

### Current state (2026-05-20)

**PRJ change:** `micro_water_content_swelling_slope = 0` set in all three PRJ
files (Pi-path default). This means the swelling stress comes entirely from the
augmentation disjoining pressure `Pi = rho_lR * K * exp(−xi)`.

**Status: NEEDS VERIFICATION with current binary.**
All five Pi-path PRJ tags are set (`slope = 0`,
`accumulate_swelling_contributions = true`,
`use_micro_liquid_density_for_pi = true`,
`use_micro_liquid_density_for_micro_pressure = true`). The historical
early-return bug that silenced `sigma_sw` when `slope = 0` is resolved by the
`accumulate_swelling_contributions` flag. A rebuild and verify run is needed to
confirm MAE < 1 % with this configuration.

### What needs to happen (in order)

1. **Rebuild OGS:**
   ```bash
   cd /Users/vinaykumar/git/build/release-omp-mfront && ninja RichardsMechanics ogs
   ```

2. **Verify K values:** K was calibrated when the Pi-path was active (old
   `slope = 0.1` caused the augmentation block to run and early-return before
   slope). With `accumulate_swelling_contributions = true` and `slope = 0`, the
   same Pi block runs — K values should hold.
   ```bash
   cd /Users/vinaykumar/git/ogs/Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIV
   python ms33_calibrate_K.py --verify
   # If MAE ≥ 1%: python ms33_calibrate_K.py  (bisection)
   ```
   See `../ANCHORS_MS33_ModelIV/AGENTS.md` Action 5.

3. **Recheck transport_porosity ≥ 0** node-wise in output VTUs.

### K values in the PRJ files

| PRJ file | K [J/kg] | Villar target |
|---|---:|---|
| `ms33_modelI_dd1400.prj` | 4981.81 | 1.504 MPa |
| `ms33_modelI_dd1600.prj` | 23423.8 | 5.824 MPa |
| `ms33_modelI_dd1800.prj` | 105429.7 | 22.556 MPa |

All use `hamaker_constant = 5.1e-21 J` (literature), `specific_surface = 523`
(code units, gives h₀ ≈ 1.4 nm), `vdw_augmentation_decay_length = 1e-6`.
See `../ANCHORS_MS33_ModelIV/AGENTS.md` Action 5 for the literature justification
and the λ unit-convention note.

### Gate to restore PASS status

Run `python ms33_calibrate_K.py --verify` after the Action 3 fix. All three
densities must show < 1 % error. Then regenerate `ms33_modelI_ps_path.pdf` and
`ms33_modelI_ks_path.pdf` (postprocess script: `ms33_postprocess_modelI.py`).
Update `../ANCHORS_MS33_RUN_SUMMARY_2026-05-20_anchors_dd.md` with the new run date.

## Binding rule — record output and result paths

Every figure / number produced here must carry its full provenance chain
(result artefact, postprocess script, `.prj` / `.pvd` / `.vtu` / `.log`) in the
run summary and the postprocess docstring. The status deck must show a source
line on every slide that plots or pulls this data. See
`../ANCHORS_MS33_ModelIV/AGENTS.md` and
`materialmodels/src/TPM/VK_SB_EURAD_DSM/AGENTS.md`.

Model I artefacts: `ms33_modelI_dd{1400,1600,1800}.prj` -> `*.pvd`;
`ms33_postprocess_modelI.py` -> `ms33_modelI_ps_path.pdf`,
`ms33_modelI_ks_path.pdf`, `ms33_modelI_villar_benchmark_calibration.csv/.json`;
logs `anchors_dd*_full.log`.
