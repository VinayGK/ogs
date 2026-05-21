# Agent instructions — ANCHORS MS33 Model III

Part of the **ANCHORS MS33** benchmark suite (Models I, III, IV, VII), DSM
native hierarchical model in `RichardsMechanics`.

**Authoritative MS33 instructions:** `../ANCHORS_MS33_ModelIV/AGENTS.md` — read
it for the suite-wide context, the negative-macro-porosity defect, and the
output-path provenance rule.

## This model

Model III (`ms33_modelIII_gap2mm.prj`): gap-closure case, 50 mm × 70 mm
axisymmetric cylinder with 2 mm lateral gap. Reference dry density 1600 kg/m³.

Previous results: final gap aperture 1.933 mm (from 2.000 mm); mean effective
stress (bottom/centre/top) 3.726 / 7.340 / 7.545 MPa.

### Current state (2026-05-20)

**PRJ change:** `micro_water_content_swelling_slope = 0` in the clay medium
(Pi-path default). K = 23423.8 J/kg (reference calibration at ρ_d = 1600 kg/m³).

**Status: NEEDS VERIFICATION with current binary.**
All five Pi-path PRJ tags are set (`slope = 0`,
`accumulate_swelling_contributions = true`,
`use_micro_liquid_density_for_pi = true`,
`use_micro_liquid_density_for_micro_pressure = true`). The historical
early-return bug that silenced `sigma_sw` when `slope = 0` is resolved by the
`accumulate_swelling_contributions` flag. A rebuild and verify run is needed to
confirm outputs are consistent with the Pi-path configuration.

### What needs to happen (in order)

1. **Rebuild OGS:**
   ```bash
   cd /Users/vinaykumar/git/build/release-omp-mfront && ninja RichardsMechanics ogs
   ```

2. **Verify K(dd1600) from Model I first:**
   ```bash
   cd /Users/vinaykumar/git/ogs/Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIV
   python ms33_calibrate_K.py --verify
   # If MAE ≥ 1%: python ms33_calibrate_K.py  (bisection)
   ```
   Model III uses K = 23423.8 J/kg (same as Model I dd1600 calibration). If
   K(dd1600) is recalibrated in Model I, update this PRJ to match before rerunning.

3. **Rerun Model III:**
   ```bash
   OGS=/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs
   OUTDIR=/Users/vinaykumar/git/ogs/Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIII
   $OGS -o $OUTDIR -l warn $OUTDIR/ms33_modelIII_gap2mm.prj \
       > $OUTDIR/ms33_modelIII_run_$(date +%Y%m%d).log 2>&1
   ```

4. **Check outputs:** gap aperture and mean effective stress must be in the
   expected ranges (see previous results above). Check
   `transport_porosity >= 0` node-wise.

### Restore PASS status gate

After rerun: gap aperture plausible (0.9–2.0 mm closure range), mean stress
finite. Check `transport_porosity >= 0`. Regenerate PDFs via
`ms33_postprocess_modelIII.py`. No strict numerical gate (unlike Model I);
plausibility check against the previous values is sufficient.

Note: Model III also carries a second medium (gap zone, medium id=1) with near-zero
n_s ≈ 0.015. The Pi contribution from the gap zone is negligible (~2.6 % of clay
value) but the DSM exchange path is still active — monitor for numerical issues
in the gap zone after the first rerun with the new flags.

## Binding rule — record output and result paths

Every figure / number produced here must carry its full provenance chain
(result artefact, postprocess script, `.prj` / `.pvd` / `.vtu` / `.log`) in the
run summary and the postprocess docstring. The status deck must show a source
line on every slide that plots or pulls this data. See
`../ANCHORS_MS33_ModelIV/AGENTS.md` and
`materialmodels/src/TPM/VK_SB_EURAD_DSM/AGENTS.md`.

Model III artefacts: `ms33_modelIII_gap2mm.prj` -> `ms33_modelIII_gap2mm.pvd`;
`ms33_postprocess_modelIII.py` -> `ms33_modelIII_gap_aperture.pdf`,
`ms33_modelIII_mean_stress.pdf`; log `anchors_modelIII_full.log`.
