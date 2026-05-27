# ANCHORS MS33 — native DSM + MFront MCC variant

Native branch (`dsm_native_hierarchical`). Combines the native `<potential_exchange>` (DSM disjoining-pressure source)
with the MFront `ModCamClay_semiExpl_constE` mechanical constitutive.

## Status

| Model | dd | K (J/kg) | Status | sigma_swell at convergence |
|---|---|---|---|---|
| I (dd1400) | 1400 | 5500 | ✓ runs | 1.120 MPa (Dixon median 1.12) |
| I (dd1600) | 1600 | 13050 | ✓ runs | 2.610 MPa (Dixon median 2.61) |
| I (dd1800) | 1800 | 31280 | ✓ runs | 6.093 MPa (Dixon median 6.09) |
| III (gap2mm) | 1600 | 13050 | ✗ FAILS at TS~64 (mid-ramp) | n/a |
| IV (pellets) | 1600 (clay), 900 (pellet) | 13050 | ✗ FAILS at TS 1 (initial state, pellet Vr0=3.09) | n/a |
| VII (free-swell) | 1600 | 13050 | ✗ FAILS at TS~78 (end of suction ramp, Bishop's-cutoff transition) | n/a |

For Model I (isotropic, single-element axisymmetric), MCC stays elastic throughout
(pc_char = 1e10 Pa is never reached under isotropic swelling), so the MCC result is
numerically identical to the LE_NATIVE variant — but the MCC constitutive machinery
is wired up and a non-zero PreConsolidationPressure / VolumeRatio state is tracked.

For Models III / IV / VII, the OGS-MFront-MCC integrator fails on the trial state
handed to it by the macro-Newton (failure persists down to dt=0.01 s). Both
ConstE and absP MFront-MCC variants exhibit this. This is an upstream-MFront
integrator robustness issue, not a PRJ-side fix.

## How to run

```bash
NATIVE_OGS=/path/to/build/bin/ogs
$NATIVE_OGS ANCHORS_MS33_MCC_NATIVE/ModelI_dd1600/ms33_modelI_dd1600_mcc_native.prj
```

K calibration anchored to Dixon et al. (2023) MX-80 median at each dd
(per-dd calibration as agreed for Model I; clay zone of III/IV/VII inherits the
dd=1.6 anchor K=13050 J/kg).
