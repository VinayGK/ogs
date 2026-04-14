# DSM Leftover Cleanup (2026-04-14)

Scope: DSM-related assets in `ProcessLib/RichardsMechanics`, `MaterialLib/SolidModels/MFront`, and `Tests/Data/RichardsMechanics`.

## What was cleaned

1. Removed stale XML key from BEACON report deck:
   - File: `BEACON_report_figures/1a01/1a01_native_reportplot.prj`
   - Change: removed `<potential_role_mapping>micro_macro_potential_role_mapping_mode</potential_role_mapping>`.
   - Reason: this key is no longer consumed by the current DSM parameter parser.

2. Removed stale mention from ANCHORS notes:
   - File: `ANCHORS_MS33_ModelI/README.md`
   - Change: dropped `potential_role_mapping=...` from native-vs-MFront notes.
   - Reason: documentation now matches the active parameter set.

3. Removed outdated wording in CIEMAT overlay script:
   - File: `ANCHORS_MS33_ModelI/plot_ms33_vs_ciemat_swelling_pressure_dd.py`
   - Change: replaced "scaffold" wording with "reduced model/setup".
   - Reason: avoids temporary/prototype wording for now-versioned benchmark assets.

## Intentionally kept (not leftovers)

1. `Tests/Data/RichardsMechanics/Mathematica_DSM_Driver/*`
   - Kept because these files explicitly document and run the Mathematica-side DSM driver.

2. Generic OGS `notebook` mentions outside DSM process logic
   - Example: `ProcessLib/Tests.cmake`, `ProcessLib/HeatConduction/Tests.cmake`.
   - Kept because these are global OGS test infrastructure references, not DSM transfer leftovers.

## Verification after cleanup

- Re-scan for `potential_role_mapping` in DSM code/tests: only historical/binary outputs were excluded and no active DSM source/deck usage remains.
- No constitutive or solver logic changed in this cleanup; this pass is naming/config/documentation alignment only.
