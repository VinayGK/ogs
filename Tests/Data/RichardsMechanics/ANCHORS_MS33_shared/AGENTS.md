# ANCHORS_MS33_shared — auditable extraction framework

Worklog for the harmonized post-processing layer used by all four EURAD-2 MS33
models (I, III, IV, VII) on both the MFront bridge branch (current data) and
the native dsm_micromacro branch (data pending).

## What is here

```
audit_extract_ms33.py    Main extraction tool (PVD → CSV + JSON + figures).
probes_modelI.yaml       Probe spec for Model I (single-element confined Villar swell).
probes_modelIII.yaml     Probe spec for Model III (clay + 2 mm soft-elastic gap).
probes_modelIV.yaml      Probe spec for Model IV (clay/pellet stack).
probes_modelVII.yaml     Probe spec for Model VII (free-swelling cylinder).
AGENTS.md                This file.
```

All scientific anchors (zone definitions, sampling locations, marker times,
triptych field selection, units) live in the YAML probe specs. The script
treats them as data — auditable through diff, no code edit required to change
what a plot reports.

## Provenance summary

For every model the script emits `model<X>_<tag>_summary.json` next to the
outputs, containing:
- OGS repo git HEAD
- absolute PVD path + md5
- probes YAML path + md5
- list of CSV columns
- list of fields requested but missing from the VTUs
- list of generated figures
- recap of the probe spec relevant to interpretation (zones, marker days,
  stress-column doc)

This means any figure can be traced back to the exact (code, data, probe spec)
triplet that produced it.

## Universal rules in this framework

(Derived from the design Q&A captured in the parent conversation,
2026-05-28; cf. probes_*.yaml headers for per-model citations.)

| Rule | Choice | Notes |
| :--- | :----- | :---- |
| Zone averaging | Volume-avg plotted, centroid stored as audit column | r-weighted axisymmetric volume; centroid = cell closest to zone centroid. |
| Time-axis | Continuous curve + markers at fixed_output_times | markers_days lives per-YAML; default {20, 100, 200} d. |
| Stress label | `mean_stress_MPa` (neutral) | Under BishopsSaturationCutoff cutoff=1 (all PRJs) sigma_eff ≡ sigma_total numerically. Rename if cutoff changes. |
| Triptych fields | pressure / displacement magnitude / mean compressive swelling stress | Matches "Pressure, Displacement, Swelling-Stress Fields" frame titles in VK_SB_EURAD_DSM_output.tex (frames 184, 428, 472, 509). |
| Model IV interface | Line probe at z=0.035 m, r ∈ [0, 0.025], 21 samples | Mesh-convergent; single scalar per timestep. |
| Model III gap aperture | 2 mm − max(u_r) at clay/gap interface (r=0.025) | Mesh-robust outer-free-surface rule applied to the clay zone. |
| Native scope | Branch-agnostic; --tag {mfront,native} | Picks up native PVDs the moment they exist; no special-casing in the script. |

## Stress flavor — source-code citation

The VTU field `sigma` is `sigma_eff` as defined by
`ProcessLib/ConstitutiveRelations/EffectiveStressData.h:28`
(reflection: `reflectWithName("sigma", &Self::sigma_eff)`). Under
`MaterialLib/MPL/Properties/BishopsSaturationCutoff.cpp:38`
(`return S_L < S_L_max_ ? 0. : 1.`) with `cutoff_value=1` in every MS33 PRJ,
χ = 0 on the unsaturated branch and χ = 1 only at S = 1. Combined with p_L = 0
at the saturated endpoint, this makes sigma_eff and sigma_total *numerically
identical* throughout the simulation in these particular PRJs. The label
`mean_stress_MPa` carries no claim about which interpretation applies — both
are correct for these data.

## Validation runs (smoke tests)

| Date       | Model | PVD                                                                   | Final-state value                | Reference / check |
| :--------- | :---- | :-------------------------------------------------------------------- | :------------------------------- | :---------------- |
| 2026-05-28 | I dd1600 | tex/.../ANCHORS_MS33_ModelI/ms33_modelI_dd1600.pvd                  | mean_stress = 5.823 MPa at 200 d | Villar target exp(6.77·1.6 − 9.07) ≈ 5.824 MPa — match < 0.05%. |
| 2026-05-28 | III   | tex/.../ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.pvd                | clay σ̄ = 5.295 MPa; gap_aperture = 1.71 mm at 200 d | Consistent with 0.29 mm of the 2 mm gap closed; clay slightly relaxed below Villar 1600 target because some swelling went into the gap. |
| 2026-05-28 | IV    | tex/.../ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.pvd                 | clay σ̄ = 7.57 MPa; pellet σ̄ = 0.99 MPa; interface line = 3.76 MPa at 200 d | Interface value sits between zone bulks as expected; deck cites 3.42 MPa (different sampling rule or older snapshot, ~10% drift). |
| 2026-05-28 | VII   | tex/.../ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.pvd          | σ̄ = 4.78 MPa, void ratio = 2.58 at 240 d | Free-swelling expanded mesh; ρ_d = 777 kg/m³ at final consistent with e=2.58 and ρ_s=2780. |

## Known limitations / open TODOs

- The deck triptych PNGs (`ms33_model{1,3,4,7}_fields_triptych.png`) are
  ParaView screenshots; this framework produces matplotlib triptychs with the
  same three fields but a different visual style. If pixel-identical ParaView
  output is needed, a `.pvsm` state file would have to be added separately —
  not a blocker for auditability since the matplotlib output is reproducible
  from the script.
- Native-branch runs for Models I/III/IV/VII do not exist on disk; only
  Model V has native PVDs in `ogs-worktrees/dsm_native_hierarchical_wt/`. The
  framework will run against native PVDs as soon as they appear; no code
  change needed.
- Model VII canonical PVD (`ms33_modelVII_freeswelling.pvd`) only — the
  exploratory v0..v7 variants in the LE_RERUN directory are not part of the
  audit. The user runs those by hand for investigation.
- PVDs for Models III/IV/VII output only at fixed_output_times {0, 20 d,
  100 d, 200 d}, so the evolution PDFs show piecewise-linear interpolation
  between four points. For denser curves the PRJs would need `<output>` to
  declare every-N-step output rather than only fixed_output_times.

## Side change — ANCHORS_MS33_ModelI/summarize_ms33_model_i.py

The existing summarizer's column `mean_total_stress_MPa` is misnamed (the
underlying field is sigma_eff per the source-code finding above). The column
has been renamed to `mean_stress_MPa`, with `mean_total_stress_MPa` kept as a
legacy alias for backward compatibility with downstream tooling. The original
numerical value is unchanged; only the column name and an explanatory comment
were edited.

Cf. GUARDRAIL §5.1: user approved the relabel as part of the audit-framework
design loop on 2026-05-28.

## Usage examples

Single-case Model III:
```
python audit_extract_ms33.py \
    --model III \
    --pvd /path/to/ms33_modelIII_gap2mm.pvd \
    --tag mfront --probes probes_modelIII.yaml \
    -o out_III/
```

Model I three-density sweep (gives the Villar overlay):
```
python audit_extract_ms33.py --model I \
    --pvd .../ms33_modelI_dd1400.pvd --rho-dry 1400 \
    --pvd .../ms33_modelI_dd1600.pvd --rho-dry 1600 \
    --pvd .../ms33_modelI_dd1800.pvd --rho-dry 1800 \
    --tag mfront --probes probes_modelI.yaml \
    -o out_I/
```

Native-branch run (when native PVDs exist):
```
python audit_extract_ms33.py \
    --model IV --tag native \
    --pvd .../native/ms33_modelIV_pellets.pvd \
    --probes probes_modelIV.yaml -o out_IV_native/
```

## DONE / pending

- DONE (2026-05-28): probes_model{I,III,IV,VII}.yaml committed.
- DONE (2026-05-28): audit_extract_ms33.py committed; smoke-tested on all four MFront PVDs.
- DONE (2026-05-28): summarize_ms33_model_i.py column relabel to `mean_stress_MPa` (+ legacy alias).
- TODO: rerun Models I/III/IV/VII with output_every_n_step to densify the evolution curves between snapshots — out of scope here; needs OGS execution.
- TODO: produce native-branch PVDs for Models I/III/IV/VII to populate the `--tag native` half of every figure — needs OGS execution.
- TODO (optional): add a ParaView state file (`.pvsm`) if pixel-identical-to-deck triptychs are wanted; matplotlib triptychs from this framework are auditable but visually different.
