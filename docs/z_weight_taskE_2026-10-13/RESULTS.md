# Task E results: Model VII, OFF / level / level+Z (2026-10-13)

Task: `docs/CLOUD_PROMPT_Z_weight_taskE_2026-10-13.md` item E. Design: `docs/dsm_z_weight_design_2026-10.md`.
Status words: computed (a run or a script), derived (subtraction of computed values), open.
No comparison with other teams; no interpretation beyond the differences stated.

## Provenance (computed)

- Source: branch `dsm_z_weight_2026-10`, commit cc5e71e12a (`git describe`: submission-config-anchors-T4.1-T4.2-2026-10-09-1-gcc5e71e12a).
  Results branch `dsm_z_weight_2026-10_taskE` cut from cc5e71e12a; no deck and no reference frame changed.
- Machine: Vinay's MacBook Pro (Apple Silicon, 18 cores, Darwin 27.0.0), Apple clang (`/usr/bin/c++`), not a cloud VM.
- Build: Release, Ninja, `OGS_BUILD_TESTING=ON`, `OGS_USE_MFRONT=OFF`, plus the local recipe flags
  (`-D_LIBCPP_TEMPLATE_VIS=`, shared CPM cache, homebrew libomp OpenMP, no PETSc/LIS/NetCDF); configure log in `logs/cmake_configure.log`.
  `bin/ogs` md5 01d24ddc645891c7441fb4e1d657128f, `lib/libRichardsMechanics.dylib` md5 cea8804c4cc2a11537b9a284b2508c7f.
- Unit tests: `testrunner --gtest_filter='RichardsMechanics.DSMZWeighted*:RichardsMechanics.DSMMaxwellPair*'`:
  2 tests, both PASSED (`DSMMaxwellPairIntegrabilityIdentity`, `DSMZWeightedPartnerMaxwellIdentity`); `logs/testrunner.log`.
- Runs: `run_taskE.sh` (six processes concurrently, `OMP_NUM_THREADS=1`, `OGS_ASM_THREADS` unset = 1 assembly thread),
  each from `Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelVII/` with its own `-o` folder (outside the repo; VTU series not committed).
  Logs show the switches as parsed: OFF `swelling_stress_K_level = true, swelling_stress_form = step`; level and level+Z
  `swelling_stress_form = level`; only the two `_levelz` logs carry the INFO line `DIAGNOSTIC, NOT ADOPTED (Z-weighted partner, 2026-10-13)`.

## Step counts and wall times (computed)

| deck | accepted | rejected | OGS "Simulation completed" [s] | process wall [s] |
|---|---|---|---|---|
| ms33_modelVII_freeswelling (OFF) | 881 | 0 | 205.1 | 212 |
| ms33_modelVII_freeswelling_level | 837 | 0 | 211.1 | 218 |
| ms33_modelVII_freeswelling_levelz | 833 | 0 | 215.5 | 223 |
| ms33_modelVII_freeswelling_ladder550 (OFF) | 1150 | 0 | 250.7 | 258 |
| ms33_modelVII_freeswelling_ladder550_level | 1097 | 0 | 253.3 | 260 |
| ms33_modelVII_freeswelling_ladder550_levelz | 1106 | 0 | 261.4 | 269 |

The six ran concurrently on one machine; wall times are not isolated timings.

## OFF reproduces the record (computed)

- Step counts equal the record in the deck headers: 881 / 0 (240-d deck) and 1150 / 0 (ladder550).
- `vtkdiff` of the final frame against the committed reference beside the deck, every field of the deck's
  `<test_definition>` at its tolerances (displacement 1e-9/1e-7; saturation 1e-8/1e-8; porosity, transport_porosity,
  micro_porosity, micro_water_content 1e-8/1e-7; dry_density_solid 1e-5/1e-8; sigma, swelling_stress, pressure,
  micro_pressure 1e3/1e-2; abs/rel): all 11 fields PASS for both decks (`logs/vtkdiff_OFF.log`).
  `ms33_modelVII_freeswelling_ts_881_t_20736000.000000.vtu` and `ms33_modelVII_freeswelling_ladder550_ts_1150_t_47520000.000000.vtu`:
  every vtkdiff norm is 0, and a direct comparison of all 15 point-data arrays gives max |difference| = 0 (bitwise equal).

## Definitions of the reported quantities (computed by `extract_probes.py` from written point data)

- Probes: nearest mesh node to r = 0 and z = 0.070 (top, node 1040), 0.04025 (centre, node 598), 0.0105 m (bottom, node 156);
  all three are exact nodes of `ms33_cylinder_r25_h70.vtu` (coordinates as written in the VTU, i.e. the undeformed mesh).
- rho_d = `dry_density_solid`; e = 2780/rho_d - 1 (rho_s = 2780 kg/m3 from the task prompt; the deck states phi0 = 1 - 1600/2780);
  phi_M = `transport_porosity`; n_l/phi = `micro_water_content` / `porosity`;
  swelling stress = -tr(`swelling_stress`)/3, mean total stress = -tr(`sigma`)/3 (compression positive; trace = the first
  three of the four written components rr, zz, thetatheta).
- "phi_M <= 1e-6 first at": time of the first OUTPUT frame with `transport_porosity` <= 1e-6 at the probe; resolution is the
  output grid (5-d frames from 5 d to 200 d in all six decks), not the step at which it happened.
- Histories at every output frame: `csv/<prefix>_<top|centre|bottom>.csv`; per-run summaries `csv/<prefix>_summary.json`.
- Tables below: `make_tables.py csv/` (regenerates this section exactly).

## Tables (computed; difference table derived)

#### 240-d decks, t = 200 d

| run | probe | rho_d [kg/m3] | e [-] | phi_M [-] | n_l/phi [-] | -tr(sigma_sw)/3 [MPa] | -tr(sigma)/3 [MPa] |
|---|---|---|---|---|---|---|---|
| OFF | top | 1366.31 | 1.0347 | 0 | 1.0000 | 6.7501 | 0.0567 |
| OFF | centre | 1365.96 | 1.0352 | 0 | 1.0000 | 6.7742 | 0.0698 |
| OFF | bottom | 1366.71 | 1.0341 | 0 | 1.0000 | 6.7400 | 0.0591 |
| level | top | 1275.51 | 1.1795 | 0 | 1.0000 | 9.7426 | 0.0667 |
| level | centre | 1275.51 | 1.1795 | 0 | 1.0000 | 9.7426 | 0.0666 |
| level | bottom | 1275.52 | 1.1795 | 0 | 1.0000 | 9.7428 | 0.0666 |
| level+Z | top | 1275.50 | 1.1795 | 0 | 1.0000 | 9.7425 | 0.0667 |
| level+Z | centre | 1275.49 | 1.1795 | 0 | 1.0000 | 9.7424 | 0.0666 |
| level+Z | bottom | 1275.52 | 1.1795 | 0 | 1.0000 | 9.7427 | 0.0666 |

#### 240-d decks, t = 240 d

| run | probe | rho_d [kg/m3] | e [-] | phi_M [-] | n_l/phi [-] | -tr(sigma_sw)/3 [MPa] | -tr(sigma)/3 [MPa] |
|---|---|---|---|---|---|---|---|
| OFF | top | 1367.84 | 1.0324 | 0 | 1.0000 | 6.7716 | 0.1239 |
| OFF | centre | 1367.50 | 1.0329 | 0 | 1.0000 | 6.7961 | 0.1367 |
| OFF | bottom | 1368.26 | 1.0318 | 0 | 1.0000 | 6.7645 | 0.1259 |
| level | top | 1277.95 | 1.1754 | 0 | 1.0000 | 9.7346 | 0.1347 |
| level | centre | 1277.97 | 1.1753 | 0 | 1.0000 | 9.7350 | 0.1343 |
| level | bottom | 1278.15 | 1.1750 | 0 | 1.0000 | 9.7375 | 0.1336 |
| level+Z | top | 1277.94 | 1.1754 | 0 | 1.0000 | 9.7346 | 0.1344 |
| level+Z | centre | 1277.96 | 1.1753 | 0 | 1.0000 | 9.7350 | 0.1340 |
| level+Z | bottom | 1278.15 | 1.1750 | 0 | 1.0000 | 9.7375 | 0.1336 |

#### 240-d decks, end of run and first frame with phi_M <= 1e-6

| run | probe | node (r, z) [m] | t_end [d] | rho_d end [kg/m3] | e end [-] | phi_M end [-] | phi_M <= 1e-6 first at [d] |
|---|---|---|---|---|---|---|---|
| OFF | top | 1040 (0, 0.07) | 240 | 1367.84 | 1.0324 | 0 | 65 |
| OFF | centre | 598 (0, 0.04025) | 240 | 1367.50 | 1.0329 | 0 | 55 |
| OFF | bottom | 156 (0, 0.0105) | 240 | 1368.26 | 1.0318 | 0 | 25 |
| level | top | 1040 (0, 0.07) | 240 | 1277.95 | 1.1754 | 0 | 75 |
| level | centre | 598 (0, 0.04025) | 240 | 1277.97 | 1.1753 | 0 | 60 |
| level | bottom | 156 (0, 0.0105) | 240 | 1278.15 | 1.1750 | 0 | 25 |
| level+Z | top | 1040 (0, 0.07) | 240 | 1277.94 | 1.1754 | 0 | 80 |
| level+Z | centre | 598 (0, 0.04025) | 240 | 1277.96 | 1.1753 | 0 | 65 |
| level+Z | bottom | 156 (0, 0.0105) | 240 | 1278.15 | 1.1750 | 0 | 30 |

#### ladder550 decks, t = 200 d

| run | probe | rho_d [kg/m3] | e [-] | phi_M [-] | n_l/phi [-] | -tr(sigma_sw)/3 [MPa] | -tr(sigma)/3 [MPa] |
|---|---|---|---|---|---|---|---|
| OFF | top | 1366.31 | 1.0347 | 0 | 1.0000 | 6.7501 | 0.0567 |
| OFF | centre | 1365.96 | 1.0352 | 0 | 1.0000 | 6.7742 | 0.0698 |
| OFF | bottom | 1366.71 | 1.0341 | 0 | 1.0000 | 6.7400 | 0.0591 |
| level | top | 1275.51 | 1.1795 | 0 | 1.0000 | 9.7426 | 0.0667 |
| level | centre | 1275.51 | 1.1795 | 0 | 1.0000 | 9.7426 | 0.0666 |
| level | bottom | 1275.52 | 1.1795 | 0 | 1.0000 | 9.7428 | 0.0666 |
| level+Z | top | 1275.50 | 1.1795 | 0 | 1.0000 | 9.7425 | 0.0667 |
| level+Z | centre | 1275.49 | 1.1795 | 0 | 1.0000 | 9.7424 | 0.0666 |
| level+Z | bottom | 1275.52 | 1.1795 | 0 | 1.0000 | 9.7427 | 0.0666 |

#### ladder550 decks, t = 240 d

| run | probe | rho_d [kg/m3] | e [-] | phi_M [-] | n_l/phi [-] | -tr(sigma_sw)/3 [MPa] | -tr(sigma)/3 [MPa] |
|---|---|---|---|---|---|---|---|
| OFF | top | 1367.81 | 1.0324 | 0 | 1.0000 | 6.7691 | 0.1234 |
| OFF | centre | 1367.47 | 1.0330 | 0 | 1.0000 | 6.7932 | 0.1364 |
| OFF | bottom | 1368.21 | 1.0318 | 0 | 1.0000 | 6.7590 | 0.1258 |
| level | top | 1277.76 | 1.1757 | 0 | 1.0000 | 9.7330 | 0.1334 |
| level | centre | 1277.76 | 1.1757 | 0 | 1.0000 | 9.7330 | 0.1333 |
| level | bottom | 1277.77 | 1.1757 | 0 | 1.0000 | 9.7332 | 0.1333 |
| level+Z | top | 1277.75 | 1.1757 | 0 | 1.0000 | 9.7329 | 0.1333 |
| level+Z | centre | 1277.74 | 1.1757 | 0 | 1.0000 | 9.7328 | 0.1333 |
| level+Z | bottom | 1277.76 | 1.1757 | 0 | 1.0000 | 9.7331 | 0.1333 |

#### ladder550 decks, t = 550 d

| run | probe | rho_d [kg/m3] | e [-] | phi_M [-] | n_l/phi [-] | -tr(sigma_sw)/3 [MPa] | -tr(sigma)/3 [MPa] |
|---|---|---|---|---|---|---|---|
| OFF | top | 1367.82 | 1.0324 | 0 | 1.0000 | 6.7699 | 0.1234 |
| OFF | centre | 1367.47 | 1.0329 | 0 | 1.0000 | 6.7940 | 0.1364 |
| OFF | bottom | 1368.22 | 1.0318 | 0 | 1.0000 | 6.7598 | 0.1258 |
| level | top | 1277.81 | 1.1756 | 0 | 1.0000 | 9.7337 | 0.1334 |
| level | centre | 1277.81 | 1.1756 | 0 | 1.0000 | 9.7337 | 0.1333 |
| level | bottom | 1277.83 | 1.1756 | 0 | 1.0000 | 9.7339 | 0.1333 |
| level+Z | top | 1277.80 | 1.1756 | 0 | 1.0000 | 9.7335 | 0.1333 |
| level+Z | centre | 1277.80 | 1.1756 | 0 | 1.0000 | 9.7335 | 0.1333 |
| level+Z | bottom | 1277.82 | 1.1756 | 0 | 1.0000 | 9.7338 | 0.1333 |

#### ladder550 decks, end of run and first frame with phi_M <= 1e-6

| run | probe | node (r, z) [m] | t_end [d] | rho_d end [kg/m3] | e end [-] | phi_M end [-] | phi_M <= 1e-6 first at [d] |
|---|---|---|---|---|---|---|---|
| OFF | top | 1040 (0, 0.07) | 550 | 1367.82 | 1.0324 | 0 | 65 |
| OFF | centre | 598 (0, 0.04025) | 550 | 1367.47 | 1.0329 | 0 | 55 |
| OFF | bottom | 156 (0, 0.0105) | 550 | 1368.22 | 1.0318 | 0 | 25 |
| level | top | 1040 (0, 0.07) | 550 | 1277.81 | 1.1756 | 0 | 75 |
| level | centre | 598 (0, 0.04025) | 550 | 1277.81 | 1.1756 | 0 | 60 |
| level | bottom | 156 (0, 0.0105) | 550 | 1277.83 | 1.1756 | 0 | 25 |
| level+Z | top | 1040 (0, 0.07) | 550 | 1277.80 | 1.1756 | 0 | 80 |
| level+Z | centre | 598 (0, 0.04025) | 550 | 1277.80 | 1.1756 | 0 | 65 |
| level+Z | bottom | 156 (0, 0.0105) | 550 | 1277.82 | 1.1756 | 0 | 30 |

#### Differences (derived: subtraction of the computed rows)

| deck set | t [d] | probe | pair | rho_d [kg/m3] | e [-] | n_l/phi [-] | -tr(sigma_sw)/3 [kPa] | -tr(sigma)/3 [kPa] |
|---|---|---|---|---|---|---|---|---|
| 240-d | 200 | top | level - OFF | -90.796 | 0.14484 | 0.00e+00 | 2992.50 | 9.938 |
| 240-d | 200 | top | level+Z - level | -0.014 | 0.00002 | 0.00e+00 | -0.17 | -0.006 |
| 240-d | 200 | centre | level - OFF | -90.452 | 0.14432 | 0.00e+00 | 2968.41 | -3.130 |
| 240-d | 200 | centre | level+Z - level | -0.016 | 0.00003 | 0.00e+00 | -0.19 | 0.002 |
| 240-d | 200 | bottom | level - OFF | -91.186 | 0.14541 | 0.00e+00 | 3002.79 | 7.536 |
| 240-d | 200 | bottom | level+Z - level | -0.005 | 0.00001 | 0.00e+00 | -0.07 | 0.020 |
| 240-d | 240 | top | level - OFF | -89.891 | 0.14296 | 0.00e+00 | 2962.95 | 10.764 |
| 240-d | 240 | top | level+Z - level | -0.011 | 0.00002 | 0.00e+00 | 0.05 | -0.337 |
| 240-d | 240 | centre | level - OFF | -89.529 | 0.14242 | 0.00e+00 | 2938.88 | -2.481 |
| 240-d | 240 | centre | level+Z - level | -0.007 | 0.00001 | 0.00e+00 | 0.05 | -0.242 |
| 240-d | 240 | bottom | level - OFF | -90.111 | 0.14324 | 0.00e+00 | 2973.01 | 7.736 |
| 240-d | 240 | bottom | level+Z - level | -0.002 | 0.00000 | 0.00e+00 | 0.01 | -0.059 |
| ladder550 | 200 | top | level - OFF | -90.796 | 0.14484 | 0.00e+00 | 2992.50 | 9.938 |
| ladder550 | 200 | top | level+Z - level | -0.014 | 0.00002 | 0.00e+00 | -0.17 | -0.006 |
| ladder550 | 200 | centre | level - OFF | -90.452 | 0.14432 | 0.00e+00 | 2968.41 | -3.130 |
| ladder550 | 200 | centre | level+Z - level | -0.016 | 0.00003 | 0.00e+00 | -0.19 | 0.002 |
| ladder550 | 200 | bottom | level - OFF | -91.186 | 0.14541 | 0.00e+00 | 3002.79 | 7.536 |
| ladder550 | 200 | bottom | level+Z - level | -0.005 | 0.00001 | 0.00e+00 | -0.07 | 0.020 |
| ladder550 | 240 | top | level - OFF | -90.054 | 0.14324 | 0.00e+00 | 2963.93 | 9.940 |
| ladder550 | 240 | top | level+Z - level | -0.015 | 0.00002 | 0.00e+00 | -0.17 | -0.006 |
| ladder550 | 240 | centre | level - OFF | -89.709 | 0.14273 | 0.00e+00 | 2939.84 | -3.129 |
| ladder550 | 240 | centre | level+Z - level | -0.016 | 0.00003 | 0.00e+00 | -0.19 | 0.002 |
| ladder550 | 240 | bottom | level - OFF | -90.444 | 0.14382 | 0.00e+00 | 2974.21 | 7.533 |
| ladder550 | 240 | bottom | level+Z - level | -0.005 | 0.00001 | 0.00e+00 | -0.07 | 0.020 |
| ladder550 | 550 | top | level - OFF | -90.005 | 0.14316 | 0.00e+00 | 2963.77 | 9.941 |
| ladder550 | 550 | top | level+Z - level | -0.015 | 0.00002 | 0.00e+00 | -0.17 | -0.006 |
| ladder550 | 550 | centre | level - OFF | -89.660 | 0.14265 | 0.00e+00 | 2939.69 | -3.130 |
| ladder550 | 550 | centre | level+Z - level | -0.016 | 0.00003 | 0.00e+00 | -0.19 | 0.002 |
| ladder550 | 550 | bottom | level - OFF | -90.391 | 0.14373 | 0.00e+00 | 2974.07 | 7.535 |
| ladder550 | 550 | bottom | level+Z - level | -0.005 | 0.00001 | 0.00e+00 | -0.07 | 0.020 |

## Differences stated (derived from the tables; no interpretation)

- level vs OFF: at 200 d rho_d is lower by 90.5 to 91.2 kg/m3 (e higher by 0.144 to 0.145) and the mean swelling stress is
  higher by 2.97 to 3.00 MPa at all three probes; the mean total stress differs by -3.1 to +9.9 kPa. phi_M closes at the probe
  frames 75/60/25 d (level) vs 65/55/25 d (OFF) for top/centre/bottom. At 240 d and 550 d (both deck sets) rho_d is lower by
  89.5 to 90.4 kg/m3 (e higher by 0.142 to 0.144), the mean swelling stress higher by 2.94 to 2.97 MPa, and the mean total
  stress differs by -3.1 to +10.8 kPa.
- level+Z vs level: rho_d differs by at most 0.016 kg/m3, e by at most 3e-5, mean swelling stress by at most 0.19 kPa and
  mean total stress by at most 0.34 kPa at the reported times; phi_M closes one output frame (5 d) later at every probe
  (80/65/30 d vs 75/60/25 d). Accepted steps: 833 vs 837 (240-d), 1106 vs 1097 (ladder550); 0 rejected in all runs.
- At every reported time (200, 240, 550 d) and every probe, phi_M = 0 exactly and n_l/phi = 1 in all six runs.

## Open

- Whether the level form and/or the Z partner are adopted for the ANCHORS submission is Vinay's call (design note section 6).
- Effects inside the 5-d output intervals (exact closure time of phi_M) are not resolved by these frames.
