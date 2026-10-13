# Cloud session prompt: task E of the Z-weight brief (2026-10-13)

You start from github.com/VinayGK/ogs, branch `dsm_z_weight_2026-10` (tag 8bbcace35c + the Z-weighted partner,
see docs/dsm_z_weight_design_2026-10.md on the branch). Tasks A to D of the brief are DONE on that branch and
recorded; do not redo them. Your task is E only.

1. Build (Release; MFront may be OFF; OGS_BUILD_TESTING on so `bin/testrunner` exists):
   cmake -G Ninja -S . -B build -DCMAKE_BUILD_TYPE=Release -DOGS_BUILD_TESTING=ON -DOGS_USE_MFRONT=OFF
   ninja -C build ogs vtkdiff testrunner
   Then: build/bin/testrunner --gtest_filter='RichardsMechanics.DSMZWeighted*' must PASS.
2. Run, one thread each (OMP_NUM_THREADS=1), from Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelVII/:
   OFF  : ms33_modelVII_freeswelling.prj          -> expect 881 accepted / 0 rejected steps (record on the MBP).
          NOTE: the submitted deck runs the STEP form with K-level (swelling_stress_k_level = true, no
          swelling_stress_form tag). The Z partner belongs to the LEVEL form and is refused otherwise, so
          "ON" is two changes; the middle deck separates them.
   LEVEL: ms33_modelVII_freeswelling_level.prj    (+ swelling_stress_form = level only)
   ON   : ms33_modelVII_freeswelling_levelz.prj   (+ level form + micro_potential_partner_z_weight = true)
   and the 550 d load variant: ms33_modelVII_freeswelling_ladder550{,_level,_levelz}.prj.
3. OFF must reproduce the record: vtkdiff of the OFF final frame against the reference frame committed
   beside the deck (ms33_modelVII_freeswelling_ts_881_t_20736000.000000.vtu; ladder550: ts_1150) at the
   tolerances of the deck's test_definition. If it does not, STOP and report the differences; do not
   continue to the comparison.
4. Report per run (probes: top, centre, bottom on the axis r = 0: z = 0.070, 0.04025, 0.0105 m; nearest
   node): end dry density, void ratio e = rho_s/rho_d - 1 (rho_s = 2780), macro porosity, the time the
   macro porosity reaches zero, micro water content against the pore space (n_l / phi), swelling stress and
   mean total stress (compression positive) at 200 d (and 240 d, and the ladder end), plus CSV histories
   of these at every output frame. Values as computed; no comparison with other teams (that is done locally).
5. Rules: no value, sign or parameter without a source (the deck, the code on the branch, this note);
   predicted is not verified; one writer per file; no Co-Authored-By trailer; commit your results on a
   branch `dsm_z_weight_2026-10_taskE`; push only that branch to github; never force-push.
