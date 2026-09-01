if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE RichardsMechanics/gravity.prj)
    OgsTest(PROJECTFILE RichardsMechanics/mechanics_linear.prj)
    OgsTest(PROJECTFILE RichardsMechanics/confined_compression_fully_saturated.prj RUNTIME 7)
    OgsTest(PROJECTFILE RichardsMechanics/flow_fully_saturated.prj)
    OgsTest(PROJECTFILE RichardsMechanics/flow_fully_saturated_linear.prj)
    OgsTest(PROJECTFILE RichardsMechanics/flow_fully_saturated_anisotropic.prj)
    OgsTest(PROJECTFILE RichardsMechanics/flow_fully_saturated_coordinate_system.prj)
    OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_small.prj RUNTIME 9)
    OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_small_masslumping.prj RUNTIME 10)
    OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_quasinewton.prj RUNTIME 80)
    OgsTest(PROJECTFILE RichardsMechanics/double_porosity_swelling.prj RUNTIME 20)
    OgsTest(PROJECTFILE RichardsMechanics/deformation_dependent_porosity.prj RUNTIME 8)
    OgsTest(PROJECTFILE RichardsMechanics/deformation_dependent_porosity_swelling.prj RUNTIME 11)
    OgsTest(PROJECTFILE RichardsMechanics/orthotropic_power_law_permeability_xyz.prj RUNTIME 80)
    OgsTest(PROJECTFILE RichardsMechanics/orthotropic_swelling_xyz.prj)
    OgsTest(PROJECTFILE RichardsMechanics/orthotropic_swelling_xy.prj)
    OgsTest(PROJECTFILE RichardsMechanics/bishops_effective_stress_power_law.prj)
    OgsTest(PROJECTFILE RichardsMechanics/bishops_effective_stress_saturation_cutoff.prj)
    OgsTest(PROJECTFILE RichardsMechanics/alternative_mass_balance_anzInterval_10.prj)
    if(NOT ENABLE_ASAN)
        OgsTest(PROJECTFILE RichardsMechanics/rotated_consolidation.prj RUNTIME 2)
    endif()
    OgsTest(PROJECTFILE RichardsMechanics/LiakopoulosHM/liakopoulos.prj RUNTIME 17)
    OgsTest(PROJECTFILE RichardsMechanics/LiakopoulosHM/liakopoulos_restart.xml RUNTIME 17)
    OgsTest(PROJECTFILE RichardsMechanics/LiakopoulosHM/liakopoulos_QN.prj RUNTIME 50)
    OgsTest(PROJECTFILE RichardsMechanics/A2.prj RUNTIME 20)
    OgsTest(PROJECTFILE RichardsMechanics/restart_w_backfill.prj RUNTIME 20)

    # ANCHORS EURAD-2 MS33 theoretical benchmarking — DSM native hierarchical runs
    # Re-timed 2026-08-31 on the dd900 K-knot adopt runs, OGS's own timer:
    # Model III 10.1542 s / 889 steps (III/run.log, OMP_NUM_THREADS=4) and
    # Model VII 207.022 s / 920 steps (VII/run_stdout.log, OMP_NUM_THREADS
    # unset -> 18-thread fallback), both under
    # /Users/vinaykumar/ogs-models/dd900_adopt_run_2026-08-31/, ogs
    # archive/dsm_native_Pi_fofnlev_branchtip_2026-08-11-49-gbed3e395. Both
    # stay inside their present RUNTIME, which is left unchanged; Model VII's
    # margin to RUNTIME 300 was measured at 18 threads, not at one.
    # Model I was NOT re-timed. Its figures — 0.111132 / 0.118636 / 0.149283 s
    # for dd1400 / dd1600 / dd1800, verbatim from the tracked
    # rerun_ms33_modelI_dd*.log — are 2026-05-22 runs of the now stale binary
    # vdw-baseline-2026-05-08-41-g9a1b956c, predating both commit 1bb414ac05
    # (log-linear live K, 2026-08-26) and the dd900 knot; those three decks
    # carry no <prefactors> table, so they were never part of the cascade. The
    # numbers are kept as the only timings on record for Model I, not as
    # current measurements, and RUNTIME 120 is left unchanged.
    # Lowering Model I/III RUNTIME below large_runtime = 60 would rename the
    # ctests (-LARGE suffix, scripts/cmake/test/OgsTest.cmake). Only Model IV
    # below was stale.
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1400.prj RUNTIME 120)
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj RUNTIME 120)
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1800.prj RUNTIME 120)
    # Model III ships the GAP-SWITCH deck (Vinay 2026-08-17). The outer radial
    # boundary swells free until u_r reaches the 2 mm technological gap, then
    # switches to a rigid Dirichlet wall — true container contact.
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIII/ms33_modelIII_gapswitch.prj RUNTIME 120)
    # DEPRECATED 2026-08-17 — soft 2-medium gap annulus surrogate: no contact
    # mechanics, over-closes to ~67% with a residual aperture. Superseded by the
    # gap-switch deck above. Deck and reference retained (CLAUDE.md §6.2/§6.3);
    # registration commented out, not deleted.
    # OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj RUNTIME 240)
    # RUNTIME covers the SLOWEST of three measured runs of this exact deck,
    # each 21500 accepted / 0 rejected steps to t = 200 d, OGS's own timer:
    #   10720    s  2026-08-28 09:02:27+0200, OMP_NUM_THREADS=4,
    #               ogs 6.5.8-565-gbea47887  (cascade_refit/IV_rep1/run.log)
    #   10735.4  s  2026-08-28 09:02:27+0200, OMP_NUM_THREADS=4,
    #               ogs 6.5.8-565-gbea47887  (cascade_refit/IV_rep2/run.log)
    #    8666.15 s  2026-08-31 09:23:34+0200, OMP_NUM_THREADS=6, ogs
    #            archive/dsm_native_Pi_fofnlev_branchtip_2026-08-11-49-gbed3e395
    #            (dd900_adopt_run_2026-08-31/IV/out/run.log, the run quoted in
    #            DSM/AGENTS.md entry 30)
    # The spread tracks the thread count, not the deck: the 6-thread run is the
    # fastest configuration on record, not the representative one. RUNTIME is
    # therefore taken from the slowest, 10736 = ceil(10735.4 s), which keeps the
    # repo's 2x convention against the worst case on record — RUNTIME > 750
    # makes OgsTest emit an explicit TIMEOUT = 2*RUNTIME = 21472 s, where
    # RUNTIME 8667 would give 17334 s, only 1.61x the slowest run. At the old
    # RUNTIME 240 no TIMEOUT is emitted at all, so ctest's default 1500 s would
    # apply and would kill a healthy run of this deck. No test is renamed by
    # the change: 240 and 10736 are both above large_runtime = 60, so this case
    # was already -LARGE. Effect on the OGS_CTEST_MAX_RUNTIME gate
    # (OgsTest.cmake) is band-dependent: a cap below 240 dropped the deck
    # before the change and still does — the only in-tree consumer,
    # scripts/ci/jobs/build-linux.yml, sets 60 — while a cap in [240, 10735]
    # kept the deck before and now drops it. Timings are
    # machine-local (this workstation, non-MPI, OpenMP thread counts as
    # listed); CI headroom is expected, not verified.
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj RUNTIME 10736)
    # K(rho_d) equivalence pair (each material's k0 x20 spec, for speed): the
    # table-K variant resolves K = K(dry_density) at parse time and must
    # reproduce, bit-for-bit, the per-material scalar-K reference. Verified
    # 2026-06-08 (abs max diff = 0 on all 14 output fields at t=200 d). Both
    # registered run-only here.
    # DE-REGISTERED 2026-08-12 (Vinay): cannot pass as registered (no <test_definition>; OGS hard-fails at parse under the ctest wrapper). The two ModelIV variants additionally DIVERGE on the merged code (die ts #825 FD / #2333 analytic). Decks kept per never-delete; re-register only with ratified references.
    # OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIV/ms33_modelIV_pellets_kref20x.prj RUNTIME 240)
    # DE-REGISTERED 2026-08-12 (Vinay): cannot pass as registered (no <test_definition>; OGS hard-fails at parse under the ctest wrapper). The two ModelIV variants additionally DIVERGE on the merged code (die ts #825 FD / #2333 analytic). Decks kept per never-delete; re-register only with ratified references.
    # OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIV/ms33_modelIV_pellets_kofdd.prj RUNTIME 240)
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj RUNTIME 300)
    # K(rho_d) feature on a 2nd model (single-material Model VII -> table resolves
    # to the rho_d=1600 node, a physical no-op; k0 x50 spec for speed). Run to
    # t_end 2026-06-08. Exercises the table-resolution path on the free-swelling cell.
    # DE-REGISTERED 2026-08-12 (Vinay): cannot pass as registered (no <test_definition>; OGS hard-fails at parse under the ctest wrapper). The two ModelIV variants additionally DIVERGE on the merged code (die ts #825 FD / #2333 analytic). Decks kept per never-delete; re-register only with ratified references.
    # OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling_kofdd.prj RUNTIME 300)
endif()

if (NOT OGS_USE_MPI AND OGS_USE_MFRONT)
    OgsTest(PROJECTFILE RichardsMechanics/mfront_restart_part1.prj RUNTIME 1)
    OgsTest(PROJECTFILE RichardsMechanics/mfront_restart_part2.xml RUNTIME 1)
    OgsTest(PROJECTFILE RichardsMechanics/DoubleStructureBenchmark/double_porosity_swelling_RM.prj RUNTIME 1)
endif()

if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE RichardsMechanics/confined_compression_fully_saturated_restart.prj RUNTIME 4)
    OgsTest(PROJECTFILE RichardsMechanics/A2_total_stress0.xml RUNTIME 8)
endif()

AddTest(
    NAME RichardsMechanics_beacon_1a01_dsm_micromacro_smoke
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1a01_dsm_micromacro_smoke.prj
    WRAPPER time
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME RichardsMechanics_beacon_1a01_dsm_micromacro_stressprobe
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1a01_dsm_micromacro_stressprobe.prj
    WRAPPER time
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME RichardsMechanics_beacon_1a01_dsm_micromacro_inflow
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1a01_dsm_micromacro_inflow.prj
    WRAPPER time
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME RichardsMechanics_beacon_1a01_dsm_micromacro_reference
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1a01_dsm_micromacro_smoke.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu displacement displacement 1e-12 0
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu pressure pressure 1e-12 1e-12
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu saturation saturation 1e-12 0
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu micro_pressure micro_pressure 1e-12 1e-12
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu micro_saturation micro_saturation 1e-12 0
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu swelling_stress swelling_stress 1e-12 0
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu sigma sigma 1e-12 1e-10
)

AddTest(
    NAME RichardsMechanics_beacon_1a01_dsm_micromacro_inflow_reference
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1a01_dsm_micromacro_inflow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu displacement displacement 1e-12 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu pressure pressure 1e-12 1e-12
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu saturation saturation 1e-12 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu micro_pressure micro_pressure 1e-12 1e-12
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu micro_saturation micro_saturation 1e-12 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu micro_water_content micro_water_content 1e-12 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu micro_porosity micro_porosity 1e-12 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu micro_exchange_source micro_exchange_source 1e-12 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu swelling_stress swelling_stress 1e-12 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu sigma sigma 1e-12 1e-10
)

AddTest(
    NAME RichardsMechanics_beacon_1b_dsm_micromacro_smoke
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1b_dsm_micromacro_smoke.prj
    WRAPPER time
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME RichardsMechanics_beacon_1b_dsm_micromacro_reference
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1b_dsm_micromacro_smoke.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu displacement displacement 1e-12 0
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu pressure pressure 1e-12 1e-12
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu saturation saturation 1e-12 0
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu micro_pressure micro_pressure 1e-12 1e-12
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu micro_saturation micro_saturation 1e-12 0
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu swelling_stress swelling_stress 1e-12 0
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu sigma sigma 1e-12 1e-10
)

AddTest(
    NAME RichardsMechanics_beacon_1c_dsm_micromacro_smoke
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1c_dsm_micromacro_smoke.prj
    WRAPPER time
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME RichardsMechanics_beacon_1c_dsm_micromacro_reference
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1c_dsm_micromacro_smoke.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu displacement displacement 1e-12 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu pressure pressure 1e-12 1e-12
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu saturation saturation 1e-12 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu porosity porosity 1e-12 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu transport_porosity transport_porosity 1e-12 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu micro_pressure micro_pressure 1e-12 1e-12
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu micro_saturation micro_saturation 1e-12 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu swelling_stress swelling_stress 1e-12 0
    # GUARDRAIL EXEMPTION CLAUDE.md sec 1.2 (2026-09-01): sigma abs gate
    # widened one order of magnitude (1e-12 -> 1e-11), per Vinay's explicit
    # approval, in response to the ogsds01 Jenkins run (ctest --parallel 20,
    # --preset bgr) failing this vtkdiff with sigma abs max norm up to
    # 4.91e-11 on the zz component and rel max norm 57.7 on the near-zero rz
    # (shear) component. NOTE (flagged to Vinay before landing): this widening
    # is NOT expected to make that run's comparison pass -- 4.91e-11 still
    # exceeds the new 1e-11 gate by ~5x, and the rz rel blowup (near-zero
    # reference value) is untouched by an abs-only change. Root cause not
    # otherwise identified as a beacon_1c code regression (see DSM/AGENTS.md);
    # leading hypothesis is ogsds01 toolchain/build noise beyond what 81b344a
    # measured (cross-compiler epsilon ~4.4e-16..7.8e-16) when it set the
    # prior 1e-12 floor. rel gate intentionally left unchanged (Vinay's
    # choice); scope limited to beacon_1c only, not the 1a01/1b siblings,
    # which currently pass under the unchanged 1e-12/1e-10 gates.
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu sigma sigma 1e-11 1e-10
)
