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
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1400.prj RUNTIME 120)
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj RUNTIME 120)
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1800.prj RUNTIME 120)
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj RUNTIME 240)
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj RUNTIME 240)
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
    NAME RichardsMechanics_double_porosity_swelling_dsm_micromacro_constbc_reference
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS double_porosity_swelling_dsm_micromacro_constbc.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    dsm_micromacro_constbc_reference_t_1000.000000.vtu double_porosity_swelling_dsm_micromacro_constbc_t_1000.000000.vtu pressure pressure 1e-16 1e-12
    dsm_micromacro_constbc_reference_t_1000.000000.vtu double_porosity_swelling_dsm_micromacro_constbc_t_1000.000000.vtu saturation saturation 4e-15 0
    dsm_micromacro_constbc_reference_t_1000.000000.vtu double_porosity_swelling_dsm_micromacro_constbc_t_1000.000000.vtu porosity porosity 1e-16 0
    dsm_micromacro_constbc_reference_t_1000.000000.vtu double_porosity_swelling_dsm_micromacro_constbc_t_1000.000000.vtu transport_porosity transport_porosity 1e-16 0
    dsm_micromacro_constbc_reference_t_1000.000000.vtu double_porosity_swelling_dsm_micromacro_constbc_t_1000.000000.vtu micro_pressure micro_pressure 1e-16 1e-12
    dsm_micromacro_constbc_reference_t_1000.000000.vtu double_porosity_swelling_dsm_micromacro_constbc_t_1000.000000.vtu micro_saturation micro_saturation 4e-15 0
    dsm_micromacro_constbc_reference_t_1000.000000.vtu double_porosity_swelling_dsm_micromacro_constbc_t_1000.000000.vtu swelling_stress swelling_stress 5e-14 0
)

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
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu displacement displacement 1e-16 0
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu pressure pressure 1e-16 1e-12
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu saturation saturation 1e-14 0
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu micro_pressure micro_pressure 1e-16 1e-12
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu micro_saturation micro_saturation 1e-14 0
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu swelling_stress swelling_stress 1e-16 0
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu sigma sigma 1e-16 1e-10
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
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu displacement displacement 1e-16 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu pressure pressure 1e-16 1e-12
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu saturation saturation 1e-14 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu micro_pressure micro_pressure 1e-16 1e-12
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu micro_saturation micro_saturation 1e-14 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu micro_water_content micro_water_content 1e-16 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu micro_porosity micro_porosity 1e-16 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu micro_exchange_source micro_exchange_source 1e-16 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu swelling_stress swelling_stress 1e-16 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu sigma sigma 1e-16 1e-10
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
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu displacement displacement 1e-16 0
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu pressure pressure 1e-16 1e-12
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu saturation saturation 1e-16 0
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu micro_pressure micro_pressure 1e-16 1e-12
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu micro_saturation micro_saturation 1e-16 0
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu swelling_stress swelling_stress 1e-16 0
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu sigma sigma 1e-16 1e-10
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
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu displacement displacement 1e-16 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu pressure pressure 1e-16 1e-12
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu saturation saturation 1e-14 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu porosity porosity 1e-16 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu transport_porosity transport_porosity 1e-16 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu micro_pressure micro_pressure 1e-16 1e-12
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu micro_saturation micro_saturation 1e-14 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu swelling_stress swelling_stress 1e-16 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu sigma sigma 1e-16 1e-10
)
