if(NOT DEFINED OGS_EXE)
    message(FATAL_ERROR "OGS_EXE is required.")
endif()

if(NOT DEFINED SOURCE_PATH)
    message(FATAL_ERROR "SOURCE_PATH is required.")
endif()

if(NOT DEFINED VTKDIFF_EXE)
    message(FATAL_ERROR "VTKDIFF_EXE is required.")
endif()

if(NOT DEFINED BINARY_PATH)
    message(FATAL_ERROR "BINARY_PATH is required.")
endif()

file(MAKE_DIRECTORY "${BINARY_PATH}")

function(run_ogs project_file)
    execute_process(
        COMMAND "${OGS_EXE}" -o "${BINARY_PATH}" "${project_file}"
        WORKING_DIRECTORY "${SOURCE_PATH}"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE stdout
        ERROR_VARIABLE stderr
    )

    if(NOT result EQUAL 0)
        message(
            FATAL_ERROR
                "ogs failed for '${project_file}' in '${SOURCE_PATH}'.\n"
                "stdout:\n${stdout}\n"
                "stderr:\n${stderr}"
        )
    endif()
endfunction()

function(run_vtkdiff file_a file_b array_a array_b abs_tol rel_tol)
    execute_process(
        COMMAND
            "${VTKDIFF_EXE}"
            -a "${array_a}"
            -b "${array_b}"
            --abs "${abs_tol}"
            --rel "${rel_tol}"
            "${file_a}"
            "${file_b}"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE stdout
        ERROR_VARIABLE stderr
    )

    if(NOT result EQUAL 0)
        message(
            FATAL_ERROR
                "vtkdiff failed for '${array_a}' comparing\n"
                "  ${file_a}\n"
                "  ${file_b}\n"
                "stdout:\n${stdout}\n"
                "stderr:\n${stderr}"
        )
    endif()
endfunction()

run_ogs(mfront_restart_part1_notebook_mcc_tuller_native.prj)
run_ogs(mfront_restart_part1_notebook_mcc_tuller_bridge.prj)

set(step_indices 0 1000)
set(step_times 0.000000 1000.000000)

foreach(list_index RANGE 0 1)
    list(GET step_indices ${list_index} ts)
    list(GET step_times ${list_index} time)
    set(
        native_file
        "${BINARY_PATH}/mfront_restart_part1_notebook_mcc_tuller_native_ts_${ts}_t_${time}.vtu"
    )
    set(
        bridge_file
        "${BINARY_PATH}/mfront_restart_part1_notebook_mcc_tuller_bridge_ts_${ts}_t_${time}.vtu"
    )

    run_vtkdiff("${native_file}" "${bridge_file}" displacement displacement 1e-12 0)
    run_vtkdiff("${native_file}" "${bridge_file}" pressure pressure 1e-6 0)
    run_vtkdiff("${native_file}" "${bridge_file}" sigma sigma 1e-9 0)
    run_vtkdiff("${native_file}" "${bridge_file}" epsilon epsilon 1e-12 0)
    run_vtkdiff("${native_file}" "${bridge_file}" saturation saturation 1e-12 0)
    run_vtkdiff("${native_file}" "${bridge_file}" velocity velocity 1e-15 0)
    run_vtkdiff("${native_file}" "${bridge_file}" ElasticStrain ElasticStrain 1e-15 0)
    run_vtkdiff(
        "${native_file}" "${bridge_file}" EquivalentPlasticStrain
        EquivalentPlasticStrain 1e-15 0
    )
    run_vtkdiff(
        "${native_file}" "${bridge_file}" PreConsolidationPressure
        PreConsolidationPressure 1e-15 0
    )
    run_vtkdiff(
        "${native_file}" "${bridge_file}" PlasticVolumetricStrain
        PlasticVolumetricStrain 1e-12 0
    )
    run_vtkdiff("${native_file}" "${bridge_file}" VolumeRatio VolumeRatio 1e-15 0)
    run_vtkdiff("${native_file}" "${bridge_file}" swelling_stress swelling_stress 1e-15 0)
    run_vtkdiff(
        "${native_file}" "${bridge_file}" transport_porosity
        transport_porosity 1e-15 0
    )
    run_vtkdiff(
        "${native_file}" "${bridge_file}" dry_density_solid dry_density_solid
        1e-9 0
    )
endforeach()
