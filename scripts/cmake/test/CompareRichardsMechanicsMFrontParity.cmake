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

run_ogs(mfront_parity_1element_native.prj)
run_ogs(mfront_parity_1element_bridge.prj)

foreach(ts RANGE 0 4)
    set(
        native_file
        "${BINARY_PATH}/mfront_parity_1element_native_ts_${ts}_t_${ts}.000000.vtu"
    )
    set(
        bridge_file
        "${BINARY_PATH}/mfront_parity_1element_bridge_ts_${ts}_t_${ts}.000000.vtu"
    )

    run_vtkdiff("${native_file}" "${bridge_file}" displacement displacement 1e-15 0)
    run_vtkdiff("${native_file}" "${bridge_file}" pressure pressure 1e-15 0)
    # sigma still differs at machine precision because the native and bridge
    # paths do not round identically.
    run_vtkdiff("${native_file}" "${bridge_file}" sigma sigma 5e-12 0)
    run_vtkdiff("${native_file}" "${bridge_file}" epsilon epsilon 1e-15 0)
    run_vtkdiff("${native_file}" "${bridge_file}" saturation saturation 1e-15 0)
    run_vtkdiff("${native_file}" "${bridge_file}" swelling_stress swelling_stress 1e-15 0)
endforeach()
