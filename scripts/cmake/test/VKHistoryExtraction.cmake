if(NOT DEFINED PYTHON_EXECUTABLE)
    message(FATAL_ERROR "PYTHON_EXECUTABLE is required.")
endif()
if(NOT DEFINED PROJECT_SOURCE_DIR)
    message(FATAL_ERROR "PROJECT_SOURCE_DIR is required.")
endif()
if(NOT DEFINED BINARY_DIR)
    message(FATAL_ERROR "BINARY_DIR is required.")
endif()
if(NOT DEFINED SOURCE_CSV)
    message(FATAL_ERROR "SOURCE_CSV is required.")
endif()
if(NOT DEFINED OUTPUT_CSV)
    message(FATAL_ERROR "OUTPUT_CSV is required.")
endif()

file(MAKE_DIRECTORY "${BINARY_DIR}")

execute_process(
    COMMAND
        "${PYTHON_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/scripts/vk_extract_observables.py"
        --glob "${BINARY_DIR}/beacon_1a01_vk_notebook_roles_t_*.vtu"
        --output "${OUTPUT_CSV}"
        --arrays
        pressure
        micro_pressure
        micro_saturation
        micro_water_content
        micro_porosity
        micro_exchange_source
        swelling_stress
        sigma
    RESULT_VARIABLE extract_result
)

if(NOT extract_result EQUAL 0)
    message(
        FATAL_ERROR
            "Notebook-role history extraction failed with exit code ${extract_result}."
    )
endif()

if(NOT EXISTS "${OUTPUT_CSV}")
    message(FATAL_ERROR "Notebook-role history CSV was not written.")
endif()

file(SIZE "${OUTPUT_CSV}" output_size)
if(output_size EQUAL 0)
    message(FATAL_ERROR "Notebook-role history CSV is empty.")
endif()
