
# Compile test to check intel compiler version
add_executable(test_intel_compiler_version test_intel_compiler_version.cpp)

add_subdirectory(_temp)
add_subdirectory(containers)
add_subdirectory(inout)
add_subdirectory(green)
add_subdirectory(mpi)
add_subdirectory(tools)

# Add test to check the hydrostatics module
set(
        test_hydrostatics_src
        test_hydrostatics.cpp
        ${containers_module_files}
        ${hydrostatic_module_files}
        ${mesh_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable(
                    test_hydrostatics_cmd
                    ${test_hydrostatics_src}
                )
set_mkl_lib_target( test_hydrostatics_cmd )

set(
        hydrostatics_input_files
        "E:/sergio/developments/SeaMotions/tests/tests_data/hydrostatics/dike_properties.dat"
        "E:/sergio/developments/SeaMotions/tests/tests_data/hydrostatics/dike_mesh.poly"
        "E:/sergio/developments/SeaMotions/tests/tests_data/hydrostatics/semisub_properties.dat"
        "E:/sergio/developments/SeaMotions/tests/tests_data/hydrostatics/semisub_mesh.poly"
    )
add_test(
            NAME test_hydrostatics
            COMMAND test_hydrostatics_cmd
            ${hydrostatics_input_files}
        )

# Add test to check the dispersion relation
set(test_dispersion_src
    test_dispersion_relation.cpp
    ../src/config.hpp
    ../src/math/math_tools.hpp
    ../src/math/math_tools.cpp
    ../src/math/math_tools.txx
    ../src/tools.hpp
    ../src/tools.cpp
    ../src/tools.txx
    ${waves_module_files}
)

add_executable(test_dispersion_cmd ${test_dispersion_src})
set(
    files_path_dispersion
    "${CMAKE_SOURCE_DIR}/tests/tests_data/solution_dispersion_real.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/solution_dispersion_imag.dat"
    )
add_test(NAME test_dispersion
        COMMAND test_dispersion_cmd
        ${files_path_dispersion}
        )