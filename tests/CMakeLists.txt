
# Compile test to check intel compiler version
add_executable(test_intel_compiler_version test_intel_compiler_version.cpp)

add_subdirectory(_temp)
add_subdirectory(containers)
add_subdirectory(green)
add_subdirectory(mpi)
add_subdirectory(tools)

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
    ../src/waves.hpp
    ../src/waves.cpp
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