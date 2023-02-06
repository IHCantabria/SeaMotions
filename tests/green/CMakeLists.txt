
# Define common files to compile green function tests
set(green_common_src
    ../../src/math/chebyshev.hpp
    ../../src/math/chebyshev.cpp
    ../../src/config.hpp
    ../../src/containers.hpp
    ../../src/green/common.hpp
    ../../src/green/common.cpp
    ../../src/green/dipole.hpp
    ../../src/green/dipole.cpp
    ../../src/green/pulsating.hpp
    ../../src/green/pulsating.cpp
    ../../src/green/source.hpp
    ../../src/green/source.cpp
    ../../src/math/math_interface.hpp
    ../../src/math/math_tools.hpp
    ../../src/math/math_tools.cpp
    ../../src/math/special_math.hpp
    ../../src/math/special_math.cpp
    ../../src/tools.hpp
    ../../src/tools.cpp
    test_panel_geometries.hpp
    test_panel_geometries.cpp
)

# Add test for dipole potential function
set(test_dipole_potential_src
    test_dipole_potential.cpp
    ${green_common_src}
    )
add_executable(test_dipole_potential ${test_dipole_potential_src})
set_mkl_lib_target(test_dipole_potential)

# Add test for source potential function
set(test_source_potential_src
    test_source_potential.cpp
    ${green_common_src}
)

add_executable(test_source_potential_cmd ${test_source_potential_src})
set_mkl_lib_target(test_source_potential_cmd)
add_test(NAME test_source_potential COMMAND ${test_source_potential_src})


# Add test for source velocity function
set(test_source_velocity_src
    test_source_velocity.cpp
    ${green_common_src}
)

add_executable(test_source_velocity_cmd ${test_source_velocity_src})
set_mkl_lib_target(test_source_velocity_cmd)
add_test(NAME test_source_velocity COMMAND test_source_velocity_cmd)


# Configure source performance test Hess vs Newman
set(performance_source_newman_hess_src 
    performance_hess_vs_newman.cpp
    ${green_common_src}
)

add_executable(performance_source_newman_hess ${performance_source_newman_hess_src})
set_mkl_lib_target(performance_source_newman_hess)


# Add test to check the romberg integration methods with the
# analytical expressions derived in the infinite water depth
# formulation of the pulsating source.
set(test_pulsating_inf_depth_romberg_src
    test_pulsating_inf_depth_romberg.cpp
    ${green_common_src}
    )

add_executable(test_pulsating_inf_depth_romberg_cmd ${test_pulsating_inf_depth_romberg_src})
set_mkl_lib_target(test_pulsating_inf_depth_romberg_cmd)

set(expint_file_path "${CMAKE_SOURCE_DIR}/tests/tests_data/green_inf_depth_tables/table_1_int.dat")
set(expint_dx_file_path "${CMAKE_SOURCE_DIR}/tests/tests_data/green_inf_depth_tables/table_2_int.dat")
set(wave_term_file_path "${CMAKE_SOURCE_DIR}/tests/tests_data/green_inf_depth_tables/table_3_int.dat")
set(wave_term_dx_file_path "${CMAKE_SOURCE_DIR}/tests/tests_data/green_inf_depth_tables/table_4_int.dat")
add_test(NAME test_pulsating_inf_depth_romberg
    COMMAND test_pulsating_inf_depth_romberg_cmd
    ${expint_file_path}
    ${expint_dx_file_path}
    ${wave_term_file_path}
    ${wave_term_dx_file_path}
)


# Add test to check the wave term Newman series
# expansion for infinite water depth
set(test_wave_term_inf_depth_src
    test_wave_term_inf_depth.cpp
    ${green_common_src}
)

add_executable(test_wave_term_inf_depth_cmd ${test_wave_term_inf_depth_src})
set_mkl_lib_target(test_wave_term_inf_depth_cmd)

add_test(NAME test_wave_term_inf_depth
    COMMAND test_wave_term_inf_depth_cmd)


# Add test to check the finite water depth 
# wave term
set(test_wave_term_fin_depth_src
    test_wave_term_fin_depth.cpp
    ../../src/config.hpp
    ../../src/math/math_tools.hpp
    ../../src/math/math_tools.cpp
    ../../src/math/math_tools.txx
    ../../src/green/pulsating_fin_depth.hpp
    ../../src/green/pulsating_fin_depth.cpp
    ../../src/math/chebyshev.hpp
    ../../src/math/chebyshev.cpp
    ../../src/math/special_math.hpp
    ../../src/math/special_math.cpp
    ../../src/tools.hpp
    ../../src/tools.cpp
    ../../src/tools.txx
    ../../src/waves.hpp
    ../../src/waves.cpp
)

add_executable(test_wave_term_fin_depth_cmd ${test_wave_term_fin_depth_src})

set(
    test_wave_term_fin_dpeth_files
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/john_series.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/john_series_dr.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/john_series_dz.dat"
    )
add_test(
        NAME test_wave_term_fin_depth
        COMMAND test_wave_term_fin_depth_cmd
        ${test_wave_term_fin_dpeth_files}
        )

# Add test to check the finite water depth 
# wave term (John series) integration using
# quadrature rules.
set(test_wave_term_john_quad_src
    test_wave_term_john_quad.cpp
    ../../src/config.hpp
    ../../src/math/math_tools.hpp
    ../../src/math/math_tools.cpp
    ../../src/math/math_tools.txx
    ../../src/green/pulsating_fin_depth.hpp
    ../../src/green/pulsating_fin_depth.cpp
    ../../src/math/chebyshev.hpp
    ../../src/math/chebyshev.cpp
    ../../src/math/special_math.hpp
    ../../src/math/gauss.hpp
    ../../src/math/gauss.cpp
    ../../src/math/special_math.cpp
    ../../src/tools.hpp
    ../../src/tools.cpp
    ../../src/tools.txx
    ../../src/waves.hpp
    ../../src/waves.cpp
)

add_executable(
                test_wave_term_john_quad_cmd
                ${test_wave_term_john_quad_src}
                )   