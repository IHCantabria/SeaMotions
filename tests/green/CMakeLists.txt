
# Define file chunks to compile green functions
set(green_common_inf_depth_src
    ../../src/math/chebyshev.hpp
    ../../src/math/chebyshev.cpp
    ../../src/config.hpp
    ../../src/containers.hpp
    ../../src/green/common.hpp
    ../../src/green/common.cpp
    ../../src/green/dipole.hpp
    ../../src/green/dipole.cpp
    ../../src/green/pulsating_inf_depth.hpp
    ../../src/green/pulsating_inf_depth.cpp
    ../../src/green/pulsating_inf_depth_series.hpp
    ../../src/green/pulsating_inf_depth_series.cpp
    ../../src/green/pulsating_inf_depth_cheby.hpp
    ../../src/green/pulsating_inf_depth_cheby.cpp
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

set(
    green_common_fin_depth_src
    ../../src/green/integrals_db.hpp
    ../../src/green/integrals_db.cpp
    ../../src/green/pulsating_fin_depth.hpp
    ../../src/green/pulsating_fin_depth.cpp
    ../../src/green/pulsating_fin_depth_cheby.hpp
    ../../src/green/pulsating_fin_depth_cheby.cpp
    ../../src/waves.hpp
    ../../src/waves.cpp
    ${green_common_inf_depth_src}
)

set(
    inf_depth_integrals_src
    ../../src/green/inf_depth_coeffs/R11.hpp
    ../../src/green/inf_depth_coeffs/R11.cpp
    ../../src/green/inf_depth_coeffs/R11A_dX.hpp
    ../../src/green/inf_depth_coeffs/R11A_dX.cpp
    ../../src/green/inf_depth_coeffs/R11B_dX.hpp
    ../../src/green/inf_depth_coeffs/R11B_dX.cpp
    ../../src/green/inf_depth_coeffs/R12.hpp
    ../../src/green/inf_depth_coeffs/R12.cpp
    ../../src/green/inf_depth_coeffs/R12_dX.hpp
    ../../src/green/inf_depth_coeffs/R12_dX.cpp
    ../../src/green/inf_depth_coeffs/R21.hpp
    ../../src/green/inf_depth_coeffs/R21.cpp
    ../../src/green/inf_depth_coeffs/R21_dX.hpp
    ../../src/green/inf_depth_coeffs/R21_dX.cpp
    ../../src/green/inf_depth_coeffs/R22.hpp
    ../../src/green/inf_depth_coeffs/R22.cpp
    ../../src/green/inf_depth_coeffs/R22_dX.hpp
    ../../src/green/inf_depth_coeffs/R22_dX.cpp
)

set(
    fin_depth_integrals_src
    ../../src/green/fin_depth_coeffs/L1.hpp
    ../../src/green/fin_depth_coeffs/L1.cpp
    ../../src/green/fin_depth_coeffs/L1_dA.hpp
    ../../src/green/fin_depth_coeffs/L1_dA.cpp
    ../../src/green/fin_depth_coeffs/L1_dB.hpp
    ../../src/green/fin_depth_coeffs/L1_dB.cpp
    ../../src/green/fin_depth_coeffs/L2.hpp
    ../../src/green/fin_depth_coeffs/L2.cpp
    ../../src/green/fin_depth_coeffs/L3.hpp
    ../../src/green/fin_depth_coeffs/L3.cpp
    ../../src/green/fin_depth_coeffs/L3_dA.hpp
    ../../src/green/fin_depth_coeffs/L3_dA.cpp
    ../../src/green/fin_depth_coeffs/L3_dB.hpp
    ../../src/green/fin_depth_coeffs/L3_dB.cpp
    ../../src/green/fin_depth_coeffs/M1.hpp
    ../../src/green/fin_depth_coeffs/M1.cpp
    ../../src/green/fin_depth_coeffs/M1_dA.hpp
    ../../src/green/fin_depth_coeffs/M1_dA.cpp
    ../../src/green/fin_depth_coeffs/M1_dB.hpp
    ../../src/green/fin_depth_coeffs/M1_dB.cpp
    ../../src/green/fin_depth_coeffs/M2.hpp
    ../../src/green/fin_depth_coeffs/M2.cpp
    ../../src/green/fin_depth_coeffs/M3.hpp
    ../../src/green/fin_depth_coeffs/M3.cpp
    ../../src/green/fin_depth_coeffs/M3_dA.hpp
    ../../src/green/fin_depth_coeffs/M3_dA.cpp
    ../../src/green/fin_depth_coeffs/M3_dB.hpp
    ../../src/green/fin_depth_coeffs/M3_dB.cpp
    )

# Add test for dipole potential function
set(test_dipole_potential_src
    test_dipole_potential.cpp
    ${green_common_inf_depth_src}
    )
add_executable(test_dipole_potential ${test_dipole_potential_src})
set_mkl_lib_target(test_dipole_potential)

# Add test for source potential function
set(test_source_potential_src
    test_source_potential.cpp
    ${green_common_inf_depth_src}
)

add_executable(test_source_potential_cmd ${test_source_potential_src})
set_mkl_lib_target(test_source_potential_cmd)
add_test(NAME test_source_potential COMMAND ${test_source_potential_src})


# Add test for source velocity function
set(test_source_velocity_src
    test_source_velocity.cpp
    ${green_common_inf_depth_src}
)

add_executable(test_source_velocity_cmd ${test_source_velocity_src})
set_mkl_lib_target(test_source_velocity_cmd)
add_test(NAME test_source_velocity COMMAND test_source_velocity_cmd)


# Configure source performance test Hess vs Newman
set(performance_source_newman_hess_src 
    performance_hess_vs_newman.cpp
    ${green_common_inf_depth_src}
)

add_executable(performance_source_newman_hess ${performance_source_newman_hess_src})
set_mkl_lib_target(performance_source_newman_hess)


# Add test to check the romberg integration methods with the
# analytical expressions derived in the infinite water depth
# formulation of the pulsating source.
set(test_pulsating_inf_depth_romberg_src
    test_pulsating_inf_depth_romberg.cpp
    ${green_common_inf_depth_src}
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
# expansion for infinite water depth wave term
set(test_wave_term_inf_depth_series_src
    test_wave_term_inf_depth_series.cpp
    ${green_common_inf_depth_src}
    ${inf_depth_integrals_src}
)

add_executable(test_wave_term_inf_depth_series_cmd ${test_wave_term_inf_depth_series_src})
set_mkl_lib_target(test_wave_term_inf_depth_series_cmd)

add_test(NAME test_wave_term_inf_depth_series
    COMMAND test_wave_term_inf_depth_series_cmd)


# Add test to check the chebyshev polynomials 
# evaluation of the infinite water depth wave term
set(test_wave_term_inf_depth_src
    test_wave_term_inf_depth.cpp
    ${green_common_fin_depth_src}
    ${inf_depth_integrals_src}
    ${fin_depth_integrals_src}
)

add_executable(test_wave_term_inf_depth_cmd ${test_wave_term_inf_depth_src})
set_mkl_lib_target(test_wave_term_inf_depth_cmd)

add_test(NAME test_wave_term_inf_depth
    COMMAND test_wave_term_inf_depth_cmd)


# Add test to check the finite water depth 
# wave term
set(test_wave_term_fin_depth_src
    test_wave_term_fin_depth.cpp
    ${green_common_fin_depth_src}
    ${fin_depth_integrals_src}
    ${inf_depth_integrals_src}
)

add_executable(test_wave_term_fin_depth_cmd ${test_wave_term_fin_depth_src})

set(
    test_wave_term_fin_depth_files
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/john_series.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/john_series_dr.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/john_series_dz.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/fin_depth_integral.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/fin_depth_integral_dr.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/fin_depth_integral_dz.dat"
    )
add_test(
        NAME test_wave_term_fin_depth
        COMMAND test_wave_term_fin_depth_cmd
        ${test_wave_term_fin_depth_files}
        )

# Add test to check the finite water depth 
# wave term (John series) integration using
# quadrature rules.
set(test_wave_term_john_quad_src
    test_wave_term_john_quad.cpp
    ${green_common_fin_depth_src}
    ${fin_depth_integrals_src}
)

add_executable(
                test_wave_term_john_quad_cmd
                ${test_wave_term_john_quad_src}
                )


# Configure test to check the value of the wave term integrals
# separately
set(
    test_wave_term_fin_depth_integrals_src
    test_wave_term_fin_depth_integrals.cpp
    ${green_common_fin_depth_src}
    ${fin_depth_integrals_src}
)

add_executable(
                test_wave_term_fin_depth_integrals_cmd
                ${test_wave_term_fin_depth_integrals_src}
)

set(
    test_wave_term_fin_depth_integrals_args
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/L1.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/L1_dA.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/L1_dB.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/L2.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/L3.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/L3_dA.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/L3_dB.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/M1.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/M1_dA.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/M1_dB.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/M2.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/M3.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/M3_dA.dat"
    "${CMAKE_SOURCE_DIR}/tests/tests_data/green_fin_depth_tables/M3_dB.dat"
)
add_test(
        NAME test_wave_term_fin_depth_integrals
        COMMAND test_wave_term_fin_depth_integrals_cmd
        ${test_wave_term_fin_depth_integrals_args}
)