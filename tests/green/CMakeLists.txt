set( 
        panel_geometry_src
        test_panel_geometries.hpp
        test_panel_geometries.cpp
    )

# Add test for dipole potential function
set(
        test_dipole_potential_src
        test_dipole_potential.cpp
        ${containers_module_files}
        ${panel_geometry_src}
        ${green_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable( 
                    test_dipole_potential 
                    ${test_dipole_potential_src}
                )
set_mkl_lib_target( test_dipole_potential )

# Add test for source potential function
set(
        test_source_potential_src
        test_source_potential.cpp
        ${containers_module_files}
        ${panel_geometry_src}
        ${green_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable(
                    test_source_potential_cmd 
                    ${test_source_potential_src}
                )
set_mkl_lib_target( test_source_potential_cmd )
add_test(
            NAME test_source_potential 
            COMMAND test_source_potential_cmd
        )


# Add test for source velocity function
set(    
        test_source_velocity_src
        test_source_velocity.cpp
        ${containers_module_files}
        ${panel_geometry_src}
        ${green_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable(
                    test_source_velocity_cmd 
                    ${test_source_velocity_src}
                )
set_mkl_lib_target( test_source_velocity_cmd )
add_test(
            NAME test_source_velocity 
            COMMAND test_source_velocity_cmd
        )


# Configure source performance test Hess vs Newman
set(    
        performance_source_newman_hess_src  
        performance_hess_vs_newman.cpp
        ${containers_module_files}
        ${panel_geometry_src}
        ${green_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable(
                    performance_source_newman_hess 
                    ${performance_source_newman_hess_src}
                )
set_mkl_lib_target(performance_source_newman_hess)


# Add test to check the romberg integration methods with the
# analytical expressions derived in the infinite water depth
# formulation of the pulsating source.
set(
        test_pulsating_inf_depth_romberg_src
        test_pulsating_inf_depth_romberg.cpp
        ${containers_module_files}
        ${panel_geometry_src}
        ${green_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable(
                    test_pulsating_inf_depth_romberg_cmd 
                    ${test_pulsating_inf_depth_romberg_src}
                )
set_mkl_lib_target(test_pulsating_inf_depth_romberg_cmd)

set(expint_file_path "${CMAKE_SOURCE_DIR}/tests/tests_data/green_inf_depth_tables/table_1_int.dat")
set(expint_dx_file_path "${CMAKE_SOURCE_DIR}/tests/tests_data/green_inf_depth_tables/table_2_int.dat")
set(wave_term_file_path "${CMAKE_SOURCE_DIR}/tests/tests_data/green_inf_depth_tables/table_3_int.dat")
set(wave_term_dx_file_path "${CMAKE_SOURCE_DIR}/tests/tests_data/green_inf_depth_tables/table_4_int.dat")
add_test(
            NAME test_pulsating_inf_depth_romberg
            COMMAND test_pulsating_inf_depth_romberg_cmd
            ${expint_file_path}
            ${expint_dx_file_path}
            ${wave_term_file_path}
            ${wave_term_dx_file_path}
        )


# Add test to check the wave term Newman series
# expansion for infinite water depth wave term
set(
        test_wave_term_inf_depth_series_src
        test_wave_term_inf_depth_series.cpp
        ${containers_module_files}
        ${panel_geometry_src}
        ${green_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable(
                    test_wave_term_inf_depth_series_cmd 
                    ${test_wave_term_inf_depth_series_src}
                )
set_mkl_lib_target(test_wave_term_inf_depth_series_cmd)

add_test(
            NAME test_wave_term_inf_depth_series
            COMMAND test_wave_term_inf_depth_series_cmd
        )


# Add test to check the chebyshev polynomials 
# evaluation of the infinite water depth wave term
set(
        test_wave_term_inf_depth_src
        test_wave_term_inf_depth.cpp
        ${containers_module_files}
        ${panel_geometry_src}
        ${green_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable(
                    test_wave_term_inf_depth_cmd
                    ${test_wave_term_inf_depth_src}
                )
set_mkl_lib_target(test_wave_term_inf_depth_cmd)

add_test(
            NAME test_wave_term_inf_depth
            COMMAND test_wave_term_inf_depth_cmd
        )


# Add test to check the finite water depth 
# wave term
set(
        test_wave_term_fin_depth_src
        test_wave_term_fin_depth.cpp
        ${containers_module_files}
        ${panel_geometry_src}
        ${green_module_files}
        ${math_module_files}
        ${tools_module_files}
        ${waves_module_files}
    )

add_executable(
                    test_wave_term_fin_depth_cmd 
                    ${test_wave_term_fin_depth_src}
                )
set_mkl_lib_target( test_wave_term_fin_depth_cmd )

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
set(
        test_quad_green_element_src
        test_quad_green_element.cpp
        ${containers_module_files}
        ${panel_geometry_src}
        ${green_module_files}
        ${math_module_files}
        ${mesh_module_files}
        ${tools_module_files}
        ${waves_module_files}
    )

add_executable(
                    test_quad_green_element_cmd
                    ${test_quad_green_element_src}
                )
set_mkl_lib_target( test_quad_green_element_cmd )


# Configure test to check the value of the wave term integrals
# separately
set(
        test_wave_term_fin_depth_integrals_src
        test_wave_term_fin_depth_integrals.cpp
        ${containers_module_files}
        ${panel_geometry_src}
        ${green_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable(
                    test_wave_term_fin_depth_integrals_cmd
                    ${test_wave_term_fin_depth_integrals_src}
                )
set_mkl_lib_target( test_wave_term_fin_depth_integrals_cmd )

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


# Add test to check the evaluation of the chebyshev polynomials with the input
# format of the integral databases
set(
        test_chebyshev_integraldb_src
        test_chebyshev_integraldb.cpp
        ${containers_module_files}
        ${green_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable( test_chebyshev_integraldb_cmd ${test_chebyshev_integraldb_src} )
set_mkl_lib_target( test_chebyshev_integraldb_cmd )

set(
        test_chebyshev_dbeval_args
        "E:/sergio/developments/SeaMotions/tests/tests_data/chebyshev_integral_db.dat"
    )
add_test(
            NAME    test_chebyshev_integraldb
            COMMAND test_chebyshev_integraldb_cmd
            ${test_chebyshev_dbeval_args}
        )