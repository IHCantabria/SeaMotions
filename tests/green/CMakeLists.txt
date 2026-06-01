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


# Add test to compare Newman integration vs Gauss rule over quadrilateral
set(
        test_gauss_source_quad_src
        test_gauss_source_quadrilateral.cpp
        "../../src/config.hpp"
        ${containers_module_files}
        ${green_module_files}
        ${interfaces_module_files}
        ${hydrostatic_module_files}
        ${inout_module_files}
        ${math_module_files}
        ${mesh_module_files}
        ${raos_module_files}
        ${tools_module_files}
        ${waves_module_files}
        ${version_module_files}
    )

add_executable(
                    test_gauss_source_quadrilateral_cmd
                    ${test_gauss_source_quad_src}
                )
set_mkl_lib_target( test_gauss_source_quadrilateral_cmd )
set_mpi_target( test_gauss_source_quadrilateral_cmd )
set_hdf5_target( test_gauss_source_quadrilateral_cmd )
add_test(
            NAME test_gauss_source_quadrilateral
            COMMAND test_gauss_source_quadrilateral_cmd
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
        ${interfaces_module_files}
        ${math_module_files}
        ${mesh_module_files}
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
        "${CMAKE_SOURCE_DIR}/tests/tests_data/green_depth_tables/G_Agt1.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/green_depth_tables/G_Agt1_dR.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/green_depth_tables/G_Agt1_dz.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/green_depth_tables/G_Agt1_dzeta.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/green_depth_tables/G_Alt1.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/green_depth_tables/G_Alt1_dR.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/green_depth_tables/G_Alt1_dz.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/green_depth_tables/G_Alt1_dzeta.dat"
    )

add_test(
            NAME test_wave_term_fin_depth
            COMMAND test_wave_term_fin_depth_cmd
            ${test_wave_term_fin_depth_files}
        )



# Add test to compare full finite-depth wave integration against split wave plus exact steady terms
set(
        test_wave_term_panel_full_vs_split_src
        test_wave_term_panel_full_vs_split.cpp
    ../../src/containers/panel_geom.cpp
    ../../src/green/common.cpp
    ../../src/green/dipole.cpp
    ../../src/green/source.cpp
    ../../src/math/chebyshev.cpp
    ../../src/math/gauss.cpp
    ../../src/math/integration.cpp
    ../../src/math/math_tools.cpp
    ../../src/math/shape_functions.cpp
    ../../src/math/special_math.cpp
    ../../src/math/topology.cpp
    ../../src/tools.cpp
    ../../src/waves/wave_dispersion_base_fo.cpp
    ../../src/waves/wave_dispersion_so.cpp
    ../../src/waves/waves_common.cpp
    )

add_executable(
                    test_wave_term_panel_full_vs_split_cmd
                    ${test_wave_term_panel_full_vs_split_src}
                )
set_mkl_lib_target( test_wave_term_panel_full_vs_split_cmd )
add_test(
            NAME test_wave_term_panel_full_vs_split
            COMMAND test_wave_term_panel_full_vs_split_cmd
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
        ${interfaces_module_files}
        ${panel_geometry_src}
        ${green_module_files}
        ${math_module_files}
        ${mesh_module_files}
        ${tools_module_files}
        ${waves_module_files}
    )

add_executable(
                    test_wave_term_fin_depth_integrals_cmd
                    ${test_wave_term_fin_depth_integrals_src}
                )
set_mkl_lib_target( test_wave_term_fin_depth_integrals_cmd )

set(
        test_wave_term_fin_depth_integrals_args
        "${CMAKE_SOURCE_DIR}/tests/tests_data/green_depth_tables/L1.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/green_depth_tables/L2.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/green_depth_tables/L3.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/green_depth_tables/M1.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/green_depth_tables/M2.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/green_depth_tables/M3.dat"
    )

add_test(
            NAME    test_wave_term_fin_depth_integrals
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


# Add test to check the integration of the potential along the panel
set(
        test_potential_panel_quad_src
        potential_panel_quad_src.cpp
        ${containers_module_files}
        ${interfaces_module_files}
        ${green_module_files}
        ${math_module_files}
        ${mesh_module_files}
        ${tools_module_files}
        ${waves_module_files}
)

add_executable(
                    test_potential_panel_quad_cmd
                    ${test_potential_panel_quad_src}
                )
set_mkl_lib_target( test_potential_panel_quad_cmd )


# Add test to check the integrals database performance
set(
        test_integrals_performance_src
        test_integrals_performance.cpp
        ${containers_module_files}
        ${green_module_files}
        ${math_module_files}
        ${mesh_module_files}
        ${tools_module_files}
)

add_executable(
                    test_integrals_performance_cmd
                    ${test_integrals_performance_src}
                )
set_mkl_lib_target( test_integrals_performance_cmd )


set(
        test_steady_analytic_panel_src
        test_steady_analytic_panels.cpp
        ${containers_module_files}
        ${green_module_files}
        ${math_module_files}
        ${mesh_module_files}
        ${tools_module_files}
)

add_executable(
                    test_steady_analytic_panels_cmd
                    ${test_steady_analytic_panel_src}
                )
set_mkl_lib_target( test_steady_analytic_panels_cmd )


# Define test to check the integration of the velocity
# field along a panel
set(
        test_velocity_integration_src
        test_velocity_integration.cpp
        ../../src/config.hpp
        ${containers_module_files}
        ${green_module_files}
        ${interfaces_module_files}
        ${math_module_files}
        ${mesh_module_files}
        ${tools_module_files}
        ${waves_module_files}
)


add_executable(
                    test_velocity_integration_cmd
                    ${test_velocity_integration_src}
                )
set_mkl_lib_target( test_velocity_integration_cmd )

# Define Green R=0 test
set(
        test_green_r_zero_src
        test_green_r_zero.cpp
        ../../src/config.hpp
        ${containers_module_files}
        ${green_module_files}
        ${math_module_files}
        ${mesh_module_files}
        ${tools_module_files}
        ${waves_module_files}
)

add_executable(
                    test_green_r_zero_cmd
                    ${test_green_r_zero_src}
                )
set_mkl_lib_target( test_green_r_zero_cmd )


# Add test for source potential function
set(
        test_multipole_expansions_src
        test_multipole_expansions.cpp
        ${containers_module_files}
        ${inout_module_files}
        ${interfaces_module_files}
        ${panel_geometry_src}
        ${green_module_files}
        ${math_module_files}
        ${mesh_module_files}
        ${tools_module_files}
        ${version_module_files}
        ${waves_module_files}
    )

add_executable(
                    test_multipole_expansions_cmd 
                    ${test_multipole_expansions_src}
                )
set_mkl_lib_target( test_multipole_expansions_cmd )
set_mpi_target( test_multipole_expansions_cmd )
set_hdf5_target( test_multipole_expansions_cmd )
add_test(
            NAME test_multipole_expansions 
            COMMAND test_multipole_expansions_cmd
        )


# Add test to validate rotated-contour triple Hankel integration against
# the reference database generated by aux_tools/generate_triple_hankel_reference_database.py.
set(
        test_triple_hankel_mod_database_src
        test_triple_hankel_mod_database.cpp
        "../../src/config.hpp"
        "../../src/green/triple_hankel_far_field.hpp"
        "../../src/green/triple_hankel_far_field.cpp"
        ${math_module_files}
        ${tools_module_files}
    )

add_executable(
                    test_triple_hankel_mod_database_cmd
                    ${test_triple_hankel_mod_database_src}
                )
set_hdf5_target( test_triple_hankel_mod_database_cmd )
set_mkl_lib_target( test_triple_hankel_mod_database_cmd )


set(
        test_triple_hankel_mod_database_file
        "${CMAKE_SOURCE_DIR}/aux_data/2_check_database_files/0_frequency_domain/triple_hankel_reference.h5"
    )

add_test(
            NAME test_triple_hankel_mod_database
            COMMAND test_triple_hankel_mod_database_cmd
            ${test_triple_hankel_mod_database_file}
        )
set_tests_properties(
                        test_triple_hankel_mod_database
                        PROPERTIES FIXTURES_REQUIRED triple_hankel_db
                    )


# Add test to verify the 1D, 2D, and 3D Chebyshev database folding algorithms.
# No external data files needed – the test uses the built-in coefficient headers
# and a small synthetic polynomial for the 2D case.
set(
        test_fold_database_src
        test_fold_database.cpp
        "../../src/config.hpp"
        ${containers_module_files}
        ${inout_module_files}
        ${interfaces_module_files}
        ${panel_geometry_src}
        ${green_module_files}
        ${math_module_files}
        ${mesh_module_files}
        ${raos_module_files}
        ${tools_module_files}
        ${version_module_files}
        ${waves_module_files}
    )

add_executable(
                    test_fold_database_cmd
                    ${test_fold_database_src}
                )
set_mkl_lib_target( test_fold_database_cmd )
set_hdf5_target( test_fold_database_cmd )
add_test(
            NAME    test_fold_database
            COMMAND test_fold_database_cmd
        )


# Add test for the time-domain Green's function evaluators
set(
        test_time_domain_evaluator_src
        test_time_domain_evaluator.cpp
        "../../src/config.hpp"
        ${containers_module_files}
        ${inout_module_files}
        ${interfaces_module_files}
        ${panel_geometry_src}
        ${green_module_files}
        ${math_module_files}
        ${mesh_module_files}
        ${raos_module_files}
        ${tools_module_files}
        ${version_module_files}
        ${waves_module_files}
    )

add_executable(
                    test_time_domain_evaluator_cmd
                    ${test_time_domain_evaluator_src}
                )
set_mkl_lib_target( test_time_domain_evaluator_cmd )
set_hdf5_target( test_time_domain_evaluator_cmd )

set(
        test_time_domain_evaluator_args
        # Total function reference  (argv[1..7])
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdt.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtx.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtxx.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtt.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttx.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttxx.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttt.dat"
        # G0-only reference  (argv[8..14])
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdt_G0.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtx_G0.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtxx_G0.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtt_G0.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttx_G0.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttxx_G0.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttt_G0.dat"
        # Residual-only reference  (argv[15..21])
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdt_residual.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtx_residual.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtxx_residual.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtt_residual.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttx_residual.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttxx_residual.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttt_residual.dat"
    )

add_test(
            NAME    test_time_domain_evaluator
            COMMAND test_time_domain_evaluator_cmd
            ${test_time_domain_evaluator_args}
        )


# Add test for the Gauss-point-optimised time-domain evaluators
set(
        test_time_domain_evaluator_gauss_src
        test_time_domain_evaluator_gauss.cpp
        "../../src/config.hpp"
        ${containers_module_files}
        ${inout_module_files}
        ${interfaces_module_files}
        ${panel_geometry_src}
        ${green_module_files}
        ${math_module_files}
        ${mesh_module_files}
        ${raos_module_files}
        ${tools_module_files}
        ${version_module_files}
        ${waves_module_files}
    )

add_executable(
                    test_time_domain_evaluator_gauss_cmd
                    ${test_time_domain_evaluator_gauss_src}
                )
set_mkl_lib_target( test_time_domain_evaluator_gauss_cmd )
set_hdf5_target( test_time_domain_evaluator_gauss_cmd )

# Reuse the same reference data files as test_time_domain_evaluator
set(
        test_time_domain_evaluator_gauss_args
        # Total function reference  (argv[1..7])
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdt.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtx.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtxx.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtt.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttx.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttxx.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttt.dat"
        # G0-only reference  (argv[8..14])
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdt_G0.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtx_G0.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtxx_G0.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtt_G0.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttx_G0.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttxx_G0.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttt_G0.dat"
        # Residual-only reference  (argv[15..21])
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdt_residual.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtx_residual.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtxx_residual.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtt_residual.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttx_residual.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttxx_residual.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttt_residual.dat"
    )

add_test(
            NAME    test_time_domain_evaluator_gauss
            COMMAND test_time_domain_evaluator_gauss_cmd
            ${test_time_domain_evaluator_gauss_args}
        )

# ============================================================================
# Test suite for the 2D interpolators (Bilinear, Catmull-Rom, NCS) applied
# to all seven time-domain integral tables: Gt, Gtt, Gttt, Gtx, Gttx, Gtxx, Gttxx.
# ============================================================================
set( td_data_dir "${CMAKE_SOURCE_DIR}/src/green/td_data" )
set(
    test_interp2d_gt_src
    test_interp2d_gt.cpp
    "../../src/config.hpp"
    "../../src/math/data_fit/interp2d.hpp"
    "${td_data_dir}/td_database.hpp"
    "${td_data_dir}/td_common_axes.cpp"
    "${td_data_dir}/td_Gt_fcn.cpp"
    "${td_data_dir}/td_Gtt_fcn.cpp"
    "${td_data_dir}/td_Gttt_fcn.cpp"
    "${td_data_dir}/td_Gtx_fcn.cpp"
    "${td_data_dir}/td_Gttx_fcn.cpp"
    "${td_data_dir}/td_Gtxx_fcn.cpp"
    "${td_data_dir}/td_Gttxx_fcn.cpp"
)
add_executable( test_interp2d_gt_cmd ${test_interp2d_gt_src} )
set_hdf5_target( test_interp2d_gt_cmd )

set( td_db_dir "${CMAKE_SOURCE_DIR}/aux_data/0_integrals_database/1_time_domain" )
set( td_ref_dir "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain" )
set(
    test_interp2d_gt_args
    "${td_db_dir}/Gt.h5"
    "${td_db_dir}/Gtt.h5"
    "${td_db_dir}/Gttt.h5"
    "${td_db_dir}/Gtx.h5"
    "${td_db_dir}/Gttx.h5"
    "${td_db_dir}/Gtxx.h5"
    "${td_db_dir}/Gttxx.h5"
    "${td_ref_dir}/dGdt.dat"
    "${td_ref_dir}/dGdtt.dat"
    "${td_ref_dir}/dGdttt.dat"
    "${td_ref_dir}/dGdtx.dat"
    "${td_ref_dir}/dGdttx.dat"
    "${td_ref_dir}/dGdtxx.dat"
    "${td_ref_dir}/dGdttxx.dat"
)
add_test(
    NAME    test_interp2d_gt
    COMMAND test_interp2d_gt_cmd ${test_interp2d_gt_args}
)


# ============================================================================
# Test for the time-domain asymptotic Green function routines
# (src/green/time_domain_asymptotic.hpp / .cpp).
#
# Test 1 – spot_checks: C++ output vs Python float64 reference at a 4×3 grid.
# Test 2 – finiteness: all 600 (beta, mu) points in each .dat file produce a
#           finite, non-NaN result.
# ============================================================================
set(
        test_time_domain_asymptotic_src
        test_time_domain_asymptotic.cpp
        "../../src/config.hpp"
        "../../src/math/math_constants.hpp"
        "../../src/green/time_domain_asymptotic.hpp"
        "../../src/green/time_domain_asymptotic.cpp"
    )

add_executable(
                    test_time_domain_asymptotic_cmd
                    ${test_time_domain_asymptotic_src}
                )

set(
        test_time_domain_asymptotic_args
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdt.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtt.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttt.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdtx.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/time_domain/dGdttx.dat"
    )

add_test(
            NAME    test_time_domain_asymptotic
            COMMAND test_time_domain_asymptotic_cmd
            ${test_time_domain_asymptotic_args}
        )