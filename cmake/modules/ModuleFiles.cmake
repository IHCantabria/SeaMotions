# Pack containers class definition
set(
        containers_module_files
        ${CMAKE_SOURCE_DIR}/src/containers/body_def.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/body_def.cpp
        ${CMAKE_SOURCE_DIR}/src/containers/containers.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/load_condition.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/load_condition.cpp
        ${CMAKE_SOURCE_DIR}/src/containers/logger.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/logger.cpp
        ${CMAKE_SOURCE_DIR}/src/containers/matlin_group.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/mpi_config.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/mpi_config.cpp
        ${CMAKE_SOURCE_DIR}/src/containers/mpi_timer.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/panel_geom.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/panel_geom.cpp
        ${CMAKE_SOURCE_DIR}/src/containers/panel_geom_t.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/panel_geom_t.txx
        ${CMAKE_SOURCE_DIR}/src/containers/panel_geom_list.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/performance_stats.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/performance_stats.cpp
        ${CMAKE_SOURCE_DIR}/src/containers/rad_diff_data.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/rad_diff_data.cpp
        ${CMAKE_SOURCE_DIR}/src/containers/simulation_data.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/simulation_data.cpp
        ${CMAKE_SOURCE_DIR}/src/containers/source_node.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/source_node.cpp
    )

# Pack green finite water depth integral coefficients
set(
        green_findepth_coeffs_src
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/L0_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/L0_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/L0.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/L1_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/L1_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/L1.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/L2.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/L3_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/L3_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/L3.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/Linf_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/Linf_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/Linf.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/M1_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/M1_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/M1.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/M2.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/M3_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/M3_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/low/M3.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/L0_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/L0_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/L0.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/L1_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/L1_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/L1.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/L2.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/L3_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/L3_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/L3.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/Linf_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/Linf_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/Linf.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/M1_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/M1_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/M1.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/M2.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/M3_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/M3_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/high/M3.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/L0_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/L0_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/L0.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/L1_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/L1_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/L1.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/L2.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/L3_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/L3_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/L3.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/Linf_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/Linf_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/Linf.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/M1_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/M1_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/M1.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/M2.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/M3_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/M3_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/medium/M3.hpp
    )

# Pack green infinite water depth integral coefficients
set(
        green_infdepth_coeffs_src
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/low/R00_dX.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/low/R11_dX.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/low/R11.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/high/R00_dX.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/high/R11_dX.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/high/R11.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/medium/R00_dX.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/medium/R11_dX.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/medium/R11.hpp
    )

# Pack green module files
set(
        green_module_files
        ${CMAKE_SOURCE_DIR}/src/green/chebyshev_evaluator_base.hpp
        ${CMAKE_SOURCE_DIR}/src/green/chebyshev_evaluator_interface.hpp
        ${CMAKE_SOURCE_DIR}/src/green/chebyshev_traits_macros.hpp
        ${CMAKE_SOURCE_DIR}/src/green/chebyshev_traits.hpp
        ${CMAKE_SOURCE_DIR}/src/green/common.hpp
        ${CMAKE_SOURCE_DIR}/src/green/common.cpp
        ${CMAKE_SOURCE_DIR}/src/green/dipole.hpp
        ${CMAKE_SOURCE_DIR}/src/green/dipole.cpp
        ${CMAKE_SOURCE_DIR}/src/green/kochin.hpp
        ${CMAKE_SOURCE_DIR}/src/green/kochin.cpp
        # ${CMAKE_SOURCE_DIR}/src/green/pulsating_fin_depth.hpp
        # ${CMAKE_SOURCE_DIR}/src/green/pulsating_fin_depth.cpp
        # ${CMAKE_SOURCE_DIR}/src/green/pulsating_fin_depth_cheby.hpp
        # ${CMAKE_SOURCE_DIR}/src/green/pulsating_fin_depth_cheby.cpp
        # ${CMAKE_SOURCE_DIR}/src/green/pulsating_inf_depth_cheby.hpp
        # ${CMAKE_SOURCE_DIR}/src/green/pulsating_inf_depth_cheby.cpp
        # ${CMAKE_SOURCE_DIR}/src/green/pulsating_inf_depth_cheby_v2.hpp
        # ${CMAKE_SOURCE_DIR}/src/green/pulsating_inf_depth_cheby_v2.cpp
        # ${CMAKE_SOURCE_DIR}/src/green/pulsating_inf_depth_series.hpp
        # ${CMAKE_SOURCE_DIR}/src/green/pulsating_inf_depth_series.cpp
        # ${CMAKE_SOURCE_DIR}/src/green/pulsating_inf_depth.hpp
        # ${CMAKE_SOURCE_DIR}/src/green/pulsating_inf_depth.cpp
        # ${CMAKE_SOURCE_DIR}/src/green/pulsating_inf_depth_v2.hpp
        # ${CMAKE_SOURCE_DIR}/src/green/pulsating_inf_depth_v2.cpp
        ${CMAKE_SOURCE_DIR}/src/green/source.hpp
        ${CMAKE_SOURCE_DIR}/src/green/source.cpp
        ${CMAKE_SOURCE_DIR}/src/green/source.txx
        ${green_findepth_coeffs_src}
        ${green_infdepth_coeffs_src}
    )

# Pack hydrostatic module files
set(
        hydrostatic_module_files
        ${CMAKE_SOURCE_DIR}/src/hydrostatics.hpp
        ${CMAKE_SOURCE_DIR}/src/hydrostatics.cpp
        ${CMAKE_SOURCE_DIR}/src/initial_stability.hpp
        ${CMAKE_SOURCE_DIR}/src/initial_stability.txx
    )

# Pack inout module files
set(
        inout_module_files
        ${CMAKE_SOURCE_DIR}/src/inout/input.hpp
        ${CMAKE_SOURCE_DIR}/src/inout/input.cpp
        ${CMAKE_SOURCE_DIR}/src/inout/output.hpp
        ${CMAKE_SOURCE_DIR}/src/inout/output.cpp
        ${CMAKE_SOURCE_DIR}/src/inout/reader.hpp
        ${CMAKE_SOURCE_DIR}/src/inout/reader.cpp
        ${CMAKE_SOURCE_DIR}/src/inout/reader.txx
        ${CMAKE_SOURCE_DIR}/src/inout/vtu.hpp
        ${CMAKE_SOURCE_DIR}/src/inout/vtu.cpp
)
# Pack Green function interfaces module files
set(
        interfaces_module_files
        # ${CMAKE_SOURCE_DIR}/src/interfaces/grf_interface.hpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/grf_interface.cpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/grfdn_interface.hpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/grfdn_interface.cpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/grfdx_interface.hpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/grfdx_interface.cpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/grfdy_interface.hpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/grfdy_interface.cpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/grfdz_interface.hpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/grfdz_interface.cpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/gwf_interface.hpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/gwf_interface.cpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/gwfdn_interface.hpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/gwfdn_interface.cpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/gwfdx_interface.hpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/gwfdx_interface.cpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/gwfdy_interface.hpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/gwfdy_interface.cpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/gwfdz_interface.hpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/gwfdz_interface.cpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/hmf_interface.hpp
        # ${CMAKE_SOURCE_DIR}/src/interfaces/hmf_interface.cpp
        ${CMAKE_SOURCE_DIR}/src/interfaces/gwfcns_interface_t.hpp
        ${CMAKE_SOURCE_DIR}/src/interfaces/gwfcns_interface_t.txx
        ${CMAKE_SOURCE_DIR}/src/interfaces/hydrostatic_force_interface.hpp
        ${CMAKE_SOURCE_DIR}/src/interfaces/kochin_interface.hpp
        ${CMAKE_SOURCE_DIR}/src/interfaces/kochin_interface.cpp
        ${CMAKE_SOURCE_DIR}/src/interfaces/stab_interface_t.hpp
        ${CMAKE_SOURCE_DIR}/src/interfaces/stab_interface_t.txx
    )

# Pack math module files
set(
        math_module_files
        ${CMAKE_SOURCE_DIR}/src/math/chebyshev.hpp
        ${CMAKE_SOURCE_DIR}/src/math/chebyshev.cpp
        ${CMAKE_SOURCE_DIR}/src/math/chebyshev.txx
        ${CMAKE_SOURCE_DIR}/src/math/chebyshev_evaluation.hpp
        ${CMAKE_SOURCE_DIR}/src/math/custensor/base.hpp
        ${CMAKE_SOURCE_DIR}/src/math/custensor/binary_expression.hpp
        ${CMAKE_SOURCE_DIR}/src/math/custensor/custensor_ext_expression.hpp
        ${CMAKE_SOURCE_DIR}/src/math/custensor/custensor_impl.hpp
        ${CMAKE_SOURCE_DIR}/src/math/custensor/custensor.hpp
        ${CMAKE_SOURCE_DIR}/src/math/custensor/operators.hpp
        ${CMAKE_SOURCE_DIR}/src/math/custensor/scalar_expression.hpp
        ${CMAKE_SOURCE_DIR}/src/math/custensor/unary_expression.hpp
        ${CMAKE_SOURCE_DIR}/src/math/bessel_factory.hpp
        ${CMAKE_SOURCE_DIR}/src/math/euler_transforms.hpp
        ${CMAKE_SOURCE_DIR}/src/math/euler_transforms.cpp
        ${CMAKE_SOURCE_DIR}/src/math/euler_transforms.txx
        ${CMAKE_SOURCE_DIR}/src/math/gauss.hpp
        ${CMAKE_SOURCE_DIR}/src/math/gauss.cpp
        ${CMAKE_SOURCE_DIR}/src/math/integration.hpp
        ${CMAKE_SOURCE_DIR}/src/math/integration.cpp
        ${CMAKE_SOURCE_DIR}/src/math/integration.txx
        ${CMAKE_SOURCE_DIR}/src/math/math_interface.hpp
        ${CMAKE_SOURCE_DIR}/src/math/math_tools.hpp
        ${CMAKE_SOURCE_DIR}/src/math/math_tools.cpp
        ${CMAKE_SOURCE_DIR}/src/math/math_tools.txx
        ${CMAKE_SOURCE_DIR}/src/math/nonlinear_solvers_math.hpp
        ${CMAKE_SOURCE_DIR}/src/math/random.hpp
        ${CMAKE_SOURCE_DIR}/src/math/random.cpp
        ${CMAKE_SOURCE_DIR}/src/math/scalapack_solver.hpp
        ${CMAKE_SOURCE_DIR}/src/math/sparse/sparse_containers.hpp
        ${CMAKE_SOURCE_DIR}/src/math/sparse/sparse_math.hpp
        ${CMAKE_SOURCE_DIR}/src/math/sparse/sparse_solver.hpp
        ${CMAKE_SOURCE_DIR}/src/math/sparse/sparse_solver.cpp
        ${CMAKE_SOURCE_DIR}/src/math/sparse/sparse_tools.hpp
        ${CMAKE_SOURCE_DIR}/src/math/sparse/sparse_tools.cpp
        ${CMAKE_SOURCE_DIR}/src/math/shape_functions.hpp
        ${CMAKE_SOURCE_DIR}/src/math/shape_functions.cpp
        ${CMAKE_SOURCE_DIR}/src/math/special_math.hpp
        ${CMAKE_SOURCE_DIR}/src/math/special_math.cpp
        ${CMAKE_SOURCE_DIR}/src/math/topology.hpp
        ${CMAKE_SOURCE_DIR}/src/math/topology.cpp
        ${CMAKE_SOURCE_DIR}/src/math/topology.txx
    )

# Pack mesh module files
set(
        mesh_module_files
        ${CMAKE_SOURCE_DIR}/src/mesh/gmsh_reader/gmsh_binary_tools.hpp
        ${CMAKE_SOURCE_DIR}/src/mesh/gmsh_reader/gmsh_element.hpp
        ${CMAKE_SOURCE_DIR}/src/mesh/gmsh_reader/gmsh_node.hpp
        ${CMAKE_SOURCE_DIR}/src/mesh/gmsh_reader/read_gmsh.hpp
        ${CMAKE_SOURCE_DIR}/src/mesh/gmsh_reader/read_gmsh.cpp
        ${CMAKE_SOURCE_DIR}/src/mesh/mesh.hpp
        ${CMAKE_SOURCE_DIR}/src/mesh/mesh.cpp
        ${CMAKE_SOURCE_DIR}/src/mesh/mesh_group.hpp
        ${CMAKE_SOURCE_DIR}/src/mesh/mesh_group.cpp
        ${CMAKE_SOURCE_DIR}/src/mesh/mesh_refinement.hpp
        ${CMAKE_SOURCE_DIR}/src/mesh/mesh_refinement.cpp
        ${CMAKE_SOURCE_DIR}/src/mesh/mesh_operations.hpp
        ${CMAKE_SOURCE_DIR}/src/mesh/mesh_operations.cpp
        ${CMAKE_SOURCE_DIR}/src/mesh/rigid_body_mesh.hpp
        ${CMAKE_SOURCE_DIR}/src/mesh/rigid_body_mesh.cpp
        ${CMAKE_SOURCE_DIR}/src/mesh/tools.hpp
        ${CMAKE_SOURCE_DIR}/src/mesh/tools.cpp
    )


# Pack mpi interface files
set(
        mpiint_module_files
        ${CMAKE_SOURCE_DIR}/src/mpi_interface.hpp
    )

# Pack solvers files
set(
        sm_freq_module_files
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/seamotions_freq.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/diffraction.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/diffraction.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/frequency_solver.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/frequency_solver.txx
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/froude_krylov.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/froude_krylov.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/formulation_kernel_backend.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/formulation_kernel_backend.txx
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/global_static_matrix.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/global_static_matrix.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/hydromechanics.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/hydromechanics.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/panel_fields.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/panel_fields.txx
        # ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/potential.hpp
        # ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/potential.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/qtf.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/qtf.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/qtf_indirect_method.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/qtf_indirect_method.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/raos.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/raos.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/tools.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/tools.cpp
        # ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/velocities.hpp
        # ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/velocities.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/wave_elevation.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/wave_elevation.cpp
    )

# Pack SeaMotions Stability Solver files
set(
        sm_stab_module_files
        ${CMAKE_SOURCE_DIR}/src/solvers/stability/seamotions_stab.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/stability/gz_point.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/stability/gz_point.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/stability/hydrostatic_force_nlin.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/stability/hydrostatic_force_nlin.txx
        ${CMAKE_SOURCE_DIR}/src/solvers/stability/stab_input.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/stability/stab_input.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/stability/stab_output.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/stability/stab_output.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/stability/stab_output.txx
        ${CMAKE_SOURCE_DIR}/src/solvers/stability/stab_solver.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/stability/stab_solver.cpp
    )

# Pack tools module files
set(
        tools_module_files
        ${CMAKE_SOURCE_DIR}/src/tools.hpp
        ${CMAKE_SOURCE_DIR}/src/tools.cpp
        ${CMAKE_SOURCE_DIR}/src/tools.txx
        ${CMAKE_SOURCE_DIR}/src/static_tools.hpp
    )

# Pack waves module files
set(
        waves_module_files
        ${CMAKE_SOURCE_DIR}/src/waves/wave_dispersion_base_fo.hpp
        ${CMAKE_SOURCE_DIR}/src/waves/wave_dispersion_base_fo.cpp
        ${CMAKE_SOURCE_DIR}/src/waves/wave_dispersion_fo.hpp
        ${CMAKE_SOURCE_DIR}/src/waves/wave_dispersion_fo.txx
        ${CMAKE_SOURCE_DIR}/src/waves/wave_dispersion_so.hpp
        ${CMAKE_SOURCE_DIR}/src/waves/wave_dispersion_so.cpp
        ${CMAKE_SOURCE_DIR}/src/waves/waves_common.hpp
        ${CMAKE_SOURCE_DIR}/src/waves/waves_common.cpp
    )

# Pack version files
set(
        version_module_files
        ${CMAKE_SOURCE_DIR}/src/version.hpp
        ${CMAKE_SOURCE_DIR}/src/version.cpp
)