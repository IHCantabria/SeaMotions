# Pack containers class definition
set(
        containers_module_files
        ${CMAKE_SOURCE_DIR}/src/containers/body_def.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/body_def.cpp
        ${CMAKE_SOURCE_DIR}/src/containers/containers.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/mpi_config.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/mpi_config.cpp
        ${CMAKE_SOURCE_DIR}/src/containers/panel_geom.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/panel_geom.cpp
        ${CMAKE_SOURCE_DIR}/src/containers/panel_geom_list.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/performance_stats.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/performance_stats.cpp
        ${CMAKE_SOURCE_DIR}/src/containers/source_node.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/source_node.cpp
    )

# Pack green finite water depth integral coefficients
set(
        green_findepth_coeffs_src
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/L1_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/L1_dA.cpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/L1_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/L1_dB.cpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/L1.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/L1.cpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/L2.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/L2.cpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/L3_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/L3_dA.cpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/L3_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/L3_dB.cpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/L3.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/L3.cpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/M1_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/M1_dA.cpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/M1_dB.cpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/M1_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/M1.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/M1.cpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/M2.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/M2.cpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/M3_dA.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/M3_dA.cpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/M3_dB.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/M3_dB.cpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/M3.hpp
        ${CMAKE_SOURCE_DIR}/src/green/fin_depth_coeffs/M3.cpp
    )

# Pack green infinite water depth integral coefficients
set(
        green_infdepth_coeffs_src
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R11_dX.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R11_dX.cpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R11.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R11.cpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R11A_dX.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R11A_dX.cpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R11B_dX.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R11B_dX.cpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R12_dX.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R12_dX.cpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R12.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R12.cpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R21_dX.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R21_dX.cpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R21.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R21.cpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R22_dX.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R22_dX.cpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R22.hpp
        ${CMAKE_SOURCE_DIR}/src/green/inf_depth_coeffs/R22.cpp
    )

# Pack green module files
set(
        green_module_files
        ${CMAKE_SOURCE_DIR}/src/green/chebyshev_inf_depth.hpp
        ${CMAKE_SOURCE_DIR}/src/green/common.hpp
        ${CMAKE_SOURCE_DIR}/src/green/common.cpp
        ${CMAKE_SOURCE_DIR}/src/green/dipole.hpp
        ${CMAKE_SOURCE_DIR}/src/green/dipole.cpp
        ${CMAKE_SOURCE_DIR}/src/green/integrals_db.hpp
        ${CMAKE_SOURCE_DIR}/src/green/integrals_db.cpp
        ${CMAKE_SOURCE_DIR}/src/green/pulsating_fin_depth.hpp
        ${CMAKE_SOURCE_DIR}/src/green/pulsating_fin_depth.cpp
        ${CMAKE_SOURCE_DIR}/src/green/pulsating_fin_depth_cheby.hpp
        ${CMAKE_SOURCE_DIR}/src/green/pulsating_fin_depth_cheby.cpp
        ${CMAKE_SOURCE_DIR}/src/green/pulsating_inf_depth_cheby.hpp
        ${CMAKE_SOURCE_DIR}/src/green/pulsating_inf_depth_cheby.cpp
        ${CMAKE_SOURCE_DIR}/src/green/pulsating_inf_depth_series.hpp
        ${CMAKE_SOURCE_DIR}/src/green/pulsating_inf_depth_series.cpp
        ${CMAKE_SOURCE_DIR}/src/green/pulsating_inf_depth.hpp
        ${CMAKE_SOURCE_DIR}/src/green/pulsating_inf_depth.cpp
        ${CMAKE_SOURCE_DIR}/src/green/source.hpp
        ${CMAKE_SOURCE_DIR}/src/green/source.cpp
        ${green_findepth_coeffs_src}
        ${green_infdepth_coeffs_src}
    )

# Pack hydrostatic module files
set(
        hydrostatic_module_files
        ${CMAKE_SOURCE_DIR}/src/hydrostatics.hpp
        ${CMAKE_SOURCE_DIR}/src/hydrostatics.cpp
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
)


# Pack Green function interfaces module files
set(
        interfaces_module_files
        ${CMAKE_SOURCE_DIR}/src/interfaces/grf_interface.hpp
        ${CMAKE_SOURCE_DIR}/src/interfaces/grf_interface.cpp
        ${CMAKE_SOURCE_DIR}/src/interfaces/grfdn_interface.hpp
        ${CMAKE_SOURCE_DIR}/src/interfaces/grfdn_interface.cpp
        ${CMAKE_SOURCE_DIR}/src/interfaces/gwf_interface.hpp
        ${CMAKE_SOURCE_DIR}/src/interfaces/gwf_interface.cpp
        ${CMAKE_SOURCE_DIR}/src/interfaces/gwfdn_interface.hpp
        ${CMAKE_SOURCE_DIR}/src/interfaces/gwfdn_interface.cpp
        ${CMAKE_SOURCE_DIR}/src/interfaces/hmf_interface.hpp
        ${CMAKE_SOURCE_DIR}/src/interfaces/hmf_interface.cpp
    )

# Pack math module files
set(
        math_module_files
        ${CMAKE_SOURCE_DIR}/src/math/chebyshev.hpp
        ${CMAKE_SOURCE_DIR}/src/math/chebyshev.cpp
        ${CMAKE_SOURCE_DIR}/src/math/integration.hpp
        ${CMAKE_SOURCE_DIR}/src/math/integration.cpp
        ${CMAKE_SOURCE_DIR}/src/math/integration.txx
        ${CMAKE_SOURCE_DIR}/src/math/gauss.hpp
        ${CMAKE_SOURCE_DIR}/src/math/gauss.cpp
        ${CMAKE_SOURCE_DIR}/src/math/math_interface.hpp
        ${CMAKE_SOURCE_DIR}/src/math/math_tools.hpp
        ${CMAKE_SOURCE_DIR}/src/math/math_tools.cpp
        ${CMAKE_SOURCE_DIR}/src/math/math_tools.txx
        ${CMAKE_SOURCE_DIR}/src/math/nonlinear_solvers_math.hpp
        ${CMAKE_SOURCE_DIR}/src/math/scalapack_solver.hpp
        ${CMAKE_SOURCE_DIR}/src/math/shape_functions.hpp
        ${CMAKE_SOURCE_DIR}/src/math/shape_functions.cpp
        ${CMAKE_SOURCE_DIR}/src/math/special_math.hpp
        ${CMAKE_SOURCE_DIR}/src/math/special_math.cpp
        ${CMAKE_SOURCE_DIR}/src/math/topology.hpp
        ${CMAKE_SOURCE_DIR}/src/math/topology.cpp
    )

# Pack mesh module files
set(
        mesh_module_files
        ${CMAKE_SOURCE_DIR}/src/mesh/mesh.hpp
        ${CMAKE_SOURCE_DIR}/src/mesh/mesh.cpp
        ${CMAKE_SOURCE_DIR}/src/mesh/mesh_group.hpp
        ${CMAKE_SOURCE_DIR}/src/mesh/mesh_group.cpp
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
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/freq_solver_tools.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/freq_solver_tools.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/gf_intensities.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/gf_intensities.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/hydromechanics.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/hydromechanics.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/potential.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/potential.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/raos.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/raos.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/wave_elevation.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/wave_elevation.cpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/froude_krylov.hpp
        ${CMAKE_SOURCE_DIR}/src/solvers/freq_domain/froude_krylov.cpp
    )

# Pack tools module files
set(
        tools_module_files
        ${CMAKE_SOURCE_DIR}/src/tools.hpp
        ${CMAKE_SOURCE_DIR}/src/tools.cpp
        ${CMAKE_SOURCE_DIR}/src/tools.txx
    )

# Pack waves module files
set(
        waves_module_files
        ${CMAKE_SOURCE_DIR}/src/waves.hpp
        ${CMAKE_SOURCE_DIR}/src/waves.cpp
    )

# Pack version files
set(
        version_module_files
        ${CMAKE_SOURCE_DIR}/src/version.hpp
        ${CMAKE_SOURCE_DIR}/src/version.cpp
)