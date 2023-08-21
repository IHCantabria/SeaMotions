# Pack containers class definition
set(
        containers_module_files
        ${CMAKE_SOURCE_DIR}/src/containers/containers.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/panel_geom.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/panel_geom.cpp
        ${CMAKE_SOURCE_DIR}/src/containers/panel_geom_list.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/performance_stats.hpp
        ${CMAKE_SOURCE_DIR}/src/containers/performance_stats.cpp
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
        ${CMAKE_SOURCE_DIR}/src/math/special_math.hpp
        ${CMAKE_SOURCE_DIR}/src/math/special_math.cpp
        ${CMAKE_SOURCE_DIR}/src/math/topology.hpp
        ${CMAKE_SOURCE_DIR}/src/math/topology.cpp
    )

# Pack mesh module files
set(
        mesh_module_files
        ${CMAKE_SOURCE_DIR}/src/mesh/mesh.hpp
        ${CMAKE_SOURCE_DIR}/src/mesh/tools.hpp
        ${CMAKE_SOURCE_DIR}/src/mesh/tools.cpp
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