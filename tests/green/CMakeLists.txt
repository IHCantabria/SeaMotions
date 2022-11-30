
# Add test for dipole potential function
set(test_dipole_potential_src
        test_dipole_potential.cpp
        ../../src/config.hpp
        ../../src/config.cpp
        ../../src/green/common.hpp
        ../../src/green/common.cpp
        ../../src/containers.hpp
        ../../src/green/dipole.hpp
        ../../src/green/dipole.cpp
        ../../src/math_interface.hpp
        ../../src/math_tools.hpp
        ../../src/math_tools.cpp
        ../../src/tools.hpp
        ../../src/tools.cpp
        )
add_executable(test_dipole_potential ${test_dipole_potential_src})
set_mkl_lib_target(test_dipole_potential)


# Add test for source velocity function
set(test_source_velocity_src
        test_source_velocity.cpp
        ../../src/config.hpp
        ../../src/config.cpp
        ../../src/green/common.hpp
        ../../src/green/common.cpp
        ../../src/containers.hpp
        ../../src/green/dipole.hpp
        ../../src/green/dipole.cpp
        ../../src/math_interface.hpp
        ../../src/math_tools.hpp
        ../../src/math_tools.cpp
        ../../src/green/source.hpp
        ../../src/green/source.cpp
        ../../src/tools.hpp
        ../../src/tools.cpp
)

add_executable(test_source_velocity ${test_source_velocity_src})
set_mkl_lib_target(test_source_velocity)


# Configure source performance test Hess vs Newman
set(performance_source_newman_hess_src 
        performance_hess_vs_newman.cpp
        ../../src/config.hpp
        ../../src/config.cpp
        ../../src/green/common.hpp
        ../../src/green/common.cpp
        ../../src/containers.hpp
        ../../src/green/dipole.hpp
        ../../src/green/dipole.cpp
        ../../src/math_interface.hpp
        ../../src/math_tools.hpp
        ../../src/math_tools.cpp
        ../../src/green/source.hpp
        ../../src/green/source.cpp
        test_panel_geometries.hpp
        test_panel_geometries.cpp
        ../../src/tools.hpp
        ../../src/tools.cpp
)

add_executable(performance_source_newman_hess ${performance_source_newman_hess_src})
set_mkl_lib_target(performance_source_newman_hess)