
# Define common files to compile green function tests
set(green_common_src
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