
# Add test for dipole potential function
set(test_dipole_potential_src
        test_dipole_potential.cpp
        ../../src/config.hpp
        ../../src/config.cpp
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