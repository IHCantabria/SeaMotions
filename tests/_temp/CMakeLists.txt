
# Test MKL norm and its usage. It also tests the math_interface module
add_executable(test_norm
test_norm.cpp
../../src/config.hpp
../../src/config.cpp
../../src/math_interface.hpp
)
set_mkl_lib_target(test_norm)

# Add test to asset the performance of interl compiler vectors and the VM MKL implementation
add_executable(vector_performance
    vector_performance.cpp
    ../../src/tools.hpp
    ../../src/tools.cpp
)
set_mkl_lib_target(vector_performance)