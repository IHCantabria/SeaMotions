
# Add test for the cross product
set(test_cross_product_src
    test_cross_product.cpp
    ../../src/config.hpp
    ../../src/config.cpp
    ../../src/math_tools.hpp
    ../../src/math_tools.cpp
    ../../src/math_tools.txx
    )
    
add_executable(test_cross_product_cmd ${test_cross_product_src})
set_mkl_lib_target(test_cross_product_cmd)
add_test(NAME test_cross_product COMMAND test_cross_product_cmd)


# Add test for the vector arithmetics
set(test_vector_arithmetic_src
    test_vector_arithmetics.cpp
    ../../src/config.hpp
    ../../src/config.cpp
    ../../src/math_tools.hpp
    ../../src/math_tools.cpp
)

add_executable(test_vector_arithmetic_cmd ${test_vector_arithmetic_src})
set_mkl_lib_target(test_vector_arithmetic_cmd)
add_test(NAME test_vector_arithmetic COMMAND test_vector_arithmetic_cmd)


# Add test to check the romberg algorithm implemented
set(test_romberg_src
    test_romberg.cpp
    ../../src/config.hpp
    ../../src/config.cpp
    ../../src/math_tools.hpp
    ../../src/math_tools.cpp
    ../../src/math_tools.txx
    ../../src/tools.hpp
    ../../src/tools.cpp
    ../../src/tools.txx
)

add_executable(test_romberg_cmd ${test_romberg_src})
set_mkl_lib_target(test_romberg_cmd)
add_test(NAME test_romberg COMMAND test_romberg_cmd)


# Add test to check the factorial function
set(test_factorial_src
    test_factorial.cpp
    ../../src/config.hpp
    ../../src/config.cpp
    ../../src/math_tools.hpp
    ../../src/math_tools.cpp
    ../../src/math_tools.txx
    ../../src/tools.hpp
    ../../src/tools.cpp
    ../../src/tools.txx
)

add_executable(test_factorial_cmd ${test_factorial_src})
set_mkl_lib_target(test_factorial_cmd)
add_test(NAME test_factorial COMMAND test_factorial_cmd)


# Configure performance test for auxiliary bessel functions
set(performance_aux_bessel_src
        test_performance_aux_bessel.cpp
        ../../src/config.hpp
        ../../src/config.cpp
        ../../src/containers.hpp
        ../../src/special_math.hpp
        ../../src/special_math.cpp
        ../../src/tools.hpp
        ../../src/tools.cpp
        ../../src/tools.txx
)
add_executable(test_performance_aux_bessel_cmd ${performance_aux_bessel_src})

# Added test to check the psi function
set(test_psi_function_src
    test_psi_function.cpp
    ../../src/config.hpp
    ../../src/config.cpp
    ../../src/math_tools.hpp
    ../../src/math_tools.cpp
    ../../src/math_tools.txx
    ../../src/special_math.hpp
    ../../src/special_math.cpp
    ../../src/tools.hpp
    ../../src/tools.cpp
    ../../src/tools.txx
)

add_executable(test_psi_function_cmd ${test_psi_function_src})
set_mkl_lib_target(test_psi_function_cmd)
