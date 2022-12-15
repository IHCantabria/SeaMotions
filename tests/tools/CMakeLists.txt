
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
    ../containers.hpp
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
set(psi_refdata_file_path "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/psi.dat")
add_test(NAME test_psi_function COMMAND test_psi_function_cmd ${psi_refdata_file_path})

# Configure test to check the accuracy of the bessel functions
set(test_bessel_functions_src
    test_bessel_functions.cpp
    ../containers.hpp
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

add_executable(test_bessel_functions_cmd ${test_bessel_functions_src})
set_mkl_lib_target(test_bessel_functions_cmd)
set(file_path_besselj0 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/besselj0.dat")
set(file_path_besselj1 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/besselj1.dat")
set(file_path_bessely0 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/bessely0.dat")
set(file_path_bessely1 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/bessely1.dat")
set(file_path_struve0 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/struve0.dat")
set(file_path_struve1 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/struve1.dat")
set(file_path_besseli0 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/besseli0.dat")
set(file_path_besseli1 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/besseli1.dat")
set(file_path_besselk0 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/besselk0.dat")
set(file_path_besselk1 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/besselk1.dat")
add_test(NAME test_bessel_functions 
        COMMAND test_bessel_functions_cmd 
                ${file_path_besselj0}
                ${file_path_besselj1}
                ${file_path_bessely0}
                ${file_path_bessely1}
                ${file_path_struve0}
                ${file_path_struve1}
                ${file_path_besseli0}
                ${file_path_besseli1}
                ${file_path_besselk0}
                ${file_path_besselk1}
                )