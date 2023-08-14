
# Add test for the cross product
set(
        test_cross_product_src
        test_cross_product.cpp
        ../../src/config.hpp
        ../../src/math/math_tools.hpp
        ../../src/math/math_tools.cpp
        ../../src/math/math_tools.txx
    )
    
add_executable( test_cross_product_cmd ${test_cross_product_src} )
set_mkl_lib_target( test_cross_product_cmd )

add_test(
            NAME test_cross_product 
            COMMAND test_cross_product_cmd
        )


# Add test for the vector arithmetics
set(
        test_vector_arithmetic_src
        test_vector_arithmetics.cpp
        ../../src/config.hpp
        ../../src/math/math_tools.hpp
        ../../src/math/math_tools.cpp
    )

add_executable( test_vector_arithmetic_cmd ${test_vector_arithmetic_src} )
set_mkl_lib_target( test_vector_arithmetic_cmd )

add_test(
            NAME test_vector_arithmetic 
            COMMAND test_vector_arithmetic_cmd
        )


# Add test to check the romberg algorithm implemented
set(
        test_romberg_src
        test_romberg.cpp
        ../../src/config.hpp
        ../../src/math/math_tools.hpp
        ../../src/math/math_tools.cpp
        ../../src/math/math_tools.txx
        ../../src/tools.hpp
        ../../src/tools.cpp
        ../../src/tools.txx
    )

add_executable( test_romberg_cmd ${test_romberg_src} )
set_mkl_lib_target( test_romberg_cmd )

add_test(
            NAME test_romberg 
            COMMAND test_romberg_cmd
        )


# Add test to check the factorial function
set(
        test_factorial_src
        test_factorial.cpp
        ../../src/config.hpp
        ../../src/math/math_tools.hpp
        ../../src/math/math_tools.cpp
        ../../src/math/math_tools.txx
        ../../src/tools.hpp
        ../../src/tools.cpp
        ../../src/tools.txx
    )

add_executable( test_factorial_cmd ${test_factorial_src} )
set_mkl_lib_target( test_factorial_cmd )

add_test(
            NAME test_factorial 
            COMMAND test_factorial_cmd
        )


# Configure performance test for auxiliary bessel functions
set(
        test_performance_aux_bessel_src
        test_performance_aux_bessel.cpp
        ../../src/config.hpp
        ${containers_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable( test_performance_aux_bessel_cmd ${test_performance_aux_bessel_src} )
set_mkl_lib_target( test_performance_aux_bessel_cmd )

# Added test to check the psi function
set(    
        test_psi_function_src
        test_psi_function.cpp
        ../containers.hpp
        ${containers_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable( test_psi_function_cmd ${test_psi_function_src} )
set_mkl_lib_target( test_psi_function_cmd )

set( psi_refdata_file_path "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/psi.dat" )
add_test(
            NAME test_psi_function 
            COMMAND test_psi_function_cmd 
            ${psi_refdata_file_path}
        )

# Configure test to check the accuracy of the bessel functions
set(
        test_bessel_functions_src
        test_bessel_functions.cpp
        ../../src/config.hpp
        ../containers.hpp
        ${containers_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable( test_bessel_functions_cmd ${test_bessel_functions_src} )
set_mkl_lib_target( test_bessel_functions_cmd )

set( file_path_besselj0 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/besselj0.dat" )
set( file_path_besselj1 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/besselj1.dat" )
set( file_path_bessely0 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/bessely0.dat" )
set( file_path_bessely1 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/bessely1.dat" )
set( file_path_struve0 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/struve0.dat" )
set( file_path_struve1 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/struve1.dat" )
set( file_path_besseli0 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/besseli0.dat" )
set( file_path_besseli1 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/besseli1.dat" )
set( file_path_besselk0 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/besselk0.dat" )
set( file_path_besselk1 "${CMAKE_SOURCE_DIR}/tests/tests_data/special_math/besselk1.dat" )
add_test(
            NAME test_bessel_functions 
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


# Configure test to check the accuracy of the legendre polynomials implementation
set(
        test_legendre_poly_src
        test_legendre_poly.cpp
        ../../src/config.hpp
        ../containers.hpp
        ${containers_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable( test_legendre_poly_cmd ${test_legendre_poly_src} )
set_mkl_lib_target( test_legendre_poly_cmd )

set(
        legendre_file_path 
        "${CMAKE_SOURCE_DIR}/tests/tests_data/legendre_poly.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/legendre_poly_der.dat"
    )

add_test(
            NAME test_legendre_poly
            COMMAND test_legendre_poly_cmd
            ${legendre_file_path}
        )


# Configure test to check the accuracy of the chebyshev polynomials implementation
set(
        test_chebyshev_poly_src
        test_chebyshev_poly.cpp
        ../../src/config.hpp
        ../containers.hpp
        ${containers_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable( test_chebyshev_poly_cmd ${test_chebyshev_poly_src} )
set_mkl_lib_target( test_chebyshev_poly_cmd )
set(
        chebyshev_file_path 
        "${CMAKE_SOURCE_DIR}/tests/tests_data/chebyshev_poly.dat"
        "${CMAKE_SOURCE_DIR}/tests/tests_data/chebyshev_poly_der.dat"
    )
add_test(
            NAME test_chebyshev_poly
            COMMAND test_chebyshev_poly_cmd
            ${chebyshev_file_path}
        )

# Configure test to check exponential integral i
set(
        test_expint_src
        test_expint.cpp
        ../../src/config.hpp
        ../containers.hpp
        ${containers_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable( test_expint_cmd ${test_expint_src} )
set_mkl_lib_target( test_expint_cmd )

set( file_path_expint_i "${CMAKE_SOURCE_DIR}/tests/tests_data/expint_i.dat" )
add_test(
            NAME test_exp_int
            COMMAND test_expint_cmd
            ${file_path_expint_i}
        )


# Add test to check the generation of the Chebyshev roots
set(
        test_chebyshev_roots_src
        test_chebyshev_roots.cpp
        ../../src/config.hpp
        ${containers_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable( test_chebyshev_roots_cmd ${test_chebyshev_roots_src} )
set_mkl_lib_target( test_chebyshev_roots_cmd )

set( file_path_cheby_roots "${CMAKE_SOURCE_DIR}/tests/tests_data/cheby_roots.dat" )
add_test(
            NAME test_chebyshev_roots
            COMMAND test_chebyshev_roots_cmd
            ${file_path_cheby_roots}
        )


# Add test to check the bisection method
set(
        test_nonlinear_1d_solvers_src
        test_nonlinear_1d_solvers.cpp
        ../../src/config.hpp
        ${containers_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable( test_nonlinear_1d_solvers_cmd ${test_nonlinear_1d_solvers_src} )
set_mkl_lib_target( test_nonlinear_1d_solvers_cmd )

add_test(
            NAME test_nonlinear_1d_solvers
            COMMAND test_nonlinear_1d_solvers_cmd
        )


# Add test to check the complex integration functions
set(
        test_complex_integration_src
        test_complex_integration.cpp
        ../../src/config.hpp
        ${containers_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable( test_complex_integration_cmd ${test_complex_integration_src} )
set_mkl_lib_target( test_complex_integration_cmd )

add_test(NAME test_complex_integration
        COMMAND test_complex_integration_cmd
        )


# Add test to check the sign function
set(
        test_sign_src
        test_sign.cpp
        ../../src/config.hpp
        ${containers_module_files}
        ${math_module_files}
        ${tools_module_files}
    )

add_executable(test_sign_cmd ${test_sign_src})
set_mkl_lib_target( test_sign_cmd )

add_test(
        NAME test_sign
        COMMAND test_sign_cmd
        )


# Add test to check Haar wavelets integration method
set(
        test_haar_quad_src
        test_haar_quad.cpp
        ../../src/config.hpp
        ${math_module_files}
    )

add_executable( test_haar_quad_cmd ${test_haar_quad_src} )

# Add test to check quadrature integration for 2D shapes
set(
        test_quadrature_2d
        test_quadrature_2d.cpp
        ../../src/config.hpp
        ${math_module_files}
        ${tools_module_files}
    )

add_executable( test_quadrature_2d_cmd ${test_quadrature_2d} )
set_mkl_lib_target( test_quadrature_2d_cmd )

set( 
        test_quadrature_2d_files
        E:/sergio/developments/SeaMotions/tests/tests_data/triangle_elems_mesh.poly 
        E:/sergio/developments/SeaMotions/tests/tests_data/square_elems_mesh.poly
    )

add_test(
            NAME    test_quadrature_2d
            COMMAND test_quadrature_2d_cmd
            ${test_quadrature_2d_files}
        )