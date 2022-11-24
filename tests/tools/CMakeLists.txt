
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
