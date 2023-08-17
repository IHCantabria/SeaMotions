
set(containers_common_src
    ../../src/config.hpp
    ../../src/tools.hpp
    ../../src/tools.cpp
    ../../src/tools.txx
)

# Add test to check if the points are inside a polygon
set(test_inside_polygon_src
    test_inside_polygon.cpp
    ${containers_module_files}
    ${containers_common_src}
    ${math_module_files}
)

add_executable( test_inside_polygon_cmd ${test_inside_polygon_src} )
set_mkl_lib_target( test_inside_polygon_cmd )
set(path_data_1 "${CMAKE_SOURCE_DIR}/tests/tests_data/is_inside/is_inside_1.dat")
set(path_data_2 "${CMAKE_SOURCE_DIR}/tests/tests_data/is_inside/is_inside_2.dat")
set(path_data_3 "${CMAKE_SOURCE_DIR}/tests/tests_data/is_inside/is_inside_3.dat")
add_test(NAME test_inside_polygon COMMAND test_inside_polygon_cmd ${path_data_1} ${path_data_2} ${path_data_3})