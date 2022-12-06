
set(containers_common_src

)

# Add test to check if the points are inside a polygon
set(test_inside_polygon_src
    test_inside_polygon.cpp
    ${containers_common_src}
)

add_executable(test_inside_polygon_cmd ${test_inside_polygon_src})