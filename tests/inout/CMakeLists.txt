
# Add test to check the input system
set(
        test_input_src
        test_input.cpp
        ${containers_module_files}
        ${inout_module_files}
        ${math_module_files}
        ${tools_module_files}
        ${version_module_files}
    )

add_executable( test_input_cmd ${test_input_src} )
set_mkl_lib_target( test_input_cmd )