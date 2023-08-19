
# Add test to check the preprocessor flag for MPI builds is 
# configured correctly
add_executable( test_mpi_build_0_cmd test_mpi_build.cpp )

add_test(
            NAME test_mpi_build_0
            COMMAND test_mpi_build_0_cmd
            0
        )
    
add_executable( test_mpi_build_1_cmd test_mpi_build.cpp )
set_mpi_target( test_mpi_build_1_cmd )

add_test(
            NAME test_mpi_build_1
            COMMAND test_mpi_build_1_cmd
            1
        )

# Add second example on how to send and receive MPI signals
add_executable( test_mpi_send_recv_cmd test_send_recv.cpp )
set_mpi_target( test_mpi_send_recv_cmd )

# Add program to send and receive variable sized data in
# between processes
set(
        test_mpi_status_probe_src
        test_status_probe.cpp
        ../../src/config.hpp
        ${math_module_files}
)
add_executable( test_mpi_status_probe_cmd ${test_mpi_status_probe_src} )
set_mpi_target( test_mpi_status_probe_cmd )
set_mkl_lib_target( test_mpi_status_probe_cmd )
