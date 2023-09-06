cmake_minimum_required(VERSION 3.9)
project(SEAMOTIONS CXX)

###########################################################
################# Include CMake Modules ###################
###########################################################
# Appends the cmake/modules path to MAKE_MODULE_PATH variable.
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})

# Include modules
include( AuxTools )
include( GetGitRevisionDescription )
include( ModuleFiles )


###########################################################
#################### Configure project ####################
###########################################################
# Configure project to work with ctest
enable_testing()

# Configure project ot be compiled with a specific standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Include  Get Git Version
git_describe(VERSION --tags --dirty=-d)

# Parse the version information into pieces
string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+(.*)" "\\1" VERSION_SHA1 "${VERSION}")
set(VERSION_SHORT "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

# Configure version file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/version.cpp.in
               ${CMAKE_CURRENT_SOURCE_DIR}/src/version.cpp)
set(version_file "${CMAKE_CURRENT_BINARY_DIR}/version.cpp")

# Define precision to compile the program: simple or double
if ("${PRECISION}" STREQUAL "SIMPLE")
    set(PRECISION_FLAG "-DSIMPLE_PREC")
elseif("${PRECISION}" STREQUAL "DOUBLE")
    set(PRECISION_FLAG "-DDOUBLE_PREC")
else()
    message(FATAL_ERROR "PRECISION VALUE: ${PRECISION} not valid. Accepted values: SIMPLE | DOUBLE")
endif()

###########################################################
########## Define Compilation configuration ###############
###########################################################
# Add precision flag to the compilation options
set( compile_option_flags ${PRECISION_FLAG} )

# Add warning options to the compiler
if ( WIN32 )
    set( compile_option_flags /Wall /Wextra /EHsc )
elseif ( UNIX )
    set( compile_option_flags -Wall -Wextra -EHsc )
endif()

# Check for Debug or release build
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set( compile_option_flags ${compile_option_flags} "-DDEBUG_BUILD")
endif()

# Apply compile options
add_compile_options( ${compile_option_flags} )

###########################################################
#################### Find Dependencies ####################
###########################################################
# Find MKL
find_package(MKL CONFIG REQUIRED)

function(set_mkl_lib_target target)
    # Add target include directories to work with external libraries
    target_include_directories(${target} PUBLIC ${MKL_INCLUDE})

    # Add target libraries to link with
    target_link_libraries(${target} PUBLIC mkl_intel_lp64_dll mkl_intel_thread_dll mkl_core_dll mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64)
endfunction()


# Find MPI
find_package( MPI REQUIRED )

function( set_mpi_target target )
    target_compile_options( ${target} PUBLIC "-DMPI_BUILD" )
    target_include_directories( ${target} PRIVATE ${mpi_include_path} )
    target_link_libraries( ${target} PRIVATE MPI::MPI_CXX )
endfunction( )

# Find HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX HL)

function( set_hdf5_target target )
    target_compile_options( ${target} PUBLIC "-D_HDF5_BUILD" )
    target_include_directories( ${target} PRIVATE ${HDF5_INCLUDE_DIRS} )
    target_link_libraries( ${target} PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} )
endfunction( )

###########################################################
################### Configure Tests #######################
###########################################################
add_subdirectory(tests)


###########################################################
################ SeaMotions Executables ###################
###########################################################

# Add SeaMotions frequency solver
set(
        seamotions_freq_src
        ./src/seamotions_freq.cpp
        ${containers_module_files}
        ${hydrostatic_module_files}
        ${inout_module_files}
        ${math_module_files}
        ${mesh_module_files}
        ${tools_module_files}
        ${version_module_files}
    )

add_executable(
                    seamotions_freq
                    ${seamotions_freq_src}
                )
set_mkl_lib_target( seamotions_freq )
set_mpi_target( seamotions_freq )
set_hdf5_target( seamotions_freq )