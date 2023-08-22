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
    set( compile_option_flags /Wall /Wextra /EHsc)
elseif ( UNIX )
    set( compile_option_flags -Wall -Wextra -EHsc)
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
    target_link_libraries(${target} PUBLIC mkl_intel_lp64_dll mkl_intel_thread_dll mkl_core_dll)
endfunction()


# Find MPI
find_package( MPI REQUIRED )

function( set_mpi_target target )
    target_compile_options( ${target} PUBLIC "-D_MPI_BUILD" )
    target_include_directories( ${target} PRIVATE ${mpi_include_path} )
    target_link_libraries( ${target} PRIVATE MPI::MPI_CXX )
endfunction( )

###########################################################
################### Configure Tests #######################
###########################################################
add_subdirectory(tests)


###########################################################
################ SeaMotions Executables ###################
###########################################################

# Add executable to calculate chebyshev coefficients for the 
# infinite water depth expansion series
set(chebyshev_inf_depth_src
    ./src/config.hpp
    ./src/green/fit_chebyshev_inf_depth.cpp
    ${containers_module_files}
    ${green_module_files}
    ${math_module_files}
    ${tools_module_files}
)

add_executable(chebyshev_inf_depth_cmd ${chebyshev_inf_depth_src})
set_mkl_lib_target(chebyshev_inf_depth_cmd)

# Add library for finite water depth green function
set(
    lib_finite_water_depth_src
    ./src/config.hpp
    ./src/math/math_tools.hpp
    ./src/math/math_tools.cpp
    ./src/math/math_tools.txx
    ./src/green/pulsating_fin_depth.hpp
    ./src/green/pulsating_fin_depth.cpp
    ./src/math/chebyshev.hpp
    ./src/math/chebyshev.cpp
    ./src/math/special_math.hpp
    ./src/math/gauss.hpp
    ./src/math/gauss.cpp
    ./src/math/special_math.cpp
    ./src/tools.hpp
    ./src/tools.cpp
    ./src/tools.txx
    ./src/waves.hpp
    ./src/waves.cpp
)

add_library(lib_finite_water_depth ${lib_finite_water_depth_src})
install(
        TARGETS lib_finite_water_depth
        DESTINATION "${CMAKE_SOURCE_DIR}/libs"
        )

# Add test executable to check the version passes to the software
set(SOURCE_FILES src/main.cpp src/version.hpp src/version.cpp)
add_executable(version_exec ${SOURCE_FILES})