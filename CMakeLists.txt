cmake_minimum_required(VERSION 3.9)
project(SEAMOTIONS CXX)

###########################################################
#################### Configure project ####################
###########################################################
# Appends the cmake/modules path to MAKE_MODULE_PATH variable.
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})

# Include  Get Git Version
include(GetGitRevisionDescription)
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

message(STATUS "Floating precision: ${PRECISION}")

###########################################################
#################### Find Dependencies ####################
###########################################################
# Find MKL
find_package(MKL CONFIG REQUIRED)

function(set_mkl_lib_target target)
    target_compile_options(${target} PUBLIC ${PRECISION_FLAG} -Wall -Wextra)
    target_include_directories(${target} PUBLIC ${MKL_INCLUDE})
    target_link_libraries(${target} PUBLIC mkl_intel_lp64 mkl_intel_thread mkl_core)
endfunction()

function(kkk)
    message(STATUS "function kkk working")
    foreach (arg IN LISTS ARGN)
        message(STATUS "Arg " ${arg})
    endforeach()
    # message(STATUS arg1)
endfunction()

# Add tests
add_subdirectory(tests)

set(SOURCE_FILES src/main.cpp src/version.hpp src/version.cpp)
add_executable(test ${SOURCE_FILES})