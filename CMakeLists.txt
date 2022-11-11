cmake_minimum_required(VERSION 3.9)
project(SEAMOTIONS CXX)

#################### Configure project ####################
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
set(VERSION_SHORT "${version_major}.${version_minor}.${version_patch}")

# Configure version file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/version.cpp.in
               ${CMAKE_CURRENT_SOURCE_DIR}/src/version.cpp)
set(version_file "${CMAKE_CURRENT_BINARY_DIR}/version.cpp")

message(STATUS ${VERSION_SHORT})
message(STATUS ${CMAKE_CURRENT_SOURCE_DIR})
message(STATUS ${CMAKE_CURRENT_BINARY_DIR})

set(SOURCE_FILES src/main.cpp src/version.hpp src/version.cpp)
add_executable(test ${SOURCE_FILES})