# Parallel HDF5 Build Guide

SeaMotions reads and writes large hydrodynamic databases through HDF5. To avoid serialization bottlenecks you need an HDF5 build compiled with MPI collective IO enabled. The recipes below create relocatable installations that you can reference via `HDF5_ROOT` when configuring SeaMotions.

As SeaMotions is designed to work with Intel OneAPI platform it would be necessary to compile HDF5 using the Intel compilers and environment.

SeaMotions has been testest by using HDF5 1.14 version <https://github.com/HDFGroup/hdf5/releases>. It is recommendable to use this version to avoid any incompatiblity errors with HDF5 API.


## CMake Presets
CMake Presets are the easiest way to pre-configure the build including toolchains and specific configurations. HDF5 does not come with a specific preset for Intel compilers in combination with parallel HDF5 compilation. It will be necessary to create a new item in `configurationPresets` section inside file CMakePresets.json:
```bash
{
    "name": "ci-StdShar-IntelMPI",
    "description": "Intel MPI Standard Config for x64 (Release)",
    "inherits": "ci-StdShar-Intel", 
    "displayName": "Windows Intel MPI Parallel",
    
    "cacheVariables": 
    {
    "HDF5_ENABLE_PARALLEL": "ON",
    "HDF5_BUILD_CPP_LIB": "OFF"
    }
}
```
a new one inside `buildPresets` section inside file CMakePresets.json
```bash
{
    "name": "ci-StdShar-IntelMPI",
    "description": "Intel MPI Standard Build for x64 (Release)",
    "configurePreset": "ci-StdShar-IntelMPI",
    "verbose": true,
    "inherits": [
    "ci-x64-Release-Intel"
    ]
}
```
in the case of linux compilation it may be necessary to force the name of the compilers by adding these entries in `cacheVariables` subdictionary:
```
"CMAKE_C_COMPILER": "icx",
"CMAKE_CXX_COMPILER": "icpx",
"CMAKE_Fortran_COMPILER": "ifort"
```
a new one inside `testPresets` in CMakePresets.json file:
```bash
{
    "name": "ci-StdShar-IntelMPI",
    "configurePreset": "ci-StdShar-IntelMPI",
    "inherits": [
    "ci-x64-Release-Intel"
    ]
}
```
a new one inside `packagePresets` in CMakePresets.json file:
```bash
{
    "name": "ci-StdShar-IntelMPI",
    "configurePreset": "ci-StdShar-IntelMPI",
    "inherits": "ci-x64-Release-Intel"
}
```
a new one inside `workflowPresets` in CMakePresets.json file:
```bash
{
    "name": "ci-StdShar-IntelMPI",
    "steps": [
    {"type": "configure", "name": "ci-StdShar-IntelMPI"},
    {"type": "build", "name": "ci-StdShar-IntelMPI"},
    {"type": "test", "name": "ci-StdShar-IntelMPI"},
    {"type": "package", "name": "ci-StdShar-IntelMPI"}
    ]
}
```

## Configuration (Linux/macOS)
1. Open a terminal and move to the root folder of HDF5 where CMakePresets.json file is located.
2. Load Intel OneAPI environment in order to have all the libraries and dependencies available.
```bash
setenvars.sh
```
3. Configure build by using CMake presets:
```bash
cmake -G "Unix Makefiles" --presets ci-StdShar-IntelMPI --fresh
```
4. Move one level above with the terminal where it should be located a folder with name `build114` if you are trying to compile `HDF5-1.14.X`
5. Move into build folder and compile
```bash
cmake --build .
```
or using GNU Make system
```bash
make
```

## Configuration (Windows)
1. Open an "Intel OneAPI command prompt for Intel 64 for Visual Studio"
2. Load Intel OneAPI environment in order to have all the libraries and dependencies available if they do not trigger automatically.
```bash
setenvars.sh
```
3. Configure build by using CMake presets:
```bash
cmake --presets ci-StdShar-IntelMPI --fresh
```
4. Move one level above with the terminal where it should be located a folder with name `build114` if you are trying to compile `HDF5-1.14.X`
5. Move into build folder and compile
```bash
cmake --build .
```
or using Ninja system
```bash
ninja
``` 

## Installation
By using CMake and from the build folder use:
```bash
cmake --install .
```
or in linux
```bash
make install
```
or in Windows
```bash
ninja install
```
This process should create a folder at the same level than `build114` that should be called 
`install114` where binaries, includes, libs and cmake auxiliary scripts should located.
