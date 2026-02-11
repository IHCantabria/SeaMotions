# SeaMotions Installation Guide

This document collects the end-to-end steps to configure, build, test, and install SeaMotions across Linux, macOS, and Windows environments.

SeaMotions is based on Intel OneApi platform + CMake in order to have full cross platform compatibility and to ensure robustness and efficiency. So, first it will be necessary to download Intel OneAPI Base toolkit <https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html> and Intel OneAPI HPC toolkit (for the MPI) <https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html>

## Prerequisites
- **CMake:** 3.24 or newer.
- **C++ compiler:** with C++17 support
- **MPI implementation:** Intel MPI
- **Parallel HDF5:** 1.14+ with MPI IO enabled. See [PARALLEL_HDF5.md](PARALLEL_HDF5.md) for build instructions.
- **Python 3.9+:** with NumPy and h5py for optional database tooling.
- **Ninja:** (optional) for faster multi-platform builds.

## Clone the Repository

```bash
git clone https://github.com/SeaMotions/SeaMotions.git
cd SeaMotions
```

## Configure and Build (Linux/macOS)
1. Open a terminal and move to the root folder of SeaMotions project
2. Load Intel OneAPI environment in order to have all the libraries and dependencies available.
```bash
setenvars.sh
```
3. Ensure HDF5 library is available in you computer and in the system `PATH`
4. Configure build by using CMake
```bash
mkdir build
cmake -S . -B build \ 
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=icpx \ 
      -DPRECISION=DOUBLE
```
5. Build solution. Change directory to SeaMotions/build and build
solution by using CMake:
```bash
cmake --build build \ 
      --config Release \
      --target seamotions_freq seamotions_stab
```
or by using GNU Make system:
```bash
make
```

## Configure and Build (Windows)
1. Open an "Intel OneAPI command prompt for Intel 64 for Visual Studio"
2. Load Intel OneAPI environment in order to have all the libraries and dependencies available if they do not trigger automatically.
```bash
setenvars.sh
```
3. Ensure HDF5 library is available in you computer and in the system `PATH`
3. Configure build by using CMake:

```bash
cmake -S . -B build -G "Ninja" \ 
      -DCMAKE_BUILD_TYPE=Release \ 
      -DCMAKE_CXX_COMPILER=icpx \ 
      -DPRECISION=DOUBLE
```

```bash
cmake --build build --config Release --target seamotions_freq seamotions_stab
```

> **Note**: in any case HDF5 libraries are not detected during CMake configuration it might be usefull to use `-DHDF5_DIR=/path/to/your/hdf5/cmake` as CMake command line input argument during configuration stage.

## Install Binaries and Data
To install SeaMotions software move to `/path/to/SeaMotions/build` and use:
```bash
cmake --install . --prefix /path/to/bin/folder
```
Optionally, is it possible to do through the Make system in Linux
```bash
make install
```
or in Windows:
```bash
ninja install
```

The install step copies executables to `prefix/bin`, shared assets to `prefix/share/seamotions`, and headers for downstream integration. Update your `PATH` (and on Windows, the DLL search path) so the binaries and MPI/HDF5 runtimes are visible.

## Run the Test Suite

```bash
cmake --build build --target test
# or
cd build && ctest --output-on-failure
```