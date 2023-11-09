
#ifndef __config_hpp
#define __config_hpp

#include <complex>

// Check if the program has been build in debug mode
#ifdef DEBUG_BUILD
constexpr bool      _DEBUG_BUILD        = true;
#else
constexpr bool      _DEBUG_BUILD        = false;
#endif

#ifdef SIMPLE_PREC
typedef float cusfloat;
typedef std::complex<float> cuscomplex;
#define MKL_Complex8 std::complex<float>

    #ifdef MPI_BUILD
    #include "mpi.h"
    #define mpi_cusfloat MPI_FLOAT
    #define mpi_cuscomplex MPI_COMPLEX
    #endif

    #ifdef _HDF5_BUILD
    #include "H5Cpp.h"
    #define cusfloat_h5 H5::PredType::NATIVE_FLOAT
    #endif

constexpr int FLOATING_PRECISION = 32;
constexpr cusfloat EPS_PRECISION = 1e-6;
constexpr cusfloat EPS_PRECISION_ORDER = -6;
#else
typedef double cusfloat;
typedef std::complex<double> cuscomplex;
#define MKL_Complex16 std::complex<double>

    #ifdef MPI_BUILD
    #include "mpi.h"
    #define mpi_cusfloat MPI_DOUBLE
    #define mpi_cuscomplex MPI_DOUBLE_COMPLEX
    #endif

    #ifdef _HDF5_BUILD
    #include "H5Cpp.h"
    #define cusfloat_h5 H5::PredType::NATIVE_DOUBLE
    #endif

constexpr int FLOATING_PRECISION = 32;
constexpr cusfloat EPS_PRECISION = 1e-14;
constexpr cusfloat EPS_PRECISION_ORDER = -14;
#endif

constexpr cusfloat ZEROTH_EPS = 1e-302;

#endif