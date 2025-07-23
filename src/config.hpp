
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
    #define int_h5 H5::PredType::NATIVE_INT
    #endif

constexpr int       FLOATING_PRECISION  = 32;
constexpr cusfloat  EPS_PRECISION       = 1e-6;
constexpr cusfloat  EPS_PRECISION_ORDER = -6;
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
    #define int_h5 H5::PredType::NATIVE_INT
    #endif

constexpr int       FLOATING_PRECISION  = 64;
constexpr cusfloat  EPS_PRECISION       = 1e-14;
constexpr cusfloat  EPS_PRECISION_ORDER = -14;
#endif

#define MEMALINGR alignas(32)

constexpr cusfloat  ZEROTH_EPS              = 1E-14;
constexpr cusfloat  FIELD_POINT_LOCAL_TOL   = 1E-2;
constexpr int       NUM_GP                  = 4;                // Number of Gauss Points used for numerical integration
constexpr int       NUM_GP2                 = NUM_GP*NUM_GP;    // Squared number of gauss points ( just for convenience along the code )
constexpr int       NUM_GP3                 = NUM_GP2*NUM_GP;   // Cubic number of gauss points ( just for convenience along the code )
constexpr int       NUM_KN                  = 30;               // Maximum number of imaginary wave number roots used in the john series expansion

#endif