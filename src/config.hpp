
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

    #ifdef _MPI_BUILD
    #include "mpi.h"
    typedef MPI_FLOAT mpi_float;
    #endif

constexpr int FLOATING_PRECISION = 32;
constexpr cusfloat EPS_PRECISION = 1e-6;
constexpr cusfloat EPS_PRECISION_ORDER = -6;
#else
typedef double cusfloat;
typedef std::complex<double> cuscomplex;

    #ifdef _MPI_BUILD
    #include "mpi.h"
    typedef MPI_DOUBLE mpi_double;
    #endif

constexpr int FLOATING_PRECISION = 32;
constexpr cusfloat EPS_PRECISION = 1e-14;
constexpr cusfloat EPS_PRECISION_ORDER = -14;
#endif

#endif