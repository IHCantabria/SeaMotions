
#ifndef __config_hpp
#define __config_hpp

#include <complex>
#include <vector>

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

typedef std::vector<size_t> vector_st;

#define MEMALINGR alignas(32)

constexpr int       G_ON                    = 1;                // Flag used as template argument to SET the calculation of the potential green function
constexpr int       G_OFF                   = 0;                // Flag used as template argument to NOT SET the calculation of the potential green function
constexpr int       DGDR_ON                 = 1;                // Flag used as template argument to SET the calculation of the derivative of the potential green function with respect to the horizontal radius
constexpr int       DGDR_OFF                = 0;                // Flag used as template argument to NOT SET the calculation of the derivative of the potential green function with respect to the horizontal radius
constexpr int       DGDZ_ON                 = 1;                // Flag used as template argument to SET the calculation of the derivative of the potential green function with respect to the vertical coordinate
constexpr int       DGDZ_OFF                = 0;                // Flag used as template argument to NOT SET the calculation of the derivative of the potential green function with respect to the vertical coordinate
constexpr int       DIFFRAC_PANEL_CODE      = 0;                // Flag used to mark the panel as diffraction
constexpr cusfloat  FIELD_POINT_LOCAL_TOL   = 1E-2;             // Tolerance distance to assume a field point to be in the source point.
constexpr cusfloat  FS_SEL_THR              = 1e-3;             // Free surface selection threshold. All the points which z is below this threshold are considered to be in the free surface
constexpr int       FSLID_ON                = 1;                // Flag used as template argument to REMOVE log singularity for free surface panels
constexpr int       FSLID_OFF               = 0;                // Flag used as template argument to NOT REMOVE log singularity for free surface panels
constexpr int       LID_PANEL_CODE          = 1;                // Flag used to mark the panel as lid
constexpr int       NUM_GP                  = 2;                // Number of Gauss Points used for numerical integration
constexpr int       NUM_GP2                 = NUM_GP*NUM_GP;    // Squared number of gauss points ( just for convenience along the code )
constexpr int       NUM_GP3                 = NUM_GP2*NUM_GP;   // Cubic number of gauss points ( just for convenience along the code )
constexpr int       NUM_KN                  = 30;               // Maximum number of imaginary wave number roots used in the john series expansion
constexpr int       PF_ON                   = 1;                // Flag used as template argument to SET the software to use the potential formulation to calculate the potential
constexpr int       PF_OFF                  = 0;                // Flag used as template argument to NOT SET the software to use the potential formulation to calculate the potential
constexpr int       MPI_ROOT_PROC_ID        = 0;                // ID for the root process when using the MPI environment
constexpr int       STATIC_LOOP_ON          = 1;                // Flag used as template argument to SET that the bounds of the internal loops of the calling function are known at compile time
constexpr int       STATIC_LOOP_OFF         = 0;                // Flag used as template argument to NOT SET that the bounds of the internal loops of the calling function are not known at compile time
constexpr cusfloat  ZEROTH_EPS              = 1E-14;

#endif