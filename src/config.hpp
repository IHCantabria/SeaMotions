
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

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
    #include "hdf5.h"
    #define cusfloat_h5 H5T_NATIVE_FLOAT
    #define int_h5 H5T_NATIVE_INT
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
    #include "hdf5.h"
    #define cusfloat_h5 H5T_NATIVE_DOUBLE
    #define int_h5 H5T_NATIVE_INT
    #endif

constexpr int       FLOATING_PRECISION  = 64;
constexpr cusfloat  EPS_PRECISION       = 1e-14;
constexpr cusfloat  EPS_PRECISION_ORDER = -14;
#endif

typedef std::vector<size_t> vector_st;

#define MEMALINGR alignas(32)

constexpr int           G_ON                    = 1;                // Flag used as template argument to SET the calculation of the potential green function
constexpr int           G_OFF                   = 0;                // Flag used as template argument to NOT SET the calculation of the potential green function
constexpr int           DGDC_ON                 = 1;                // Flag used as template argument to SET the calculation of the gradient of the potential green function
constexpr int           DGDC_OFF                = 0;                // Flag used as template argument to NOT SET the calculation of the gradient of the potential green function
constexpr int           DGDN_ON                 = 1;                // Flag used as template argument to SET the calculation of the normal derivative of the potential green function
constexpr int           DGDN_OFF                = 0;                // Flag used as template argument to NOT SET the calculation of the normal derivative of the potential green function
constexpr int           DGDR_ON                 = 1;                // Flag used as template argument to SET the calculation of the derivative of the potential green function with respect to the horizontal radius
constexpr int           DGDR_OFF                = 0;                // Flag used as template argument to NOT SET the calculation of the derivative of the potential green function with respect to the horizontal radius
constexpr int           DGDZ_ON                 = 1;                // Flag used as template argument to SET the calculation of the derivative of the potential green function with respect to the vertical coordinate
constexpr int           DGDZ_OFF                = 0;                // Flag used as template argument to NOT SET the calculation of the derivative of the potential green function with respect to the vertical coordinate
constexpr cusfloat      FIELD_POINT_LOCAL_TOL   = 1E-2;             // Tolerance distance to assume a field point to be in the source point.
constexpr cusfloat      FS_SEL_THR              = 1e-3;             // Free surface selection threshold. All the points which z is below this threshold are considered to be in the free surface
constexpr int           FSLID_ON                = 1;                // Flag used as template argument to REMOVE log singularity for free surface panels
constexpr int           FSLID_OFF               = 0;                // Flag used as template argument to NOT REMOVE log singularity for free surface panels
constexpr cusfloat      MIN_PANEL_AREA          = 1e-5;             // Minimum panel area to be considered a workable panel during mesh refinement process
constexpr int           NUM_GP                  = 2;                // Number of Gauss Points used for numerical integration
constexpr int           NUM_GP2                 = NUM_GP*NUM_GP;    // Squared number of gauss points ( just for convenience along the code )
constexpr int           NUM_GP3                 = NUM_GP2*NUM_GP;   // Cubic number of gauss points ( just for convenience along the code )
constexpr int           NUM_KN                  = 30;               // Maximum number of imaginary wave number roots used in the john series expansion
constexpr int           ON                      = 1;                // Generic flag to SET an option
constexpr int           OFF                     = 0;                // Generic flag to NOT SET an option
constexpr int           PANEL_LOC_AW            = 1;                // Flag to mark that panel location is above water
constexpr int           PANEL_LOC_FS            = 0;                // Flag to mark that panel location is on the free surface
constexpr int           PANEL_LOC_UW            = -1;               // Flag to mark that panel location is under water
constexpr int           PF_ON                   = 1;                // Flag used as template argument to SET the software to use the potential formulation to calculate the potential
constexpr int           PF_OFF                  = 0;                // Flag used as template argument to NOT SET the software to use the potential formulation to calculate the potential
constexpr std::size_t   QTF_FAR_N               = 30;               // Number of coefficients for the expansion series of the far field second order potential used in the QTF calculation
constexpr int           MPI_ROOT_PROC_ID        = 0;                // ID for the root process when using the MPI environment
constexpr int           STATIC_LOOP_ON          = 1;                // Flag used as template argument to SET that the bounds of the internal loops of the calling function are known at compile time
constexpr int           STATIC_LOOP_OFF         = 0;                // Flag used as template argument to NOT SET that the bounds of the internal loops of the calling function are not known at compile time
constexpr cusfloat      W_ASYMPT_HIGH           = 100.0;            // High frequency asymptotic regime threshold
constexpr cusfloat      W_ASYMPT_LOW            = 0.0001;           // Low frequency asymptotic regime threshold
constexpr cusfloat      ZEROTH_EPS              = 1E-14;

// Common folder names used across the codebase
constexpr const char* MESH_FOLDER_NAME                  = "mesh";
constexpr const char* RESULTS_FOLDER_NAME               = "1_results";
constexpr const char* RESULTS_MESH_FOLDER_NAME          = "0_mesh";
constexpr const char* RESULTS_FIELD_POINTS_FOLDER_NAME  = "1_field_points";
constexpr const char* RESULTS_PARAVIEW_FOLDER_NAME      = "2_paraview";


/************************************************************/
/****** Define class enums to be used along the code  *******/
/************************************************************/

// Data type enum. Used to classify the boundary subtypes to be calculated
enum class BoundarySubtypeE: int
{
    DIFFRAC,    // Diffraction: referring to the active zones for the analysis that are submerged and not in the free surface
    WL,         // Waterline: referring to the zone of the floating body that is on the intersection between the submerged part and the free surface
    PC          // Partition Circle: referring to the zone of the floating body that is the outside boundary of the free surface circle partition used for second order potential calculation.
};

// Data type enum. Used to classify the data type to be calculated.
enum class ComplexDataTypeE: int
{
    COMPLEX,        // Complex data type
    IMAGINARY,      // Imaginary part of the data type
    MAGNITUDE,      // Magnitude of the data type
    PHASE,          // Phase of the data type
    REAL            // Real part of the data type
};

// Free regime enum. Used to classify the frequency regime
// in frequency solver.
enum class FreqRegimeE: int
{
    REGULAR,        // Regular frequency regime
    ASYMPT_LOW,     // Low frequency asymptotic regime
    ASYMPT_HIGH     // High frequency asymptotic regime
};


// Field types enum. Used to classify the physical field type to be calculated.
enum class FieldTypeE: int
{
    POTENTIAL,                  // Potential field type
    PRESSURE,                   // Pressure field type
    RELATIVE_WAVE_ELEVATION,    // Relative wave elevation field type
    VELOCITY,                   // Velocity field type
    VELOCITY_DN,                // Velocity normal derivative field type
    VELOCITY_X,                 // Velocity X field type
    VELOCITY_Y,                 // Velocity Y field type
    VELOCITY_Z,                 // Velocity Z field type
    WAVE_ELEVATION              // Wave elevation field type
};

// Field component enum. Used to classify the field component to be calculated.
enum class FieldComponentE: int
{
    DIFFRACTED,     // Diffracted field component
    INCIDENT,       // Incident field component
    RADIATED,       // Radiated field component
    SCATTERED,      // Scattered field component
    TOTAL           // Total field component
};


// Panel Type enum. Used to classify panel types according to its 
// functionality
enum class PanelTypeE: int
{
    NONE,           // No type. Used for default empty states.
    DIFFRAC,        // Diffraction-Radation panel. Used as BC type for floating bodies
    INT_LID,        // Internal lid panel. Used as BC type for internal free surface in floating bodies to suppress spurious irregular frequencies.
    EXT_LID,        // External lid panel. Used as BC type for external free surface zones where it is necessary damp stationary waves that arise due to the lack of viscosity in the potential flow model.
    QTF_LID,        // QTF lid panel. Used as BC type to impose additional equations in the near field region of the floating bodies so it is possible to solve the second order potential.
    FIELD_POINT     // Field point panel. Used as type for field points where the potential and its derivatives are calculated.
};


// QTF type enum. Used to classify the type of second order potential to be used
enum class QTFSOModelE: int
{
    NONE        = -1,   // No second order potential model. Used for default empty states.
    PINKSTER    = 0,    // Second order potential model based on Pinkster's theory.
    INDIRECT    = 1,    // Second order potential model based on indirect method proposed by X.Chen.
    DIRECT      = 2     // Second order potential model based on direct method used by J.N.Newman
};


// QTF type enum. Used to classify the type of QTF
// term to be calculated.
enum class QTFTypeE: int
{
    QTF_DIFF_CODE      = 0,    // QTF difference term
    QTF_SUM_CODE       = 1     // QTF sum term
};

// Recalculation type enum. Switch used to recalculate steady part during the frequency loop or not.
enum class RecalcSteadyE: int
{
    OFF   = 0,    // No recalculation of the steady part during the frequency loop
    ON    = 1     // Recalculation of the steady part during the frequency loop
};