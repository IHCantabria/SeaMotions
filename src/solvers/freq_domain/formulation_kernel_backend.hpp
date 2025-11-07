
#pragma once

// Include local modules
#include "../../containers/matlin_group.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../containers/simulation_data.hpp"
#include "../../green/source.hpp"
#include "../../inout/input.hpp"
#include "../../interfaces/gwfcns_interface_t.hpp"
#include "../../math/integration.hpp"
#include "../../math/scalapack_solver.hpp"
#include "../../mesh/mesh_group.hpp"
#include "../../static_tools.hpp"
#include "../../waves/waves_common.hpp"

// Declare auxiliary macros
#define COL_MAJOR_INDEX( index, row_count, col_count, num_rows_local ) index_cm = col_count * num_rows_local + row_count;
#define ROW_MAJOR_INDEX( index, row_count, col_count, num_cols_local ) index_rm = row_count * num_cols_local + col_count;

// Declare function
template<std::size_t N, int mode_pf>
struct FormulationKernelBackend
{
private:
    // Declare private variables
    GWFcnsInterfaceT<N*N>   _gwfcns_interf;                     // Wave part functor interface used for the integration over the panel by using Gauss Points
    Input*                  _input                  = nullptr;  // Input system to have access to the case configuration
    int                     _is_condition_number    = false;    // Switch to enable or disable the computation of the Condition number of the system matrixes for all the available formulations
    MeshGroup*              _mesh_gp                = nullptr;  // Mesh group describing the target case topology
    MpiConfig*              _mpi_config             = nullptr;  // Pointer to MpiConfig to have access to MPI configuration
    MLGCmpx*                _pf_gp                  = nullptr;  // Group of matrix data to storage Potential formulation data (Column-Major arranged data)
    MLGCmpx*                _pot_gp                 = nullptr;  // Group of matrix data to storage Potential matrix data (This is required for Source and Potential formulations) (Row-Major arranged data)
    cuscomplex*             _ppf_rhs                = nullptr;  // Partial potential formulation RHS. It storage the information at each vertical chunk of the potential matrix (see PF-RHS) depending of the processor
    SclCmpx*                _solver                 = nullptr;  // ScaLapack solver for complex numbers
    MLGCmpx*                _sf_gp                  = nullptr;  // Group of matrix data to storage Source formulation data (Colum-Major arranged data)

    // Declare private methods
    void _build_steady_matrixes( 
                                    void 
                                );

    void _build_rhs( 
                                    cusfloat w
                    );

    void _build_wave_matrixes( 
                                    cusfloat w
                            );

    void _initialize(                
                                    void 
                    );
    

public:
    // Declare public attributes
    cusfloat    exec_time_build_steady  = 0.0;
    cusfloat    exec_time_build_wave    = 0.0;
    cusfloat    exec_time_solve_pf      = 0.0;
    cusfloat    exec_time_solve_sf      = 0.0;
    int         ipm_cols_np             = 0;
    int         ipm_sc                  = 0;
    int         ipm_ed                  = 0;

    // Define class constructors
    FormulationKernelBackend( Input* input, MpiConfig* mpi_config, MeshGroup* mesh_gp );

    ~FormulationKernelBackend( );

    // Define class public methods
    void solve(             
                            cusfloat w 
                );

    void update_results( 
                            SimulationData* sim_data 
                        );

};

// Include class method definitions
#include "formulation_kernel_backend.txx"