
#pragma once

// Include local modules
#include "../../config.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../containers/simulation_data.hpp"
#include "formulation_kernel_backend.hpp"
#include "../../hydrostatics.hpp"
#include "../../inout/input.hpp"
#include "../../inout/output.hpp"
#include "../../mesh/mesh_group.hpp"
#include "../../tools.hpp"


/*******************************************************/
/************** Declare Module Macros ******************/
/*******************************************************/
#define REDUCE_FO_ROOT( field_name, field_type, data_type )     \
    MPI_Reduce(                                                 \
                    this->sim_data->field_name,                 \
                    this->sim_data->field_name##_p0,            \
                    this->sim_data->field_type##_np,            \
                    mpi_##data_type,                            \
                    MPI_SUM,                                    \
                    this->mpi_config->proc_root,                \
                    MPI_COMM_WORLD                              \
                );                                              \

/*******************************************************/
/*********** Declare FrequencySolver Class *************/
/*******************************************************/
template<std::size_t N, int mode_pf>
class FrequencySolver
{
public:
    // Declare public attributes
    FormulationKernelBackend<NUM_GP, PF_OFF>*   kernel          = nullptr;
    Hydrostatics**                              hydrostatics    = nullptr;
    Input*                                      input           = nullptr;
    MeshGroup*                                  mesh_gp         = nullptr;
    MeshGroup*                                  mesh_fs_qtf_gp  = nullptr;
    MpiConfig*                                  mpi_config      = nullptr;
    Output*                                     output          = nullptr;
    SimulationData*                             sim_data        = nullptr;
    

    /**** Define constructors ****/
    FrequencySolver( ) = default;

    FrequencySolver( Input* input );

    /**** Declare destructor ****/
    ~FrequencySolver( );

    /**** Declare class methods ****/
    void calculate_first_order( );

    void calculate_hydrostatics( );

    void initialize_mesh_groups( );

    void initialize_output_system( );

};

// Include class definitions
#include "frequency_solver.txx"