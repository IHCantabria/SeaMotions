
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


// Define maximum pipe out message width to align 
// elapsed times and to have a more clearer view 
// of the output
constexpr int FSOL_MSG_WIDTH = 50;


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
private:
    /**** Declare class private methods ****/
    void _calculate_global_static_matrixes( );

    void _calculate_hydrostatics( );

    void _generate_formulation_kernel( );

    void _initialize_mesh_groups( );

    void _initialize_output_system( );

public:
    // Declare public attributes
    FormulationKernelBackend<NUM_GP, PF_OFF>*   kernel          = nullptr;  // Kernel of the formulation. It solves the raddiation diffraction problem. It can be CPU or GPU based.
    Hydrostatics**                              hydrostatics    = nullptr;  // 
    Input*                                      input           = nullptr;  // 
    MeshGroup*                                  mesh_gp         = nullptr;  // 
    MeshGroup*                                  mesh_fs_qtf_gp  = nullptr;  // 
    MpiConfig*                                  mpi_config      = nullptr;  // 
    Output*                                     output          = nullptr;  // 
    SimulationData*                             sim_data        = nullptr;  // 

    /**** Define constructors ****/

    // Default Constructor
    FrequencySolver( ) = default;

    // Normal constructor. Builds case configuration and initialization
    // based on the contents of the input system.
    FrequencySolver( Input* input, MpiConfig* mpi_config_in );

    /**** Declare destructor ****/
    ~FrequencySolver( );

    /* Declare class public methods */
    void calculate_first_order( );

};

// Include class definitions
#include "frequency_solver.txx"