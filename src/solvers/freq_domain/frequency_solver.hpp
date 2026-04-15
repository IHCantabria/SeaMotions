
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
#include "../../containers/field_mesh_data.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../containers/morison_element.hpp"
#include "../../containers/rad_diff_data.hpp"
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
    /*** Declare local type aliases ***/
    using RDDFC     = RadDiffData<cuscomplex, RDDFDConfig>;
    using RDDQC     = RadDiffData<cuscomplex, RDDQTFConfig>;
    using FMD       = FieldMeshData<cuscomplex, RDDFDConfig>;
    using vec_fmd   = std::vector<FMD*>;
    using vec_med   = std::vector<MorisonElement*>;

    /**** Declare class private attributes ****/
    cut::CusTensor<cuscomplex>          _drag_force_aux     ;           // Auxiliary storage for the calculation of drag forces in slender elements
    vec_fmd                             _field_points       ;           // Storage of field points for radiation and diffraction calculations
    cut::CusTensor<cuscomplex>          _inertial_force_aux ;           // Auxiliary storage for the calculation of inertial forces in slender elements
    std::vector<vec_med>                _morison_elements   ;           // Storage of Morison elements for the calculation of inertial and drag forces in slender elements
    RDDQC*                              _qtf_wl_fields      = nullptr;  // Storage of waterline field points for QTF calculations
    RDDQC*                              _qtf_bern_fields    = nullptr;  // Storage of velocity field for the calculation of bernoulli term in QTFs
    RDDQC*                              _qtf_fs_fields      = nullptr;  // Storage of field points over the QTF free surface lid
    RDDQC*                              _qtf_fs_wl_fields   = nullptr;  // Storage of field points over the QTF free surface lid perimeter
    std::vector<Mesh*>                  _qtf_annuli_meshes  ;           // Stack of annulus meshes for QTF field sampling
    std::vector<RDDQC*>                 _qtf_annuli_fields  ;           // RadDiffData stack for annulus field points
    std::vector<std::vector<cusfloat>>  _qtf_annuli_weights;            // Weights per annulus field point set

    /**** Declare class private methods ****/
    void _calculate_field_points_values( 
                                                    std::size_t freq_index,
                                                    cusfloat    ang_freq 
                                        );

    template<FreqRegimeE freq_regime>
    void _calculate_first_order_coeffs( 
                                                    std::size_t freq_index,
                                                    cusfloat    ang_freq
                                        );

    void _calculate_mean_drift( 
                                                    std::size_t freq_index,
                                                    bool        is_multi_head 
                                );

    void _calculate_morison_forces( 
                                                    cusfloat    ang_freq 
                                    );

    void _calculate_global_static_matrixes( );

    void _calculate_hydrostatics( );

    template<QTFTypeE qtf_type>
    void _calculate_quadratic_terms( 
                                                    std::size_t freq_index_i,
                                                    std::size_t freq_index_j,
                                                    bool        is_multi_head
                                    );

    template<QTFTypeE qtf_type>
    void _calculate_and_distribute_qtf(
                                                    std::size_t freq_index_i,
                                                    std::size_t freq_index_j
                                        );


    void _check_system_mesh(  
                                                    void 
                                    );

    void _generate_formulation_kernel( );

    void _initialize_field_data( );

    void _initialize_mesh_groups( );

    void _initialize_morison_elements( );

    void _initialize_output_system( );

    void _prepare_results_folder( );

    void _save_qtf( void ) const;

public:
    // Declare public attributes
    FormulationKernelBackend<NUM_GP, mode_pf>*  kernel          = nullptr;  // Kernel of the formulation. It solves the raddiation diffraction problem. It can be CPU or GPU based.
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

    void calculate_second_order( );

};

// Include class definitions
#include "frequency_solver.txx"