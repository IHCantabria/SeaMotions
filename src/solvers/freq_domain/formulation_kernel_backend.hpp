
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
#include "../../containers/matlin_group.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../containers/rad_diff_data.hpp"
#include "../../containers/simulation_data.hpp"
#include "../../green/farfield_integrator.hpp"
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
    /*** Declare local type aliases ***/
    using RDDQC = RadDiffData<cuscomplex, RDDQTFConfig>;

    /* Declare private variables */
    FarFieldIntegrator<QTFTypeE::QTF_DIFF_CODE> _qtf_diff_code_integrator;          // Far-field integrator for the QTF difference code term
    FarFieldIntegrator<QTFTypeE::QTF_SUM_CODE>  _qtf_sum_code_integrator;           // Far-field integrator for the QTF sum code term
    GWFcnsInterfaceT<N*N>                       _gwfcns_interf;                     // Wave part functor interface used for the integration over the panel by using Gauss Points (first order potential evaluation)
    GWFcnsInterfaceT<N*N>                       _gwfcns_interf_qb;                  // Wave part functor interface used for the integration over the panel by using Gauss Points (second order potential evaluation, Qb term)
    GWFcnsInterfaceT<N*N>                       _gwfcns_interf_qf;                  // Wave part functor interface used for the integration over the panel by using Gauss Points (second order potential evaluation, Qf term)
    GWFcnsInterfaceT<1>                         _gwfcns_interf_qfa;                 // Wave part functor interface used for the integration over the panel by using Gauss Points (second order potential evaluation, Qf term annuli)
    GWFcnsInterfaceT<N>                         _gwfcns_interf_wl;                  // Wave part functor interface used for the integration over the panel by using Gauss Points (second order potential evaluation for waterline points)
    Input*                                      _input                  = nullptr;  // Input system to have access to the case configuration
    int                                         _is_condition_number    = false;    // Switch to enable or disable the computation of the Condition number of the system matrixes for all the available formulations
    MeshGroup*                                  _mesh_gp                = nullptr;  // Mesh group describing the target case topology
    MpiConfig*                                  _mpi_config             = nullptr;  // Pointer to MpiConfig to have access to MPI configuration
    MLGCmpx*                                    _pf_gp                  = nullptr;  // Group of matrix data to storage Potential formulation data (Column-Major arranged data)
    MLGCmpx*                                    _pot_gp                 = nullptr;  // Group of matrix data to storage Potential matrix data (This is required for Source and Potential formulations) (Row-Major arranged data)
    cuscomplex*                                 _ppf_rhs                = nullptr;  // Partial potential formulation RHS. It storage the information at each vertical chunk of the potential matrix (see PF-RHS) depending of the processor
    SclCmpx*                                    _solver                 = nullptr;  // ScaLapack solver for complex numbers
    int                                         _steady_mat_type        = 0;        // Steady matrix type to be used in the formulation ( 0: REGULAR+Lf, 1: HF )
    MLGCmpx*                                    _sf_gp                  = nullptr;  // Group of matrix data to storage Source formulation data (Colum-Major arranged data)
    const std::vector<RDDQC*>*                  _qtf_annuli_fields      = nullptr;  // Optional annulus fields for QTF integration
    WaveDispersionSO                            _wdso                   ;           // Wave dispersion object for second order calculations

    /* Declare private methods */
    template<FreqRegimeE freq_regime>
    void _build_steady_matrixes( 
                                    void 
                                );

    void _build_rhs( 
                                    cusfloat w
                    );

    void _build_rhs_so( 
                                    QTFTypeE                        qtf_type,
                                    std::size_t                     freq_i,
                                    std::size_t                     freq_j,
                                    cuscomplex*                     raos_hist,
                                    RDDQC*                          qtf_body_fields,
                                    RDDQC*                          qtf_body_fields_wl,
                                    RDDQC*                          qtf_fs_fields,
                                    RDDQC*                          qtf_fs_fields_wl
                    );

    void _build_rhs_so( 
                                    QTFTypeE                        qtf_type,
                                    std::size_t                     freq_i,
                                    std::size_t                     freq_j,
                                    cuscomplex*                     raos_hist,
                                    RDDQC*                          qtf_body_fields,
                                    RDDQC*                          qtf_body_fields_wl,
                                    RDDQC*                          qtf_fs_fields,
                                    RDDQC*                          qtf_fs_fields_wl,
                                    const std::vector<RDDQC*>*      qtf_annuli_fields
                    );

    template<FreqRegimeE freq_regime>
    void _build_wave_matrixes( 
                                    cusfloat w
                            );

    // template<FreqRegimeE freq_regime>
    // void _build_wave_matrixes_2( 
    //                                 cusfloat w
    //                             );

    template<FreqRegimeE freq_regime>
    void _build_wave_matrixes_so( 
                                    QTFTypeE                qtf_type,
                                    std::size_t             freq_i,
                                    std::size_t             freq_j,
                                    cuscomplex              raos_hist,
                                    RDDQC*                  qtf_body_fields,
                                    RDDQC*                  qtf_body_fields_wl,
                                    RDDQC*                  qtf_fs_fields,
                                    RDDQC*                  qtf_fs_fields_wl,
                                    WaveDispersionSO*       wd
                                );

    void _process_far_field(
                                    QTFTypeE                qtf_type,
                                    std::size_t             freq_i,
                                    std::size_t             freq_j,
                                    cuscomplex*             qb_rhs
                            );

    template<int GP, int mode_fslid, typename GwfInterf, typename RhsFunc, typename IntegrateFunc>
    void _process_qtf_rhs_panels(
                                    RDDQC*                  panel_fields,
                                    GwfInterf&              gwf_interf_local,
                                    RhsFunc&&               rhs_func,
                                    IntegrateFunc&&         integrate,
                                    bool                    mode_wl,
                                    cusfloat*               field_points_d,
                                    cusfloat*               field_points_x,
                                    cusfloat*               field_points_y,
                                    cusfloat*               field_points_z,
                                    cuscomplex*             qb_rhs
                                );

    template<typename RhsFunc>
    void _process_qtf_rhs_annuli(
                                    RDDQC*                  annulus_fields,
                                    const cusfloat*         field_weights,
                                    GWFcnsInterfaceT<1>&    gwf_interf_local,
                                    RhsFunc&&               rhs_func,
                                    cusfloat*               field_points_d,
                                    cusfloat*               field_points_x,
                                    cusfloat*               field_points_y,
                                    cusfloat*               field_points_z,
                                    cuscomplex*             qb_rhs
                                );

    void _initialize(                
                                    void 
                    );
    

public:
    /* Declare public attributes */
    cusfloat    exec_time_build_steady  = 0.0;
    cusfloat    exec_time_build_wave    = 0.0;
    cusfloat    exec_time_solve_pf      = 0.0;
    cusfloat    exec_time_solve_sf      = 0.0;
    int         ipm_cols_np             = 0;
    int         ipm_sc                  = 0;
    int         ipm_ed                  = 0;

    /* Define class constructors */
    FormulationKernelBackend( 
                                        Input*                              input, 
                                        MpiConfig*                          mpi_config, 
                                        MeshGroup*                          mesh_gp 
                            );

    ~FormulationKernelBackend( );

    /* Declare class public methods */
    template<typename RDDConfig>
    void    compute_fields(
                                        std::size_t                         freq_index,
                                        cusfloat                            ang_freq,
                                        cuscomplex*                         raos,
                                        RadDiffData<cuscomplex, RDDConfig>* rad_diff_data
                            );

    void    configure_second_order( 
                                        cusfloat                            partition_circle    
                                    );

    int     size(
                                        void
                );
    
    template<FreqRegimeE freq_regime>
    void    solve(             
                                        cusfloat w 
                );

    template<QTFTypeE qtf_type>
    void    solve_so( 
                                        std::size_t                 freq_i,
                                        std::size_t                 freq_j,
                                        cuscomplex*                 raos_hist,
                                        RDDQC*                      qtf_body_fields,
                                        RDDQC*                      qtf_body_fields_wl,
                                        RDDQC*                      qtf_fs_fields,
                                        RDDQC*                      qtf_fs_fields_wl,
                                        const std::vector<RDDQC*>*  qtf_annuli_fields
                    );

    void    update_results( 
                                        SimulationData* sim_data 
                        );

    template<QTFTypeE qtf_type>
    void    update_results_so( 
                                        SimulationData* sim_data 
                                );

    void    set_qtf_annuli_fields(
                                        const std::vector<RDDQC*>* annuli_fields
                                );

};

// Include class method definitions
#include "formulation_kernel_backend.txx"