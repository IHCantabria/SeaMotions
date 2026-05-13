
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

// Include general usage libraries
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>

// Include local module
#include "../../containers/logger.hpp"
#include "../../containers/mpi_timer.hpp"
#include "diffraction.hpp"
#include "frequency_solver.hpp"
#include "froude_krylov.hpp"
#include "global_static_matrix.hpp"
#include "../../green/kochin.hpp"
#include "hydromechanics.hpp"
#include "../../hydrostatics.hpp"
#include "../../math/integration.hpp"
#include "../../math/math_tools.hpp"
#include "../../mesh/mesh.hpp"
#include "qtf_indirect_method.hpp"
#include "qtf.hpp"
#include "raos.hpp"


/*****************************************************************/
/************** Define Module auxiliar functions  ****************/
/*****************************************************************/
inline constexpr void _field_point_index( 
                                            bool            is_multi_head,
                                            std::size_t     heads_np,
                                            std::size_t     fp_np,
                                            std::size_t     freq_index_i,
                                            std::size_t     freq_index_j,
                                            std::size_t     ih1,
                                            std::size_t     ih2,
                                            std::size_t     k,
                                            std::size_t&    idx1_i,
                                            std::size_t&    idx1_j
                                        )
{
    idx1_i =  freq_index_i * ( heads_np * fp_np ) + ih1 * fp_np + k;
    if ( is_multi_head )
    {
        idx1_j = freq_index_j * ( heads_np * fp_np ) + ih2 * fp_np + k;
    }
    else
    {
        idx1_j = freq_index_j * ( heads_np * fp_np ) + ih1 * fp_np + k;
    }
}


inline constexpr std::size_t _qtf_index_offset(
                                                    bool        is_multi_head,
                                                    std::size_t bodies_np,
                                                    std::size_t heads_np,
                                                    std::size_t dofs_np,
                                                    std::size_t ih1,
                                                    std::size_t ih2,
                                                    std::size_t body_index
                                                )
{
    std::size_t offset = 0;
    if ( is_multi_head )
    {
        offset  = (
                        ih1 * ( dofs_np * bodies_np * heads_np )
                        +
                        ih2 * ( dofs_np * bodies_np )
                        +
                        body_index * dofs_np
                    );
    }
    else
    {
        offset  = ih1 * ( dofs_np * bodies_np ) + body_index * dofs_np;
    }

    return offset;
}


inline constexpr void _rao_offset( 
                                    bool            is_multi_head,
                                    std::size_t     bodies_np,
                                    std::size_t     dofs_np,
                                    std::size_t     ih1,
                                    std::size_t     ih2,
                                    std::size_t     j,
                                    std::size_t&    idx1_i,
                                    std::size_t&    idx1_j
                                )
{
    idx1_i =  ih1 * ( dofs_np * bodies_np ) + dofs_np * j;
    if ( is_multi_head )
    {
        idx1_j = ih2 * ( dofs_np * bodies_np ) + dofs_np * j;
    }
    else
    {
        idx1_j = ih1 * ( dofs_np * bodies_np ) + dofs_np * j;
    }
}


inline void _store_indirect_kochin_data(
                                            Input*          input,
                                            MeshGroup*      mesh_gp,
                                            std::size_t     freq_index,
                                            cusfloat        ang_freq,
                                            cuscomplex*     intensities,
                                            cuscomplex*     raos,
                                            SimulationData* sim_data
                                        )
{
    const std::size_t source_nodes_np      = static_cast<std::size_t>( mesh_gp->source_nodes_tnp );
    const std::size_t bodies_np            = static_cast<std::size_t>( input->bodies_np );
    const std::size_t dofs_np              = static_cast<std::size_t>( input->dofs_np );
    const std::size_t heads_np             = static_cast<std::size_t>( input->heads_np );
    const std::size_t kochin_heads_stride  = bodies_np * static_cast<std::size_t>( input->kochin_np );
    const std::size_t kochin_rad_stride    = bodies_np * static_cast<std::size_t>( input->kochin_np );

    cusfloat            wave_num    = w2k( ang_freq, input->water_depth, input->grav_acc );
    KochinInterface     kochin(
                                    mesh_gp->source_nodes[0],
                                    wave_num,
                                    input->water_depth,
                                    0,
                                    false
                                );
    std::vector<cuscomplex> pert_src( heads_np * source_nodes_np, cuscomplex( 0.0, 0.0 ) );

    for ( std::size_t ih = 0; ih < heads_np; ++ih )
    {
        for ( std::size_t body_id = 0; body_id < bodies_np; ++body_id )
        {
            for ( int source_id = mesh_gp->source_nodes_cnp[body_id]; source_id < mesh_gp->source_nodes_cnp[body_id+1]; ++source_id )
            {
                std::size_t pert_idx = ih * source_nodes_np + static_cast<std::size_t>( source_id );
                for ( std::size_t dof_id = 0; dof_id < dofs_np; ++dof_id )
                {
                    std::size_t rao_idx = ih * ( bodies_np * dofs_np ) + body_id * dofs_np + dof_id;
                     // Radiation intensities are now per (body, DOF): column = body_id * dofs_np + dof_id
                    std::size_t src_idx = ( body_id * dofs_np + dof_id ) * source_nodes_np + static_cast<std::size_t>( source_id );
                    pert_src[pert_idx] += cuscomplex( 0.0, -ang_freq ) * raos[rao_idx] * intensities[src_idx];
                }

                // Diffraction intensities start at column dofs_np * bodies_np
                pert_src[pert_idx] += intensities[( dofs_np * bodies_np + ih ) * source_nodes_np + static_cast<std::size_t>( source_id )];
            }
        }
    }

    std::size_t kochin_heads_offset = freq_index * static_cast<std::size_t>( sim_data->qtf_kochin_heads_np );
    for ( std::size_t ih = 0; ih < heads_np; ++ih )
    {
        calculate_kochin_coefficients(
                                        input,
                                        mesh_gp,
                                        &kochin,
                                        &( pert_src[ih * source_nodes_np] ),
                                        &( sim_data->qtf_kochin_pert_cos_freqs[kochin_heads_offset + ih * kochin_heads_stride] ),
                                        &( sim_data->qtf_kochin_pert_sin_freqs[kochin_heads_offset + ih * kochin_heads_stride] )
                                    );
    }

    std::size_t kochin_rad_offset = freq_index * static_cast<std::size_t>( sim_data->qtf_kochin_rad_np );
    // Loop over all (body, DOF) radiation mode pairs
    for ( std::size_t ib = 0; ib < bodies_np; ++ib )
    {
        for ( std::size_t dof_id = 0; dof_id < dofs_np; ++dof_id )
        {
            const std::size_t mode_id = ib * dofs_np + dof_id;
            calculate_kochin_coefficients(
                                            input,
                                            mesh_gp,
                                            &kochin,
                                            &( intensities[mode_id * source_nodes_np] ),
                                            &( sim_data->qtf_kochin_rad_cos_freqs[kochin_rad_offset + mode_id * kochin_rad_stride] ),
                                            &( sim_data->qtf_kochin_rad_sin_freqs[kochin_rad_offset + mode_id * kochin_rad_stride] )
                                        );
        }
    }
}


/*****************************************************************/
/*************** Define Frequency Solver class  ******************/
/*****************************************************************/
template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_calculate_field_points_values( 
                                                                    std::size_t freq_index,
                                                                    cusfloat    ang_freq 
                                                                )
{
    if ( this->input->is_calc_mdrift || this->input->out_qtf_so_model == QTFSOModelE::INDIRECT )
    {
        // Calculate velocity field at Bernoulli points for QTF calculations
        this->kernel->template compute_fields<RDDQTFConfig>(
                                                                freq_index, 
                                                                ang_freq,
                                                                this->sim_data->raos,
                                                                this->_qtf_bern_fields
                                                            );
        
        // Calculate waterline field points values for QTF calculations
        this->kernel->template compute_fields<RDDQTFConfig>( 
                                                                freq_index,
                                                                ang_freq,
                                                                this->sim_data->raos,
                                                                this->_qtf_wl_fields
                                                            );

    }

    if ( this->input->out_qtf_so_model == QTFSOModelE::DIRECT || this->input->out_qtf_so_model == QTFSOModelE::INDIRECT )
    {
        if ( this->_qtf_fs_fields != nullptr )
        {
            this->kernel->template compute_fields<RDDQTFConfig>(
                                                                    freq_index,
                                                                    ang_freq,
                                                                    this->sim_data->raos,
                                                                    this->_qtf_fs_fields
                                                                );
        }

        if ( this->input->out_qtf_so_model == QTFSOModelE::DIRECT && this->_qtf_fs_wl_fields != nullptr )
        {
            this->kernel->template compute_fields<RDDQTFConfig>(
                                                                    freq_index,
                                                                    ang_freq,
                                                                    this->sim_data->raos,
                                                                    this->_qtf_fs_wl_fields
                                                                );
        }

        for ( auto* annulus_field : this->_qtf_annuli_fields )
        {
            this->kernel->template compute_fields<RDDQTFConfig>(
                                                                    freq_index,
                                                                    ang_freq,
                                                                    this->sim_data->raos,
                                                                    annulus_field
                                                                );
        }
    }

    for ( std::size_t i=0; i<this->input->field_points_np; i++ )
    {
        // Calculate field points values for radiation and diffraction calculations
        this->kernel->template compute_fields<RDDFDConfig>( 
                                                                freq_index,
                                                                ang_freq,
                                                                this->sim_data->raos,
                                                                this->_field_points[i]->get_rdd( )
                                                            );

        // Append new step data to exporter interface
        this->_field_points[i]->add_step( ang_freq );

    }

    // Loop over bodies
    for ( std::size_t i=0; i<static_cast<std::size_t>( this->input->bodies_np ); i++ )
    {
        if ( this->input->bodies[i]->use_morison_elements )
        {
            // Loop over morison elements of the current body
            for ( std::size_t j=0; j<static_cast<std::size_t>( this->input->bodies[i]->morison_elements_names.size( ) ); j++ )
            {
                // Calculate field values over morison element field points
                this->kernel->template compute_fields<RDDMorisonConfig>( 
                                                                            freq_index,
                                                                            ang_freq,
                                                                            this->sim_data->raos,
                                                                            this->_morison_elements[i][j]->get_rdd( )
                                                                        );
            }
        }
    }
}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::calculate_first_order( void )
{
    /*******************************************************/
    /*   Calculate regular frequency band coefficients     */
    /*******************************************************/

    if ( this->input != nullptr && this->input->out_memory_report )
    {
        std::filesystem::path report_path = std::filesystem::path( this->input->case_fopath ) / "memory_report_first_order.csv";
        this->write_memory_report( report_path.string( ) );
    }

    // Loop over frequencies
    for ( std::size_t i=0; i<static_cast<std::size_t>(this->input->angfreqs_np); i++ )
    {
        LOG_TASK_SS( freq, "FO - Calculating period: " << std::setw( 10 ) << std::fixed << std::setprecision( 3 ) << angfreq_to_period( this->input->angfreqs[i] ) << " s" )

        // Start time measurement for the current frequency
        MpiTimer freq_timer;

        // Solve radiation-diffraction problem for the current frequency
        this->kernel->template solve<FreqRegimeE::REGULAR>( input->angfreqs[i] );
        this->kernel->update_results( sim_data );

        // Calculate first order coefficients
        this->_calculate_first_order_coeffs<FreqRegimeE::REGULAR>( i, this->input->angfreqs[i] );

        // Calculate field points values at the positions required
        this->_calculate_field_points_values( i, this->input->angfreqs[i] );

        // Calculate second order coefficients given by first order potential solution
        this->_calculate_mean_drift( i, false );

        // Calculate morison coefficients and forces for slender elements
        this->_calculate_morison_forces( input->angfreqs[i] );
        
        // Print out execution times
        LOG_TASK_TIME( freq, freq_timer )
        
    }

    /*******************************************************/
    /*   Calculate low frequency asymptotic coefficients   */
    /*******************************************************/
    LOG_TASK_SS( freq_lf, "FO - Calculating period: " << std::setw( 10 ) << std::fixed << std::setprecision( 3 ) << angfreq_to_period( W_ASYMPT_LOW ) << " s" )

    // Start time measurement
    MpiTimer freq_low_timer;

    // Solve radiation-diffraction problem for the current frequency
    this->kernel->template solve<FreqRegimeE::ASYMPT_LOW>( W_ASYMPT_LOW );
    this->kernel->update_results( sim_data );

    // Calculate first order coefficients
    this->_calculate_first_order_coeffs<FreqRegimeE::ASYMPT_LOW>( -1, W_ASYMPT_LOW );

    // Print out execution times
    LOG_TASK_TIME( freq_lf, freq_low_timer )

    /*******************************************************/
    /*   Calculate high frequency asymptotic coefficients  */
    /*******************************************************/
    LOG_TASK_SS( freq_hf, "FO - Calculating period: " << std::setw( 10 ) << std::fixed << std::setprecision( 3 ) << angfreq_to_period( W_ASYMPT_HIGH ) << " s" )

    // Start time measurement
    MpiTimer freq_high_timer;

    // Solve radiation-diffraction problem for the current frequency
    this->kernel->template solve<FreqRegimeE::ASYMPT_HIGH>( W_ASYMPT_HIGH );
    this->kernel->update_results( sim_data );

    // Calculate first order coefficients
    this->_calculate_first_order_coeffs<FreqRegimeE::ASYMPT_HIGH>( -1, W_ASYMPT_HIGH );

    // Print out execution times
    LOG_TASK_TIME( freq_hf, freq_high_timer )

}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::calculate_second_order( void )
{
    // Check if second order forces are to be calculated
    if ( !this->input->out_qtf )
    {
        return;
    }

    if ( this->input != nullptr && this->input->out_memory_report )
    {
        std::filesystem::path report_path = std::filesystem::path( this->input->case_fopath ) / "memory_report_second_order.csv";
        this->write_memory_report( report_path.string( ) );
    }

    // Recover partition circle radius from mesh if not already cached
    this->sim_data->qtf_pc_radius = this->input->bodies[0]->mesh_fs_qtf->get_fs_radius( );

    // Configure second order integrators with the current simulation data
    this->kernel->configure_second_order( 
                                            this->sim_data
                                        );

    // Loop over frequencies
    for ( std::size_t i=0; i<static_cast<std::size_t>(this->input->angfreqs_np); i++ )
    {
        for ( std::size_t j=0; j<static_cast<std::size_t>(this->input->angfreqs_np); j++ )
        {
            LOG_TASK_SS( freq, "SO - Calculating period: " << std::setw( 10 ) << std::fixed << std::setprecision( 3 ) << angfreq_to_period( this->input->angfreqs[i] ) << " s" << " - " << angfreq_to_period( this->input->angfreqs[j] ) << " s" )

            // Start time measurement for the current frequency
            MpiTimer freq_timer;

            // Calculate second order potential contribution
            if ( this->input->out_qtf_so_model == QTFSOModelE::INDIRECT )
            {
                calculate_secord_force_indirect(
                                                    this->input,
                                                    this->mpi_config,
                                                    this->mesh_gp,
                                                    i,
                                                    j,
                                                    QTFTypeE::QTF_DIFF_CODE,
                                                    this->_qtf_bern_fields,
                                                    this->_qtf_fs_fields,
                                                    this->_qtf_wl_fields,
                                                    this->sim_data
                                                );

                calculate_secord_force_indirect(
                                                    this->input,
                                                    this->mpi_config,
                                                    this->mesh_gp,
                                                    i,
                                                    j,
                                                    QTFTypeE::QTF_SUM_CODE,
                                                    this->_qtf_bern_fields,
                                                    this->_qtf_fs_fields,
                                                    this->_qtf_wl_fields,
                                                    this->sim_data
                                                );
            }
            else if ( i != j )
            {
                if ( this->input->out_qtf_so_model == QTFSOModelE::PINKSTER )
                {
                    if ( this->mpi_config->is_root( ) )
                    {
                        calculate_pinkster(
                                                this->input,
                                                this->mpi_config,
                                                this->mesh_gp,
                                                this->input->angfreqs[i],
                                                this->input->angfreqs[j],
                                                this->sim_data->qtf_diff_secord_force
                                            );
                    }
                }
                else if ( this->input->out_qtf_so_model == QTFSOModelE::DIRECT )
                {
                    // Solve radiation-diffraction problem for the current frequency
                    this->kernel->template solve_so<QTFTypeE::QTF_DIFF_CODE>( 
                                                                                i,
                                                                                j,
                                                                                this->sim_data->raos_hist,
                                                                                this->_qtf_bern_fields,
                                                                                this->_qtf_wl_fields,
                                                                                this->_qtf_fs_fields,
                                                                                this->_qtf_fs_wl_fields,
                                                                                &(this->_qtf_annuli_fields)
                                                                            );
                    this->kernel->template update_results_so<QTFTypeE::QTF_DIFF_CODE>( sim_data );

                    calculate_diffraction_forces_lin(
                                                        this->input,
                                                        this->mpi_config,
                                                        this->mesh_gp,
                                                        sim_data->panels_potential_so,
                                                        std::abs( this->input->angfreqs[i] - this->input->angfreqs[j] ),
                                                        sim_data->qtf_diff_secord_force,
                                                        sim_data->panels_pressure_so,
                                                        true
                                                    );

                    MPI_Allreduce(
                                        MPI_IN_PLACE,
                                        sim_data->qtf_diff_secord_force,
                                        sim_data->qtf_np,
                                        mpi_cuscomplex,
                                        MPI_SUM,
                                        MPI_COMM_WORLD
                                    );
                }

                MPI_Barrier( MPI_COMM_WORLD );

                // Distribute second order force
                if ( this->mpi_config->is_root( ) )
                {
                    qtf_distribute_matrix_data(
                                                    this->input,
                                                    i,
                                                    j,
                                                    this->sim_data->qtf_diff_secord_force,
                                                    this->sim_data->qtf_diff_freqs,
                                                    1,
                                                    0
                                                );
                    if ( this->input->out_qtf_comp )
                    {
                        qtf_distribute_matrix_data(
                                                    this->input,
                                                    i,
                                                    j,
                                                    this->sim_data->qtf_diff_secord_force,
                                                    this->sim_data->qtf_diff_secord_force_freqs,
                                                    1,
                                                    0
                                                );
                    }
                }
            }

            if ( this->mpi_config->is_root( ) && this->input->out_qtf_so_model == QTFSOModelE::INDIRECT && this->input->out_qtf_comp )
            {
                if ( i != j )
                {
                    qtf_distribute_matrix_data(
                                                    this->input,
                                                    i,
                                                    j,
                                                    this->sim_data->qtf_diff_secord_force,
                                                    this->sim_data->qtf_diff_secord_force_freqs,
                                                    1,
                                                    0
                                                );
                }

                qtf_distribute_matrix_data(
                                                this->input,
                                                i,
                                                j,
                                                this->sim_data->qtf_sum_secord_force,
                                                this->sim_data->qtf_sum_secord_force_freqs,
                                                1,
                                                0
                                            );
            }

            this->_calculate_and_distribute_qtf<QTFTypeE::QTF_DIFF_CODE>( i, j );
            this->_calculate_and_distribute_qtf<QTFTypeE::QTF_SUM_CODE>( i, j );

            // Print out execution times
            LOG_TASK_TIME( freq, freq_timer )
        }
    }

    // Save second order forces to disk
    this->_save_qtf<QTFTypeE::QTF_DIFF_CODE>( );
    
}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_calculate_far_field_diffraction( 
                                                                        std::size_t freq_index_i
                                                                    )
{
    // Define local auxiliar variables
    cusfloat    a               = this->input->wave_amplitude;
    cusfloat    c               = 0.0;
    cusfloat    ccf             = 0.0;
    cusfloat    ep              = 0.0;
    cusfloat    g               = this->input->grav_acc;
    cusfloat    h               = this->input->water_depth;   
    std::size_t heads_offset    = freq_index_i * this->input->heads_np * QTF_FAR_N;
    cuscomplex  int_mod         = cuscomplex( 0.0, 0.0 );
    cuscomplex  int_value_cos   = cuscomplex( 0.0, 0.0 );
    cuscomplex  int_value_sin   = cuscomplex( 0.0, 0.0 );
    cusfloat    k               = this->input->wave_numbers[freq_index_i];
    PanelGeom*  panel           = nullptr;
    cusfloat    rho             = 0.0;
    std::size_t source_offset   = 0;
    cusfloat    theta           = 0.0;
    cusfloat    w               = this->input->angfreqs[freq_index_i];
    cusfloat    nu              = pow2s( w ) / g;

    // Define scaling factor
    cuscomplex  sf              = cuscomplex( 0.0, -g * a / w );
    
    // Calculate wave factor c
    c = -pow2s( k ) / ( h * ( pow2s( k ) - pow2s( nu ) ) + nu );

    for ( std::size_t ih=0; ih<static_cast<std::size_t>(this->input->heads_np); ih++ )
    {
        ep = 1.0;
        source_offset = this->mesh_gp->panels_tnp * ih;
        for ( std::size_t n=0; n<QTF_FAR_N; n++ )
        {
            for ( std::size_t i=0; i<static_cast<std::size_t>(this->mesh_gp->panels_tnp); i++ )
            {
                panel           = this->mesh_gp->panels[i];
                int_value_cos   = cuscomplex( 0.0, 0.0 );
                int_value_sin   = cuscomplex( 0.0, 0.0 );
                if ( panel->type == PanelTypeE::DIFFRAC )
                {
                    // Calculate horizontal radius and angle of the current panel
                    rho     = sqrt( pow2s( panel->center[0] ) + pow2s( panel->center[1] ) );
                    theta   = atan2( panel->center[1], panel->center[0] );
                    
                    // Calculate finite depth factor
                    ccf = cosh_cosh_factor( k, h, panel->center[2] );

                    // Calculate integral value for the current panel
                    int_mod       = this->sim_data->intensities[source_offset + i] * ccf * std::cyl_bessel_j( static_cast<cusfloat>( n ), k * rho ) * panel->area;
                    int_value_cos += int_mod * cos( n * theta );
                    int_value_sin += int_mod * sin( n * theta );
                }

            }

            const std::size_t idx = heads_offset + ih * QTF_FAR_N + n;
            this->sim_data->qtf_b_cos[idx] = cuscomplex( 0.0, 2*PI*c*ep ) * int_value_cos;
            this->sim_data->qtf_b_sin[idx] = cuscomplex( 0.0, 2*PI*c*ep ) * int_value_sin;
            ep = 2.0;
        }
    }

}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_calculate_far_field_incident( 
                                                                    std::size_t freq_index_i
                                                                )
{
    // Define local auxiliar variables
    cusfloat    a               = this->input->wave_amplitude;
    cusfloat    ep              = 0.0;
    cusfloat    g               = this->input->grav_acc;
    std::size_t heads_offset    = freq_index_i * this->input->heads_np * QTF_FAR_N;
    cuscomplex  in              = cuscomplex( 1.0, 0.0 );
    cusfloat    w               = this->input->angfreqs[freq_index_i];

    // Define scaling factor
    cuscomplex  sf              = cuscomplex( 0.0, -g * a / w );

    for ( std::size_t ih=0; ih<static_cast<std::size_t>(this->input->heads_np); ih++ )
    {
        ep = 1.0;
        
        for ( std::size_t n=0; n<QTF_FAR_N; n++ )
        {
            const std::size_t idx = heads_offset + ih * QTF_FAR_N + n;
            this->sim_data->qtf_a_cos[idx] = sf * ep * in * cos( n * this->input->heads[ih] );
            this->sim_data->qtf_a_sin[idx] = sf * ep * in * sin( n * this->input->heads[ih] );
            in *= cuscomplex( 0.0, -1 );
            ep = 2.0;
        }
    }

}


template<std::size_t N, int mode_pf>
template<FreqRegimeE freq_regime>
void FrequencySolver<N, mode_pf>::_calculate_first_order_coeffs( 
                                                                    std::size_t freq_index, 
                                                                    cusfloat    ang_freq 
                                                                )
{
    // Calculate added mass and damping
    calculate_hydromechanic_coeffs_lin( 
                                            this->input,
                                            this->mpi_config,
                                            this->mesh_gp,
                                            this->sim_data->panels_potential,
                                            ang_freq,
                                            this->sim_data->added_mass,
                                            this->sim_data->damping_rad,
                                            this->sim_data->panels_pressure
                                        );

    // Calculate wave excitation forces
    if ( freq_regime == FreqRegimeE::REGULAR || freq_regime == FreqRegimeE::ASYMPT_LOW )
    {
        // Calculate diffraction forces
        calculate_diffraction_forces_lin(
                                            this->input,
                                            this->mpi_config,
                                            this->mesh_gp,
                                            this->sim_data->panels_potential,
                                            ang_freq,
                                            this->sim_data->wave_diffrac,
                                            this->sim_data->panels_pressure
                                    );
        
        // Calculate Froude-Krylov first order forces
        calculate_froude_krylov_fo(
                                                this->input,
                                                this->mpi_config,
                                                this->mesh_gp,
                                                ang_freq,
                                                this->sim_data->froude_krylov
                                    );
    }
    else
    {
        clear_vector( this->sim_data->wave_exc_np, this->sim_data->wave_diffrac );
        clear_vector( this->sim_data->wave_exc_np, this->sim_data->froude_krylov );
    }
    
    // Reduce data into the root processor variablesç
    REDUCE_FO_ROOT( added_mass,    hydmech,  cusfloat   )
    REDUCE_FO_ROOT( damping_rad,   hydmech,  cusfloat   )
    REDUCE_FO_ROOT( froude_krylov, wave_exc, cuscomplex )
    REDUCE_FO_ROOT( wave_diffrac,  wave_exc, cuscomplex )
    
    // Generate total wave excitation forces
    if ( this->mpi_config->is_root( ) )
    {
        sv_add(
                    this->sim_data->wave_exc_np,
                    this->sim_data->wave_diffrac_p0,
                    this->sim_data->froude_krylov_p0,
                    this->sim_data->wave_exc_p0
                );
    }
    
    // Calculate RAOs
    if ( this->mpi_config->is_root( )  && freq_regime == FreqRegimeE::REGULAR )
    {
        calculate_raos(
                            this->input,
                            this->sim_data->structural_mass_p0,
                            this->sim_data->added_mass_p0,
                            this->sim_data->damping_rad_p0,
                            this->sim_data->hydrostiff_p0,
                            this->sim_data->wave_diffrac_p0,
                            this->sim_data->froude_krylov_p0,
                            this->sim_data->morison_drag_force_p0,
                            this->sim_data->morison_inertial_force_p0,
                            ang_freq,
                            this->sim_data->raos
                        );
    }

    // Distribute RAOs to all processors
    MPI_Bcast(
                    this->sim_data->raos,
                    this->sim_data->wave_exc_np,
                    mpi_cuscomplex,
                    this->mpi_config->proc_root,
                    MPI_COMM_WORLD
                );

    // Copy RAOs history only for regular frequencies (asymptotic uses a sentinel index).
    if constexpr( freq_regime == FreqRegimeE::REGULAR )
    {
        if ( this->input->is_calc_mdrift || this->input->out_qtf )
        {
            copy_vector(
                            this->sim_data->wave_exc_np,
                            this->sim_data->raos,
                            &(this->sim_data->raos_hist[freq_index * this->sim_data->wave_exc_np])
                        );
        }
    }

    
    // Storage results
    // MpiTimer storage_timer;
    const std::size_t field_count = this->kernel->size( ) * ( this->input->dofs_np * this->input->bodies_np + this->input->heads_np );
    auto reduce_field = [this, field_count](cuscomplex* src, cuscomplex* dst)
    {
        MPI_Reduce(
                        src,
                        dst,
                        field_count,
                        mpi_cuscomplex,
                        MPI_SUM,
                        this->mpi_config->proc_root,
                        MPI_COMM_WORLD
                    );
    };

    const bool need_source_history = ( freq_regime == FreqRegimeE::REGULAR )
                                    &&
                                    ( this->input->out_sources || this->input->out_qtf_so_model == QTFSOModelE::INDIRECT );

    if ( need_source_history )
    {
        reduce_field( this->sim_data->intensities, this->sim_data->intensities_p0 );
    }

    if ( this->mpi_config->is_root( ) && freq_regime == FreqRegimeE::REGULAR && this->input->out_qtf_so_model == QTFSOModelE::INDIRECT )
    {
        _store_indirect_kochin_data(
                                        this->input,
                                        this->mesh_gp,
                                        freq_index,
                                        ang_freq,
                                        this->sim_data->intensities_p0,
                                        this->sim_data->raos,
                                        this->sim_data
                                    );
    }

    if ( this->input->out_potential && freq_regime == FreqRegimeE::REGULAR )
    {
        reduce_field( this->sim_data->panels_potential, this->sim_data->panels_potential_p0 );
    }

    if ( this->input->out_pressure && freq_regime == FreqRegimeE::REGULAR )
    {
        reduce_field( this->sim_data->panels_pressure, this->sim_data->panels_pressure_p0 );
    }
    
    if ( this->mpi_config->is_root( ) )
    {
        if ( this->input->out_hydmech )
        {
            if constexpr( freq_regime == FreqRegimeE::REGULAR )
            {
                this->output->save_hydromechanics_format(
                                                            freq_index,
                                                            _DN_ADDED_MASS,
                                                            this->sim_data->added_mass_p0
                                                        );

                this->output->save_hydromechanics_format(
                                                            freq_index,
                                                            _DN_DAMPING_RAD,
                                                            this->sim_data->damping_rad_p0
                                                        );
            }
            else
            {
                // Select dataset names depending on the frequency regime
                std::string added_mass_dn   = _DN_ADDED_MASS_LF;
                std::string damping_dn      = _DN_DAMPING_RAD_LF;
                if ( freq_regime == FreqRegimeE::ASYMPT_HIGH )
                {
                    added_mass_dn   = _DN_ADDED_MASS_HF;
                    damping_dn      = _DN_DAMPING_RAD_HF;
                }

                // Storage asymptotic added mass and damping
                this->output->save_hydromechanics_asympt_format(
                                                                    added_mass_dn,
                                                                    this->sim_data->added_mass
                                                                );

                this->output->save_hydromechanics_asympt_format(
                                                                    damping_dn,
                                                                    this->sim_data->damping_rad
                                                                );
            }
            
        }

        if ( this->input->out_diffrac )
        {
            if constexpr( freq_regime == FreqRegimeE::REGULAR )
            {
                this->output->save_wave_exciting_format(
                                                            freq_index,
                                                            _DN_DIFFRAC,
                                                            this->sim_data->wave_diffrac_p0
                                                        );
            }
            else
            {
                std::string diffrac_dn   = _DN_DIFFRAC_HF;
                if ( freq_regime == FreqRegimeE::ASYMPT_LOW )
                {
                    diffrac_dn   = _DN_DIFFRAC_LF;
                }

                this->output->save_wave_exciting_asympt_format(
                                                                    diffrac_dn,
                                                                    this->sim_data->wave_diffrac_p0
                                                                );
            }
        }

        if ( this->input->out_fk )
        {
            if constexpr( freq_regime == FreqRegimeE::REGULAR )
            {
                this->output->save_wave_exciting_format(
                                                            freq_index,
                                                            _DN_FK,
                                                            this->sim_data->froude_krylov_p0
                                                        );
            }
            else
            {
                std::string fk_dn   = _DN_FK_HF;
                if ( freq_regime == FreqRegimeE::ASYMPT_LOW )
                {
                    fk_dn   = _DN_FK_LF;
                }

                this->output->save_wave_exciting_asympt_format(
                                                                    fk_dn,
                                                                    this->sim_data->froude_krylov_p0
                                                                );
            }
        }

        if ( this->input->out_wex )
        {
            if constexpr( freq_regime == FreqRegimeE::REGULAR )
            {
                this->output->save_wave_exciting_format(
                                                            freq_index,
                                                            _DN_WEX,
                                                            this->sim_data->wave_exc_p0
                                                        );
            }
            else
            {
                std::string wex_dn   = _DN_WEX_HF;
                if ( freq_regime == FreqRegimeE::ASYMPT_LOW )
                {
                    wex_dn   = _DN_WEX_LF;
                }

                this->output->save_wave_exciting_asympt_format(
                                                                    wex_dn,
                                                                    this->sim_data->wave_exc_p0
                                                                );
            }   
        }

        if ( this->input->out_morison )
        {
            this->output->save_wave_exciting_format(
                                                        freq_index,
                                                        _DN_MORISON_DRAG,
                                                        this->sim_data->morison_drag_force_p0
                                                    );

            this->output->save_wave_exciting_format(
                                                        freq_index,
                                                        _DN_MORISON_INERTIAL,
                                                        this->sim_data->morison_inertial_force_p0
                                                    );
        }

        if ( this->input->out_raos && freq_regime == FreqRegimeE::REGULAR )
        {
            this->output->save_wave_exciting_format(
                                                        freq_index,
                                                        _DN_RAO,
                                                        this->sim_data->raos
                                                    );
        }

        // Storage sources
        if ( this->input->out_sources && freq_regime == FreqRegimeE::REGULAR )
        {
            this->output->save_fields_data( 
                                                freq_index,
                                                _DN_SRC_INT,
                                                this->sim_data->intensities_p0
                                            );
        }

        // Storage panels potential
        if ( this->input->out_potential && freq_regime == FreqRegimeE::REGULAR )
        {
            this->output->save_fields_data( 
                                                freq_index,
                                                _DN_POT_INT,
                                                this->sim_data->panels_potential_p0
                                            );
        }

        // Storage panels pressure
        if ( this->input->out_pressure && freq_regime == FreqRegimeE::REGULAR )
        {
            this->output->save_fields_data( 
                                                freq_index,
                                                _DN_PRESS_INT,
                                                this->sim_data->panels_pressure_p0
                                            );
        }

    }
}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_calculate_mean_drift( 
                                                            std::size_t freq_index,
                                                            bool        is_multi_head
                                                        )
{
    // Calculate QTF coefficients if required
    if ( this->input->is_calc_mdrift )
    {
        this->_calculate_quadratic_terms<QTFTypeE::QTF_DIFF_CODE>( 
                                                                        freq_index,
                                                                        freq_index,
                                                                        is_multi_head 
                                                                    );

        if ( this->mpi_config->is_root( ) )
        {
            // Storage data if any
            this->output->save_wave_exciting_format(
                                                        freq_index,
                                                        _DN_MDRIFT,
                                                        this->sim_data->qtf
                                                    );
    
            if ( this->input->out_qtf_comp )
            {
                output->save_wave_exciting_format(
                                                    freq_index,
                                                    _DN_MDRIFT_WL,
                                                    sim_data->qtf_diff_wl
                                                );
                
                output->save_wave_exciting_format(
                                                    freq_index,
                                                    _DN_MDRIFT_BERN,
                                                    sim_data->qtf_diff_bern
                                                );
    
                output->save_wave_exciting_format(
                                                    freq_index,
                                                    _DN_MDRIFT_ACC,
                                                    sim_data->qtf_diff_acc
                                                );
                
                output->save_wave_exciting_format(
                                                    freq_index,
                                                    _DN_MDRIFT_MOM,
                                                    sim_data->qtf_diff_mom
                                                );
            }
        }
    }

    if ( this->input->out_qtf_so_model == QTFSOModelE::DIRECT )
    {
        this->_calculate_far_field_incident( freq_index );
        this->_calculate_far_field_diffraction( freq_index );
    }

}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_calculate_morison_forces( cusfloat ang_freq )
{
    // Clear morison forces vector buffer
    clear_vector( this->input->heads_np * this->input->dofs_np * this->input->bodies_np, this->sim_data->morison_drag_force         );
    clear_vector( this->input->heads_np * this->input->dofs_np * this->input->bodies_np, this->sim_data->morison_inertial_force     );

    // Loop over bodies
    for ( std::size_t i=0; i<static_cast<std::size_t>( this->input->bodies_np ); i++ )
    {
        // Loop over morison elements of the current body
        for ( std::size_t j=0; j<static_cast<std::size_t>( this->input->bodies[i]->morison_elements_names.size( ) ); j++ )
        {
            // Calculate hydrodynamic morison forces for the current element
            this->_morison_elements[i][j]->calculate_hydrodynamic_forces( 
                                                                                ang_freq,
                                                                                &(this->sim_data->raos[ 6*i ]),
                                                                                this->_inertial_force_aux,
                                                                                this->_drag_force_aux
                                                                            );
            
            // Add inertial and drag forces to body total morison forces
            for ( std::size_t ih=0; ih<static_cast<std::size_t>( this->input->heads_np ); ih++ )
            {
                std::size_t offset = ih * ( this->input->dofs_np * this->input->bodies_np ) + i * this->input->dofs_np;
                for ( std::size_t id=0; id<static_cast<std::size_t>( this->input->dofs_np ); id++ )
                {
                    std::size_t idx = ih * this->input->dofs_np + id;
                    this->sim_data->morison_drag_force[ offset + id ]       += this->_drag_force_aux[idx];
                    this->sim_data->morison_inertial_force[ offset + id ]   += this->_inertial_force_aux[idx];
                }
            }
        }
    }

    // Reduce Morison forces into the root processor variable
    REDUCE_FO_ROOT( morison_drag_force, wave_exc, cuscomplex )
    REDUCE_FO_ROOT( morison_inertial_force, wave_exc, cuscomplex )
}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_calculate_global_static_matrixes( void )
{
    // Calculate global structural mass
    if ( this->mpi_config->is_root( ) )
    {
        calculate_global_structural_mass(
                                            input,
                                            sim_data->structural_mass_p0
                                        );
    }


    // Calculate global stiffness matrix
    if ( this->mpi_config->is_root( ) )
    {
        calculate_global_hydstiffness(
                                            input,
                                            hydrostatics,
                                            sim_data->hydrostiff_p0
                                        );
    }
}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_calculate_hydrostatics( void )
{
    LOG_TASK_SS( hydro, "Calculating hydrostatics..." )
    MpiTimer hydro_timer;

    this->hydrostatics = new Hydrostatics*[this->input->bodies_np];
    for ( int i=0; i<input->bodies_np; i++ )
    {
        this->hydrostatics[i]   =   new Hydrostatics( 
                                                        this->input->bodies[i]->mesh,
                                                        this->input->water_density,
                                                        this->input->grav_acc,
                                                        this->input->bodies[i]->mass,
                                                        this->input->bodies[i]->cog,
                                                        this->input->bodies[i]->rad_inertia,
                                                        this->mpi_config
                                                    );
    }

    LOG_TASK_TIME( hydro, hydro_timer )

}


template<std::size_t N, int mode_pf>
template<QTFTypeE qtf_type>
void FrequencySolver<N, mode_pf>::_calculate_quadratic_terms( 
                                                                    std::size_t freq_index_i,
                                                                    std::size_t freq_index_j,
                                                                    bool        is_multi_head
                                                            )
{
    // Get local variables from structures for easier access
    cusfloat                                ang_freq_i      = this->input->angfreqs[freq_index_i];
    cusfloat                                ang_freq_j      = this->input->angfreqs[freq_index_j];
    std::size_t                             bodies_np       = static_cast<std::size_t>( this->mesh_gp->meshes_np );
    std::size_t                             dofs_np         = static_cast<std::size_t>( this->input->dofs_np );
    cusfloat                                grav_acc        = this->input->grav_acc;
    std::size_t                             heads_np        = static_cast<std::size_t>( this->input->heads_np );
    Input*                                  input           = this->input;
    PanelData<cuscomplex, RDDQTFConfig>*    paneld          = nullptr;
    cuscomplex*                             qtf_values      = nullptr;
    cuscomplex*                             qtf_wl          = nullptr;
    cuscomplex*                             qtf_bern        = nullptr;
    cuscomplex*                             qtf_acc         = nullptr;
    cuscomplex*                             qtf_mom         = nullptr;
    cuscomplex*                             raos_i          = &( this->sim_data->raos_hist[ freq_index_i * this->sim_data->wave_exc_np ] );
    cuscomplex*                             raos_j          = &( this->sim_data->raos_hist[ freq_index_j * this->sim_data->wave_exc_np ] );
    RadDiffData<cuscomplex, RDDQTFConfig>*  rdd_bern        = this->_qtf_bern_fields;
    RadDiffData<cuscomplex, RDDQTFConfig>*  rdd_rwel        = this->_qtf_wl_fields;
    cusfloat                                rhow            = this->input->water_density;
    cuscomplex*                             vel_x           = nullptr;
    cuscomplex*                             vel_y           = nullptr;
    cuscomplex*                             vel_z           = nullptr;
    cusfloat                                wave_amplitude  = this->input->wave_amplitude;

    if constexpr( qtf_type == QTFTypeE::QTF_DIFF_CODE )
    {
        qtf_values   = this->sim_data->qtf;
        qtf_wl      = this->sim_data->qtf_diff_wl;
        qtf_bern    = this->sim_data->qtf_diff_bern;
        qtf_acc     = this->sim_data->qtf_diff_acc;
        qtf_mom     = this->sim_data->qtf_diff_mom;
    }
    else if constexpr ( qtf_type == QTFTypeE::QTF_SUM_CODE )
    {
        qtf_values  = this->sim_data->qtf;
        qtf_wl      = this->sim_data->qtf_sum_wl;
        qtf_bern    = this->sim_data->qtf_sum_bern;
        qtf_acc     = this->sim_data->qtf_sum_acc;
        qtf_mom     = this->sim_data->qtf_sum_mom;
    }


    // Define aux variables to be used along the function
    static constexpr std::size_t ngp  = PanelGeom::gauss_points_np;
    static constexpr std::size_t ngp2 = PanelGeom::gauss_points_np * PanelGeom::gauss_points_np;

    cusfloat*   body_cog            = nullptr;
    cuscomplex  daux                = 0.0;
    std::size_t fp_np               = 0;
    std::size_t idx0                = 0;
    std::size_t idx1_i              = 0;
    std::size_t idx1_j              = 0;
    std::size_t ih2_end             = 0;
    std::size_t ih2_start           = 0;
    cuscomplex  int_mod_1d[ngp]     ;
    cuscomplex  int_mod_2d[ngp2]    ;
    cuscomplex  int_val             = 0.0;
    cusfloat*   normal_vec          = nullptr;

    // Set second heading loop bounds according with the multi-heading option
    if ( is_multi_head )
    {
        ih2_start   = 0;
        ih2_end     = heads_np;
    }
    else
    {
        ih2_start   = 0;
        ih2_end     = 1;
    }

    // Clear QTF input vector to ensure that previous data will not be storaged
    // erroneously
    std::size_t heads_factor_np = ( is_multi_head ) ? pow2s( heads_np ) : heads_np;
    std::size_t qtf_size        = heads_factor_np * bodies_np * dofs_np;

    clear_vector( qtf_size, qtf_values  );
    clear_vector( qtf_size, qtf_wl      );
    clear_vector( qtf_size, qtf_bern    ); 
    clear_vector( qtf_size, qtf_acc     );
    clear_vector( qtf_size, qtf_mom     );

    // Calculate second order force due to relative wave
    // elevation at the WL
    for ( std::size_t ih1=0; ih1<heads_np; ih1++ )
    {
        for ( std::size_t ih2=ih2_start; ih2<ih2_end; ih2++ )
        {
            for ( std::size_t i=0; i<rdd_rwel->get_size_local( ); i++  )
            {
                // Get panel
                paneld      = &(rdd_rwel->panel_data[i]);
                fp_np       = paneld->field_points_np;
                normal_vec  = paneld->panel_geom->normal_vec;

                // Loop over field points of the panel
                for ( std::size_t k=0; k<fp_np; k++ )
                {
                    // Calculate indexes for the two field points
                    _field_point_index( 
                                            is_multi_head,
                                            heads_np,
                                            fp_np,
                                            freq_index_i,
                                            freq_index_j,
                                            ih1,
                                            ih2,
                                            k,
                                            idx1_i,
                                            idx1_j
                                        );

                    // Calculate integrand value depending on the QTF type
                    if constexpr( qtf_type == QTFTypeE::QTF_DIFF_CODE )
                    {
                        int_mod_1d[k] = paneld->wev_rel_total[idx1_i] * std::conj( paneld->wev_rel_total[idx1_j] );
                    }
                    else if constexpr( qtf_type == QTFTypeE::QTF_SUM_CODE )
                    {
                        int_mod_1d[k] = paneld->wev_rel_total[idx1_i] * paneld->wev_rel_total[idx1_j];
                    }
                }

                // Integrate over the panel
                int_val = 0.0;
                gauss1d_loop<ngp>( int_val, int_mod_1d, paneld->panel_geom->len_wl );

                // Get global QTF index offset
                idx0    =   _qtf_index_offset(
                                                is_multi_head,
                                                bodies_np,
                                                heads_np,
                                                dofs_np,
                                                ih1,
                                                ih2,
                                                paneld->body_id
                                            );

                // Calculate final qtf value and store appropriately
                for ( std::size_t r=0; r<dofs_np; r++ )
                {
                    daux                = - 0.25 * grav_acc * rhow * int_val * normal_vec[r];
                    qtf_values[idx0+r]  += daux;
                    qtf_wl[idx0+r]      += daux;
                }
                
            }
        }
    }

    // Calculate second order force due to the bernouilly contribution
    for ( std::size_t ih1=0; ih1<heads_np; ih1++ )
    {
        for ( std::size_t ih2=ih2_start; ih2<ih2_end; ih2++ )
        {
            for ( std::size_t i=0; i<rdd_bern->get_size_local( ); i++  )
            {
                // Get panel
                paneld      = &(rdd_bern->panel_data[i]);
                fp_np       = paneld->field_points_np;
                normal_vec  = paneld->panel_geom->normal_vec;

                // Get field data
                vel_x       = paneld->vel_x_total.data();
                vel_y       = paneld->vel_y_total.data();
                vel_z       = paneld->vel_z_total.data();

                // Loop over field points of the panel
                for ( std::size_t k=0; k<fp_np; k++ )
                {
                    // Calculate indexes for the two field points
                    _field_point_index( 
                                            is_multi_head,
                                            heads_np,
                                            fp_np,
                                            freq_index_i,
                                            freq_index_j,
                                            ih1,
                                            ih2,
                                            k,
                                            idx1_i,
                                            idx1_j
                                        );

                    if constexpr( qtf_type == QTFTypeE::QTF_DIFF_CODE )
                    {
                        int_mod_2d[k]   =   (
                                                vel_x[idx1_i] * std::conj( vel_x[idx1_j] )
                                                +
                                                vel_y[idx1_i] * std::conj( vel_y[idx1_j] )
                                                +
                                                vel_z[idx1_i] * std::conj( vel_z[idx1_j] )
                                            );
                    }
                    else if constexpr( qtf_type == QTFTypeE::QTF_SUM_CODE )
                    {
                        int_mod_2d[k]   =   (
                                                vel_x[idx1_i] * vel_x[idx1_j]
                                                +
                                                vel_y[idx1_i] * vel_y[idx1_j]
                                                +
                                                vel_z[idx1_i] * vel_z[idx1_j]
                                            );
                    }
                }

                // Integrate over the panel
                int_val = 0.0;
                gauss2d_loop<ngp>( int_val, int_mod_2d, paneld->panel_geom );

                // Get global QTF index offset
                idx0    =   _qtf_index_offset(
                                                is_multi_head,
                                                bodies_np,
                                                heads_np,
                                                dofs_np,
                                                ih1,
                                                ih2,
                                                paneld->body_id
                                            );

                // Calculate final qtf value and store appropriately
                for ( std::size_t r=0; r<dofs_np; r++ )
                {
                    daux                = 0.25 * rhow * int_val * normal_vec[r];
                    qtf_values[idx0+r]  += daux;
                    qtf_bern[idx0+r]    += daux;
                }

            }
        }
    }

    // Calculate second order force due to acceleration term
    cusfloat    cog_to_fp[3];               clear_vector( 3, cog_to_fp );
    cuscomplex  cog_to_fp_c[3];             clear_vector( 3, cog_to_fp );
    cuscomplex  point_disp_i[3];            clear_vector( 3, point_disp_i );
    cuscomplex  point_disp_j[3];            clear_vector( 3, point_disp_j );
    cuscomplex  rao_rot_i[3];               clear_vector( 3, rao_rot_i );
    cuscomplex  rao_rot_j[3];               clear_vector( 3, rao_rot_j );
    cuscomplex  rao_trans_i[3];             clear_vector( 3, rao_trans_i );
    cuscomplex  rao_trans_j[3];             clear_vector( 3, rao_trans_j );
    cuscomplex  vel_x_acc_i, vel_x_acc_j;
    cuscomplex  vel_y_acc_i, vel_y_acc_j;
    cuscomplex  vel_z_acc_i, vel_z_acc_j;

    for ( std::size_t ih1=0; ih1<heads_np; ih1++ )
    {
        for ( std::size_t ih2=ih2_start; ih2<ih2_end; ih2++ )
        {
            for ( std::size_t i=0; i<rdd_bern->get_size_local( ); i++  )
            {
                // Get panel
                paneld      = &(rdd_bern->panel_data[i]);
                body_cog    = paneld->panel_geom->body_cog;
                fp_np       = paneld->field_points_np;
                normal_vec  = paneld->panel_geom->normal_vec;

                // Get field data
                vel_x       = paneld->vel_x_total.data( );
                vel_y       = paneld->vel_y_total.data( );
                vel_z       = paneld->vel_z_total.data( );

                // Get RAO values
                _rao_offset( 
                                        is_multi_head,
                                        bodies_np,
                                        dofs_np,
                                        ih1,
                                        ih2,
                                        paneld->body_id,
                                        idx1_i,
                                        idx1_j
                                    );
                
                for ( std::size_t r=0; r<3; r++ )
                {
                    rao_rot_i[r]    = raos_i[idx1_i+3+r] * wave_amplitude;
                    rao_rot_j[r]    = raos_j[idx1_j+3+r] * wave_amplitude;
                    rao_trans_i[r]  = raos_i[idx1_i+r]   * wave_amplitude;
                    rao_trans_j[r]  = raos_j[idx1_j+r]   * wave_amplitude;
                }
                
                // Loop over field points of the panel
                for ( std::size_t k=0; k<fp_np; k++ )
                {
                    // Calculate indexes for the two field points
                    _field_point_index( 
                                            is_multi_head,
                                            heads_np,
                                            fp_np,
                                            freq_index_i,
                                            freq_index_j,
                                            ih1,
                                            ih2,
                                            k,
                                            idx1_i,
                                            idx1_j
                                        );

                    // Define vector from cog to field point
                    sv_sub( 3, &(paneld->field_points[3*k]), body_cog, cog_to_fp );
                    for ( std::size_t r=0; r<3; r++ )
                    {
                        cog_to_fp_c[r]  = cuscomplex( cog_to_fp[r], 0.0 );
                    }

                    // Calculate first order displacement of the panel centre
                    clear_vector( 3, point_disp_i );

                    cross(
                                rao_rot_i,
                                cog_to_fp_c,
                                point_disp_i
                        );
                    sv_add(
                                3,
                                point_disp_i,
                                rao_trans_i,
                                point_disp_i
                            );

                    clear_vector( 3, point_disp_j );

                    cross(
                                rao_rot_j,
                                cog_to_fp_c,
                                point_disp_j
                        );
                    sv_add(
                                3,
                                point_disp_j,
                                rao_trans_j,
                                point_disp_j
                            );

                    // Get velocity pressure term
                    vel_x_acc_i = rhow * cuscomplex( 0.0, -ang_freq_i ) * vel_x[idx1_i];
                    vel_y_acc_i = rhow * cuscomplex( 0.0, -ang_freq_i ) * vel_y[idx1_i];
                    vel_z_acc_i = rhow * cuscomplex( 0.0, -ang_freq_i ) * vel_z[idx1_i];

                    vel_x_acc_j = rhow * cuscomplex( 0.0, -ang_freq_j ) * vel_x[idx1_j];
                    vel_y_acc_j = rhow * cuscomplex( 0.0, -ang_freq_j ) * vel_y[idx1_j];
                    vel_z_acc_j = rhow * cuscomplex( 0.0, -ang_freq_j ) * vel_z[idx1_j];

                    // Calculate point displacement
                    if constexpr( qtf_type == QTFTypeE::QTF_DIFF_CODE )
                    {
                        int_mod_2d[k]   = 0.25 * (
                                                    point_disp_i[0] * std::conj( vel_x_acc_j )
                                                    +
                                                    point_disp_i[1] * std::conj( vel_y_acc_j )
                                                    +
                                                    point_disp_i[2] * std::conj( vel_z_acc_j )
                                                    +
                                                    std::conj( point_disp_j[0] ) * vel_x_acc_i
                                                    +
                                                    std::conj( point_disp_j[1] ) * vel_y_acc_i
                                                    +
                                                    std::conj( point_disp_j[2] ) * vel_z_acc_i
                                                );
                    }
                    else if constexpr( qtf_type == QTFTypeE::QTF_SUM_CODE )
                    {
                        int_mod_2d[k]   = 0.25 * (
                                                    point_disp_i[0] * vel_x_acc_j
                                                    +
                                                    point_disp_i[1] * vel_y_acc_j
                                                    +
                                                    point_disp_i[2] * vel_z_acc_j
                                                    +
                                                    point_disp_j[0] * vel_x_acc_i
                                                    +
                                                    point_disp_j[1] * vel_y_acc_i
                                                    +
                                                    point_disp_j[2] * vel_z_acc_i
                                                );
                    }
                }

                // Integrate over the panel
                int_val = 0.0;
                gauss2d_loop<ngp>( int_val, int_mod_2d, paneld->panel_geom );

                // Get global QTF index offset
                idx0    =   _qtf_index_offset(
                                                is_multi_head,
                                                bodies_np,
                                                heads_np,
                                                dofs_np,
                                                ih1,
                                                ih2,
                                                paneld->body_id
                                            );

                // Calculate final qtf value and store appropriately
                for ( std::size_t r=0; r<dofs_np; r++ )
                {
                    qtf_values[idx0+r]  += int_val * normal_vec[r];
                    qtf_acc[idx0+r]     += int_val * normal_vec[r];
                }
            }
        }
    }

    // Calculate second order force due to momentum
    cusfloat    ang_i_2                         = pow2s( ang_freq_i );
    cusfloat    ang_j_2                         = pow2s( ang_freq_j );
    cuscomplex  conj_vec[3];                    clear_vector( 3, conj_vec );
    cuscomplex  hydro_force_i[input->dofs_np];  clear_vector( input->dofs_np, hydro_force_i );
    cuscomplex  hydro_force_j[input->dofs_np];  clear_vector( input->dofs_np, hydro_force_j );
    cuscomplex  mom_i[3];                       clear_vector( 3, mom_i );

    for ( std::size_t ih1=0; ih1<heads_np; ih1++ )
    {
        for ( std::size_t ih2=ih2_start; ih2<ih2_end; ih2++ )
        {
            for ( std::size_t j=0; j<bodies_np; j++ )
            {
                // Define chunk index
                idx0    = ih1 * ( dofs_np * bodies_np ) + j * dofs_np;

                if ( is_multi_head )
                {
                    idx1_i  = ih1 * ( dofs_np * bodies_np ) + j * dofs_np;
                    idx1_j  = ih2 * ( dofs_np * bodies_np ) + j * dofs_np;
                }
                else
                {
                    idx1_i  = ih1 * ( dofs_np * bodies_np ) + j * dofs_np;
                    idx1_j  = ih1 * ( dofs_np * bodies_np ) + j * dofs_np;
                }

                // Calculate total hydrodynamic force
                calculate_mass_acceleration(
                                                input,
                                                raos_i,
                                                ang_i_2,
                                                j,
                                                idx1_i,
                                                hydro_force_i
                                            );

                calculate_mass_acceleration(
                                                input,
                                                raos_j,
                                                ang_j_2,
                                                j,
                                                idx1_j,
                                                hydro_force_j
                                            );

                // Add moment due to translational forces
                if constexpr( qtf_type == QTFTypeE::QTF_DIFF_CODE )
                {
                    cuscomplex      scale_f( 0.25, 0.0 );

                    clear_vector(   3,                  mom_i                                                       );
                    conj_vector(    3,                  hydro_force_j,          conj_vec                            );
                    cross(          &(raos_i[idx1_i+3]),conj_vec,               mom_i                               );
                    svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
                    sv_add(         3,                  &(qtf_values[idx0]),    mom_i,      &(qtf_values[idx0])     );
                    sv_add(         3,                  &(qtf_mom[idx0]),       mom_i,      &(qtf_mom[idx0])        );

                    clear_vector(   3,                  mom_i                                                       );
                    conj_vector(    3,                  &(raos_j[idx1_j+3]),    conj_vec                            );
                    cross(          conj_vec,           hydro_force_i,          mom_i                               );
                    svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
                    sv_add(         3,                  &(qtf_values[idx0]),    mom_i,      &(qtf_values[idx0])     );
                    sv_add(         3,                  &(qtf_mom[idx0]),       mom_i,      &(qtf_mom[idx0])        );


                    // Add moment due to rotational force
                    clear_vector(   3,                  mom_i                                                       );
                    conj_vector(    3,                  &(hydro_force_j[3]),    conj_vec                            );
                    cross(          &(raos_i[idx1_i+3]),conj_vec,               mom_i                               );
                    svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
                    sv_add(         3,                  &(qtf_values[idx0+3]),  mom_i,      &(qtf_values[idx0+3])   );
                    sv_add(         3,                  &(qtf_mom[idx0+3]),     mom_i,      &(qtf_mom[idx0+3])      );

                    clear_vector(   3,                  mom_i                                                       );
                    conj_vector(    3,                  &(raos_j[idx1_j+3]),    conj_vec                            );
                    cross(          conj_vec,           &(hydro_force_i[3]),    mom_i                               );
                    svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
                    sv_add(         3,                  &(qtf_values[idx0+3]),  mom_i,      &(qtf_values[idx0+3])   );
                    sv_add(         3,                  &(qtf_mom[idx0+3]),     mom_i,      &(qtf_mom[idx0+3])      );

                }
                else if constexpr( qtf_type == QTFTypeE::QTF_SUM_CODE )
                {
                    cuscomplex      scale_f( 0.25, 0.0 );

                    clear_vector(   3,                  mom_i                                                       );
                    cross(          &(raos_i[idx1_i+3]),hydro_force_j,          mom_i                               );
                    svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
                    sv_add(         3,                  &(qtf_values[idx0]),    mom_i,      &(qtf_values[idx0])     );
                    sv_add(         3,                  &(qtf_mom[idx0]),       mom_i,      &(qtf_mom[idx0])        );

                    clear_vector(   3,                  mom_i                                                       );
                    cross(          &(raos_j[idx1_j+3]),hydro_force_i,          mom_i                               );
                    svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
                    sv_add(         3,                  &(qtf_values[idx0]),    mom_i,      &(qtf_values[idx0])     );
                    sv_add(         3,                  &(qtf_mom[idx0]),       mom_i,      &(qtf_mom[idx0])        );

                    // Add moment due to rotational force
                    clear_vector(   3,                  mom_i                                                       );
                    cross(          &(raos_i[idx1_i+3]),&(hydro_force_j[3]),    mom_i                               );
                    svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
                    sv_add(         3,                  &(qtf_values[idx0+3]),  mom_i,      &(qtf_values[idx0+3])   );
                    sv_add(         3,                  &(qtf_mom[idx0+3]),     mom_i,      &(qtf_mom[idx0+3])      );

                    clear_vector(   3,                  mom_i                                                       );
                    cross(          &(raos_j[idx1_j+3]),&(hydro_force_i[3]),    mom_i                               );
                    svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
                    sv_add(         3,                  &(qtf_values[idx0+3]),  mom_i,      &(qtf_values[idx0+3])   );
                    sv_add(         3,                  &(qtf_mom[idx0+3]),     mom_i,      &(qtf_mom[idx0+3])      );

                }
                
            }
        }
    }

    // Sum-up contributions from all processes
    auto reduce_qtf = [qtf_size](cuscomplex* buffer)
    {
        MPI_Allreduce(
                            MPI_IN_PLACE,
                            buffer,
                            qtf_size,
                            mpi_cuscomplex,
                            MPI_SUM,
                            MPI_COMM_WORLD
                        );
    };

    reduce_qtf( qtf_values );
    reduce_qtf( qtf_wl     );
    reduce_qtf( qtf_bern   );
    reduce_qtf( qtf_acc    );
    reduce_qtf( qtf_mom    );

}


template<std::size_t N, int mode_pf>
template<QTFTypeE qtf_type>
void FrequencySolver<N, mode_pf>::_calculate_and_distribute_qtf(
                                                                    std::size_t freq_index_i,
                                                                    std::size_t freq_index_j
                                                                )
{
    this->_calculate_quadratic_terms<qtf_type>( freq_index_i, freq_index_j, true );

    if ( this->mpi_config->is_root( ) )
    {
        auto distribute = [&](cuscomplex* src, cuscomplex* dst, int is_sop, int is_qtf)
        {
            qtf_distribute_matrix_data( this->input, freq_index_i, freq_index_j, src, dst, is_sop, is_qtf );
        };

        if constexpr ( qtf_type == QTFTypeE::QTF_DIFF_CODE )
        {
            distribute( this->sim_data->qtf, this->sim_data->qtf_diff_freqs, 0, 1 );
            if ( input->out_qtf_comp )
            {
                distribute( this->sim_data->qtf_diff_wl,   this->sim_data->qtf_diff_wl_freqs,   0, 0 );
                distribute( this->sim_data->qtf_diff_bern, this->sim_data->qtf_diff_bern_freqs, 0, 0 );
                distribute( this->sim_data->qtf_diff_acc,  this->sim_data->qtf_diff_acc_freqs,  0, 0 );
                distribute( this->sim_data->qtf_diff_mom,  this->sim_data->qtf_diff_mom_freqs,  0, 0 );
            }
        }
        else
        {
            distribute( this->sim_data->qtf, this->sim_data->qtf_sum_freqs, 0, 1 );
            if ( this->input->out_qtf_comp )
            {
                distribute( this->sim_data->qtf_sum_wl,   this->sim_data->qtf_sum_wl_freqs,   0, 0 );
                distribute( this->sim_data->qtf_sum_bern, this->sim_data->qtf_sum_bern_freqs, 0, 0 );
                distribute( this->sim_data->qtf_sum_acc,  this->sim_data->qtf_sum_acc_freqs,  0, 0 );
                distribute( this->sim_data->qtf_sum_mom,  this->sim_data->qtf_sum_mom_freqs,  0, 0 );
            }
        }
    }

    MPI_Barrier( MPI_COMM_WORLD );
}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_check_system_mesh( void )
{
    if ( this->mpi_config->is_root( ) )
    {
        namespace fs = std::filesystem;
        fs::path report_path = fs::path( this->input->case_fopath ) / "mesh_quality_report.txt";
        std::ofstream report_file( report_path.string( ) );
        if ( report_file.is_open( ) )
        {
            cusfloat min_lambda = 0.0;
            cusfloat max_k = 0.0;
            this->input->get_wave_quality_params( min_lambda, max_k );
            report_file << "SeaMotions mesh quality report\n\n";
            for ( int i=0; i<this->input->bodies_np; i++ )
            {
                if ( this->input->bodies[i]->mesh_total != nullptr )
                {
                    this->input->bodies[i]->mesh_total->check_quality( min_lambda, max_k, report_file );
                }
            }
            report_file.close( );
        }
        else
        {
            Logger warn_logger( this->mpi_config->is_root( ) );
            warn_logger.warn( "Could not write mesh quality report to: " + report_path.string( ) );
        }
    }
}


template<std::size_t N, int mode_pf>
FrequencySolver<N, mode_pf>::FrequencySolver( Input* input_in, MpiConfig* mpi_config_in )
{
    // Storage input sytem pointer
    this->input = input_in;

    // Create MPI environment configuration
    this->mpi_config = mpi_config_in;

    // Prepare folder for output files
    this->_prepare_results_folder( );

    // Calculate hydrostatics
    this->_calculate_hydrostatics( );

    // Create output system
    this->_initialize_output_system( );
    
    // Create mesh group
    this->_initialize_mesh_groups( );

    // Check mesh quality and compatibility with the system
    this->_check_system_mesh( );

    // Allocate space for the simulation data
    this->sim_data  = new SimulationData(
                                            this->input,
                                            this->mpi_config,
                                            this->mesh_gp->meshes_np,
                                            this->input->dofs_np,
                                            this->input->heads_np,
                                            this->mesh_gp->source_nodes_tnp
                                        );

    // Generate formulation kernel
    this->_generate_formulation_kernel( );

    // Calculate mass and hydrostatic stiffness matrixes
    // used for the calculation of the RAOs.
    this->_calculate_global_static_matrixes( );

    // Initialize field data containers if any
    this->_initialize_field_data( );

    // Initialize Morison elements if any
    this->_initialize_morison_elements( );

}


template<std::size_t N, int mode_pf>
FrequencySolver<N, mode_pf>::~FrequencySolver( void )
{
    // Delete Morison elements if any
    for ( std::size_t i=0; i<static_cast<std::size_t>( this->input->bodies_np ); i++ )
    {
        if ( this->input->bodies[i]->use_morison_elements )
        {
            for ( std::size_t j=0; j<this->input->bodies[i]->morison_elements_np; j++ )
            {
                delete this->_morison_elements[i][j];
            }
        }
    }

    // Delete field data if any
    if ( this->input->is_calc_mdrift )
    {
        if ( this->_qtf_wl_fields != nullptr )
        {
            delete this->_qtf_wl_fields;
        }

        if ( this->_qtf_bern_fields != nullptr )
        {
            delete this->_qtf_bern_fields;
        }

    }

    if ( this->_qtf_fs_fields != nullptr )
    {
        delete this->_qtf_fs_fields;
        this->_qtf_fs_fields = nullptr;
    }

    if ( this->_qtf_fs_wl_fields != nullptr )
    {
        delete this->_qtf_fs_wl_fields;
        this->_qtf_fs_wl_fields = nullptr;
    }

    for ( auto* annulus_field : this->_qtf_annuli_fields )
    {
        delete annulus_field;
    }
    this->_qtf_annuli_fields.clear();

    for ( auto* annulus_mesh : this->_qtf_annuli_meshes )
    {
        delete annulus_mesh;
    }
    this->_qtf_annuli_meshes.clear();
    this->_qtf_annuli_weights.clear();

    for ( std::size_t i=0; i<this->input->field_points_np; i++ )
    {
        delete this->_field_points[i];
    }

    // Delete simulation data
    delete this->sim_data;

    // Delete mesh groups
    delete this->mesh_gp;
    delete this->mesh_fs_qtf_gp;

    // Delete hydrostatic objects
    for ( int i=0; i<input->bodies_np; i++ )
    {
        delete this->hydrostatics[i];
    }
    delete [] this->hydrostatics;

    // Delete output system
    delete this->output;

}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_generate_formulation_kernel( void )
{
    LOG_TASK_SS( kernel, "Generating kernel..." )
    MpiTimer kernel_timer;
    this->kernel = new FormulationKernelBackend<NUM_GP, mode_pf, RecalcSteadyE::ON>( this->input, this->mpi_config, this->mesh_gp );
    LOG_TASK_TIME( kernel, kernel_timer )
}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_initialize_field_data( void )
{
    LOG_TASK_SS( fields, "Initialize field points data..." )
    MpiTimer fields_timer;

    // Initialize QTF waterline field points data container
    if ( 
            this->input->is_calc_mdrift 
            || 
            this->input->out_qtf_so_model == QTFSOModelE::INDIRECT
            || 
            this->input->out_qtf_so_model == QTFSOModelE::DIRECT
        )
    {
        this->_qtf_bern_fields  = new RadDiffData<cuscomplex, RDDQTFConfig>(
                                                                                this->mpi_config,
                                                                                this->mesh_gp,
                                                                                this->input->angfreqs_np,
                                                                                this->input->heads_np,
                                                                                this->input->dofs_np,
                                                                                false
                                                                            );

        this->_qtf_wl_fields    = new RadDiffData<cuscomplex, RDDQTFConfig>(
                                                                                    this->mpi_config,
                                                                                    this->mesh_gp,
                                                                                    this->input->angfreqs_np,
                                                                                    this->input->heads_np,
                                                                                    this->input->dofs_np,
                                                                                    true
                                                                            );
        
        
    }

    if ( this->input->out_qtf_so_model == QTFSOModelE::DIRECT || this->input->out_qtf_so_model == QTFSOModelE::INDIRECT )
    {
        if ( this->_qtf_fs_fields == nullptr && this->mesh_fs_qtf_gp != nullptr )
        {
            this->_qtf_fs_fields = new RadDiffData<cuscomplex, RDDQTFConfig>(
                                                                                this->mpi_config,
                                                                                this->mesh_fs_qtf_gp,
                                                                                this->input->angfreqs_np,
                                                                                this->input->heads_np,
                                                                                this->input->dofs_np,
                                                                                false
                                                                            );
        }

        if ( this->input->out_qtf_so_model == QTFSOModelE::DIRECT && this->_qtf_fs_wl_fields == nullptr && this->mesh_fs_qtf_gp != nullptr )
        {
            this->_qtf_fs_wl_fields = new RadDiffData<cuscomplex, RDDQTFConfig>(
                                                                                    this->mpi_config,
                                                                                    this->mesh_fs_qtf_gp,
                                                                                    this->input->angfreqs_np,
                                                                                    this->input->heads_np,
                                                                                    this->input->dofs_np,
                                                                                    true
                                                                                );
        }

        if ( this->input->out_qtf_so_model == QTFSOModelE::DIRECT )
        {
            this->_qtf_annuli_fields.clear( );
            this->_qtf_annuli_fields.reserve( this->_qtf_annuli_meshes.size( ) );
        
            std::size_t fp_np_max = 0;
            for ( Mesh* annulus_mesh : this->_qtf_annuli_meshes )
            {
                fp_np_max = std::max( fp_np_max, static_cast<std::size_t>( annulus_mesh->nodes_np ) );
            }

            cusfloat* field_points = nullptr;
            if ( fp_np_max > 0 )
            {
                field_points = generate_empty_vector<cusfloat>( fp_np_max * 3 );
            }

            // Initialize a shared panel geometry for annulus fields using the first body COG.
            this->_qtf_annulus_geom.normal_vec[0] = 0.0;
            this->_qtf_annulus_geom.normal_vec[1] = 0.0;
            this->_qtf_annulus_geom.normal_vec[2] = 1.0;
            if ( this->input->bodies_np > 0 && this->input->bodies[0] != nullptr )
            {
                for ( int r=0; r<3; ++r )
                {
                    this->_qtf_annulus_geom.body_cog[r] = this->input->bodies[0]->cog[r];
                }
            }

            for ( std::size_t annulus_idx = 0; annulus_idx < this->_qtf_annuli_meshes.size( ); ++annulus_idx )
            {
                Mesh* annulus_mesh = this->_qtf_annuli_meshes[annulus_idx];
                const std::size_t fp_np = static_cast<std::size_t>( annulus_mesh->nodes_np );
                for ( std::size_t i=0; i<fp_np; ++i )
                {
                    field_points[ 3*i + 0 ] = annulus_mesh->x[ i ];
                    field_points[ 3*i + 1 ] = annulus_mesh->y[ i ];
                    field_points[ 3*i + 2 ] = annulus_mesh->z[ i ];
                }

                cusfloat* annulus_weights = nullptr;
                if ( annulus_idx < this->_qtf_annuli_weights.size( ) && !this->_qtf_annuli_weights[annulus_idx].empty( ) )
                {
                    annulus_weights = this->_qtf_annuli_weights[annulus_idx].data( );
                }

                this->_qtf_annuli_fields.emplace_back(
                                                            new RadDiffData<cuscomplex, RDDQTFConfig>(
                                                                                                        this->mpi_config,
                                                                                                        field_points,
                                                                                                        annulus_weights,
                                                                                                        fp_np,
                                                                                                        this->input->angfreqs_np,
                                                                                                        this->input->heads_np,
                                                                                                        this->input->dofs_np,
                                                                                                        false
                                                                                                    )
                                                        );

                RDDQC* annulus_field = this->_qtf_annuli_fields.back( );
                for ( std::size_t i=0; i<annulus_field->get_size_local( ); ++i )
                {
                    annulus_field->panel_data[i].panel_geom = &(this->_qtf_annulus_geom);
                }
            }

            if ( field_points != nullptr )
            {
                mkl_free( field_points );
            }

            if ( this->kernel != nullptr )
            {
                this->kernel->set_qtf_annuli_fields( &(this->_qtf_annuli_fields) );
            }
        }
    }

    if ( this->input->out_qtf_so_model == QTFSOModelE::INDIRECT && this->mpi_config->is_root( ) )
    {
        std::size_t kochin_heads_np = static_cast<std::size_t>( this->input->heads_np * this->input->bodies_np * this->input->kochin_np );
        std::size_t kochin_rad_np   = static_cast<std::size_t>( this->input->dofs_np * this->input->bodies_np * this->input->kochin_np );
        std::size_t freqs_np        = static_cast<std::size_t>( this->input->angfreqs_np );

        this->sim_data->qtf_kochin_heads_np = static_cast<int>( kochin_heads_np );
        this->sim_data->qtf_kochin_rad_np   = static_cast<int>( kochin_rad_np );

        if ( this->sim_data->qtf_kochin_pert_cos_freqs == nullptr )
        {
            this->sim_data->qtf_kochin_pert_cos_freqs = generate_empty_vector<cuscomplex>( kochin_heads_np * freqs_np );
            this->sim_data->qtf_kochin_pert_sin_freqs = generate_empty_vector<cuscomplex>( kochin_heads_np * freqs_np );
            this->sim_data->qtf_kochin_rad_cos_freqs  = generate_empty_vector<cuscomplex>( kochin_rad_np * freqs_np );
            this->sim_data->qtf_kochin_rad_sin_freqs  = generate_empty_vector<cuscomplex>( kochin_rad_np * freqs_np );
        }
    }

    // Load field points data
    namespace   fs                      = std::filesystem;
    fs::path    case_fopath_            = fs::path( this->input->case_fopath );
    fs::path    field_points_foname_    = fs::path( RESULTS_FOLDER_NAME ) / fs::path( RESULTS_FIELD_POINTS_FOLDER_NAME );

    this->_field_points.reserve( this->input->field_points_np );
    for ( std::size_t i=0; i<this->input->field_points_np; i++ )
    {
        // Define field points output folder path and create if not exist
        fs::path fd_out_fopath_ = case_fopath_ / field_points_foname_ / fs::path( this->input->field_points[i].name );
        if ( !fs::exists( fd_out_fopath_ ) )
        {
            fs::create_directories( fd_out_fopath_ );
        }

        // Create new field point object and store in the vector
        this->_field_points.emplace_back( 
                                                new FMD(
                                                            this->input,
                                                            &(this->input->field_points[i]),
                                                            fd_out_fopath_.string( ),
                                                            this->mpi_config
                                                        )
                                            );
    }

    LOG_TASK_TIME( fields, fields_timer )
}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_initialize_mesh_groups( void )
{
    LOG_TASK_SS( mesh, "Initialize mesh groups..." )
    MpiTimer mesh_timer;
    
    // Group all meshes in a vector
    Mesh** all_meshes = new Mesh*[input->bodies_np];
    for ( int i=0; i<this->input->bodies_np; i++ )
    {
        all_meshes[i] = this->input->bodies[i]->is_mesh_total
                            ? this->input->bodies[i]->mesh_total
                            : nullptr;
    }

    // Create new mesh from the meshes of all objects
    this->mesh_gp       = new MeshGroup( 
                                            all_meshes,
                                            input->bodies_np,
                                            input->is_wl_points
                                        );

    Mesh**      fs_mesh = nullptr;
    if ( this->input->is_fs_qtf )
    {
        // Create vector with the QTF FS mesh to feed the
        // MeshGroup object constructor
        fs_mesh         = new Mesh*[1];
        fs_mesh[0]      = this->input->bodies[0]->mesh_fs_qtf;

        // Create partition circle water line if required
        if ( this->input->out_qtf_so_model == QTFSOModelE::DIRECT )
        {
            fs_mesh[0]->detect_pc_points( this->input->wl_det_prec );
        }

        // Create mesh group
        this->mesh_fs_qtf_gp    = new MeshGroup(
                                                    fs_mesh,
                                                    1,
                                                    true
                                                );

        if ( this->input->out_qtf_so_model == QTFSOModelE::DIRECT )
        {
            const cusfloat r0 = this->input->bodies[0]->mesh_fs_qtf->get_fs_radius( );
            const cusfloat r1 = static_cast<cusfloat>( 2.0 ) * r0;
            int annuli_np = this->input->qtf_annuli_np;
            if ( annuli_np < 1 )
            {
                annuli_np = 1;
            }
            const cusfloat dr = ( r1 - r0 ) / static_cast<cusfloat>( annuli_np );

            this->_qtf_annuli_meshes.clear( );
            this->_qtf_annuli_meshes.reserve( annuli_np );
            this->_qtf_annuli_weights.clear( );
            this->_qtf_annuli_weights.reserve( annuli_np );
            for ( int k=0; k<annuli_np; ++k )
            {
                const cusfloat r_in = r0 + static_cast<cusfloat>( k ) * dr;
                const cusfloat r_out = r_in + dr;
                AnnulusNodeCloud cloud = create_annulus_nodes_runtime(
                                            r_in,
                                            r_out,
                                            this->input->qtf_annuli_theta_np,
                                            this->input->gauss_order,
                                            AngularRule::Trapezoidal
                                        );

                this->_qtf_annuli_weights.emplace_back(
                                                        cloud.w.begin( ),
                                                        cloud.w.end( )
                                                    );

                Mesh annulus_mesh = Mesh::from_annulus_nodes( cloud, "QTFAnnulus" );
                this->_qtf_annuli_meshes.emplace_back( new Mesh( std::move( annulus_mesh ) ) );
            }
        }
    }

    this->mesh_gp->define_mirror_panels( );

    // Delete allocated heap memory
    delete []   all_meshes;
    delete []   fs_mesh;

    LOG_TASK_TIME( mesh, mesh_timer )

}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_initialize_morison_elements( void )
{
    LOG_TASK_SS( morison, "Initialize Morison elements..." )
    MpiTimer morison_timer;

    // Initialize Morison elements
    this->_morison_elements.reserve( this->input->bodies_np );
    for ( std::size_t i=0; i<static_cast<std::size_t>( this->input->bodies_np ); i++ )
    {
        this->_morison_elements.emplace_back( );
        for ( std::size_t j=0; j<this->input->bodies[i]->morison_elements_np; j++ )
        {
            this->_morison_elements[i].emplace_back( 
                                                        new MorisonElement(
                                                                                this->mpi_config,
                                                                                this->input,
                                                                                &(this->input->bodies[i]->morison_elements[j])
                                                                    )
                                            );
        }
    }

    // Reserve memory for auxiliary Morison elements forces
    this->_drag_force_aux       = cut::CusTensor<cuscomplex>( { this->input->heads_np, this->input->dofs_np } );
    this->_inertial_force_aux   = cut::CusTensor<cuscomplex>( { this->input->heads_np, this->input->dofs_np } );

    LOG_TASK_TIME( morison, morison_timer )
}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_initialize_output_system( void )
{
    LOG_TASK_SS( output, "Creating output system..." )
    MpiTimer output_timer;

    // Create output system
    if ( this->mpi_config->is_root( ) )
    {
        this->output = new Output( this->input );
    }

    // Storage initial parameters
    if ( this->mpi_config->is_root( ) )
    {
        // Storage frequency set
        this->output->save_frequencies( this->input->freqs );

        // Storage headings set
        this->output->save_headings( this->input->heads.data( ) );

        // Storage structural mass
        if ( this->input->out_struct_mass )
        {
            this->output->save_structural_mass( );
        }

        // Storage hydrostatic stiffness matrix
        if ( this->input->out_hydstiff )
        {
            this->output->save_hydstiffness( this->hydrostatics );
        }

        // Storage mesh
        if ( this->input->out_mesh )
        {
            this->output->save_mesh( );
        }

    }

    LOG_TASK_TIME( output, output_timer )

}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_prepare_results_folder( void )
{
    if ( this->mpi_config->is_root( ) )
    {
        namespace fs = std::filesystem;

        // Define folder path
        fs::path folder_path_( std::string( this->input->case_fopath ) );

        // Check if results folder exist and create if not
        fs::path results_foname_{ std::string( RESULTS_FOLDER_NAME ) };
        fs::path results_mesh_fipath_ = folder_path_ / results_foname_;
        if ( !fs::exists( results_mesh_fipath_) )
        {
            fs::create_directory( results_mesh_fipath_ );
        }

        // Check if mesh folder exist and create if not
        fs::path mesh_foname_1_{ std::string( RESULTS_MESH_FOLDER_NAME ) };
        fs::path plot_mesh_fipath_ = results_mesh_fipath_ / mesh_foname_1_;
        if ( !fs::exists( plot_mesh_fipath_) )
        {
            fs::create_directory( plot_mesh_fipath_ );
        }

        // Check if field points folder exist and create if not
        fs::path field_points_foname_{ std::string( RESULTS_FIELD_POINTS_FOLDER_NAME ) };
        fs::path plot_field_points_fipath_ = results_mesh_fipath_ / field_points_foname_;
        if ( !fs::exists( plot_field_points_fipath_) )
        {
            fs::create_directory( plot_field_points_fipath_ );
        }
    }
}

template<std::size_t N, int mode_pf>
template<QTFTypeE qtf_type>
void FrequencySolver<N, mode_pf>::_save_qtf( void ) const
{
    if ( this->mpi_config->is_root( ) )
    {
        const char* qtf_name = nullptr;
        cuscomplex* qtf_freqs = nullptr;

        if constexpr ( qtf_type == QTFTypeE::QTF_DIFF_CODE )
        {
            qtf_name = "qtf_diff";
            qtf_freqs = this->sim_data->qtf_diff_freqs;
        }
        else
        {
            qtf_name = "qtf_sum";
            qtf_freqs = this->sim_data->qtf_sum_freqs;
        }

        // Save data into the disk file
        this->output->save_qtf_format( qtf_name, qtf_freqs );

        if ( this->input->out_qtf_comp )
        {
            if constexpr ( qtf_type == QTFTypeE::QTF_DIFF_CODE )
            {
                // Storage QTF acceleration term
                this->output->save_qtf_format(
                                                "qtf_diff_acc",
                                                this->sim_data->qtf_diff_acc_freqs
                                            );

                // Storage QTF bernoulli term
                this->output->save_qtf_format(
                                                "qtf_diff_bern",
                                                this->sim_data->qtf_diff_bern_freqs
                                            );

                // Storage QTF momentum term
                this->output->save_qtf_format(
                                                "qtf_diff_mom",
                                                this->sim_data->qtf_diff_mom_freqs
                                            );

                // Storage QTF second order potential term
                this->output->save_qtf_format(
                                                "qtf_diff_sop",
                                                this->sim_data->qtf_diff_secord_force_freqs
                                            );

                // Storage QTF wl term
                this->output->save_qtf_format(
                                                "qtf_diff_wl",
                                                this->sim_data->qtf_diff_wl_freqs
                                            );
            }
            else
            {
                // Storage QTF acceleration term
                this->output->save_qtf_format(
                                                "qtf_sum_acc",
                                                this->sim_data->qtf_sum_acc_freqs
                                            );

                // Storage QTF bernoulli term
                this->output->save_qtf_format(
                                                "qtf_sum_bern",
                                                this->sim_data->qtf_sum_bern_freqs
                                            );

                // Storage QTF momentum term
                this->output->save_qtf_format(
                                                "qtf_sum_mom",
                                                this->sim_data->qtf_sum_mom_freqs
                                            );

                // Storage QTF second order potential term
                this->output->save_qtf_format(
                                                "qtf_sum_sop",
                                                this->sim_data->qtf_sum_secord_force_freqs
                                            );

                // Storage QTF wl term
                this->output->save_qtf_format(
                                                "qtf_sum_wl",
                                                this->sim_data->qtf_sum_wl_freqs
                                            );
            }
        }
    }
}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::write_memory_report( const std::string& file_path ) const
{
    if ( this->mpi_config == nullptr || !this->mpi_config->is_root( ) )
    {
        return;
    }

    std::ofstream report_file( file_path );
    if ( !report_file.is_open( ) )
    {
        Logger warn_logger( this->mpi_config->is_root( ) );
        warn_logger.warn( "Could not write memory report to: " + file_path );
        return;
    }

    std::vector<MemoryReportEntry> entries;
    auto is_bad_ptr = [](const void* ptr)
    {
        if ( ptr == nullptr )
        {
            return true;
        }
        const std::uintptr_t val = reinterpret_cast<std::uintptr_t>( ptr );
        const std::uint32_t  low = static_cast<std::uint32_t>( val & 0xFFFFFFFFu );
        return val < 0x10000ull
            || val == 0xFFFFFFFFFFFFFFFFull
            || low == 0xFEEEFEEEu
            || low == 0xDDDDDDDDu
            || low == 0xCDCDCDCDu
            || low == 0xCCCCCCCCu;
    };

    if ( this->sim_data != nullptr )
    {
        this->sim_data->append_memory_report( entries, "SimulationData" );
    }
    if ( this->kernel != nullptr )
    {
        this->kernel->append_memory_report( entries, "FormulationKernelBackend" );
    }

    if ( this->mesh_gp != nullptr )
    {
        this->mesh_gp->append_memory_report( entries, "MeshGroup" );
        if ( this->mesh_gp->meshes != nullptr )
        {
            for ( int i = 0; i < this->mesh_gp->meshes_np; i++ )
            {
                if ( !is_bad_ptr( this->mesh_gp->meshes[i] ) )
                {
                    this->mesh_gp->meshes[i]->append_memory_report(
                                                                    entries,
                                                                    "MeshGroup.Mesh[" + std::to_string( i ) + "]"
                                                                );
                }
            }
        }
    }

    if ( this->mesh_fs_qtf_gp != nullptr )
    {
        this->mesh_fs_qtf_gp->append_memory_report( entries, "MeshFsQtfGroup" );
        if ( this->mesh_fs_qtf_gp->meshes != nullptr )
        {
            for ( int i = 0; i < this->mesh_fs_qtf_gp->meshes_np; i++ )
            {
                if ( !is_bad_ptr( this->mesh_fs_qtf_gp->meshes[i] ) )
                {
                    this->mesh_fs_qtf_gp->meshes[i]->append_memory_report(
                                                                            entries,
                                                                            "MeshFsQtfGroup.Mesh[" + std::to_string( i ) + "]"
                                                                        );
                }
            }
        }
    }

    for ( std::size_t i = 0; i < this->_field_points.size( ); i++ )
    {
        if ( this->_field_points[i] != nullptr )
        {
            this->_field_points[i]->append_memory_report(
                                                            entries,
                                                            "FieldMeshData[" + std::to_string( i ) + "]"
                                                        );
        }
    }

    if ( this->_qtf_wl_fields != nullptr )
    {
        this->_qtf_wl_fields->append_memory_report( entries, "QTFWlFields" );
    }
    if ( this->_qtf_bern_fields != nullptr )
    {
        this->_qtf_bern_fields->append_memory_report( entries, "QTFBernFields" );
    }
    if ( this->_qtf_fs_fields != nullptr )
    {
        this->_qtf_fs_fields->append_memory_report( entries, "QTFFsFields" );
    }
    if ( this->_qtf_fs_wl_fields != nullptr )
    {
        this->_qtf_fs_wl_fields->append_memory_report( entries, "QTFFsWlFields" );
    }
    for ( std::size_t i = 0; i < this->_qtf_annuli_fields.size( ); i++ )
    {
        if ( this->_qtf_annuli_fields[i] != nullptr )
        {
            this->_qtf_annuli_fields[i]->append_memory_report(
                                                                entries,
                                                                "QTFAnnulusFields[" + std::to_string( i ) + "]"
                                                            );
        }
    }

    for ( std::size_t i = 0; i < this->_morison_elements.size( ); i++ )
    {
        for ( std::size_t j = 0; j < this->_morison_elements[i].size( ); j++ )
        {
            if ( this->_morison_elements[i][j] != nullptr )
            {
                this->_morison_elements[i][j]->append_memory_report(
                                                                        entries,
                                                                        "MorisonElement[" + std::to_string( i ) + "][" + std::to_string( j ) + "]"
                                                                    );
            }
        }
    }

    write_memory_report_csv( report_file, entries );
}