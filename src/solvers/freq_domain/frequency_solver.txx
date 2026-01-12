
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
#include <iomanip>

// Include local module
#include "../../containers/logger.hpp"
#include "../../containers/mpi_timer.hpp"
#include "diffraction.hpp"
#include "frequency_solver.hpp"
#include "froude_krylov.hpp"
#include "global_static_matrix.hpp"
#include "hydromechanics.hpp"
#include "../../hydrostatics.hpp"
#include "../../math/integration.hpp"
#include "../../math/math_tools.hpp"
#include "../../mesh/mesh.hpp"
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


/*****************************************************************/
/*************** Define Frequency Solver class  ******************/
/*****************************************************************/
template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_calculate_field_points_values( 
                                                                    std::size_t freq_index,
                                                                    cusfloat    ang_freq 
                                                                )
{
    if ( this->input->is_calc_mdrift )
    {
        // Calculate waterline field points values for QTF calculations
        this->kernel->template compute_fields<RDDQTFConfig>( 
                                                                freq_index,
                                                                ang_freq,
                                                                this->sim_data->raos,
                                                                this->_qtf_wl_fields
                                                            );

        // Calculate velocity field at Bernoulli points for QTF calculations
        this->kernel->template compute_fields<RDDQTFConfig>(
                                                                freq_index, 
                                                                ang_freq,
                                                                this->sim_data->raos,
                                                                this->_qtf_bern_fields
                                                            );
    }
}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::calculate_first_order( void )
{
    /*******************************************************/
    /*   Calculate regular frequency band coefficients     */
    /*******************************************************/

    // Loop over frequencies
    for ( std::size_t i=0; i<static_cast<std::size_t>(this->input->angfreqs_np); i++ )
    {
        LOG_TASK_SS( freq, "FO - Calculating period: " << std::setw( 10 ) << std::fixed << std::setprecision( 3 ) << angfreq_to_period( this->input->angfreqs[i] ) << " s" )

        // Start time measurement for the current frequency
        MpiTimer freq_timer;

        // Solve radiation-diffraction problem for the current frequency
        this->kernel->template solve<freq_regime_t::REGULAR>( input->angfreqs[i] );
        this->kernel->update_results( sim_data );

        // Calculate first order coefficients
        this->_calculate_first_order_coeffs<freq_regime_t::REGULAR>( i, this->input->angfreqs[i] );

        // Calculate field points values at the positions required
        this->_calculate_field_points_values( i, this->input->angfreqs[i] );

        // Calculate second order coefficients given by first order potential solution
        this->_calculate_mean_drift( i, false );
        
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
    this->kernel->template solve<freq_regime_t::ASYMPT_LOW>( W_ASYMPT_LOW );
    this->kernel->update_results( sim_data );

    // Calculate first order coefficients
    this->_calculate_first_order_coeffs<freq_regime_t::ASYMPT_LOW>( -1, W_ASYMPT_LOW );

    // Print out execution times
    LOG_TASK_TIME( freq_lf, freq_low_timer )

    /*******************************************************/
    /*   Calculate high frequency asymptotic coefficients  */
    /*******************************************************/
    LOG_TASK_SS( freq_hf, "FO - Calculating period: " << std::setw( 10 ) << std::fixed << std::setprecision( 3 ) << angfreq_to_period( W_ASYMPT_HIGH ) << " s" )

    // Start time measurement
    MpiTimer freq_high_timer;

    // Solve radiation-diffraction problem for the current frequency
    this->kernel->template solve<freq_regime_t::ASYMPT_HIGH>( W_ASYMPT_HIGH );
    this->kernel->update_results( sim_data );

    // Calculate first order coefficients
    this->_calculate_first_order_coeffs<freq_regime_t::ASYMPT_HIGH>( -1, W_ASYMPT_HIGH );

    // Print out execution times
    LOG_TASK_TIME( freq_hf, freq_high_timer )

}


template<std::size_t N, int mode_pf>
template<freq_regime_t freq_regime>
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
    if ( freq_regime == freq_regime_t::REGULAR || freq_regime == freq_regime_t::ASYMPT_LOW )
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
    if ( this->mpi_config->is_root( )  && freq_regime == freq_regime_t::REGULAR )
    {
        calculate_raos(
                            this->input,
                            this->sim_data->structural_mass_p0,
                            this->sim_data->added_mass_p0,
                            this->sim_data->damping_rad_p0,
                            this->sim_data->hydrostiff_p0,
                            this->sim_data->wave_diffrac_p0,
                            this->sim_data->froude_krylov_p0,
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

    // Copy raos to raos_hist
    if ( this->input->is_calc_mdrift || this->input->out_qtf )
    {
        copy_vector(
                        this->sim_data->wave_exc_np,
                        this->sim_data->raos,
                        &(this->sim_data->raos_hist[freq_index * this->sim_data->wave_exc_np])
                    );
    }

    
    // Storage results
    // MpiTimer storage_timer;       
    if ( this->input->out_sources && freq_regime == freq_regime_t::REGULAR )
    {
        MPI_Reduce(
                        this->sim_data->intensities,
                        this->sim_data->intensities_p0,
                        this->kernel->size( ) * ( this->input->dofs_np + this->input->heads_np ),
                        mpi_cuscomplex,
                        MPI_SUM,
                        this->mpi_config->proc_root,
                        MPI_COMM_WORLD
                    );
    }

    if ( this->input->out_potential && freq_regime == freq_regime_t::REGULAR )
    {
        MPI_Reduce(
                        this->sim_data->panels_potential,
                        this->sim_data->panels_potential_p0,
                        this->kernel->size( ) * ( this->input->dofs_np + this->input->heads_np ),
                        mpi_cuscomplex,
                        MPI_SUM,
                        this->mpi_config->proc_root,
                        MPI_COMM_WORLD
                    );
    }

    if ( this->input->out_pressure && freq_regime == freq_regime_t::REGULAR )
    {
        MPI_Reduce(
                        this->sim_data->panels_pressure,
                        this->sim_data->panels_pressure_p0,
                        this->kernel->size( ) * ( this->input->dofs_np + this->input->heads_np ),
                        mpi_cuscomplex,
                        MPI_SUM,
                        this->mpi_config->proc_root,
                        MPI_COMM_WORLD
                    );
    }
    
    if ( this->mpi_config->is_root( ) )
    {
        if ( this->input->out_hydmech )
        {
            if constexpr( freq_regime == freq_regime_t::REGULAR )
            {
                this->output->save_hydromechanics_format(
                                                            freq_index,
                                                            _DN_ADDED_MASS,
                                                            this->sim_data->added_mass
                                                        );

                this->output->save_hydromechanics_format(
                                                            freq_index,
                                                            _DN_DAMPING_RAD,
                                                            this->sim_data->damping_rad
                                                        );
            }
            else
            {
                // Select dataset names depending on the frequency regime
                std::string added_mass_dn   = _DN_ADDED_MASS_LF;
                std::string damping_dn      = _DN_DAMPING_RAD_LF;
                if ( freq_regime == freq_regime_t::ASYMPT_HIGH )
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
            if constexpr( freq_regime == freq_regime_t::REGULAR )
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
                if ( freq_regime == freq_regime_t::ASYMPT_LOW )
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
            if constexpr( freq_regime == freq_regime_t::REGULAR )
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
                if ( freq_regime == freq_regime_t::ASYMPT_LOW )
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
            if constexpr( freq_regime == freq_regime_t::REGULAR )
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
                if ( freq_regime == freq_regime_t::ASYMPT_LOW )
                {
                    wex_dn   = _DN_WEX_LF;
                }

                this->output->save_wave_exciting_asympt_format(
                                                                    wex_dn,
                                                                    this->sim_data->wave_exc_p0
                                                                );
            }   
        }

        if ( this->input->out_raos && freq_regime == freq_regime_t::REGULAR )
        {
            this->output->save_wave_exciting_format(
                                                        freq_index,
                                                        _DN_RAO,
                                                        this->sim_data->raos
                                                    );
        }

        // Storage sources
        if ( this->input->out_sources && freq_regime == freq_regime_t::REGULAR )
        {
            this->output->save_fields_data( 
                                                freq_index,
                                                _DN_SRC_INT,
                                                this->sim_data->intensities_p0
                                            );
        }

        // Storage panels potential
        if ( this->input->out_potential && freq_regime == freq_regime_t::REGULAR )
        {
            this->output->save_fields_data( 
                                                freq_index,
                                                _DN_POT_INT,
                                                this->sim_data->panels_potential_p0
                                            );
        }

        // Storage panels pressure
        if ( this->input->out_pressure && freq_regime == freq_regime_t::REGULAR )
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
        this->_calculate_quadratic_terms<QTFTypeT::QTF_DIFF_CODE>( 
                                                                        freq_index,
                                                                        freq_index,
                                                                        is_multi_head 
                                                                    );

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
template<QTFTypeT qtf_type>
void FrequencySolver<N, mode_pf>::_calculate_quadratic_terms( 
                                                                    std::size_t freq_index_i,
                                                                    std::size_t freq_index_j,
                                                                    bool        is_multi_head
                                                            )
{
    // Get local variables from structures for easier access
    cusfloat                    ang_freq_i      = this->input->angfreqs[freq_index_i];
    cusfloat                    ang_freq_j      = this->input->angfreqs[freq_index_j];
    std::size_t                 bodies_np       = static_cast<std::size_t>( this->mesh_gp->meshes_np );
    std::size_t                 dofs_np         = static_cast<std::size_t>( this->input->dofs_np );
    cusfloat                    grav_acc        = this->input->grav_acc;
    std::size_t                 heads_np        = static_cast<std::size_t>( this->input->heads_np );
    Input*                      input           = this->input;
    PanelData<RDDQTFConfig>*    paneld          = nullptr;
    cuscomplex*                 qtf_values      = nullptr;
    cuscomplex*                 qtf_wl          = nullptr;
    cuscomplex*                 qtf_bern        = nullptr;
    cuscomplex*                 qtf_acc         = nullptr;
    cuscomplex*                 qtf_mom         = nullptr;
    cuscomplex*                 raos_i          = &( this->sim_data->raos_hist[ freq_index_i * this->sim_data->wave_exc_np ] );
    cuscomplex*                 raos_j          = &( this->sim_data->raos_hist[ freq_index_j * this->sim_data->wave_exc_np ] );
    RadDiffData<RDDQTFConfig>*  rdd_bern        = this->_qtf_bern_fields;
    RadDiffData<RDDQTFConfig>*  rdd_rwel        = this->_qtf_wl_fields;
    cusfloat                    rhow            = this->input->water_density;
    cuscomplex*                 vel_x           = nullptr;
    cuscomplex*                 vel_y           = nullptr;
    cuscomplex*                 vel_z           = nullptr;
    cusfloat                    wave_amplitude  = this->input->wave_amplitude;

    if constexpr( qtf_type == QTFTypeT::QTF_DIFF_CODE )
    {
        qtf_values   = this->sim_data->qtf;
        qtf_wl      = this->sim_data->qtf_diff_wl;
        qtf_bern    = this->sim_data->qtf_diff_bern;
        qtf_acc     = this->sim_data->qtf_diff_acc;
        qtf_mom     = this->sim_data->qtf_diff_mom;
    }
    else if constexpr ( qtf_type == QTFTypeT::QTF_SUM_CODE )
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
                    if constexpr( qtf_type == QTFTypeT::QTF_DIFF_CODE )
                    {
                        int_mod_1d[k] = paneld->wev_rel_total[idx1_i] * std::conj( paneld->wev_rel_total[idx1_j] );
                    }
                    else if constexpr( qtf_type == QTFTypeT::QTF_SUM_CODE )
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
                vel_x       = paneld->vel_x_total;
                vel_y       = paneld->vel_y_total;
                vel_z       = paneld->vel_z_total;

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

                    if constexpr( qtf_type == QTFTypeT::QTF_DIFF_CODE )
                    {
                        int_mod_2d[k]   =   (
                                                vel_x[idx1_i] * std::conj( vel_x[idx1_j] )
                                                +
                                                vel_y[idx1_i] * std::conj( vel_y[idx1_j] )
                                                +
                                                vel_z[idx1_i] * std::conj( vel_z[idx1_j] )
                                            );
                    }
                    else if constexpr( qtf_type == QTFTypeT::QTF_SUM_CODE )
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
                vel_x       = paneld->vel_x_total;
                vel_y       = paneld->vel_y_total;
                vel_z       = paneld->vel_z_total;

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
                    if constexpr( qtf_type == QTFTypeT::QTF_DIFF_CODE )
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
                    else if constexpr( qtf_type == QTFTypeT::QTF_SUM_CODE )
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
                if constexpr( qtf_type == QTFTypeT::QTF_DIFF_CODE )
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
                else if constexpr( qtf_type == QTFTypeT::QTF_SUM_CODE )
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
    MPI_Allreduce(
                        MPI_IN_PLACE,
                        qtf_values,
                        qtf_size,
                        mpi_cuscomplex,
                        MPI_SUM,
                        MPI_COMM_WORLD
                    );
    
    MPI_Allreduce(
                        MPI_IN_PLACE,
                        qtf_wl,
                        qtf_size,
                        mpi_cuscomplex,
                        MPI_SUM,
                        MPI_COMM_WORLD
                    );
                    
    MPI_Allreduce(
                        MPI_IN_PLACE,
                        qtf_bern,
                        qtf_size,
                        mpi_cuscomplex,
                        MPI_SUM,
                        MPI_COMM_WORLD
                    ); 

    MPI_Allreduce(
                        MPI_IN_PLACE,
                        qtf_acc,
                        qtf_size,
                        mpi_cuscomplex,
                        MPI_SUM,
                        MPI_COMM_WORLD
                    );

    MPI_Allreduce(
                        MPI_IN_PLACE,
                        qtf_mom,
                        qtf_size,
                        mpi_cuscomplex,
                        MPI_SUM,
                        MPI_COMM_WORLD
                    ); 

}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::calculate_second_order( void )
{
    // TO DO: Second order implementation
}


template<std::size_t N, int mode_pf>
FrequencySolver<N, mode_pf>::FrequencySolver( Input* input_in, MpiConfig* mpi_config_in )
{
    // Storage input sytem pointer
    this->input = input_in;

    // Create MPI environment configuration
    this->mpi_config = mpi_config_in;

    // Calculate hydrostatics
    this->_calculate_hydrostatics( );

    // Create output system
    this->_initialize_output_system( );
    
    // Create mesh group
    this->_initialize_mesh_groups( );

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

}


template<std::size_t N, int mode_pf>
FrequencySolver<N, mode_pf>::~FrequencySolver( void )
{
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
    this->kernel = new FormulationKernelBackend<NUM_GP, PF_OFF>( this->input, this->mpi_config, this->mesh_gp );
    LOG_TASK_TIME( kernel, kernel_timer )
}


template<std::size_t N, int mode_pf>
void FrequencySolver<N, mode_pf>::_initialize_field_data( void )
{
    LOG_TASK_SS( fields, "Initialize field points data..." )
    MpiTimer fields_timer;

    // Initialize QTF waterline field points data container
    if ( this->input->is_calc_mdrift )
    {
        this->_qtf_wl_fields    = new RadDiffData<RDDQTFConfig>(
                                                                        this->mpi_config,
                                                                        this->mesh_gp,
                                                                        this->input->angfreqs_np,
                                                                        this->input->heads_np,
                                                                        this->input->dofs_np,
                                                                        true
                                                                );
        
        this->_qtf_bern_fields  = new RadDiffData<RDDQTFConfig>(
                                                                        this->mpi_config,
                                                                        this->mesh_gp,
                                                                        this->input->angfreqs_np,
                                                                        this->input->heads_np,
                                                                        this->input->dofs_np,
                                                                        false
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
    for ( int i=0; i<input->bodies_np; i++ )
    {
        all_meshes[i] = this->input->bodies[i]->mesh;
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
        if ( this->input->out_qtf_so_model == 2 )
        {
            fs_mesh[0]->detect_pc_points( this->input->wl_det_prec );
        }

        // Create mesh group
        this->mesh_fs_qtf_gp    = new MeshGroup(
                                                    fs_mesh,
                                                    1,
                                                    false
                                                );
    }

    this->mesh_gp->define_mirror_panels( );

    // Delete allocated heap memory
    delete []   all_meshes;
    delete []   fs_mesh;

    LOG_TASK_TIME( mesh, mesh_timer )

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