
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

// Include local modules
#include "formulation_kernel_backend.hpp"
#include "../../green/integrals_database.hpp"
#include "panel_fields.hpp"


/********************************************************/
/************** Define module macros ********************/
/********************************************************/
#define _QUADRATURE_PANEL_T( TA0, TA1, TA2, TA3, TA4 )                                  \
{                                                                                       \
    if constexpr( freq_regime == freq_regime_t::REGULAR )                               \
    {                                                                                   \
        quadrature_panel_t<                                                             \
                                PanelGeom,                                              \
                                GWFcnsInterfaceT<NUM_GP2>,                              \
                                wave_term_integral<TA0, TA1, TA2, TA3, TA4>,            \
                                NUM_GP                                                  \
                            >(                                                          \
                                source_i->panel,                                        \
                                this->_gwfcns_interf,                                   \
                                wave_fcn_value,                                         \
                                wave_fcn_dn_sf_value,                                   \
                                wave_fcn_dn_pf_value,                                   \
                                wave_fcn_dx_value,                                      \
                                wave_fcn_dy_value,                                      \
                                wave_fcn_dz_value                                       \
                            );                                                          \
    }                                                                                   \
    else if constexpr( freq_regime == freq_regime_t::ASYMPT_LOW )                       \
    {                                                                                   \
        quadrature_panel_t<                                                             \
                                PanelGeom,                                              \
                                GWFcnsInterfaceT<NUM_GP2>,                              \
                                wave_term_integral_zero_freq<TA0, TA1, TA2, TA3, TA4>,  \
                                NUM_GP                                                  \
                            >(                                                          \
                                source_i->panel,                                        \
                                this->_gwfcns_interf,                                   \
                                wave_fcn_value,                                         \
                                wave_fcn_dn_sf_value,                                   \
                                wave_fcn_dn_pf_value,                                   \
                                wave_fcn_dx_value,                                      \
                                wave_fcn_dy_value,                                      \
                                wave_fcn_dz_value                                       \
                            );                                                          \
    }                                                                                   \
    else                                                                                \
    {                                                                                   \
        quadrature_panel_t<                                                             \
                                PanelGeom,                                              \
                                GWFcnsInterfaceT<NUM_GP2>,                              \
                                wave_term_integral_inf_freq<TA0, TA1, TA2, TA3, TA4>,   \
                                NUM_GP                                                  \
                            >(                                                          \
                                source_i->panel,                                        \
                                this->_gwfcns_interf,                                   \
                                wave_fcn_value,                                         \
                                wave_fcn_dn_sf_value,                                   \
                                wave_fcn_dn_pf_value,                                   \
                                wave_fcn_dx_value,                                      \
                                wave_fcn_dy_value,                                      \
                                wave_fcn_dz_value                                       \
                            );                                                          \
    }                                                                                   \
}                                                                                       \


/*************************************************************************/
/****************** Define auxiliar module functions *********************/
/*************************************************************************/
template<int mode_pf, freq_regime_t freq_regime>
void _formulation_kernel_steady(
                                    bool        is_diag,
                                    PanelGeom*  panel_i,
                                    cusfloat*   field_point,
                                    cusfloat    water_depth,
                                    cusfloat&   pot_term,
                                    cusfloat*   int_dn_pf_value,
                                    cusfloat*   int_dn_sf_value,
                                    cusfloat*   vel_total
                                )
{
    // Define local variables to work with the fast solver
    const int   ndim                    = 3;
    cusfloat    field_point_i[ndim]     = { 0.0, 0.0, 0.0 };
    PanelGeom*  panel_i                 = nullptr;
    PanelGeom*  panel_mirror_i          = nullptr;
    cusfloat    vel_0[ndim]             = { 0.0, 0.0, 0.0 };
    cusfloat    vel_1[ndim]             = { 0.0, 0.0, 0.0 };
    cusfloat    vel_2[ndim]             = { 0.0, 0.0, 0.0 };
    cusfloat    vel_3[ndim]             = { 0.0, 0.0, 0.0 };
    cusfloat    vel_4[ndim]             = { 0.0, 0.0, 0.0 };
    cusfloat    vel_5[ndim]             = { 0.0, 0.0, 0.0 };

    cusfloat    pot_0                   = 0.0;
    cusfloat    pot_1                   = 0.0;
    cusfloat    pot_2                   = 0.0;
    cusfloat    pot_3                   = 0.0;
    cusfloat    pot_4                   = 0.0;
    cusfloat    pot_5                   = 0.0;
    cuscomplex  pot_term                = 0.0;
    
    // Calcualte velocity corresponding to the r0 source
    calculate_source_newman(
                                panel_i,
                                field_point,
                                0,
                                0, 
                                vel_0,
                                pot_0
                            );

    // Calculate velocity corresponding to the r1 source
    field_point_i[0]    =   field_point[0];
    field_point_i[1]    =   field_point[1];
    field_point_i[2]    =   field_point[2] + 2 * water_depth;
    calculate_source_newman(
                                panel_mirror_i,
                                field_point_i, 
                                0,
                                0, 
                                vel_1,
                                pot_1
                            );
    
    // Calculate velocity corresponding to the r2 source
    field_point_i[0]    =   field_point[0];
    field_point_i[1]    =   field_point[1];
    field_point_i[2]    =   field_point[2];
    calculate_source_newman(
                                panel_mirror_i,
                                field_point_i, 
                                0,
                                0, 
                                vel_2,
                                pot_2
                            );

    // Calculate velocity corresponding to the r3 source
    field_point_i[0]    =   field_point[0];
    field_point_i[1]    =   field_point[1];
    field_point_i[2]    =   field_point[2] + 2.0 * water_depth;
    calculate_source_newman(
                                panel_i,
                                field_point_i, 
                                0,
                                0, 
                                vel_3,
                                pot_3
                            );

    // Calculate velocity corresponding to the r4 source
    field_point_i[0]    =   field_point[0];
    field_point_i[1]    =   field_point[1];
    field_point_i[2]    =   -field_point[2] + 2.0 * water_depth;
    calculate_source_newman(
                                panel_mirror_i,
                                field_point_i, 
                                0,
                                0, 
                                vel_4,
                                pot_4
                            );

    // Calculate velocity corresponding to the r5 source
    field_point_i[0]    =   field_point[0];
    field_point_i[1]    =   field_point[1];
    field_point_i[2]    =   field_point[2] + 4.0 * water_depth;
    calculate_source_newman(
                                panel_mirror_i,
                                field_point_i, 
                                0,
                                0, 
                                vel_5,
                                pot_5
                            );

    // Compose total velocity vector
    if constexpr( freq_regime == freq_regime_t::REGULAR || freq_regime == freq_regime_t::ASYMPT_LOW )
    {
        vel_total_sf[0] = - ( vel_0[0] + vel_1[0] + vel_2[0] + vel_3[0] + vel_4[0] + vel_5[0] );
        vel_total_sf[1] = - ( vel_0[1] + vel_1[1] + vel_2[1] + vel_3[1] + vel_4[1] + vel_5[1] );
        vel_total_sf[2] = - ( vel_0[2] + vel_1[2] + vel_2[2] + vel_3[2] - vel_4[2] + vel_5[2] );

        STATIC_COND( ONLY_PF, vel_total_pf[2] = - ( - vel_0[2] + vel_1[2] + vel_2[2] - vel_3[2] + vel_4[2] + vel_5[2] );    )

        pot_term        = ( pot_0 + pot_1 + pot_2 + pot_3 + pot_4 + pot_5 ) / 4.0 / PI;
    }
    else // freq_regime == FREQ_REGIME_INF_FREQ
    {
        vel_total_sf[0] = - ( vel_0[0] + vel_1[0] - vel_2[0] - vel_3[0] - vel_4[0] - vel_5[0] );
        vel_total_sf[1] = - ( vel_0[1] + vel_1[1] - vel_2[1] - vel_3[1] - vel_4[1] - vel_5[1] );
        vel_total_sf[2] = - ( vel_0[2] + vel_1[2] - vel_2[2] - vel_3[2] + vel_4[2] - vel_5[2] );

        STATIC_COND( ONLY_PF, vel_total_pf[2] = - ( - vel_0[2] + vel_1[2] - vel_2[2] + vel_3[2] - vel_4[2] - vel_5[2] );    )

        pot_term        = ( pot_0 + pot_1 - pot_2 - pot_3 - pot_4 - pot_5 ) / 4.0 / PI;
    }

    STATIC_COND( ONLY_PF, vel_total_pf[0] = - vel_total_sf[0];                                                          )
    STATIC_COND( ONLY_PF, vel_total_pf[1] = - vel_total_sf[1];                                                          )
    
                            
    int_dn_sf_value    = (
                            this->_mesh_gp->source_nodes[j]->normal_vec[0] * vel_total_sf[0]
                            +
                            this->_mesh_gp->source_nodes[j]->normal_vec[1] * vel_total_sf[1]
                            +
                            this->_mesh_gp->source_nodes[j]->normal_vec[2] * vel_total_sf[2]
                        ) / 4.0 / PI;
                    
    

    STATIC_COND(
                    ONLY_PF,
                    int_dn_pf_value    = (
                                            this->_mesh_gp->source_nodes[i]->normal_vec[0] * vel_total_pf[0]
                                            +
                                            this->_mesh_gp->source_nodes[i]->normal_vec[1] * vel_total_pf[1]
                                            +
                                            this->_mesh_gp->source_nodes[i]->normal_vec[2] * vel_total_pf[2]
                                        ) / 4.0 / PI;
                )

    // Discard spurious values for normal derivatives
    if ( is_diag )
    {
                                int_dn_sf_value = cuscomplex( 0.0, 0.0 );
        STATIC_COND( ONLY_PF,   int_dn_pf_value = cuscomplex( 0.0, 0.0 ); )
                                vel_total[0]    = cuscomplex( 0.5, 0.0 ) * panel_i->normal_vec[0]; // Used 0.5 but need to be verified because of term: exp( -k * V3 )
                                vel_total[1]    = cuscomplex( 0.5, 0.0 ) * panel_i->normal_vec[1]; // Used 0.5 but need to be verified because of term: exp( -k * V3 )
                                vel_total[2]    = cuscomplex( 0.5, 0.0 ) * panel_i->normal_vec[2]; // Used 0.5 but need to be verified because of term: exp( -k * V3 )
    }
    else
    {
        vel_total[0] = vel_total_sf[0] / 4.0 / PI;
        vel_total[1] = vel_total_sf[1] / 4.0 / PI;
        vel_total[2] = vel_total_sf[2] / 4.0 / PI;
    }

    
}


template<typename T>
void _formulation_kernel_wave( 
                                int         is_diag,
                                SourceNode* source_i,
                                SourceNode* source_j,
                                T*          gwf_interf,
                                cuscomplex& pot_term,
                                cuscomplex& int_dn_pf_value,
                                cuscomplex& int_dn_sf_value,
                                cusfloat*   vel_total
                            )
{
    // Get memory address of the panel jth
    PanelGeom*      panel_j = source_j->panel;
    gwf_interf->set_source_j( source_j );
    
    // Calculate distance in between field point and source
    cusfloat    distn   =  std::sqrt( 
                                        pow2s( source_j->position[0] - source_i->panel->center[0] )
                                        +
                                        pow2s( source_j->position[1] - source_i->panel->center[1] )
                                    );
    cusfloat    dist    =  distn / this->_input->water_depth;
    bool        is_john = dist > 1.0;
    

    // Declare local auxiliar variables
    cusfloat    log_sing_val            = 0.0;
    cuscomplex  wave_fcn_value          = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dn_sf_value    = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dn_pf_value    = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dx_value       = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dy_value       = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dz_value       = cuscomplex( 0.0, 0.0 );

    // Reset to cero variables in other to avoid spurious data
    pot_term           = 0.0;
    int_dn_sf_value     = 0.0;
    int_dn_pf_value     = 0.0;

    // Integrate green function normal derivative along the current panel
    if ( is_diag )
    {
        if ( panel_j->type == DIFFRAC_PANEL_CODE )
        {
            int_dn_sf_value     = cuscomplex( 0.5, 0.0 );
            int_dn_pf_value     = cuscomplex( 0.5, 0.0 );

            _QUADRATURE_PANEL_T( NUM_GP2, G_ON, DGDR_OFF, DGDZ_OFF, FSLID_OFF )
            
            pot_term   =   wave_fcn_value / 4.0 / PI;
        }
        else if ( panel_j->type == LID_PANEL_CODE )
        {
            int_dn_sf_value     = -cuscomplex( 1.0, 0.0 );
            int_dn_pf_value     = -cuscomplex( 1.0, 0.0 );
            
            _QUADRATURE_PANEL_T( NUM_GP2, G_ON, DGDR_OFF, DGDZ_OFF, FSLID_ON )
            
            if constexpr( freq_regime == freq_regime_t::REGULAR )
            {
                log_sing_val        = 2.0 * ( LOG2_GAMMA - std::log( nu ) - source_i->panel->free_surface_log_int ) * source_i->panel->area;
                wave_fcn_value      += cuscomplex( nu * log_sing_val, 0.0 );
            }
            pot_term           = wave_fcn_value / 4.0 / PI;
            
        }
        
    }
    else
    {
        wave_fcn_value          = 0.0;
        wave_fcn_dn_sf_value    = 0.0;
        wave_fcn_dn_pf_value    = 0.0;

        if ( is_john && freq_regime == freq_regime_t::REGULAR )
        {
            quadrature_panel_t<
                                    PanelGeom, 
                                    GWFcnsInterfaceT<NUM_GP2>, 
                                    john_series<NUM_GP2, G_ON, DGDR_ON, DGDZ_ON>, 
                                    NUM_GP
                                >( 
                                    source_i->panel, 
                                    this->_gwfcns_interf, 
                                    wave_fcn_value,
                                    wave_fcn_dn_sf_value,
                                    wave_fcn_dn_pf_value,
                                    wave_fcn_dx_value,
                                    wave_fcn_dy_value,
                                    wave_fcn_dz_value
                                );
        }
        else
        {
            _QUADRATURE_PANEL_T( NUM_GP2, G_ON, DGDR_ON, DGDZ_ON, FSLID_OFF )
        }

        pot_term        =   wave_fcn_value / 4.0 / PI;
        int_dn_sf_value =   wave_fcn_dn_sf_value / 4.0 / PI;
        int_dn_pf_value =   wave_fcn_dn_pf_value / 4.0 / PI;
        vel_total[0]    =   wave_fcn_dx_value / 4.0 / PI;
        vel_total[1]    =   wave_fcn_dy_value / 4.0 / PI;
        vel_total[2]    =   wave_fcn_dz_value / 4.0 / PI;

    }
}


template<std::size_t N, int mode_pf>
template<freq_regime_t freq_regime>
void FormulationKernelBackend<N, mode_pf>::_build_steady_matrixes( void )
{
    /***************************************/
    /******** Fill system matrix  **********/
    /***************************************/
    // Loop over panels to integrate value
    int         col_count       = 0;
    cuscomplex  int_value_sf    = cuscomplex( 0.0, 0.0 );
    cuscomplex  int_value_pf    = cuscomplex( 0.0, 0.0 );
    int         row_count       = 0;
    int         index_cm        = 0;
    int         index_rm        = 0;

    // Define local variables to work with the fast solver
    const int   ndim                    = 3;
    cusfloat    field_point_i[ndim];    clear_vector( ndim, field_point_i );
    PanelGeom*  panel_i                 = nullptr;
    PanelGeom*  panel_mirror_i          = nullptr;
    cusfloat    vel_0[ndim];            clear_vector( ndim, vel_0 );
    cusfloat    vel_1[ndim];            clear_vector( ndim, vel_1 );
    cusfloat    vel_2[ndim];            clear_vector( ndim, vel_2 );
    cusfloat    vel_3[ndim];            clear_vector( ndim, vel_3 );
    cusfloat    vel_4[ndim];            clear_vector( ndim, vel_4 );
    cusfloat    vel_5[ndim];            clear_vector( ndim, vel_5 );
    cusfloat    vel_total_pf[ndim];     clear_vector( ndim, vel_total_pf );
    cusfloat    vel_total_sf[ndim];     clear_vector( ndim, vel_total_sf );

    cusfloat    pot_0                   = 0.0;
    cusfloat    pot_1                   = 0.0;
    cusfloat    pot_2                   = 0.0;
    cusfloat    pot_3                   = 0.0;
    cusfloat    pot_4                   = 0.0;
    cusfloat    pot_5                   = 0.0;
    cuscomplex  pot_term                = 0.0;

    // Define field points to calculate the source influence matrix
    int         field_points_np   = this->_mesh_gp->panels_tnp;
    cusfloat*   field_points      = generate_empty_vector<cusfloat>( 3 * field_points_np );

    for ( int i=0; i<this->_mesh_gp->panels_tnp; i++ )
    {
        copy_vector( 3, this->_mesh_gp->panels[i]->center, &(field_points[3*i]) );
    }

    // Loop over panels and field points to create the steady source matrix
    for ( int i=this->_solver->start_col_0; i<this->_solver->end_col_0; i++ )
    {
        // Get pointer to ith panel
        panel_i         = this->_mesh_gp->panels[i];
        panel_mirror_i  = this->_mesh_gp->panels_mirror[i];

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=0; j<field_points_np; j++ )
        {
            // Reset velocity values
            clear_vector( ndim, vel_0 );
            clear_vector( ndim, vel_1 );
            clear_vector( ndim, vel_2 );
            clear_vector( ndim, vel_3 );
            clear_vector( ndim, vel_4 );
            clear_vector( ndim, vel_5 );
            clear_vector( ndim, vel_total_pf );
            clear_vector( ndim, vel_total_sf );

            // Reset potential values
            pot_0       = 0.0;
            pot_1       = 0.0;
            pot_2       = 0.0;
            pot_3       = 0.0;
            pot_4       = 0.0;
            pot_5       = 0.0;
            pot_term    = 0.0;
            
            // Calcualte velocity corresponding to the r0 source
            calculate_source_newman(
                                        panel_i,
                                        &(field_points[3*j]), 
                                        0,
                                        0, 
                                        vel_0,
                                        pot_0
                                    );

            // Calculate velocity corresponding to the r1 source
            field_point_i[0]    =   field_points[3*j];
            field_point_i[1]    =   field_points[3*j+1];
            field_point_i[2]    =   field_points[3*j+2] + 2 * this->_input->water_depth;
            calculate_source_newman(
                                        panel_mirror_i,
                                        field_point_i, 
                                        0,
                                        0, 
                                        vel_1,
                                        pot_1
                                    );
            
            // Calculate velocity corresponding to the r2 source
            field_point_i[0]    =   field_points[3*j];
            field_point_i[1]    =   field_points[3*j+1];
            field_point_i[2]    =   field_points[3*j+2];
            calculate_source_newman(
                                        panel_mirror_i,
                                        field_point_i, 
                                        0,
                                        0, 
                                        vel_2,
                                        pot_2
                                    );

            // Calculate velocity corresponding to the r3 source
            field_point_i[0]    =   field_points[3*j];
            field_point_i[1]    =   field_points[3*j+1];
            field_point_i[2]    =   field_points[3*j+2] + 2.0 * this->_input->water_depth;
            calculate_source_newman(
                                        panel_i,
                                        field_point_i, 
                                        0,
                                        0, 
                                        vel_3,
                                        pot_3
                                    );

            // Calculate velocity corresponding to the r4 source
            field_point_i[0]    =   field_points[3*j];
            field_point_i[1]    =   field_points[3*j+1];
            field_point_i[2]    =   -field_points[3*j+2] + 2.0 * this->_input->water_depth;
            calculate_source_newman(
                                        panel_mirror_i,
                                        field_point_i, 
                                        0,
                                        0, 
                                        vel_4,
                                        pot_4
                                    );

            // Calculate velocity corresponding to the r5 source
            field_point_i[0]    =   field_points[3*j];
            field_point_i[1]    =   field_points[3*j+1];
            field_point_i[2]    =   field_points[3*j+2] + 4.0 * this->_input->water_depth;
            calculate_source_newman(
                                        panel_mirror_i,
                                        field_point_i, 
                                        0,
                                        0, 
                                        vel_5,
                                        pot_5
                                    );

            // Compose total velocity vector
            if constexpr( freq_regime == freq_regime_t::REGULAR || freq_regime == freq_regime_t::ASYMPT_LOW )
            {
                vel_total_sf[0] = - ( vel_0[0] + vel_1[0] + vel_2[0] + vel_3[0] + vel_4[0] + vel_5[0] );
                vel_total_sf[1] = - ( vel_0[1] + vel_1[1] + vel_2[1] + vel_3[1] + vel_4[1] + vel_5[1] );
                vel_total_sf[2] = - ( vel_0[2] + vel_1[2] + vel_2[2] + vel_3[2] - vel_4[2] + vel_5[2] );

                STATIC_COND( ONLY_PF, vel_total_pf[2] = - ( - vel_0[2] + vel_1[2] + vel_2[2] - vel_3[2] + vel_4[2] + vel_5[2] );    )

                pot_term        = ( pot_0 + pot_1 + pot_2 + pot_3 + pot_4 + pot_5 ) / 4.0 / PI;
            }
            else // freq_regime == FREQ_REGIME_INF_FREQ
            {
                vel_total_sf[0] = - ( vel_0[0] + vel_1[0] - vel_2[0] - vel_3[0] - vel_4[0] - vel_5[0] );
                vel_total_sf[1] = - ( vel_0[1] + vel_1[1] - vel_2[1] - vel_3[1] - vel_4[1] - vel_5[1] );
                vel_total_sf[2] = - ( vel_0[2] + vel_1[2] - vel_2[2] - vel_3[2] + vel_4[2] - vel_5[2] );

                STATIC_COND( ONLY_PF, vel_total_pf[2] = - ( - vel_0[2] + vel_1[2] - vel_2[2] + vel_3[2] - vel_4[2] - vel_5[2] );    )

                pot_term        = ( pot_0 + pot_1 - pot_2 - pot_3 - pot_4 - pot_5 ) / 4.0 / PI;
            }

            STATIC_COND( ONLY_PF, vel_total_pf[0] = - vel_total_sf[0];                                                          )
            STATIC_COND( ONLY_PF, vel_total_pf[1] = - vel_total_sf[1];                                                          )
            
                                    
            int_value_sf    = (
                                    this->_mesh_gp->source_nodes[j]->normal_vec[0] * vel_total_sf[0]
                                    +
                                    this->_mesh_gp->source_nodes[j]->normal_vec[1] * vel_total_sf[1]
                                    +
                                    this->_mesh_gp->source_nodes[j]->normal_vec[2] * vel_total_sf[2]
                                ) / 4.0 / PI;
                            
            

            STATIC_COND(
                            ONLY_PF,
                            int_value_pf    = (
                                                    this->_mesh_gp->source_nodes[i]->normal_vec[0] * vel_total_pf[0]
                                                    +
                                                    this->_mesh_gp->source_nodes[i]->normal_vec[1] * vel_total_pf[1]
                                                    +
                                                    this->_mesh_gp->source_nodes[i]->normal_vec[2] * vel_total_pf[2]
                                                ) / 4.0 / PI;
                        )

            // Discard spurious values for normal derivatives
            if ( i == j )
            {
                                        int_value_sf = cuscomplex( 0.0, 0.0 );
                STATIC_COND( ONLY_PF,   int_value_pf = cuscomplex( 0.0, 0.0 ); )
            }

            COL_MAJOR_INDEX( index_cm, row_count, col_count, this->_solver->num_rows_local )
            ROW_MAJOR_INDEX( index_rm, row_count, col_count, this->_solver->num_cols_local )

                                    this->_sf_gp->sysmat_steady[index_cm]  = int_value_sf;
                                    this->_pot_gp->sysmat_steady[index_rm] = pot_term;
            STATIC_COND( ONLY_PF,   this->_pf_gp->sysmat_steady[index_cm]  = int_value_pf; )

            // Advance row count
            row_count++;
        }

        // Advance column count
        col_count++;
    }

    // Delete heap memory associated to this block of code
    mkl_free( field_points );

    // Synchronize processes progress status
    MPI_Barrier( MPI_COMM_WORLD );
}


template<std::size_t N, int mode_pf>
void FormulationKernelBackend<N, mode_pf>::_build_rhs( 
                                                        cusfloat w
                                                    )
{
    // Declare local variables
    int         col_count   = 0;
    int         dofs_offset = this->_input->dofs_np * this->_pot_gp->sysmat_nrows;
    int         index       = 0;
    int         index_rm    = 0;
    PanelGeom*  panel_j     = nullptr;
    int         row_count   = 0;
    SourceNode* source_i    = nullptr;
    cuscomplex  wave_dx     = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_dy     = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_dz     = cuscomplex( 0.0, 0.0 );

    // Calculate wave dependent parameters
    cusfloat    k           = w2k( w, this->_input->water_depth, this->_input->grav_acc );

    // Clear potential rhs to avoid spurious valures
    STATIC_COND( ONLY_PF, clear_vector( this->_pf_gp->sysmat_nrows * ( this->_input->dofs_np + this->_input->heads_np ), this->_ppf_rhs ); )

    // Calculate potential formulation rhs
    if constexpr( ONLY_PF )
    {
        for ( int i=this->_solver->start_col_0; i<this->_solver->end_col_0; i++ )
        {
            source_i    = this->_mesh_gp->source_nodes[i];
            row_count   = 0;
            for ( int j=this->_solver->start_row_0; j<this->_solver->end_row_0; j++ )
            {
                ROW_MAJOR_INDEX( index_rm, row_count, col_count, this->_solver->num_cols_local )

                if ( source_i->panel->type == DIFFRAC_PANEL_CODE )
                {
                    for ( int id=0; id<this->_input->dofs_np; id++ )
                    {
                        index                   = id * this->_pf_gp->sysmat_nrows + j; 
                        this->_ppf_rhs[index]   += (
                                                        this->_mesh_gp->source_nodes[i]->normal_vec[id]
                                                        *
                                                        this->_mesh_gp->source_nodes[i]->panel->is_move_f
                                                    ) * this->_pot_gp->sysmat[index_rm];
                    }
    
                    for ( int id=0; id<this->_input->heads_np; id++ )
                    {
                        // Get wave potential derivatives for the panel
                        wave_dx     =   wave_potential_fo_space_dx(
                                                                        1.0,
                                                                        w,
                                                                        k,
                                                                        this->_input->water_depth,
                                                                        this->_input->grav_acc,
                                                                        this->_mesh_gp->source_nodes[i]->panel->center[0],
                                                                        this->_mesh_gp->source_nodes[i]->panel->center[1],
                                                                        this->_mesh_gp->source_nodes[i]->panel->center[2],
                                                                        this->_input->heads[id]
                                                                    );
    
                        wave_dy     =   wave_potential_fo_space_dy(
                                                                        1.0,
                                                                        w,
                                                                        k,
                                                                        this->_input->water_depth,
                                                                        this->_input->grav_acc,
                                                                        this->_mesh_gp->source_nodes[i]->panel->center[0],
                                                                        this->_mesh_gp->source_nodes[i]->panel->center[1],
                                                                        this->_mesh_gp->source_nodes[i]->panel->center[2],
                                                                        this->_input->heads[id]
                                                                    );
    
                        wave_dz     =   wave_potential_fo_space_dz(
                                                                        1.0,
                                                                        w,
                                                                        k,
                                                                        this->_input->water_depth,
                                                                        this->_input->grav_acc,
                                                                        this->_mesh_gp->source_nodes[i]->panel->center[0],
                                                                        this->_mesh_gp->source_nodes[i]->panel->center[1],
                                                                        this->_mesh_gp->source_nodes[i]->panel->center[2],
                                                                        this->_input->heads[id]
                                                                    );
                        
                        // Calculate normal derivative of the wave flow velocities for the jth panel
                        index                   = dofs_offset + id * this->_pf_gp->sysmat_nrows + j; 
                        this->_ppf_rhs[index]   += -(
                                                        wave_dx * this->_mesh_gp->source_nodes[i]->normal_vec[0]
                                                        +
                                                        wave_dy * this->_mesh_gp->source_nodes[i]->normal_vec[1]
                                                        +
                                                        wave_dz * this->_mesh_gp->source_nodes[i]->normal_vec[2]
                                                    ) * this->_pot_gp->sysmat[index_rm];
                    }
                }
    
                // Advance row count
                row_count++;

            }
            // Advance column count
            col_count++;
    
        }

        // Sum up contributions along processors
        MPI_Allreduce(
                            this->_ppf_rhs,
                            this->_pf_gp->field_values,
                            this->_pf_gp->sysmat_nrows * ( this->_input->dofs_np + this->_input->heads_np ),
                            mpi_cuscomplex,
                            MPI_SUM,
                            MPI_COMM_WORLD
                        );
    }


    // Calculate source formulation rhs
    for ( int i=0; i<this->_input->dofs_np; i++ )
    {
        for ( int j=this->_solver->start_row_0; j<this->_solver->end_row_0; j++ )
        {
            panel_j = this->_mesh_gp->source_nodes[j]->panel;
            index   = i * this->_sf_gp->sysmat_nrows + j;
            if ( panel_j->type == DIFFRAC_PANEL_CODE )
            {
                this->_sf_gp->field_values[index] = ( 
                                                        this->_mesh_gp->source_nodes[j]->normal_vec[i]
                                                        *
                                                        this->_mesh_gp->source_nodes[j]->panel->is_move_f
                                                    );

            }
            else if ( panel_j->type == LID_PANEL_CODE )
            {
                this->_sf_gp->field_values[index] = 0.0;
            }

        }
    }

    for ( int i=0; i<this->_input->heads_np; i++ )
    {
        for ( int j=this->_solver->start_row_0; j<this->_solver->end_row_0; j++ )
        {
            panel_j = this->_mesh_gp->source_nodes[j]->panel;
            index   = dofs_offset + i * this->_sf_gp->sysmat_nrows + j; 
            if ( panel_j->type == DIFFRAC_PANEL_CODE )
            {
                // Get wave potential derivatives for the panel
                wave_dx     =   wave_potential_fo_space_dx(
                                                                1.0,
                                                                w,
                                                                k,
                                                                this->_input->water_depth,
                                                                this->_input->grav_acc,
                                                                this->_mesh_gp->source_nodes[j]->panel->center[0],
                                                                this->_mesh_gp->source_nodes[j]->panel->center[1],
                                                                this->_mesh_gp->source_nodes[j]->panel->center[2],
                                                                this->_input->heads[i]
                                                            );

                wave_dy     =   wave_potential_fo_space_dy(
                                                                1.0,
                                                                w,
                                                                k,
                                                                this->_input->water_depth,
                                                                this->_input->grav_acc,
                                                                this->_mesh_gp->source_nodes[j]->panel->center[0],
                                                                this->_mesh_gp->source_nodes[j]->panel->center[1],
                                                                this->_mesh_gp->source_nodes[j]->panel->center[2],
                                                                this->_input->heads[i]
                                                            );

                wave_dz     =   wave_potential_fo_space_dz(
                                                                1.0,
                                                                w,
                                                                k,
                                                                this->_input->water_depth,
                                                                this->_input->grav_acc,
                                                                this->_mesh_gp->source_nodes[j]->panel->center[0],
                                                                this->_mesh_gp->source_nodes[j]->panel->center[1],
                                                                this->_mesh_gp->source_nodes[j]->panel->center[2],
                                                                this->_input->heads[i]
                                                            );
                
                // Calculate normal derivative of the wave flow velocities for the jth panel
                this->_sf_gp->field_values[index]   = -(
                                                            wave_dx * this->_mesh_gp->source_nodes[j]->normal_vec[0]
                                                            +
                                                            wave_dy * this->_mesh_gp->source_nodes[j]->normal_vec[1]
                                                            +
                                                            wave_dz * this->_mesh_gp->source_nodes[j]->normal_vec[2]
                                                        );
            }
            else if ( panel_j->type == LID_PANEL_CODE )
            {
                this->_sf_gp->field_values[index]  = 0.0;
            }
            
        }
    }

    // Synchronize processes progress status
    MPI_Barrier( MPI_COMM_WORLD );

}



template<std::size_t N, int mode_pf>
template<freq_regime_t freq_regime>
void FormulationKernelBackend<N, mode_pf>::_build_wave_matrixes( 
                                                                    cusfloat w
                                                                )
{
    // Clean system matrixes
                            this->_sf_gp->clear_sysmat( );
    STATIC_COND( ONLY_PF,   this->_pf_gp->clear_sysmat( ); )
                            this->_pot_gp->clear_sysmat( );
                            this->_pot_gp->clear_field_values( );
    
    // Declare local variables
    int         col_count               = 0;
    cusfloat    dist                    = 0.0;
    cusfloat    distn                   = 0.0;
    int         index_cm                = 0;
    int         index_rm                = 0;
    bool        is_john                 = false;
    cuscomplex  int_value               = cuscomplex( 0.0, 0.0 );
    cuscomplex  int_dn_sf_value         = cuscomplex( 0.0, 0.0 );
    cuscomplex  int_dn_pf_value         = cuscomplex( 0.0, 0.0 );
    cusfloat    log_sing_val            = 0.0;
    PanelGeom*  panel_j                 = nullptr;
    int         row_count               = 0;
    SourceNode* source_i                = nullptr;
    cuscomplex  wave_fcn_value          = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dn_sf_value    = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dn_pf_value    = cuscomplex( 0.0, 0.0 );
    
    // Calculate wave dependent parameters
    cusfloat    nu                      = pow2s( w ) / this->_input->grav_acc;
    
    for ( int i=this->_solver->start_col_0; i<this->_solver->end_col_0; i++ )
    {
        // Get memory address of the ith panel
        source_i = this->_mesh_gp->source_nodes[i];
        this->_gwfcns_interf.set_source_i( source_i, 1.0 );

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=this->_solver->start_row_0; j<this->_solver->end_row_0; j++ )
        {
            // Get memory address of the panel jth
            panel_j = this->_mesh_gp->source_nodes[j]->panel;
            this->_gwfcns_interf.set_source_j( this->_mesh_gp->source_nodes[j] );
            
            // Calculate distance in between field point and source
            distn   =  std::sqrt( 
                                    pow2s( this->_mesh_gp->source_nodes[j]->position[0] - source_i->panel->center[0] )
                                    +
                                    pow2s( this->_mesh_gp->source_nodes[j]->position[1] - source_i->panel->center[1] )
                                );
            dist    =  distn / this->_input->water_depth;
            is_john = dist > 1.0;
            
            int_value       = 0.0;
            int_dn_sf_value = 0.0;
            int_dn_pf_value = 0.0;
            // Integrate green function normal derivative along the current panel
            if ( i == j )
            {
                if ( panel_j->type == DIFFRAC_PANEL_CODE )
                {
                    int_dn_sf_value     = cuscomplex( 0.5, 0.0 );
                    int_dn_pf_value     = cuscomplex( 0.5, 0.0 );

                    _QUADRATURE_PANEL_T( NUM_GP2, G_ON, DGDR_OFF, DGDZ_OFF, FSLID_OFF )
                    
                    int_value   =   wave_fcn_value / 4.0 / PI;
                }
                else if ( panel_j->type == LID_PANEL_CODE )
                {
                    int_dn_sf_value     = -cuscomplex( 1.0, 0.0 );
                    int_dn_pf_value     = -cuscomplex( 1.0, 0.0 );
                    
                    _QUADRATURE_PANEL_T( NUM_GP2, G_ON, DGDR_OFF, DGDZ_OFF, FSLID_ON )
                    
                    if constexpr( freq_regime == freq_regime_t::REGULAR )
                    {
                        log_sing_val        = 2.0 * ( LOG2_GAMMA - std::log( nu ) - source_i->panel->free_surface_log_int ) * source_i->panel->area;
                        wave_fcn_value      += cuscomplex( nu * log_sing_val, 0.0 );
                    }
                    int_value           = wave_fcn_value / 4.0 / PI;
                    
                }
                
            }
            else
            {
                wave_fcn_value          = 0.0;
                wave_fcn_dn_sf_value    = 0.0;
                wave_fcn_dn_pf_value    = 0.0;

                if ( is_john && freq_regime == freq_regime_t::REGULAR )
                {
                    quadrature_panel_t<
                                            PanelGeom, 
                                            GWFcnsInterfaceT<NUM_GP2>, 
                                            john_series<NUM_GP2, G_ON, DGDR_ON, DGDZ_ON>, 
                                            NUM_GP
                                        >( 
                                            source_i->panel, 
                                            this->_gwfcns_interf, 
                                            wave_fcn_value,
                                            wave_fcn_dn_sf_value,
                                            wave_fcn_dn_pf_value
                                        );
                }
                else
                {
                    _QUADRATURE_PANEL_T( NUM_GP2, G_ON, DGDR_ON, DGDZ_ON, FSLID_OFF )
                }

                int_value       =   wave_fcn_value / 4.0 / PI;
                int_dn_sf_value =   wave_fcn_dn_sf_value / 4.0 / PI;
                int_dn_pf_value =   wave_fcn_dn_pf_value / 4.0 / PI;

            }

            // Apply the integral value accordingly
            COL_MAJOR_INDEX( index_cm, row_count, col_count, this->_solver->num_rows_local )
            ROW_MAJOR_INDEX( index_rm, row_count, col_count, this->_solver->num_cols_local )

            if ( is_john && freq_regime == freq_regime_t::REGULAR )
            {
                this->_pot_gp->sysmat[index_rm] = int_value;
                this->_sf_gp->sysmat[index_cm]  = int_dn_sf_value;

                STATIC_COND( ONLY_PF, this->_pf_gp->sysmat[index_cm]  = int_dn_pf_value; )
            }
            else
            {
                this->_pot_gp->sysmat[index_rm] = this->_pot_gp->sysmat_steady[index_rm] + int_value;
                this->_sf_gp->sysmat[index_cm]  = this->_sf_gp->sysmat_steady[index_cm]  + int_dn_sf_value;

                STATIC_COND( ONLY_PF, this->_pf_gp->sysmat[index_cm]  = this->_pf_gp->sysmat_steady[index_cm] + int_dn_pf_value; )
            }

            // Advance row count
            row_count++;
        }
        
        // Advance column count
        col_count++;

    }
    MPI_Barrier( MPI_COMM_WORLD );

}


template<std::size_t N, int mode_pf>
void FormulationKernelBackend<N, mode_pf>::_initialize( 
                                                            void 
                                                        )
{
    // Calculate steady part integral over the panels
    MPI_TIME_EXEC( this->_build_steady_matrixes<freq_regime_t::REGULAR>( ); , this->exec_time_build_steady )
    this->_steady_mat_type  = 0;
}


template<std::size_t N, int mode_pf>
FormulationKernelBackend<N, mode_pf>::FormulationKernelBackend(
                                                                    Input*      input, 
                                                                    MpiConfig*  mpi_config, 
                                                                    MeshGroup*  mesh_gp
                                                                )
{
    // Storage input arguments
    this->_input        = input;
    this->_mesh_gp      = mesh_gp;
    this->_mpi_config   = mpi_config;

    // Instantiate Scalapack solver
    this->_solver       = new   SclCmpx( 
                                            this->_mesh_gp->source_nodes_tnp,
                                            this->_input->dofs_np + this->_input->heads_np,
                                            this->_mpi_config->procs_total,
                                            this->_mpi_config->proc_rank,
                                            this->_mpi_config->proc_root,
                                            MPI_COMM_WORLD
                                        );

    // Set data distribution variables
    this->ipm_sc        = this->_solver->start_col - 1;
    this->ipm_ed        = this->ipm_sc + this->_solver->num_cols_local - 1;
    this->ipm_cols_np   = this->_solver->num_cols_local;
    
    // Allocate space for the system matrixes
    this->_sf_gp        = new   MLGCmpx(
                                            this->_mesh_gp->panels_tnp,
                                            this->ipm_cols_np,
                                            this->_mesh_gp->meshes_np,
                                            ( this->_input->dofs_np + this->_input->heads_np ),
                                            0,
                                            this->_mesh_gp->panels_tnp-1,
                                            this->ipm_sc,
                                            this->ipm_ed,
                                            true
                                        );

    this->_pot_gp       = new   MLGCmpx(
                                            this->_mesh_gp->panels_tnp,
                                            this->ipm_cols_np,
                                            this->_mesh_gp->meshes_np,
                                            ( this->_input->dofs_np + this->_input->heads_np ),
                                            0,
                                            this->_mesh_gp->panels_tnp-1,
                                            this->ipm_sc,
                                            this->ipm_ed,
                                            true
                                        );

    STATIC_COND( 
                    ONLY_PF,
                    this->_pf_gp        = new   MLGCmpx(
                                                            this->_mesh_gp->panels_tnp,
                                                            this->ipm_cols_np,
                                                            this->_mesh_gp->meshes_np,
                                                            ( this->_input->dofs_np + this->_input->heads_np ),
                                                            0,
                                                            this->_mesh_gp->panels_tnp-1,
                                                            this->ipm_sc,
                                                            this->ipm_ed,
                                                            true
                                                        );
                )

    // Allocate space for the partial potential formulation RHS vector
    STATIC_COND( ONLY_PF, this->_ppf_rhs  = generate_empty_vector<cuscomplex>( this->_pot_gp->sysmat_nrows * ( this->_input->dofs_np + this->_input->heads_np ) ); )

    // Allocate space for the wave part integration interface functor
    this->_gwfcns_interf.initialize(
                                        this->_input->angfreqs[0],
                                        this->_input->water_depth,
                                        this->_input->grav_acc
                                    );
    
    // Launch object initialization
    this->_initialize( );
}


template<std::size_t N, int mode_pf>
FormulationKernelBackend<N, mode_pf>::~FormulationKernelBackend( )
{
    delete this->_solver;
    delete this->_sf_gp;
    delete this->_pot_gp;

    STATIC_COND( ONLY_PF, delete this->_pf_gp; );
}


template<std::size_t N, int mode_pf>
int FormulationKernelBackend<N, mode_pf>::size( void )
{
    return this->_solver->num_rows;
}


template<std::size_t N, int mode_pf>
template<freq_regime_t freq_regime>
void FormulationKernelBackend<N, mode_pf>::solve( cusfloat w )
{
    // Recalculate steady part of the system matrixes if required
    if ( 
            this->_steady_mat_type != 0 
            &&
            (   
                freq_regime == freq_regime_t::REGULAR 
                || 
                freq_regime == freq_regime_t::ASYMPT_LOW 
            )
        )
    {
        MPI_TIME_EXEC( this->_build_steady_matrixes<freq_regime>( ); , this->exec_time_build_steady )
        this->_steady_mat_type  = 0;
    }
    else if ( 
                this->_steady_mat_type != 1 
                &&
                freq_regime == freq_regime_t::ASYMPT_HIGH 
            )
    {
        MPI_TIME_EXEC( this->_build_steady_matrixes<freq_regime>( ); , this->exec_time_build_steady )
        this->_steady_mat_type  = 1;
    }
    
    // Fold integrals database if required
    if constexpr( freq_regime == freq_regime_t::REGULAR )
    {
        // Fold database for current frequency and water depth
        cusfloat H = pow2s( w ) * this->_input->water_depth / this->_input->grav_acc;

        fold_database( H );

        // Update integration functor to ne new frequency
        this->_gwfcns_interf.set_ang_freq( w );
    }
    
    // Re-calculate wave dependent system matrix term
    MPI_TIME_EXEC(  this->_build_wave_matrixes<freq_regime>( w );, this->exec_time_build_wave )
                    this->_build_rhs( w );

    // Calculate system matrixes condition number if required
    if ( this->_is_condition_number )
    {
        cusfloat potm_cond  = 0.0;
        cusfloat pf_cond    = 0.0;
        cusfloat sf_cond    = 0.0;
    
                                this->_solver->Cond( this->_pot_gp->sysmat,   potm_cond   );
        STATIC_COND( ONLY_PF,   this->_solver->Cond( this->_pf_gp->sysmat,    pf_cond     ); )
                                this->_solver->Cond( this->_sf_gp->sysmat,    sf_cond     );
    
        if ( this->_mpi_config->is_root( ) )
        {
                                    std::cout << "Condition Number -> SF: " << sf_cond;
            STATIC_COND( ONLY_PF,   std::cout << " - PF: " << pf_cond; )
                                    std::cout << " - POT: " << potm_cond << std::endl;
        }
    }

    // Calculate potential through the potential formulation
    STATIC_COND(    
                    ONLY_PF, 
                    MPI_TIME_EXEC( 
                                    ( this->_solver->Solve( this->_pf_gp->sysmat, this->_pf_gp->field_values ) );, 
                                    this->exec_time_solve_pf 
                                ) 
                )
    STATIC_COND( 
                    ONLY_PF, 
                    MPI_Bcast(
                                    this->_pf_gp->field_values,
                                    this->_pf_gp->field_values_np,
                                    mpi_cuscomplex,
                                    this->_mpi_config->proc_root,
                                    MPI_COMM_WORLD                
                                );
                )

    // Corrrect source formulation FS jump in the rhs by using
    // potential formulation if possible
    if constexpr( ONLY_PF )
    {
        // Declare local block variables
        int         dofs_offset = this->_input->dofs_np * this->_pot_gp->sysmat_nrows;
        int         index       = 0;
        PanelGeom*  panel_j     = nullptr;
        double      nu          = pow2s( w ) / this->_input->grav_acc;

        for ( int i=0; i<this->_input->dofs_np; i++ )
        {
            for ( int j=this->_solver->start_row_0; j<this->_solver->end_row_0; j++ )
            {
                panel_j = this->_mesh_gp->source_nodes[j]->panel;
                if ( panel_j->type == LID_PANEL_CODE )
                {
                    index                               = i * this->_sf_gp->sysmat_nrows + j;
                    this->_sf_gp->field_values[index]   = -nu * this->_pf_gp->field_values[index];
                }

            }
        }

        for ( int i=0; i<this->_input->heads_np; i++ )
        {
            for ( int j=this->_solver->start_row_0; j<this->_solver->end_row_0; j++ )
            {
                panel_j = this->_mesh_gp->source_nodes[j]->panel;

                if ( panel_j->type == LID_PANEL_CODE )
                {
                    index                               = dofs_offset + i * this->_sf_gp->sysmat_nrows + j; 
                    this->_sf_gp->field_values[index]   = -nu * this->_pf_gp->field_values[index];
                }
                
            }
        }
    }

    // Calculate sources intensity through the source formulation
    MPI_TIME_EXEC( ( this->_solver->Solve( this->_sf_gp->sysmat, this->_sf_gp->field_values ) );, this->exec_time_solve_sf )

    MPI_Bcast(
                this->_sf_gp->field_values,
                this->_sf_gp->field_values_np,
                mpi_cuscomplex,
                this->_mpi_config->proc_root,
                MPI_COMM_WORLD                
            );

    // Calculate potential values through the source formulation if 
    // potential formulation was not enabled
    if constexpr( !(ONLY_PF) )
    {
        // Calculate potentials for each chunk of the matrix along the processors
        calculate_fields_raddif_lin(
                                        this->_input,
                                        this->_sf_gp->field_values,
                                        this->_pot_gp
                                    );

        // Sum contributions from each processor and distribute accordingly
        MPI_Allreduce(
                        MPI_IN_PLACE,
                        this->_pot_gp->field_values,
                        this->_pot_gp->field_values_np,
                        mpi_cuscomplex,
                        MPI_SUM,
                        MPI_COMM_WORLD
                    );
        
    }

}


template<std::size_t N, int mode_pf>
void FormulationKernelBackend<N, mode_pf>::update_results( SimulationData* sim_data )
{
    // Update data to simulation results
    copy_vector( 
                    this->_sf_gp->field_values_np,
                    this->_sf_gp->field_values,
                    sim_data->intensities
                );

    STATIC_COND( 
                    ONLY_PF, 
                    copy_vector( 
                                    this->_pf_gp->field_values_np, 
                                    this->_pf_gp->field_values,
                                    sim_data->panels_potential
                                );
                )
    
    STATIC_COND( 
                    !(ONLY_PF), 
                    copy_vector( 
                                    this->_pot_gp->field_values_np, 
                                    this->_pot_gp->field_values,
                                    sim_data->panels_potential
                                );
                )
}