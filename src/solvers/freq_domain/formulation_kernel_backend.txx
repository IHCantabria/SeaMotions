
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
    if constexpr( freq_regime == FreqRegimeE::REGULAR )                               \
    {                                                                                   \
        quadrature_panel_t<                                                             \
                                PanelGeom,                                              \
                                GWFcnsInterfaceT<NUM_GP2>,                              \
                                wave_term_integral<TA0, TA1, TA2, TA3, TA4>,            \
                                NUM_GP                                                  \
                            >(                                                          \
                                source_i->panel,                                        \
                                gwf_interf,                                             \
                                wave_fcn_value,                                         \
                                wave_fcn_dn_sf_value,                                   \
                                wave_fcn_dn_pf_value,                                   \
                                wave_fcn_dx_value,                                      \
                                wave_fcn_dy_value,                                      \
                                wave_fcn_dz_value                                       \
                            );                                                          \
    }                                                                                   \
    else if constexpr( freq_regime == FreqRegimeE::ASYMPT_LOW )                       \
    {                                                                                   \
        quadrature_panel_t<                                                             \
                                PanelGeom,                                              \
                                GWFcnsInterfaceT<NUM_GP2>,                              \
                                wave_term_integral_zero_freq<TA0, TA1, TA2, TA3, TA4>,  \
                                NUM_GP                                                  \
                            >(                                                          \
                                source_i->panel,                                        \
                                gwf_interf,                                             \
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
                                gwf_interf,                                             \
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
inline cuscomplex calculate_fs_singularity( 
                                                cusfloat    nu,
                                                cusfloat    area,
                                                cusfloat    free_surface_log_int
                                            )
{
    cusfloat    log_sing_val   = 2.0 * ( LOG2_GAMMA - std::log( nu ) - free_surface_log_int ) * area;
    cuscomplex  wave_fcn_value = cuscomplex( nu * log_sing_val, 0.0 );
    return wave_fcn_value;
}

static inline cuscomplex panel_scale_factor( const PanelGeom* panel, cusfloat nu )
{
    switch ( panel->type )
    {
        case PanelTypeE::DIFFRAC:
            return cuscomplex( 1.0, 0.0 );
        case PanelTypeE::INT_LID:
            return cuscomplex( -nu, 0.0 );
        case PanelTypeE::EXT_LID:
            return cuscomplex( 0.0, +nu * panel->ext_lid_damp_f );
        default:
            return cuscomplex( 1.0, 0.0 );
    }
}


template<int mode_pf, FreqRegimeE freq_regime>
void _formulation_kernel_steady(
                                    bool        is_diag,
                                    PanelGeom*  panel_i,
                                    PanelGeom*  panel_mirror_i,
                                    PanelGeom*  panel_j,
                                    cusfloat    water_depth,
                                    cuscomplex& pot_term,
                                    cuscomplex& int_dn_pf_value,
                                    cuscomplex& int_dn_sf_value,
                                    cuscomplex* vel_total
                                )
{
    // Define local variables to work with the fast solver
    const int   ndim                    = 3;
    cusfloat    field_point[ndim]       = { 0.0, 0.0, 0.0 };
    cusfloat    field_point_i[ndim]     = { 0.0, 0.0, 0.0 };
    cusfloat    vel_0[ndim]             = { 0.0, 0.0, 0.0 };
    cusfloat    vel_1[ndim]             = { 0.0, 0.0, 0.0 };
    cusfloat    vel_2[ndim]             = { 0.0, 0.0, 0.0 };
    cusfloat    vel_3[ndim]             = { 0.0, 0.0, 0.0 };
    cusfloat    vel_4[ndim]             = { 0.0, 0.0, 0.0 };
    cusfloat    vel_5[ndim]             = { 0.0, 0.0, 0.0 };
    cusfloat    vel_total_pf[ndim];     clear_vector( ndim, vel_total_pf );

    cusfloat    pot_0                   = 0.0;
    cusfloat    pot_1                   = 0.0;
    cusfloat    pot_2                   = 0.0;
    cusfloat    pot_3                   = 0.0;
    cusfloat    pot_4                   = 0.0;
    cusfloat    pot_5                   = 0.0;

    // Copy data to field point variable
    field_point[0] = panel_j->center[0];
    field_point[1] = panel_j->center[1];
    field_point[2] = panel_j->center[2];
    
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
    if constexpr( freq_regime == FreqRegimeE::REGULAR || freq_regime == FreqRegimeE::ASYMPT_LOW )
    {
        vel_total[0]    = - ( vel_0[0] + vel_1[0] + vel_2[0] + vel_3[0] + vel_4[0] + vel_5[0] );
        vel_total[1]    = - ( vel_0[1] + vel_1[1] + vel_2[1] + vel_3[1] + vel_4[1] + vel_5[1] );
        vel_total[2]    = - ( vel_0[2] + vel_1[2] + vel_2[2] + vel_3[2] - vel_4[2] + vel_5[2] );

        STATIC_COND( ONLY_PF, vel_total_pf[2] = - ( - vel_0[2] + vel_1[2] + vel_2[2] - vel_3[2] + vel_4[2] + vel_5[2] );    )

        pot_term        = ( pot_0 + pot_1 + pot_2 + pot_3 + pot_4 + pot_5 ) / 4.0 / PI;
    }
    else // freq_regime == FREQ_REGIME_INF_FREQ
    {
        vel_total[0]    = - ( vel_0[0] + vel_1[0] - vel_2[0] - vel_3[0] - vel_4[0] - vel_5[0] );
        vel_total[1]    = - ( vel_0[1] + vel_1[1] - vel_2[1] - vel_3[1] - vel_4[1] - vel_5[1] );
        vel_total[2]    = - ( vel_0[2] + vel_1[2] - vel_2[2] - vel_3[2] + vel_4[2] - vel_5[2] );

        STATIC_COND( ONLY_PF, vel_total_pf[2] = - ( - vel_0[2] + vel_1[2] - vel_2[2] + vel_3[2] - vel_4[2] - vel_5[2] );    )

        pot_term        = ( pot_0 + pot_1 - pot_2 - pot_3 - pot_4 - pot_5 ) / 4.0 / PI;
    }


    STATIC_COND( ONLY_PF, vel_total_pf[0] = - vel_total[0].real( ); )
    STATIC_COND( ONLY_PF, vel_total_pf[1] = - vel_total[1].real( ); )

                            
    int_dn_sf_value    = (
                            panel_j->normal_vec[0] * vel_total[0]
                            +
                            panel_j->normal_vec[1] * vel_total[1]
                            +
                            panel_j->normal_vec[2] * vel_total[2]
                        ) / 4.0 / PI;

    STATIC_COND(
                    ONLY_PF,
                    int_dn_pf_value = (
                                            panel_i->normal_vec[0] * vel_total_pf[0]
                                            +
                                            panel_i->normal_vec[1] * vel_total_pf[1]
                                            +
                                            panel_i->normal_vec[2] * vel_total_pf[2]
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
        vel_total[0] = vel_total[0] / 4.0 / PI;
        vel_total[1] = vel_total[1] / 4.0 / PI;
        vel_total[2] = vel_total[2] / 4.0 / PI;
    }

    
}


template<int mode_pf, FreqRegimeE freq_regime, typename T>
void _formulation_kernel_wave( 
                                bool        is_diag,
                                SourceNode* source_i,
                                SourceNode* source_j,
                                T&          gwf_interf,
                                cusfloat    water_depth,
                                cusfloat    nu,
                                bool&       is_john,
                                cuscomplex& pot_term,
                                cuscomplex& int_dn_pf_value,
                                cuscomplex& int_dn_sf_value,
                                cuscomplex* vel_total
                            )
{
    // Get memory address of the panel jth
    PanelGeom*      panel_j = source_j->panel;

    // Reset source and field point in the interface
    gwf_interf.set_source_i( source_i, 1.0 );
    gwf_interf.set_source_j( source_j );
    
    // Calculate distance in between field point and source
    is_john = source_i->panel->is_john(  
                                            source_j->position, 
                                            water_depth
                                        );
    

    // Declare local auxiliar variables
    cuscomplex  wave_fcn_value          = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dn_sf_value    = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dn_pf_value    = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dx_value       = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dy_value       = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dz_value       = cuscomplex( 0.0, 0.0 );

    // Reset to cero variables in other to avoid spurious data
    pot_term        = 0.0;
    int_dn_sf_value = 0.0;
    int_dn_pf_value = 0.0;

    // Integrate green function normal derivative along the current panel
    if ( is_diag )
    {
        if ( panel_j->type == PanelTypeE::DIFFRAC )
        {
            int_dn_sf_value     = cuscomplex( 0.5, 0.0 );
            int_dn_pf_value     = cuscomplex( 0.5, 0.0 );

            _QUADRATURE_PANEL_T( NUM_GP2, G_ON, DGDR_OFF, DGDZ_OFF, FSLID_OFF )
            
            pot_term   =   wave_fcn_value / 4.0 / PI;
        }
        else if ( panel_j->type == PanelTypeE::INT_LID )
        {
            int_dn_sf_value = -cuscomplex( 1.0, 0.0 );
            int_dn_pf_value = -cuscomplex( 1.0, 0.0 );
            
            _QUADRATURE_PANEL_T( NUM_GP2, G_ON, DGDR_OFF, DGDZ_OFF, FSLID_ON )
            
            if constexpr( freq_regime == FreqRegimeE::REGULAR )
            {
                wave_fcn_value += calculate_fs_singularity( nu, source_i->panel->area, source_i->panel->free_surface_log_int );
            }
            pot_term = wave_fcn_value / 4.0 / PI;
            
        }
        else if ( panel_j->type == PanelTypeE::EXT_LID )
        {
            int_dn_sf_value = cuscomplex( 1.0, 0.0 );
            int_dn_pf_value = cuscomplex( 1.0, 0.0 );
            
            _QUADRATURE_PANEL_T( NUM_GP2, G_ON, DGDR_OFF, DGDZ_OFF, FSLID_ON )
            
            if constexpr( freq_regime == FreqRegimeE::REGULAR )
            {
                wave_fcn_value += calculate_fs_singularity( nu, source_i->panel->area, source_i->panel->free_surface_log_int );
            }
            pot_term = wave_fcn_value / 4.0 / PI;
            
        }
        
    }
    else
    {
        wave_fcn_value          = 0.0;
        wave_fcn_dn_sf_value    = 0.0;
        wave_fcn_dn_pf_value    = 0.0;

        if ( is_john && freq_regime == FreqRegimeE::REGULAR )
        {
            quadrature_panel_t<
                                    PanelGeom, 
                                    GWFcnsInterfaceT<NUM_GP2>, 
                                    john_series<NUM_GP2, G_ON, DGDR_ON, DGDZ_ON>, 
                                    NUM_GP
                                >( 
                                    source_i->panel, 
                                    gwf_interf, 
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

        pot_term        =   wave_fcn_value       / 4.0 / PI;
        int_dn_sf_value =   wave_fcn_dn_sf_value / 4.0 / PI;
        int_dn_pf_value =   wave_fcn_dn_pf_value / 4.0 / PI;
        vel_total[0]    =   wave_fcn_dx_value    / 4.0 / PI;
        vel_total[1]    =   wave_fcn_dy_value    / 4.0 / PI;
        vel_total[2]    =   wave_fcn_dz_value    / 4.0 / PI;

    }
}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
template<FreqRegimeE freq_regime>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::_build_steady_matrixes( void )
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
            if constexpr( freq_regime == FreqRegimeE::REGULAR || freq_regime == FreqRegimeE::ASYMPT_LOW )
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

            if ( 
                    ( std::isnan( int_value_sf.real( ) ) || std::isnan( int_value_sf.imag( ) ) )
                    ||
                    ( std::isnan( pot_term.real( ) ) || std::isnan( pot_term.imag( ) ) )
                    ||
                    ( std::isnan( int_value_pf.real( ) ) || std::isnan( int_value_pf.imag( ) ) )
                )
            {
                std::cout << "Error: NaN value found when building steady source system matrix at row " << row_count << " and column " << col_count << " for SF integral." << std::endl;
                std::cout << "I: " << i << "J : " << j << std::endl;
                MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
            }

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


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::_build_rhs( 
                                                        cusfloat w
                                                    )
{
    // Declare local variables
    cusfloat    ang_freq_2  = w * w;
    int         col_count   = 0;
    int         dofs_offset = this->_input->dofs_np * this->_mesh_gp->meshes_np * this->_pot_gp->sysmat_nrows;
    int         index       = 0;
    int         index_rm    = 0;
    PanelGeom*  panel_j     = nullptr;
    int         row_count   = 0;
    SourceNode* source_i    = nullptr;
    cuscomplex  wave_dx     = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_dy     = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_dz     = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_pot    = cuscomplex( 0.0, 0.0 );

    // Calculate wave dependent parameters
    cusfloat    k           = w2k( w, this->_input->water_depth, this->_input->grav_acc );

    // Clear potential rhs to avoid spurious valures
    STATIC_COND( ONLY_PF, clear_vector( this->_pf_gp->sysmat_nrows * ( this->_input->dofs_np * this->_mesh_gp->meshes_np + this->_input->heads_np ), this->_ppf_rhs ); )

    // Calculate potential formulation rhs
    if constexpr( ONLY_PF )
    {
        for ( int i=this->_solver->start_col_0; i<this->_solver->end_col_0; i++ )
        {
            source_i    = this->_mesh_gp->source_nodes[i];
            // Determine which body owns this source panel
            const int   ib_src  = this->_mesh_gp->get_body_id_sn( i );
            row_count   = 0;
            for ( int j=this->_solver->start_row_0; j<this->_solver->end_row_0; j++ )
            {
                ROW_MAJOR_INDEX( index_rm, row_count, col_count, this->_solver->num_cols_local )

                if ( source_i->panel->type == PanelTypeE::DIFFRAC )
                {
                    // Radiation: each body-DOF pair gets its own column
                    for ( int id=0; id<this->_input->dofs_np; id++ )
                    {
                        index                   = ( ib_src * this->_input->dofs_np + id ) * this->_pf_gp->sysmat_nrows + j; 
                        this->_ppf_rhs[index]   += (
                                                        source_i->normal_vec[id]
                                                        *
                                                        source_i->panel->is_move_f
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
                                                                        source_i->panel->center[0],
                                                                        source_i->panel->center[1],
                                                                        source_i->panel->center[2],
                                                                        this->_input->heads[id]
                                                                    );
    
                        wave_dy     =   wave_potential_fo_space_dy(
                                                                        1.0,
                                                                        w,
                                                                        k,
                                                                        this->_input->water_depth,
                                                                        this->_input->grav_acc,
                                                                        source_i->panel->center[0],
                                                                        source_i->panel->center[1],
                                                                        source_i->panel->center[2],
                                                                        this->_input->heads[id]
                                                                    );
    
                        wave_dz     =   wave_potential_fo_space_dz(
                                                                        1.0,
                                                                        w,
                                                                        k,
                                                                        this->_input->water_depth,
                                                                        this->_input->grav_acc,
                                                                        source_i->panel->center[0],
                                                                        source_i->panel->center[1],
                                                                        source_i->panel->center[2],
                                                                        this->_input->heads[id]
                                                                    );
                        
                        // Calculate normal derivative of the wave flow velocities for the jth panel
                        index                   = dofs_offset + id * this->_pf_gp->sysmat_nrows + j; 
                        this->_ppf_rhs[index]   += -(
                                                        wave_dx * source_i->normal_vec[0]
                                                        +
                                                        wave_dy * source_i->normal_vec[1]
                                                        +
                                                        wave_dz * source_i->normal_vec[2]
                                                    ) * this->_pot_gp->sysmat[index_rm];
                    }
                }
                else if ( source_i->panel->type == PanelTypeE::EXT_LID )
                {
                    for ( int id=0; id<this->_input->heads_np; id++ )
                    {
                        // Get wave potential derivatives for the panel
                        wave_pot     =   wave_potential_fo_space(
                                                                    1.0,
                                                                    w,
                                                                    k,
                                                                    this->_input->water_depth,
                                                                    this->_input->grav_acc,
                                                                    source_i->panel->center[0],
                                                                    source_i->panel->center[1],
                                                                    source_i->panel->center[2],
                                                                    this->_input->heads[id]
                                                                );
                        
                        // Calculate normal derivative of the wave flow velocities for the jth panel
                        index                   = dofs_offset + id * this->_pf_gp->sysmat_nrows + j; 
                        this->_ppf_rhs[index]   += -(
                                                        cuscomplex( 0.0, 1.0 )
                                                        *
                                                        ang_freq_2
                                                        *
                                                        source_i->panel->ext_lid_damp_f
                                                        *
                                                        wave_pot
                                                        /
                                                        this->_input->grav_acc
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
                            this->_pf_gp->sysmat_nrows * ( this->_input->dofs_np * this->_mesh_gp->meshes_np + this->_input->heads_np ),
                            mpi_cuscomplex,
                            MPI_SUM,
                            MPI_COMM_WORLD
                        );
    }


    // Calculate source formulation rhs — radiation: one column per (body, DOF) pair
    for ( int ib_col=0; ib_col<this->_mesh_gp->meshes_np; ib_col++ )
    {
        for ( int id=0; id<this->_input->dofs_np; id++ )
        {
            for ( int j=this->_solver->start_row_0; j<this->_solver->end_row_0; j++ )
            {
                panel_j         = this->_mesh_gp->source_nodes[j]->panel;
                const int ib_j  = this->_mesh_gp->get_body_id_sn( j );
                index           = ( ib_col * this->_input->dofs_np + id ) * this->_sf_gp->sysmat_nrows + j;
                if ( panel_j->type == PanelTypeE::DIFFRAC )
                {
                    // BN: only panels of body ib_col contribute to radiation mode (ib_col, id)
                    this->_sf_gp->field_values[index] = ( ib_j == ib_col )
                        ? ( this->_mesh_gp->source_nodes[j]->normal_vec[id] * this->_mesh_gp->source_nodes[j]->panel->is_move_f )
                        : cuscomplex( 0.0, 0.0 );
                }
                else if ( panel_j->type == PanelTypeE::INT_LID )
                {
                    this->_sf_gp->field_values[index] = 0.0;
                }
            }
        }
    }

    for ( int i=0; i<this->_input->heads_np; i++ )
    {
        for ( int j=this->_solver->start_row_0; j<this->_solver->end_row_0; j++ )
        {
            panel_j = this->_mesh_gp->source_nodes[j]->panel;
            index   = dofs_offset + i * this->_sf_gp->sysmat_nrows + j; 
            if ( panel_j->type == PanelTypeE::DIFFRAC )
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
            else if ( panel_j->type == PanelTypeE::INT_LID )
            {
                this->_sf_gp->field_values[index]  = 0.0;
            }
            else if ( panel_j->type == PanelTypeE::EXT_LID )
            {
                // Get wave potential derivatives for the panel
                wave_pot    =   wave_potential_fo_space(
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
                this->_sf_gp->field_values[index]  = -(
                                                            cuscomplex( 0.0, 1.0 )
                                                            *
                                                            ang_freq_2
                                                            *
                                                            panel_j->ext_lid_damp_f
                                                            *
                                                            wave_pot
                                                            /
                                                            this->_input->grav_acc
                                                        );
            }
        }
    }

    // Synchronize processes progress status
    MPI_Barrier( MPI_COMM_WORLD );

}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
template<int GP, int mode_fslid, typename GwfInterf, typename RhsFunc, typename IntegrateFunc>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::_process_qtf_rhs_panels(
                                                                        RDDQC*          panel_fields,
                                                                        GwfInterf&      gwf_interf_local,
                                                                        RhsFunc&&       rhs_func,
                                                                        IntegrateFunc&& integrate,
                                                                        bool            mode_wl,
                                                                        cusfloat*       field_points_d,
                                                                        cusfloat*       field_points_x,
                                                                        cusfloat*       field_points_y,
                                                                        cusfloat*       field_points_z,
                                                                        cuscomplex*     qb_rhs
                                                                    )
{
    std::size_t heads_offset    = 0;
    cuscomplex  int_result      = cuscomplex( 0.0, 0.0 );
    PanelGeom*  panel_j         = nullptr;
    std::size_t rhs_pos         = 0;
    SourceNode  source_i        = SourceNode( );

    for ( std::size_t i=0; i<panel_fields->panel_data.size( ); i++ )
    {
        auto* panel_data = &( panel_fields->panel_data[i] );

        source_i.normal_vec = panel_data->panel_geom->normal_vec;
        source_i.panel      = panel_data->panel_geom;

        gwf_interf_local.set_source_i( &source_i, 1.0 );

        assert( panel_data->field_points_np == GP );
        for ( std::size_t id=0; id<panel_data->field_points_np; id++ )
        {
            field_points_x[id]   = panel_data->field_points[ 3 * id + 0 ];
            field_points_y[id]   = panel_data->field_points[ 3 * id + 1 ];
            field_points_z[id]   = panel_data->field_points[ 3 * id + 2 ];
            field_points_d[id]   = 0.0;
        }

        for ( std::size_t j=static_cast<std::size_t>(this->_solver->start_row_0); j<static_cast<std::size_t>(this->_solver->end_row_0); j++ )
        {
            panel_j = this->_mesh_gp->source_nodes[j]->panel;
            gwf_interf_local.set_source_j( this->_mesh_gp->source_nodes[j] );

            gwf_interf_local.template operator()<wave_term_integral<GP, G_ON, DGDR_ON, DGDZ_ON, mode_fslid>>(
                                                                                                                field_points_d,
                                                                                                                field_points_d,
                                                                                                                field_points_x,
                                                                                                                field_points_y,
                                                                                                                field_points_z
                                                                                                            );

            if ( panel_data->panel_geom == panel_j )
            {
                for ( std::size_t id=0; id<panel_data->field_points_np; id++ )
                {
                    gwf_interf_local.dG_dx[id] = 0.0;
                    gwf_interf_local.dG_dy[id] = 0.0;
                    gwf_interf_local.dG_dz[id] = -0.5;
                }
            }

            rhs_func( panel_data, panel_j, gwf_interf_local, mode_wl, panel_data->body_id );

            for ( std::size_t ih1=0; ih1<static_cast<std::size_t>(this->_input->heads_np); ih1++ )
            {
                for ( std::size_t ih2=0; ih2<static_cast<std::size_t>(this->_input->heads_np); ih2++ )
                {
                    heads_offset = ( ih1 * this->_input->heads_np + ih2 ) * panel_data->field_points_np;
                    int_result   = 0.0;
                    integrate( int_result, heads_offset, panel_data->panel_geom );

                    rhs_pos                   = ( ih1 * this->_input->heads_np + ih2 ) * this->_pf_gp->sysmat_nrows + j;
                    this->_ppf_rhs[rhs_pos]  += int_result;
                }
            }
        }
    }
}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
template<typename RhsFunc>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::_process_qtf_rhs_annuli(
                                                                        RDDQC*                  annulus_fields,
                                                                        const cusfloat*         field_weights,
                                                                        GWFcnsInterfaceT<1>&    gwf_interf_local,
                                                                        RhsFunc&&               rhs_func,
                                                                        cusfloat*               field_points_d,
                                                                        cusfloat*               field_points_x,
                                                                        cusfloat*               field_points_y,
                                                                        cusfloat*               field_points_z,
                                                                        cuscomplex*             qb_rhs
                                                                    )
{
    if ( annulus_fields == nullptr || field_weights == nullptr )
    {
        return;
    }

    const std::size_t local_panels = annulus_fields->get_size_local( );
    if ( local_panels == 0 )
    {
        return;
    }

    const std::size_t heads_np   = this->_input->heads_np;
    const std::size_t heads_np2  = heads_np * heads_np;
    const std::size_t rows_local = static_cast<std::size_t>( this->_solver->num_rows_local );
    const std::size_t acc_size   = heads_np2 * rows_local;

    cuscomplex* annuli_acc = generate_empty_vector<cuscomplex>( acc_size );
    clear_vector( acc_size, annuli_acc );

    cusfloat annulus_normal[3] = { 0.0, 0.0, 1.0 };
    PanelGeom annulus_panel_geom;
    annulus_panel_geom.normal_vec[0] = annulus_normal[0];
    annulus_panel_geom.normal_vec[1] = annulus_normal[1];
    annulus_panel_geom.normal_vec[2] = annulus_normal[2];

    SourceNode source_i;
    source_i.normal_vec = annulus_panel_geom.normal_vec;

    gwf_interf_local.set_source_i( &source_i, 1.0 );

    for ( std::size_t i=0; i<local_panels; i++ )
    {
        PanelData<cuscomplex, RDDQTFConfig>* panel_data = &(annulus_fields->panel_data[i]);
        PanelGeom* panel_geom_prev = panel_data->panel_geom;
        panel_data->panel_geom = &annulus_panel_geom;

        field_points_x[0] = panel_data->field_points[0];
        field_points_y[0] = panel_data->field_points[1];
        field_points_z[0] = panel_data->field_points[2];
        field_points_d[0] = 0.0;

        const cusfloat weight = field_weights[i];

        for ( int j=this->_solver->start_row_0; j<this->_solver->end_row_0; j++ )
        {
            PanelGeom* panel_j = this->_mesh_gp->source_nodes[j]->panel;
            gwf_interf_local.set_source_j( this->_mesh_gp->source_nodes[j] );

            gwf_interf_local.template operator( )<wave_term_integral<1, G_ON, DGDR_ON, DGDZ_ON, ON>>(
                                                                                                                field_points_d,
                                                                                                                field_points_d,
                                                                                                                field_points_x,
                                                                                                                field_points_y,
                                                                                                                field_points_z
                                                                                                            );

            rhs_func( panel_data, panel_j, gwf_interf_local, false, 0 );

            const std::size_t row_local = static_cast<std::size_t>( j - this->_solver->start_row_0 );
            for ( std::size_t ih1=0; ih1<heads_np; ih1++ )
            {
                for ( std::size_t ih2=0; ih2<heads_np; ih2++ )
                {
                    const std::size_t heads_offset = ih1 * heads_np + ih2;
                    annuli_acc[ heads_offset * rows_local + row_local ] += qb_rhs[heads_offset] * weight;
                }
            }
        }

        panel_data->panel_geom = panel_geom_prev;
    }

    MPI_Allreduce( MPI_IN_PLACE, annuli_acc, acc_size, mpi_cuscomplex, MPI_SUM, MPI_COMM_WORLD );

    for ( std::size_t ih1=0; ih1<heads_np; ih1++ )
    {
        for ( std::size_t ih2=0; ih2<heads_np; ih2++ )
        {
            const std::size_t heads_offset = ih1 * heads_np + ih2;
            for ( std::size_t jr=0; jr<rows_local; jr++ )
            {
                const std::size_t rhs_pos = heads_offset * this->_pf_gp->sysmat_nrows
                                            + static_cast<std::size_t>( this->_solver->start_row_0 )
                                            + jr;
                this->_ppf_rhs[rhs_pos]  += annuli_acc[ heads_offset * rows_local + jr ];
            }
        }
    }

    mkl_free( annuli_acc );
}

template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::_build_rhs_so( 
                                                            QTFTypeE                        qtf_type,
                                                            std::size_t                     freq_i,
                                                            std::size_t                     freq_j,
                                                            cuscomplex*                      raos_hist,
                                                            RDDQC*                          qtf_body_fields,
                                                            RDDQC*                          qtf_body_fields_wl,
                                                            RDDQC*                          qtf_fs_fields,
                                                            RDDQC*                          qtf_fs_fields_wl
                                                        )
{
    this->_build_rhs_so(
                            qtf_type,
                            freq_i,
                            freq_j,
                            raos_hist,
                            qtf_body_fields,
                            qtf_body_fields_wl,
                            qtf_fs_fields,
                            qtf_fs_fields_wl,
                            this->_qtf_annuli_fields
                        );
}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::_build_rhs_so( 
                                                            QTFTypeE                        qtf_type,
                                                            std::size_t                     freq_i,
                                                            std::size_t                     freq_j,
                                                            cuscomplex*                      raos_hist,
                                                            RDDQC*                          qtf_body_fields,
                                                            RDDQC*                          qtf_body_fields_wl,
                                                            RDDQC*                          qtf_fs_fields,
                                                            RDDQC*                          qtf_fs_fields_wl,
                                                            const std::vector<RDDQC*>*      qtf_annuli_fields
                                                        )     
{
    // Declare local variables
    cusfloat    field_points_d[NUM_GP2];
    cusfloat    field_points_x[NUM_GP2];
    cusfloat    field_points_y[NUM_GP2];
    cusfloat    field_points_z[NUM_GP2];

    // Allocate space for intermediate variables for the body rhs calculation
    std::size_t qb_size     = this->_input->heads_np * this->_input->heads_np * qtf_body_fields->panel_data[0].field_points_np;
    cuscomplex* qb_rhs      = generate_empty_vector<cuscomplex>( qb_size );

    // Clean rhs vector
    // Clear potential rhs to avoid spurious valures
    STATIC_COND( ONLY_PF, clear_vector( this->_pf_gp->sysmat_nrows * this->_pf_gp->fields_np, this->_ppf_rhs ); )

    auto qb_rhs_func = [&](auto* panel_data, PanelGeom* panel_geom, auto& gwf_interf_local, bool mode_wl, int rhs_body_id)
    {
        so_potential_qb_rhs(
                                qtf_type,
                                this->_input,
                                freq_i,
                                freq_j,
                                raos_hist,
                                rhs_body_id,
                                panel_geom->body_cog,
                                panel_data,
                                gwf_interf_local.G,
                                gwf_interf_local.dG_dx,
                                gwf_interf_local.dG_dy,
                                gwf_interf_local.dG_dz,
                                &(this->_wdso),
                                mode_wl,
                                qb_rhs
                            );
    };

    auto qf_rhs_func = [&](auto* panel_data, PanelGeom* panel_geom, auto& gwf_interf_local, bool mode_wl, int rhs_body_id)
    {
        so_potential_qf_rhs(
                                qtf_type,
                                this->_input,
                                freq_i,
                                freq_j,
                                raos_hist,
                                rhs_body_id,
                                panel_geom->body_cog,
                                panel_data,
                                gwf_interf_local.G,
                                gwf_interf_local.dG_dx,
                                gwf_interf_local.dG_dy,
                                gwf_interf_local.dG_dz,
                                &(this->_wdso),
                                mode_wl ? BoundarySubtypeE::WL : BoundarySubtypeE::DIFFRAC,
                                qb_rhs
                            );
    };

    auto qf_rhs_annuli_func = [&](auto* panel_data, PanelGeom* panel_geom, auto& gwf_interf_local, bool , int rhs_body_id)
    {
        so_potential_qf_rhs(
                                qtf_type,
                                this->_input,
                                freq_i,
                                freq_j,
                                raos_hist,
                                rhs_body_id,
                                panel_geom->body_cog,
                                panel_data,
                                gwf_interf_local.G,
                                gwf_interf_local.dG_dx,
                                gwf_interf_local.dG_dy,
                                gwf_interf_local.dG_dz,
                                &(this->_wdso),
                                BoundarySubtypeE::PC,
                                qb_rhs
                            );
    };

    auto integrate_2d = [&](cuscomplex& int_out, std::size_t offset, PanelGeom* panel_geom)
    {
        gauss2d_loop<NUM_GP2>(
                                int_out,
                                [&](int idf){ return qb_rhs[offset + idf]; },
                                panel_geom
                            );
    };

    auto integrate_1d = [&](cuscomplex& int_out, std::size_t offset, PanelGeom* panel_geom)
    {
        gauss1d_loop<NUM_GP>(
                                int_out,
                                qb_rhs + offset,
                                panel_geom->len_wl
                            );
    };

    this->_process_qtf_rhs_panels<NUM_GP2, OFF>(
                                                    qtf_body_fields,
                                                    this->_gwfcns_interf_qb,
                                                    qb_rhs_func,
                                                    integrate_2d,
                                                    false,
                                                    field_points_d,
                                                    field_points_x,
                                                    field_points_y,
                                                    field_points_z,
                                                    qb_rhs
                                            );
    MPI_Barrier( MPI_COMM_WORLD );

    this->_process_qtf_rhs_panels<NUM_GP, ON>(
                                                qtf_body_fields_wl,
                                                this->_gwfcns_interf_wl,
                                                qb_rhs_func,
                                                integrate_1d,
                                                true,
                                                field_points_d,
                                                field_points_x,
                                                field_points_y,
                                                field_points_z,
                                                qb_rhs
                                        );
    MPI_Barrier( MPI_COMM_WORLD );

    this->_process_qtf_rhs_panels<NUM_GP2, ON>(
                                                qtf_fs_fields,
                                                this->_gwfcns_interf_qf,
                                                qf_rhs_func,
                                                integrate_2d,
                                                false,
                                                field_points_d,
                                                field_points_x,
                                                field_points_y,
                                                field_points_z,
                                                qb_rhs
                                        );
    
    this->_process_qtf_rhs_panels<NUM_GP, ON>(
                                                qtf_body_fields_wl,
                                                this->_gwfcns_interf_wl,
                                                qf_rhs_func,
                                                integrate_1d,
                                                true,
                                                field_points_d,
                                                field_points_x,
                                                field_points_y,
                                                field_points_z,
                                                qb_rhs
                                        );

    this->_process_qtf_rhs_panels<NUM_GP, ON>(
                                                qtf_fs_fields_wl,
                                                this->_gwfcns_interf_wl,
                                                qf_rhs_func,
                                                integrate_1d,
                                                true,
                                                field_points_d,
                                                field_points_x,
                                                field_points_y,
                                                field_points_z,
                                                qb_rhs
                                        );
    MPI_Barrier( MPI_COMM_WORLD );

    if ( qtf_annuli_fields != nullptr )
    {
        for ( auto* annulus_field : *qtf_annuli_fields )
        {
            if ( annulus_field == nullptr || !annulus_field->has_field_weights( ) )
            {
                continue;
            }

            this->_process_qtf_rhs_annuli(
                                                annulus_field,
                                                annulus_field->get_field_weights( ),
                                                this->_gwfcns_interf_qfa,
                                                qf_rhs_annuli_func,
                                                field_points_d,
                                                field_points_x,
                                                field_points_y,
                                                field_points_z,
                                                qb_rhs
                                        );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );

    // Add Far Field contribution to the rhs
    if ( this->_mpi_config->is_root( ) )
    {
        this->_process_far_field( qtf_type, freq_i, freq_j, this->_ppf_rhs );
    }

    // Add all processors contributions to the rhs vector
    MPI_Allreduce(
                        MPI_IN_PLACE,
                        this->_ppf_rhs,
                        this->_pf_gp->sysmat_nrows * ( this->_input->heads_np * this->_input->heads_np ),
                        mpi_cuscomplex,
                        MPI_SUM,
                        MPI_COMM_WORLD
                    );

    // Deallocate local heap variables
    mkl_free( qb_rhs );
    
}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
template<FreqRegimeE freq_regime>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::_build_wave_matrixes( 
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
    auto        gwf_interf              = this->_gwfcns_interf;
    int         index_cm                = 0;
    int         index_rm                = 0;
    bool        is_john                 = false;
    cuscomplex  int_value               = cuscomplex( 0.0, 0.0 );
    cuscomplex  int_dn_sf_value         = cuscomplex( 0.0, 0.0 );
    cuscomplex  int_dn_pf_value         = cuscomplex( 0.0, 0.0 );
    cuscomplex  int_scale_f_i           = cuscomplex( 1.0, 0.0 );
    cuscomplex  int_scale_f_j           = cuscomplex( 1.0, 0.0 );
    PanelGeom*  panel_j                 = nullptr;
    int         row_count               = 0;
    SourceNode* source_i                = nullptr;
    cuscomplex  wave_fcn_value          = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dn_sf_value    = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dn_pf_value    = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dx_value       = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dy_value       = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dz_value       = cuscomplex( 0.0, 0.0 );
    
    // Calculate wave dependent parameters
    cusfloat    nu                      = pow2s( w ) / this->_input->grav_acc;
    
    for ( int i=this->_solver->start_col_0; i<this->_solver->end_col_0; i++ )
    {
        // Get memory address of the ith panel
        source_i = this->_mesh_gp->source_nodes[i];
        gwf_interf.set_source_i( source_i, 1.0 );

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=this->_solver->start_row_0; j<this->_solver->end_row_0; j++ )
        {
            // Get memory address of the panel jth
            panel_j = this->_mesh_gp->source_nodes[j]->panel;
            gwf_interf.set_source_j( this->_mesh_gp->source_nodes[j] );
            
            // Check if the field point is in the John region
            is_john = source_i->panel->is_john( 
                                                    this->_mesh_gp->source_nodes[j]->position, 
                                                    this->_input->water_depth 
                                                );

            
            int_value       = 0.0;
            int_dn_sf_value = 0.0;
            int_dn_pf_value = 0.0;
            // Integrate green function normal derivative along the current panel
            if ( i == j )
            {
                if ( panel_j->type == PanelTypeE::DIFFRAC )
                {
                    int_dn_sf_value     = cuscomplex( 0.5, 0.0 );
                    int_dn_pf_value     = cuscomplex( 0.5, 0.0 );

                    _QUADRATURE_PANEL_T( NUM_GP2, G_ON, DGDR_OFF, DGDZ_OFF, FSLID_OFF )
                    
                    int_value   =   wave_fcn_value / 4.0 / PI;
                }
                else if ( panel_j->type == PanelTypeE::INT_LID )
                {
                    int_dn_sf_value     = -cuscomplex( 1.0, 0.0 );
                    int_dn_pf_value     = -cuscomplex( 1.0, 0.0 );
                    
                    _QUADRATURE_PANEL_T( NUM_GP2, G_ON, DGDR_OFF, DGDZ_OFF, FSLID_ON )
                    
                    if constexpr( freq_regime == FreqRegimeE::REGULAR )
                    {
                        wave_fcn_value += calculate_fs_singularity( nu, source_i->panel->area, source_i->panel->free_surface_log_int );
                    }
                    int_value           = wave_fcn_value / 4.0 / PI;
                    
                }
                else if ( panel_j->type == PanelTypeE::EXT_LID )
                {
                    int_dn_sf_value     = cuscomplex( 1.0, 0.0 );
                    int_dn_pf_value     = cuscomplex( 1.0, 0.0 );
                    
                    _QUADRATURE_PANEL_T( NUM_GP2, G_ON, DGDR_OFF, DGDZ_OFF, FSLID_ON )
                    
                    if constexpr( freq_regime == FreqRegimeE::REGULAR )
                    {
                        wave_fcn_value += calculate_fs_singularity( nu, source_i->panel->area, source_i->panel->free_surface_log_int );
                    }
                    int_value           = wave_fcn_value / 4.0 / PI;
                    
                }
                
            }
            else
            {
                wave_fcn_value          = 0.0;
                wave_fcn_dn_sf_value    = 0.0;
                wave_fcn_dn_pf_value    = 0.0;

                if ( is_john && freq_regime == FreqRegimeE::REGULAR )
                {
                    quadrature_panel_t<
                                            PanelGeom, 
                                            GWFcnsInterfaceT<NUM_GP2>, 
                                            john_series<NUM_GP2, G_ON, DGDR_ON, DGDZ_ON>, 
                                            NUM_GP
                                        >( 
                                            source_i->panel, 
                                            gwf_interf, 
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

                int_value       =   wave_fcn_value / 4.0 / PI;
                int_dn_sf_value =   wave_fcn_dn_sf_value / 4.0 / PI;
                int_dn_pf_value =   wave_fcn_dn_pf_value / 4.0 / PI;

            }

            // Apply the integral value accordingly
            COL_MAJOR_INDEX( index_cm, row_count, col_count, this->_solver->num_rows_local )
            ROW_MAJOR_INDEX( index_rm, row_count, col_count, this->_solver->num_cols_local )

            // Calculate scaling factors for the integrals according to the panel types
            if ( i != j )
            {
                int_scale_f_i = panel_scale_factor( source_i->panel, nu );
                int_scale_f_j = panel_scale_factor( panel_j, nu );
            }
            else
            {
                int_scale_f_i = cuscomplex( 1.0, 0.0 );
                int_scale_f_j = cuscomplex( 1.0, 0.0 );
            }

            const bool       is_john_regular    = ( is_john && freq_regime == FreqRegimeE::REGULAR );
            const bool       i_is_diffrac       = ( source_i->panel->type == PanelTypeE::DIFFRAC );
            const bool       j_is_diffrac       = ( panel_j->type == PanelTypeE::DIFFRAC );
            const cuscomplex pot_total          = ( this->_pot_gp->sysmat_steady[index_rm] + int_value );

            this->_pot_gp->sysmat[index_rm]     = is_john_regular ? int_value : pot_total;
            this->_sf_gp->sysmat[index_cm]      = int_scale_f_j * ( j_is_diffrac ? int_dn_sf_value : ( is_john_regular ? int_value : pot_total ) );

            STATIC_COND(
                            ONLY_PF,
                            this->_pf_gp->sysmat[index_cm]  = int_scale_f_i * (
                                                                                    i_is_diffrac
                                                                                    ? ( is_john_regular ? int_dn_pf_value : ( this->_pf_gp->sysmat_steady[index_cm] + int_dn_pf_value ) )
                                                                                    : ( is_john_regular ? int_value : pot_total )
                                                                                );
                        )


            // if ( is_john && freq_regime == FreqRegimeE::REGULAR )
            // {
            //     this->_pot_gp->sysmat[index_rm] = int_scale_f * int_value;
            //     this->_sf_gp->sysmat[index_cm]  =  int_scale_f * int_dn_sf_value;

            //     STATIC_COND( ONLY_PF, this->_pf_gp->sysmat[index_cm]  = int_scale_f * int_dn_pf_value; )
            // }
            // else
            // {
            //     this->_pot_gp->sysmat[index_rm] = int_scale_f * ( this->_pot_gp->sysmat_steady[index_rm] + int_value );
            //     this->_sf_gp->sysmat[index_cm]  = int_scale_f * ( this->_sf_gp->sysmat_steady[index_cm]  + int_dn_sf_value );

            //     STATIC_COND( ONLY_PF, this->_pf_gp->sysmat[index_cm]  = int_scale_f * ( this->_pf_gp->sysmat_steady[index_cm] + int_dn_pf_value ); )
            // }

            if ( 
                    ( std::isnan( int_dn_sf_value.real( ) ) || std::isnan( int_dn_sf_value.imag( ) ) )
                    ||
                    ( std::isnan( int_value.real( ) ) || std::isnan( int_value.imag( ) ) )
                    ||
                    ( std::isnan( int_dn_pf_value.real( ) ) || std::isnan( int_dn_pf_value.imag( ) ) )
                )
            {
                std::cout << "Error: NaN value found when building wave source system matrix at row " << row_count << " and column " << col_count << " for SF integral." << std::endl;
                std::cout << "I: " << i << "J : " << j << std::endl;
                MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
            }

            // Advance row count
            row_count++;
        }
        
        // Advance column count
        col_count++;

    }
    MPI_Barrier( MPI_COMM_WORLD );

}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
template<FreqRegimeE freq_regime>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::_build_wave_matrixes_2(
                                                                                    cusfloat w
                                                                                )
{
    // Clean system matrixes
                                this->_sf_gp->clear_sysmat( );
    STATIC_COND( ONLY_PF,       this->_pf_gp->clear_sysmat( );          )
    STATIC_COND( !(ONLY_PF),    this->_pot_gp->clear_sysmat( );         )
    STATIC_COND( !(ONLY_PF),    this->_pot_gp->clear_field_values( );   )

    STATIC_COND( ONLY_PF, clear_vector( this->_pf_gp->sysmat_nrows * ( this->_input->dofs_np * this->_mesh_gp->meshes_np + this->_input->heads_np ), this->_ppf_rhs ); )

    // Declare local variables
    cusfloat    ang_freq_2              = w * w;
    int         col_count               = 0;
    int         dofs_offset_pf          = this->_input->dofs_np * this->_mesh_gp->meshes_np * this->_pf_gp->sysmat_nrows;
    int         dofs_offset_sf          = this->_input->dofs_np * this->_mesh_gp->meshes_np * this->_sf_gp->sysmat_nrows;
    auto        gwf_interf              = this->_gwfcns_interf;
    int         index                   = 0;
    int         index_cm                = 0;
    int         index_rm                = 0;
    bool        is_john                 = false;
    cuscomplex  int_dn_sf_st            = cuscomplex( 0.0, 0.0 );
    cuscomplex  int_dn_sf_wv            = cuscomplex( 0.0, 0.0 );
    cuscomplex  int_dn_pf_st            = cuscomplex( 0.0, 0.0 );
    cuscomplex  int_dn_pf_wv            = cuscomplex( 0.0, 0.0 );
    cuscomplex  int_scale_f_i           = cuscomplex( 1.0, 0.0 );
    cuscomplex  int_scale_f_j           = cuscomplex( 1.0, 0.0 );
    PanelGeom*  panel_j                 = nullptr;
    cuscomplex  pot_term_st             = cuscomplex( 0.0, 0.0 );
    cuscomplex  pot_term_wv             = cuscomplex( 0.0, 0.0 );
    int         row_count               = 0;
    SourceNode* source_i                = nullptr;
    cuscomplex  vel_total_st[3];
    cuscomplex  vel_total_wv[3];
    cuscomplex  wave_dx                 = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_dy                 = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_dz                 = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_pot                = cuscomplex( 0.0, 0.0 );

    // Calculate wave dependent parameters
    cusfloat    k                       = w2k( w, this->_input->water_depth, this->_input->grav_acc );
    cusfloat    nu                      = pow2s( w ) / this->_input->grav_acc;

    // Clear potential rhs to avoid spurious valures
    STATIC_COND( ONLY_PF, clear_vector( this->_pf_gp->sysmat_nrows * ( this->_input->dofs_np * this->_mesh_gp->meshes_np + this->_input->heads_np ), this->_ppf_rhs ); )

    for ( int i=this->_solver->start_col_0; i<this->_solver->end_col_0; i++ )
    {
        // Get memory address of the ith panel
        source_i = this->_mesh_gp->source_nodes[i];
        // Determine which body owns this source panel
        const int ib_src = this->_mesh_gp->get_body_id_sn( i );
        gwf_interf.set_source_i( source_i, 1.0 );

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=this->_solver->start_row_0; j<this->_solver->end_row_0; j++ )
        {
            // Get memory address of the panel jth
            panel_j = this->_mesh_gp->source_nodes[j]->panel;
            gwf_interf.set_source_j( this->_mesh_gp->source_nodes[j] );

            // Define row and column major indexes
            COL_MAJOR_INDEX( index_cm, row_count, col_count, this->_solver->num_rows_local )
            ROW_MAJOR_INDEX( index_rm, row_count, col_count, this->_solver->num_cols_local )

            // Calculate steady contribution
            if constexpr ( recalc_steady == RecalcSteadyE::ON )
            {
                _formulation_kernel_steady<
                                            mode_pf,
                                            freq_regime
                                        >
                                        (
                                            i == j,
                                            this->_mesh_gp->panels[i],
                                            this->_mesh_gp->panels_mirror[i],
                                            this->_mesh_gp->panels[j],
                                            this->_input->water_depth,
                                            pot_term_st,
                                            int_dn_pf_st,
                                            int_dn_sf_st,
                                            vel_total_st
                                        );
            }
            else
            {
                pot_term_st     = this->_pot_gp->sysmat_steady[ index_rm ];
                int_dn_pf_st    = this->_pf_gp->sysmat_steady[ index_cm ];
                int_dn_sf_st    = this->_sf_gp->sysmat_steady[ index_cm ];
            }

            // Calculate wave contribution
            _formulation_kernel_wave<
                                        mode_pf,
                                        freq_regime
                                    >
                                    (
                                        i == j,
                                        this->_mesh_gp->source_nodes[i],
                                        this->_mesh_gp->source_nodes[j],
                                        this->_gwfcns_interf,
                                        this->_input->water_depth,
                                        nu,
                                        is_john,
                                        pot_term_wv,
                                        int_dn_pf_wv,
                                        int_dn_sf_wv,
                                        vel_total_wv
                                    );

            // Calculate scaling factors for the integrals according to the panel types
            if ( i != j )
            {
                int_scale_f_i = panel_scale_factor( source_i->panel, nu );
                int_scale_f_j = panel_scale_factor( panel_j, nu );
            }
            else
            {
                int_scale_f_i = cuscomplex( 1.0, 0.0 );
                int_scale_f_j = cuscomplex( 1.0, 0.0 );
            }

            const bool       is_john_regular    = ( is_john && freq_regime == FreqRegimeE::REGULAR );
            const bool       i_is_diffrac       = ( source_i->panel->type == PanelTypeE::DIFFRAC );
            const bool       j_is_diffrac       = ( panel_j->type == PanelTypeE::DIFFRAC );
            const cuscomplex pot_total          = is_john_regular ? pot_term_wv  : ( pot_term_st  + pot_term_wv  );
            const cuscomplex int_dn_sf_value    = is_john_regular ? int_dn_sf_wv : ( int_dn_sf_st + int_dn_sf_wv );
            const cuscomplex int_dn_pf_value    = is_john_regular ? int_dn_pf_wv : ( int_dn_pf_st + int_dn_pf_wv );

            this->_sf_gp->sysmat[index_cm]      = int_scale_f_j * ( j_is_diffrac ? int_dn_sf_value : pot_total );

            if constexpr( ONLY_PF )
            {
                this->_pf_gp->sysmat[index_cm]  = int_scale_f_i * ( i_is_diffrac ? int_dn_pf_value : pot_total );
            }
            else
            {
                this->_pot_gp->sysmat[index_rm] = pot_total;
            }

            if ( 
                    ( std::isnan( this->_sf_gp->sysmat[index_cm].real( ) ) || std::isnan( this->_sf_gp->sysmat[index_cm].imag( ) ) )
                    ||
                    ( std::isnan( pot_total.real( ) ) || std::isnan( pot_total.imag( ) ) )
                    ||
                    ( std::isnan( int_dn_pf_wv.real( ) ) || std::isnan( int_dn_pf_wv.imag( ) ) )
                )
            {
                std::cout << "Error: NaN value found when building wave source system matrix at row " << row_count << " and column " << col_count << " for SF integral." << std::endl;
                std::cout << "I: " << i << "J : " << j << std::endl;
                MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
            }

            // Calculate potential formulation rhs
            if constexpr( ONLY_PF )
            {
                if ( source_i->panel->type == PanelTypeE::DIFFRAC )
                {
                    // Radiation: column is (ib_src * dofs_np + id) — per body, per DOF
                    for ( int id=0; id<this->_input->dofs_np; id++ )
                    {
                        index                   = ( ib_src * this->_input->dofs_np + id ) * this->_pf_gp->sysmat_nrows + j;
                        this->_ppf_rhs[index]   += (
                                                        source_i->normal_vec[id]
                                                        *
                                                        source_i->panel->is_move_f
                                                    ) * pot_total;
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
                                                                        source_i->panel->center[0],
                                                                        source_i->panel->center[1],
                                                                        source_i->panel->center[2],
                                                                        this->_input->heads[id]
                                                                    );

                        wave_dy     =   wave_potential_fo_space_dy(
                                                                        1.0,
                                                                        w,
                                                                        k,
                                                                        this->_input->water_depth,
                                                                        this->_input->grav_acc,
                                                                        source_i->panel->center[0],
                                                                        source_i->panel->center[1],
                                                                        source_i->panel->center[2],
                                                                        this->_input->heads[id]
                                                                    );

                        wave_dz     =   wave_potential_fo_space_dz(
                                                                        1.0,
                                                                        w,
                                                                        k,
                                                                        this->_input->water_depth,
                                                                        this->_input->grav_acc,
                                                                        source_i->panel->center[0],
                                                                        source_i->panel->center[1],
                                                                        source_i->panel->center[2],
                                                                        this->_input->heads[id]
                                                                    );

                        // Calculate normal derivative of the wave flow velocities for the jth panel
                        index                   = dofs_offset_pf + id * this->_pf_gp->sysmat_nrows + j;
                        this->_ppf_rhs[index]   += -(
                                                        wave_dx * source_i->normal_vec[0]
                                                        +
                                                        wave_dy * source_i->normal_vec[1]
                                                        +
                                                        wave_dz * source_i->normal_vec[2]
                                                    ) * pot_total;
                    }
                }
                else if ( source_i->panel->type == PanelTypeE::EXT_LID )
                {
                    for ( int id=0; id<this->_input->heads_np; id++ )
                    {
                        // Get wave potential derivatives for the panel
                        wave_pot     =   wave_potential_fo_space(
                                                                    1.0,
                                                                    w,
                                                                    k,
                                                                    this->_input->water_depth,
                                                                    this->_input->grav_acc,
                                                                    source_i->panel->center[0],
                                                                    source_i->panel->center[1],
                                                                    source_i->panel->center[2],
                                                                    this->_input->heads[id]
                                                                );

                        // Calculate normal derivative of the wave flow velocities for the jth panel
                        index                   = dofs_offset_pf + id * this->_pf_gp->sysmat_nrows + j;
                        this->_ppf_rhs[index]   += -(
                                                        cuscomplex( 0.0, 1.0 )
                                                        *
                                                        ang_freq_2
                                                        *
                                                        source_i->panel->ext_lid_damp_f
                                                        *
                                                        wave_pot
                                                        /
                                                        this->_input->grav_acc
                                                    ) * pot_total;
                    }
                }
            }

            // Advance row count
            row_count++;
        }

        // Advance column count
        col_count++;

    }

    // Sum up contributions along processors
    STATIC_COND(
                    ONLY_PF,
                    MPI_Allreduce(
                                        this->_ppf_rhs,
                                        this->_pf_gp->field_values,
                                        this->_pf_gp->sysmat_nrows * ( this->_input->dofs_np * this->_mesh_gp->meshes_np + this->_input->heads_np ),
                                        mpi_cuscomplex,
                                        MPI_SUM,
                                        MPI_COMM_WORLD
                                    );
                )

    // Calculate source formulation rhs — radiation: one column per (body, DOF) pair
    for ( int ib_col=0; ib_col<this->_mesh_gp->meshes_np; ib_col++ )
    {
        for ( int id=0; id<this->_input->dofs_np; id++ )
        {
            for ( int j=this->_solver->start_row_0; j<this->_solver->end_row_0; j++ )
            {
                panel_j             = this->_mesh_gp->source_nodes[j]->panel;
                const int ib_j      = this->_mesh_gp->get_body_id_sn( j );
                index               = ( ib_col * this->_input->dofs_np + id ) * this->_sf_gp->sysmat_nrows + j;
                if ( panel_j->type == PanelTypeE::DIFFRAC )
                {
                    // BN: only panels of body ib_col contribute to radiation mode (ib_col, id)
                    this->_sf_gp->field_values[index] = ( ib_j == ib_col )
                        ? ( this->_mesh_gp->source_nodes[j]->normal_vec[id]
                            * this->_mesh_gp->source_nodes[j]->panel->is_move_f )
                        : cuscomplex( 0.0, 0.0 );
                }
                else if ( panel_j->type == PanelTypeE::INT_LID )
                {
                    this->_sf_gp->field_values[index] = 0.0;
                }
            }
        }
    }

    for ( int i=0; i<this->_input->heads_np; i++ )
    {
        for ( int j=this->_solver->start_row_0; j<this->_solver->end_row_0; j++ )
        {
            panel_j = this->_mesh_gp->source_nodes[j]->panel;
            index   = dofs_offset_sf + i * this->_sf_gp->sysmat_nrows + j;
            if ( panel_j->type == PanelTypeE::DIFFRAC )
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
            else if ( panel_j->type == PanelTypeE::INT_LID )
            {
                this->_sf_gp->field_values[index]  = 0.0;
            }
            else if ( panel_j->type == PanelTypeE::EXT_LID )
            {
                // Get wave potential derivatives for the panel
                wave_pot    =   wave_potential_fo_space(
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
                this->_sf_gp->field_values[index]  = -(
                                                            cuscomplex( 0.0, 1.0 )
                                                            *
                                                            ang_freq_2
                                                            *
                                                            panel_j->ext_lid_damp_f
                                                            *
                                                            wave_pot
                                                            /
                                                            this->_input->grav_acc
                                                        );
            }
        }
    }

    // Synchronize processes progress status
    MPI_Barrier( MPI_COMM_WORLD );

}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
template<FreqRegimeE freq_regime>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::_build_wave_matrixes_so( 
                                                                        QTFTypeE            qtf_type,
                                                                        std::size_t         freq_i,
                                                                        std::size_t         freq_j,
                                                                        cuscomplex          raos_hist,
                                                                        RDDQC*              qtf_body_fields,
                                                                        RDDQC*              qtf_body_fields_wl,
                                                                        RDDQC*              qtf_fs_fields,
                                                                        RDDQC*              qtf_fs_fields_wl,
                                                                        WaveDispersionSO*   wd
                                                                    )
{
    // Clean system matrixes
                            this->_sf_gp->clear_sysmat( );
    STATIC_COND( ONLY_PF,   this->_pf_gp->clear_sysmat( ); )
                            this->_pot_gp->clear_sysmat( );
                            this->_pot_gp->clear_field_values( );
    
    // Declare local variables
    int         body_id                 = 0;
    int         col_count               = 0;
    auto        gwf_interf              = this->_gwfcns_interf;
    int         index_cm                = 0;
    int         index_rm                = 0;
    cuscomplex  int_dn_sf_st            = cuscomplex( 0.0, 0.0 );
    cuscomplex  int_dn_sf_wv            = cuscomplex( 0.0, 0.0 );
    cuscomplex  int_dn_pf_st            = cuscomplex( 0.0, 0.0 );
    cuscomplex  int_dn_pf_wv            = cuscomplex( 0.0, 0.0 );
    cuscomplex  int_scale_f             = cuscomplex( 1.0, 0.0 );
    bool        is_john                 = false;
    cusfloat    log_sing_val            = 0.0;
    PanelGeom*  panel_j                 = nullptr;
    cuscomplex  pot_term                = cuscomplex( 0.0, 0.0 );
    cuscomplex  pot_term_st             = cuscomplex( 0.0, 0.0 );
    cuscomplex  pot_term_wv             = cuscomplex( 0.0, 0.0 );
    int         qtf_sign                = ( qtf_type == QTFTypeE::QTF_DIFF_CODE ) ? -1 : 1;
    int         row_count               = 0;
    SourceNode* source_i                = nullptr;
    cuscomplex  vel_total[3]            ;
    cuscomplex  vel_total_st[3]         ;
    cuscomplex  vel_total_wv[3]         ;
    cusfloat    wi                      = this->_input->frequencies[freq_i];
    cusfloat    wj                      = this->_input->frequencies[freq_j];
    cusfloat    wi_wj                   = wi * qtf_sign * wj;
    cuscomplex  wave_fcn_value          = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dn_sf_value    = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dn_pf_value    = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dx_value       = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dy_value       = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_fcn_dz_value       = cuscomplex( 0.0, 0.0 );
    
    // Calculate wave dependent parameters
    cusfloat    nu                      = pow2s( wi_wj ) / this->_input->grav_acc;

    // Allocate space for intermediate variables for the body rhs calculation
    std::size_t qb_size = this->_input->heads_np * this->_heads_np * qtf_body_fields->panel_data[0]->field_points_np;
    cuscomplex* qb_rhs  = generate_empty_vector<cuscomplex>( qb_size );

    // Clean rhs vector
    // Clear potential rhs to avoid spurious valures
    STATIC_COND( ONLY_PF, clear_vector( this->_pf_gp->sysmat_nrows * this->_pf_gp->fields_np, this->_ppf_rhs ); )
    
    for ( int i=this->_solver->start_col_0; i<this->_solver->end_col_0; i++ )
    {
        // Get memory address of the ith panel
        source_i = this->_mesh_gp->source_nodes[i];
        gwf_interf.set_source_i( source_i, 1.0 );

        // Get current body id for the columns panels
        body_id = this->_mesh_gp->get_body_id( i );

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=this->_solver->start_row_0; j<this->_solver->end_row_0; j++ )
        {
            // Get memory address of the panel jth
            panel_j = this->_mesh_gp->source_nodes[j]->panel;
            gwf_interf.set_source_j( this->_mesh_gp->source_nodes[j] );
            
            // Calculate steady contribution
            _formulation_kernel_steady<
                                            mode_pf,
                                            freq_regime
                                        >
                                        (
                                            i == j,
                                            this->_mesh_gp->panels[i],
                                            this->_mesh_gp->panels_mirror[i],
                                            this->_mesh_gp->panels[j],
                                            this->_input->water_depth,
                                            pot_term_st,
                                            int_dn_pf_st,
                                            int_dn_sf_st,
                                            vel_total_st
                                        );

            // Calculate wave contribution
            _formulation_kernel_wave<
                                        mode_pf,
                                        freq_regime
                                    >
                                    ( 
                                        i == j,
                                        this->_mesh_gp->source_nodes[i],
                                        this->_mesh_gp->source_nodes[j],
                                        this->_gwfcns_interf,
                                        this->_input->water_depth,
                                        nu,
                                        is_john,
                                        pot_term_wv,
                                        int_dn_pf_wv,
                                        int_dn_sf_wv,
                                        vel_total_wv
                                    );

            // Apply the integral value accordingly
            COL_MAJOR_INDEX( index_cm, row_count, col_count, this->_solver->num_rows_local )
            ROW_MAJOR_INDEX( index_rm, row_count, col_count, this->_solver->num_cols_local )

            if ( source_i->panel->type != PanelTypeE::EXT_LID )
            {
                int_scale_f = cuscomplex( 1.0, 0.0 );
            }
            else
            {
                int_scale_f = cuscomplex( 0.0, - source_i->panel->ext_lid_damp_f );
            }

            if ( is_john && freq_regime == FreqRegimeE::REGULAR )
            {
                this->_pot_gp->sysmat[index_rm] = int_scale_f * pot_term_wv;
                this->_sf_gp->sysmat[index_cm]  = int_scale_f * int_dn_sf_wv;

                STATIC_COND( ONLY_PF, this->_pf_gp->sysmat[index_cm]  = int_scale_f * int_dn_pf_wv; )
            }
            else
            {
                this->_pot_gp->sysmat[index_rm] = int_scale_f * ( pot_term_st + pot_term_wv );
                this->_sf_gp->sysmat[index_cm]  = int_scale_f * ( int_dn_sf_st + int_dn_sf_wv );

                STATIC_COND( ONLY_PF, this->_pf_gp->sysmat[index_cm]  = int_scale_f * ( int_dn_pf_st + int_dn_pf_wv ); )
            }

            // Advance row count
            row_count++;

        }
        
        // Advance column count
        col_count++;

    }
    MPI_Barrier( MPI_COMM_WORLD );

    // Calculate Qb WL RHS calculation

    // Deallocate local heap variables
    mkl_free( qb_rhs );

}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
template<typename RDDConfig>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::compute_fields(
                                                                std::size_t                         freq_index,
                                                                cusfloat                            ang_freq,
                                                                cuscomplex*                         raos,
                                                                RadDiffData<cuscomplex, RDDConfig>* rad_diff_data
                                                            )
{
    /*** Get template parameters from Config ***/
    static constexpr int mode_comp      = RDDConfig::mode_comp;
    static constexpr int mode_f         = RDDConfig::mode_f;
    static constexpr int mode_dfdn      = RDDConfig::mode_dfdn;
    static constexpr int mode_dfdc      = RDDConfig::mode_dfdc;
    static constexpr int store_freqs    = RDDConfig::store_freqs;

    // Declare local auxiliary variables
    std::size_t                         body_id             = 0;
    std::size_t                         body_np             = this->_mesh_gp->meshes_np;
    cusfloat                            center_aux[3]       = { 0.0, 0.0, 0.0 };
    std::size_t                         dofs_np             = this->_input->dofs_np;
    std::size_t                         fp_np               = 0;
    cusfloat                            grav_acc            = this->_input->grav_acc;
    std::size_t                         heads_np            = this->_input->heads_np;
    std::size_t                         index_ax            = 0;
    std::size_t                         index_sc            = 0;
    std::size_t                         index_fd            = 0;
    cuscomplex                          int_dn_sf           = cuscomplex( 0.0, 0.0 );
    cuscomplex                          int_dn_sf_st        = cuscomplex( 0.0, 0.0 );
    cuscomplex                          int_dn_sf_wv        = cuscomplex( 0.0, 0.0 );
    cuscomplex                          int_dn_pf_st        = cuscomplex( 0.0, 0.0 );
    cuscomplex                          int_dn_pf_wv        = cuscomplex( 0.0, 0.0 );
    cuscomplex                          int_scale_f         = cuscomplex( 1.0, 0.0 );
    bool                                is_john             = false;
    cusfloat*                           field_point         = nullptr;
    cusfloat                            normal_vec_aux[3]   = { 0.0, 0.0, 0.0 };
    PanelGeom                           panel_aux;
    PanelData<cuscomplex, RDDConfig>*   panel_rdd;
    cuscomplex                          point_disp[3]       ;
    cuscomplex                          pot_aux             = cuscomplex( 0.0, 0.0 );
    cuscomplex                          pot_term            = cuscomplex( 0.0, 0.0 );
    cuscomplex                          pot_term_st         = cuscomplex( 0.0, 0.0 );
    cuscomplex                          pot_term_wv         = cuscomplex( 0.0, 0.0 );
    cuscomplex                          radius[3]           ;
    cuscomplex                          rao_rot[3]          ;
    cuscomplex                          rao_trans[3]        ;
    cuscomplex                          rao_val             = cuscomplex( 0.0, 0.0 );
    cusfloat                            rhow                = this->_input->water_density;
    SourceNode                          source_aux( &panel_aux, 0, 0,  0, center_aux, normal_vec_aux );
    std::size_t                         sources_np          = this->_solver->num_rows;
    cuscomplex                          source_val          = cuscomplex( 0.0, 0.0 );
    cuscomplex                          vel_aux[3]          = { 0.0, 0.0, 0.0 };
    cuscomplex                          vel_dn_aux          = 0.0;
    cuscomplex                          vel_total[3]        = { 0.0, 0.0, 0.0 };
    cuscomplex                          vel_total_st[3]     = { 0.0, 0.0, 0.0 };
    cuscomplex                          vel_total_wv[3]     = { 0.0, 0.0, 0.0 };
    cusfloat                            wave_amplitude      = this->_input->wave_amplitude;
    cusfloat                            water_depth         = this->_input->water_depth;

    // Calculate wave dependent parameters
    cusfloat    k               =   w2k( 
                                            ang_freq,
                                            water_depth,
                                            grav_acc
                                        );
    // Calculate wave dependent parameters
    cusfloat    nu              = pow2s( ang_freq ) / grav_acc;

    std::vector<cuscomplex> pot_vals( this->_input->dofs_np * body_np + this->_input->heads_np, cuscomplex( 0.0, 0.0 ) );
                                

    // Compute diffraction and radiation fields at the field points
    for ( std::size_t i=0; i<rad_diff_data->get_size_local( ); i++ )
    {
        // Get panel data instance
        panel_rdd   = &(rad_diff_data->panel_data[i]);
        body_id     = panel_rdd->body_id;
        fp_np       = panel_rdd->field_points_np;

        // Clear fields data to avoid spurious values
        panel_rdd->clear_data( );

        // Loop over field points in the current panel
        for ( std::size_t j=0; j<fp_np; j++ )
        {
            // Get current field point
            field_point             = &(panel_rdd->field_points[3*j]);
            panel_aux.center[0]     = field_point[0];
            panel_aux.center[1]     = field_point[1];
            panel_aux.center[2]     = field_point[2];
            source_aux.position[0]  = field_point[0];
            source_aux.position[1]  = field_point[1];
            source_aux.position[2]  = field_point[2];

            // Calculate incident wave field at the field point
            for ( std::size_t idh=0; idh<heads_np; idh++ )
            {
                index_fd = store_freqs * freq_index * ( heads_np * fp_np ) + idh * fp_np + j;
                STATIC_COND( 
                                ONLY_FCN,   
                                panel_rdd->pot_total[index_fd]      = wave_potential_fo_space(
                                                                                                    wave_amplitude,
                                                                                                    ang_freq,
                                                                                                    k,
                                                                                                    water_depth,
                                                                                                    grav_acc,
                                                                                                    field_point[0],
                                                                                                    field_point[1],
                                                                                                    field_point[2],
                                                                                                    this->_input->heads[idh]
                                                                                                );                                      
                            )

                STATIC_COND( 
                                ONLY_FCN,   
                                panel_rdd->pot_total_2[index_fd]    = wave_potential_fo_space(
                                                                                                    wave_amplitude,
                                                                                                    ang_freq,
                                                                                                    k,
                                                                                                    water_depth,
                                                                                                    grav_acc,
                                                                                                    field_point[0],
                                                                                                    field_point[1],
                                                                                                    field_point[2],
                                                                                                    this->_input->heads[idh]
                                                                                                );                                      
                            )

                STATIC_COND( 
                                ONLY_FCNDC,
                                panel_rdd->vel_x_total[index_fd]    = wave_potential_fo_space_dx(
                                                                                                        wave_amplitude,
                                                                                                        ang_freq,
                                                                                                        k,
                                                                                                        water_depth,
                                                                                                        grav_acc,
                                                                                                        field_point[0],
                                                                                                        field_point[1],
                                                                                                        field_point[2],
                                                                                                        this->_input->heads[idh]
                                                                                                    );                                      
                            )
                
                STATIC_COND( 
                                ONLY_FCNDC,
                                panel_rdd->vel_y_total[index_fd]    = wave_potential_fo_space_dy(
                                                                                                        wave_amplitude,
                                                                                                        ang_freq,
                                                                                                        k,
                                                                                                        water_depth,
                                                                                                        grav_acc,
                                                                                                        field_point[0],
                                                                                                        field_point[1],
                                                                                                        field_point[2],
                                                                                                        this->_input->heads[idh]
                                                                                                    );                                      
                            )

                STATIC_COND( 
                                ONLY_FCNDC,
                                panel_rdd->vel_z_total[index_fd]    = wave_potential_fo_space_dz(
                                                                                                        wave_amplitude,
                                                                                                        ang_freq,
                                                                                                        k,
                                                                                                        water_depth,
                                                                                                        grav_acc,
                                                                                                        field_point[0],
                                                                                                        field_point[1],
                                                                                                        field_point[2],
                                                                                                        this->_input->heads[idh]
                                                                                                    );                                      
                            )

                if constexpr( USE_COMP )
                {
                    STATIC_COND( ONLY_FCN,      panel_rdd->pot_incident[index_fd]   = panel_rdd->pot_total[index_fd];   )
                    STATIC_COND( ONLY_FCNDC,    panel_rdd->vel_x_incident[index_fd] = panel_rdd->vel_x_total[index_fd]; )
                    STATIC_COND( ONLY_FCNDC,    panel_rdd->vel_y_incident[index_fd] = panel_rdd->vel_y_total[index_fd]; )
                    STATIC_COND( ONLY_FCNDC,    panel_rdd->vel_z_incident[index_fd] = panel_rdd->vel_z_total[index_fd]; )
                }

                
            }

            for ( int ip=0; ip<this->_input->dofs_np + this->_input->heads_np; ip++ )
            {
                pot_vals[ip] = cuscomplex( 0.0, 0.0 );
            }

            // Calculate raddiation and diffraction fields
            // Loop over columns (source nodes) to calculate influence coefficients for each field point
            for ( std::size_t k=0; k<static_cast<std::size_t>(this->_solver->num_rows); k++ )
            {
                // Check if the field point is located at the source position
                bool is_close = assert_vector_equality( 3, field_point, this->_mesh_gp->source_nodes[k]->position, 1E-3 );

                // Calculate steady coefficients
                _formulation_kernel_steady<
                                                mode_pf,
                                                FreqRegimeE::REGULAR
                                            >
                                            (
                                                is_close,
                                                this->_mesh_gp->panels[k],
                                                this->_mesh_gp->panels_mirror[k],
                                                &panel_aux,
                                                water_depth,
                                                pot_term_st,
                                                int_dn_pf_st,
                                                int_dn_sf_st,
                                                vel_total_st
                                            );

                // Calculate wave coefficients
                _formulation_kernel_wave<
                                            mode_pf,
                                            FreqRegimeE::REGULAR
                                        >
                                        ( 
                                            is_close,
                                            this->_mesh_gp->source_nodes[k],
                                            &(source_aux),
                                            this->_gwfcns_interf,
                                            water_depth,
                                            nu,
                                            is_john,
                                            pot_term_wv,
                                            int_dn_pf_wv,
                                            int_dn_sf_wv,
                                            vel_total_wv
                                        );

                // Add contributions from steady and wave
                if ( this->_mesh_gp->panels[k]->type != PanelTypeE::EXT_LID )
                {
                    int_scale_f = cuscomplex( 1.0, 0.0 );
                }
                else
                {
                    int_scale_f = cuscomplex( 0.0, - this->_mesh_gp->panels[k]->ext_lid_damp_f );
                }

                if ( is_john )
                {
                    STATIC_COND( ONLY_FCN,      pot_term        = int_scale_f * pot_term_wv;        )
                    STATIC_COND( ONLY_FCNDC,    vel_total[0]    = int_scale_f * vel_total_wv[0];    )
                    STATIC_COND( ONLY_FCNDC,    vel_total[1]    = int_scale_f * vel_total_wv[1];    )
                    STATIC_COND( ONLY_FCNDC,    vel_total[2]    = int_scale_f * vel_total_wv[2];    )
                    STATIC_COND( ONLY_FCNDN,    int_dn_sf       = int_scale_f * int_dn_sf_wv;       )
                }
                else
                {
                    STATIC_COND( ONLY_FCN,      pot_term        = int_scale_f * ( pot_term_st     + pot_term_wv );        )
                    STATIC_COND( ONLY_FCNDC,    vel_total[0]    = int_scale_f * ( vel_total_st[0] + vel_total_wv[0] );    )
                    STATIC_COND( ONLY_FCNDC,    vel_total[1]    = int_scale_f * ( vel_total_st[1] + vel_total_wv[1] );    )
                    STATIC_COND( ONLY_FCNDC,    vel_total[2]    = int_scale_f * ( vel_total_st[2] + vel_total_wv[2] );    )
                    STATIC_COND( ONLY_FCNDN,    int_dn_sf       = int_scale_f * ( int_dn_sf_st    + int_dn_sf_wv );       )
                }

                // Loop over headings to add radiation and diffraction contribution
                for ( std::size_t idh=0; idh<heads_np; idh++ )
                {
                    // Get indexes to locate and to storage data
                    // Diffraction column: dofs_np * bodies_np + idh
                    index_sc    = ( dofs_np * body_np + idh ) * sources_np + k;
                    index_fd    = store_freqs * freq_index * ( heads_np * fp_np ) + idh * fp_np + j;

                    // Get source value for the current position
                    source_val  = this->_sf_gp->field_values[index_sc];

                    pot_vals[this->_input->dofs_np * body_np + idh] += pot_term * source_val;

                    // Calculate diffraction field contribution
                    STATIC_COND( ONLY_FCN,   pot_aux    = pot_term     * source_val;   )
                    STATIC_COND( ONLY_FCNDN, vel_dn_aux = int_dn_sf    * source_val;   )
                    STATIC_COND( ONLY_FCNDC, vel_aux[0] = vel_total[0] * source_val;   )
                    STATIC_COND( ONLY_FCNDC, vel_aux[1] = vel_total[1] * source_val;   )
                    STATIC_COND( ONLY_FCNDC, vel_aux[2] = vel_total[2] * source_val;   )

                    STATIC_COND( ONLY_FCN,   panel_rdd->pot_total[index_fd]     += pot_aux;      )
                    STATIC_COND( ONLY_FCNDN, panel_rdd->vel_dn_total[index_fd]  += vel_dn_aux;   )
                    STATIC_COND( ONLY_FCNDC, panel_rdd->vel_x_total[index_fd]   += vel_aux[0];   )
                    STATIC_COND( ONLY_FCNDC, panel_rdd->vel_y_total[index_fd]   += vel_aux[1];   )
                    STATIC_COND( ONLY_FCNDC, panel_rdd->vel_z_total[index_fd]   += vel_aux[2];   )

                    if constexpr( USE_COMP )
                    {
                        STATIC_COND( ONLY_FCN,   panel_rdd->pot_diff[index_fd]      += pot_aux;      )
                        STATIC_COND( ONLY_FCNDN, panel_rdd->vel_dn_diff[index_fd]   += vel_dn_aux;   )
                        STATIC_COND( ONLY_FCNDC, panel_rdd->vel_x_diff[index_fd]    += vel_aux[0];   )
                        STATIC_COND( ONLY_FCNDC, panel_rdd->vel_y_diff[index_fd]    += vel_aux[1];   )
                        STATIC_COND( ONLY_FCNDC, panel_rdd->vel_z_diff[index_fd]    += vel_aux[2];   )
                    }

                    for ( std::size_t idd=0; idd<dofs_np; idd++ )
                    {
                        // Get indexes to locate and to storage data
                        index_ax    = idh * ( dofs_np * body_np ) + body_id * dofs_np + idd;
                        // Radiation column: body_id * dofs_np + idd
                        index_sc    = ( body_id * dofs_np + idd ) * sources_np + k;
                        
                        // Get source value for the current position
                        source_val  = this->_sf_gp->field_values[index_sc];

                        pot_vals[body_id * dofs_np + idd] += pot_term * source_val;

                        // Get current RAO value
                        rao_val     = raos[ index_ax ];

                        // Calculate field contributions
                        STATIC_COND( ONLY_FCN,   pot_aux    = cuscomplex(0.0, -ang_freq) * pot_term     * rao_val * source_val;   )
                        STATIC_COND( ONLY_FCNDN, vel_dn_aux = cuscomplex(0.0, -ang_freq) * int_dn_sf    * rao_val * source_val;   )
                        STATIC_COND( ONLY_FCNDC, vel_aux[0] = cuscomplex(0.0, -ang_freq) * vel_total[0] * rao_val * source_val;   )
                        STATIC_COND( ONLY_FCNDC, vel_aux[1] = cuscomplex(0.0, -ang_freq) * vel_total[1] * rao_val * source_val;   )
                        STATIC_COND( ONLY_FCNDC, vel_aux[2] = cuscomplex(0.0, -ang_freq) * vel_total[2] * rao_val * source_val;   )

                        STATIC_COND( ONLY_FCN,   panel_rdd->pot_total[index_fd]     += pot_aux;   )
                        STATIC_COND( ONLY_FCNDN, panel_rdd->vel_dn_total[index_fd]  += vel_dn_aux;   )
                        STATIC_COND( ONLY_FCNDC, panel_rdd->vel_x_total[index_fd]   += vel_aux[0];   )
                        STATIC_COND( ONLY_FCNDC, panel_rdd->vel_y_total[index_fd]   += vel_aux[1];   )
                        STATIC_COND( ONLY_FCNDC, panel_rdd->vel_z_total[index_fd]   += vel_aux[2];   )

                        if constexpr( USE_COMP )
                        {
                            STATIC_COND( ONLY_FCN,   panel_rdd->pot_rad[index_fd]    += pot_aux;   )
                            STATIC_COND( ONLY_FCNDN, panel_rdd->vel_dn_rad[index_fd] += vel_dn_aux;   )
                            STATIC_COND( ONLY_FCNDC, panel_rdd->vel_x_rad[index_fd]  += vel_aux[0];   )
                            STATIC_COND( ONLY_FCNDC, panel_rdd->vel_y_rad[index_fd]  += vel_aux[1];   )
                            STATIC_COND( ONLY_FCNDC, panel_rdd->vel_z_rad[index_fd]  += vel_aux[2];   )
                        }

                    }
                }
            }
        
            // Calculate total pressure field
            for ( std::size_t idh=0; idh<heads_np; idh++ )
            {
                index_fd = store_freqs * freq_index * ( heads_np * fp_np ) + idh * fp_np + j;

                STATIC_COND(  
                                ONLY_FCN, 
                                panel_rdd->press_total[index_fd] = cuscomplex( 0.0, rhow * ang_freq ) * panel_rdd->pot_total[index_fd];
                            )

                if constexpr( USE_COMP )
                {
                    panel_rdd->press_incident[index_fd] = cuscomplex( 0.0, rhow * ang_freq ) * panel_rdd->pot_incident[index_fd];
                    panel_rdd->press_diff[index_fd]     = cuscomplex( 0.0, rhow * ang_freq ) * panel_rdd->pot_diff[index_fd];
                    panel_rdd->press_rad[index_fd]      = cuscomplex( 0.0, rhow * ang_freq ) * panel_rdd->pot_rad[index_fd];
                }
            }

            // Calculate wave elevation field
            for ( std::size_t idh=0; idh<heads_np; idh++ )
            {
                index_fd = store_freqs * freq_index * ( heads_np * fp_np ) + idh * fp_np + j;

                STATIC_COND(  
                                ONLY_FCN, 
                                panel_rdd->wev_total[index_fd]  = panel_rdd->pot_total[index_fd] * cuscomplex( 0.0, ang_freq ) / grav_acc;
                            )

                if constexpr( ONLY_FCN )
                {
                    index_ax = idh * ( dofs_np * body_np ) + body_id * dofs_np;

                    for ( int r=0; r<3; r++ )
                    {
                        radius[r]       = cuscomplex( 
                                                        panel_rdd->field_points[3*j+r] 
                                                        - 
                                                        panel_rdd->panel_geom->body_cog[r], 
                                                        0.0 
                                                    );
                        
                        rao_trans[r]    = raos[index_ax+r]   * wave_amplitude;
                        rao_rot[r]      = raos[index_ax+3+r] * wave_amplitude;
                    }

                    cross( 
                            rao_rot,
                            radius,
                            point_disp
                    );

                    sv_add(
                                3,
                                point_disp,
                                rao_trans,
                                point_disp
                            );

                    // Calculate relative water height by substracting the 
                    // vertical motion at the required point to the wave elevation
                    panel_rdd->wev_rel_total[index_fd] = ( 
                                                                panel_rdd->wev_total[index_fd] 
                                                                - 
                                                                point_disp[2]
                                                            );

                }

            }

        }

    }

    // #include <fstream>
    // std::ofstream eval_points_file( "S:/seamotions_validation/0_seamotions/1_H50/0_Box/eval_points.txt" );
    // std::ofstream pot_total_file( "S:/seamotions_validation/0_seamotions/1_H50/0_Box/pot_total_file.txt" );
    // std::ofstream pot_total_2_file( "S:/seamotions_validation/0_seamotions/1_H50/0_Box/pot_total_2_file.txt" );
    // std::ofstream vel_x_total_file( "S:/seamotions_validation/0_seamotions/1_H50/0_Box/vel_x_total.txt" );
    // std::ofstream vel_y_total_file( "S:/seamotions_validation/0_seamotions/1_H50/0_Box/vel_y_total.txt" );
    // std::ofstream vel_z_total_file( "S:/seamotions_validation/0_seamotions/1_H50/0_Box/vel_z_total.txt" );
    // std::ofstream wev_total_file( "S:/seamotions_validation/0_seamotions/1_H50/0_Box/wev_total.txt" );
    // std::ofstream rwev_total_file( "S:/seamotions_validation/0_seamotions/1_H50/0_Box/rwev_total.txt" );

    // std::size_t count_fid = 0;
    // for ( std::size_t i=0; i<rad_diff_data->get_size_local( ); i++ )
    // {
    //     fp_np = rad_diff_data->panel_data[i].field_points_np;
    //     for ( std::size_t j=0; j<fp_np; j++ )
    //     {
    //         eval_points_file << count_fid << " ";
    //         eval_points_file << rad_diff_data->panel_data[i].field_points[3*j+0] << " ";
    //         eval_points_file << rad_diff_data->panel_data[i].field_points[3*j+1] << " ";
    //         eval_points_file << rad_diff_data->panel_data[i].field_points[3*j+2] << "\n";
    //         count_fid++;
    //     }
    // }

    // for ( std::size_t idh=0; idh<heads_np; idh++ )
    // {
    //     for ( std::size_t i=0; i<rad_diff_data->get_size_local( ); i++ )
    //     {
    //         fp_np = rad_diff_data->panel_data[i].field_points_np;
    //         for ( std::size_t j=0; j<fp_np; j++ )
    //         {
    //             index_fd = freq_index * ( heads_np * fp_np ) + idh * fp_np + j;
    //             pot_total_file << rad_diff_data->panel_data[i].pot_total[index_fd].real( ) << " " << rad_diff_data->panel_data[i].pot_total[index_fd].imag( ) << "\n";
    //             pot_total_2_file << rad_diff_data->panel_data[i].pot_total_2[index_fd].real( ) << " " << rad_diff_data->panel_data[i].pot_total_2[index_fd].imag( ) << "\n";
    //             vel_x_total_file << rad_diff_data->panel_data[i].vel_x_total[index_fd].real( ) << " " << rad_diff_data->panel_data[i].vel_x_total[index_fd].imag( ) << "\n";
    //             vel_y_total_file << rad_diff_data->panel_data[i].vel_y_total[index_fd].real( ) << " " << rad_diff_data->panel_data[i].vel_y_total[index_fd].imag( ) << "\n";
    //             vel_z_total_file << rad_diff_data->panel_data[i].vel_z_total[index_fd].real( ) << " " << rad_diff_data->panel_data[i].vel_z_total[index_fd].imag( ) << "\n";
    //             wev_total_file << rad_diff_data->panel_data[i].wev_total[index_fd].real( ) << " " << rad_diff_data->panel_data[i].wev_total[index_fd].imag( ) << "\n";
    //             rwev_total_file << rad_diff_data->panel_data[i].wev_rel_total[index_fd].real( ) << " " << rad_diff_data->panel_data[i].wev_rel_total[index_fd].imag( ) << "\n";
    //         }
    //     }
    // }

    // pot_total_file.close( );
    // pot_total_2_file.close( );
    // vel_x_total_file.close( );
    // vel_y_total_file.close( );
    // vel_z_total_file.close( );
    // wev_total_file.close( );
    // rwev_total_file.close( );

    MPI_Barrier( MPI_COMM_WORLD );

}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::configure_second_order( 
                                                                        SimulationData* sim_data
                                                                    )
{
    this->_qtf_diff_code_integrator = FarFieldIntegrator<QTFTypeE::QTF_DIFF_CODE>( 
                                                                                        this->_input,
                                                                                        0,
                                                                                        1,
                                                                                        sim_data->qtf_pc_radius,
                                                                                        sim_data->qtf_a_cos,
                                                                                        sim_data->qtf_a_sin,
                                                                                        sim_data->qtf_b_cos,
                                                                                        sim_data->qtf_b_sin,
                                                                                        QTF_FAR_N
                                                                                    );

    this->_qtf_sum_code_integrator  = FarFieldIntegrator<QTFTypeE::QTF_SUM_CODE>( 
                                                                                        this->_input,
                                                                                        0,
                                                                                        1,
                                                                                        sim_data->qtf_pc_radius,
                                                                                        sim_data->qtf_a_cos,
                                                                                        sim_data->qtf_a_sin,
                                                                                        sim_data->qtf_b_cos,
                                                                                        sim_data->qtf_b_sin,
                                                                                        QTF_FAR_N
                                                                                    );
}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::_initialize( 
                                                            void 
                                                        )
{
    if constexpr ( recalc_steady == RecalcSteadyE::OFF )
    {
        //Calculate steady part integral over the panels
        MPI_TIME_EXEC( this->_build_steady_matrixes<FreqRegimeE::REGULAR>( ); , this->exec_time_build_steady )
        this->_steady_mat_type  = 0;
    }
}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
FormulationKernelBackend<N, mode_pf, recalc_steady>::FormulationKernelBackend(
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
                                            this->_input->dofs_np * this->_mesh_gp->meshes_np + this->_input->heads_np,
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
                                            ( this->_input->dofs_np * this->_mesh_gp->meshes_np + this->_input->heads_np ),
                                            0,
                                            this->_mesh_gp->panels_tnp-1,
                                            this->ipm_sc,
                                            this->ipm_ed,
                                            true,
                                            recalc_steady
                                        );

    STATIC_COND( 
                    !(ONLY_PF),
                    this->_pot_gp       = new   MLGCmpx(
                                                            this->_mesh_gp->panels_tnp,
                                                            this->ipm_cols_np,
                                                            this->_mesh_gp->meshes_np,
                                                            ( this->_input->dofs_np * this->_mesh_gp->meshes_np + this->_input->heads_np ),
                                                            0,
                                                            this->_mesh_gp->panels_tnp-1,
                                                            this->ipm_sc,
                                                            this->ipm_ed,
                                                            true,
                                                            recalc_steady
                                                        );
                )

    std::size_t  fo_len = static_cast<std::size_t>( this->_input->dofs_np * this->_mesh_gp->meshes_np ) + this->_input->heads_np;
    std::size_t  so_len = this->_input->heads_np * this->_input->heads_np;
    std::size_t  pf_len = ( ( mode_pf > 0 ) && ( so_len > fo_len ) ) ? so_len : fo_len;
    STATIC_COND( 
                    ONLY_PF,
                    this->_pf_gp        = new   MLGCmpx(
                                                            this->_mesh_gp->panels_tnp,
                                                            this->ipm_cols_np,
                                                            this->_mesh_gp->meshes_np,
                                                            pf_len,
                                                            0,
                                                            this->_mesh_gp->panels_tnp-1,
                                                            this->ipm_sc,
                                                            this->ipm_ed,
                                                            true,
                                                            recalc_steady
                                                        );
                )

    // Allocate space for the partial potential formulation RHS vector
    STATIC_COND( ONLY_PF, this->_ppf_rhs  = generate_empty_vector<cuscomplex>( this->_pf_gp->sysmat_nrows * pf_len ); )

    // Allocate space for the wave part integration interface functor
    this->_gwfcns_interf.initialize(
                                        this->_input->angfreqs[0],
                                        this->_input->water_depth,
                                        this->_input->grav_acc,
                                        false
                                    );

    if ( this->_input->out_qtf_so_model == QTFSOModelE::DIRECT )
    {
        this->_gwfcns_interf_qb.initialize(
                                                this->_input->angfreqs[0],
                                                this->_input->water_depth,
                                                this->_input->grav_acc,
                                                true
                                            );

        this->_gwfcns_interf_qf.initialize(
                                                this->_input->angfreqs[0],
                                                this->_input->water_depth,
                                                this->_input->grav_acc,
                                                true
                                            );

        this->_gwfcns_interf_qfa.initialize(
                                                this->_input->angfreqs[0],
                                                this->_input->water_depth,
                                                this->_input->grav_acc,
                                                true
                                            );

        this->_gwfcns_interf_wl.initialize(
                                                this->_input->angfreqs[0],
                                                this->_input->water_depth,
                                                this->_input->grav_acc,
                                                true
                                            );
    }

    // Create second order wave properties object
    if ( this->_input->out_qtf_so_model == QTFSOModelE::DIRECT )
    {
        this->_wdso = WaveDispersionSO(
                                            this->_input->wave_amplitude,
                                            this->_input->wave_amplitude,
                                            this->_input->angfreqs[0],
                                            this->_input->angfreqs[0],
                                            this->_input->heads[0],
                                            this->_input->heads[0],
                                            this->_input->water_depth,
                                            this->_input->grav_acc
                                        );
    }
    
    // Launch object initialization
    this->_initialize( );
}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
FormulationKernelBackend<N, mode_pf, recalc_steady>::~FormulationKernelBackend( )
{
    if ( this->_input->out_qtf_so_model == QTFSOModelE::DIRECT )
    {
        delete this->_wdso;
    }

    delete this->_solver;
    delete this->_sf_gp;
    delete this->_pot_gp;

    STATIC_COND( ONLY_PF, delete this->_pf_gp; );
}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
int FormulationKernelBackend<N, mode_pf, recalc_steady>::size( void )
{
    return this->_solver->num_rows;
}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
template<FreqRegimeE freq_regime>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::solve( cusfloat w )
{
    // Recalculate steady part of the system matrixes if required
    if constexpr ( recalc_steady == RecalcSteadyE::OFF )
    {
        if ( 
                this->_steady_mat_type != 0 
                &&
                (   
                    freq_regime == FreqRegimeE::REGULAR 
                    || 
                    freq_regime == FreqRegimeE::ASYMPT_LOW 
                )
            )
        {
            MPI_TIME_EXEC( this->_build_steady_matrixes<freq_regime>( ); , this->exec_time_build_steady )
            this->_steady_mat_type  = 0;
        }
        else if ( 
                    this->_steady_mat_type != 1 
                    &&
                    freq_regime == FreqRegimeE::ASYMPT_HIGH 
                )
        {
            MPI_TIME_EXEC( this->_build_steady_matrixes<freq_regime>( ); , this->exec_time_build_steady )
            this->_steady_mat_type  = 1;
        }
    }
    
    // Fold integrals database if required
    if constexpr( freq_regime == FreqRegimeE::REGULAR )
    {
        // Fold database for current frequency and water depth
        cusfloat H = pow2s( w ) * this->_input->water_depth / this->_input->grav_acc;

        fold_database( H );

        // Update integration functor to ne new frequency
        this->_gwfcns_interf.set_ang_freq( w );
    }
    
    // Re-calculate wave dependent system matrix term
    MPI_TIME_EXEC(  this->_build_wave_matrixes_2<freq_regime>( w );, this->exec_time_build_wave )
                    // this->_build_rhs( w );

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
        int         dofs_offset = this->_input->dofs_np * this->_pf_gp->sysmat_nrows;
        int         index       = 0;
        PanelGeom*  panel_j     = nullptr;
        double      nu          = pow2s( w ) / this->_input->grav_acc;

        for ( int i=0; i<this->_input->dofs_np; i++ )
        {
            for ( int j=this->_solver->start_row_0; j<this->_solver->end_row_0; j++ )
            {
                panel_j = this->_mesh_gp->source_nodes[j]->panel;
                if ( panel_j->type == PanelTypeE::INT_LID )
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

                if ( panel_j->type == PanelTypeE::INT_LID )
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


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
template<QTFTypeE qtf_type>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::solve_so( 
                                                        std::size_t                 freq_i,
                                                        std::size_t                 freq_j,
                                                        cuscomplex*                 raos_hist,
                                                        RDDQC*                      qtf_body_fields,
                                                        RDDQC*                      qtf_body_fields_wl,
                                                        RDDQC*                      qtf_fs_fields,
                                                        RDDQC*                      qtf_fs_fields_wl,
                                                        const std::vector<RDDQC*>*  qtf_annuli_fields
                                                    )
{
    // Compute wave frequency depending on the QTF type
    static_assert( ( qtf_type == QTFTypeE::QTF_DIFF_CODE ) || ( qtf_type == QTFTypeE::QTF_SUM_CODE ), "Invalid QTF type" );
    cusfloat w = 0.0;
    if constexpr( qtf_type == QTFTypeE::QTF_DIFF_CODE )
    {
        w = std::abs( this->_input->angfreqs[freq_i] - this->_input->angfreqs[freq_j] );
    }
    else
    {
        w = this->_input->angfreqs[freq_i] + this->_input->angfreqs[freq_j];
    }

    if constexpr ( recalc_steady == RecalcSteadyE::OFF )
    {
        // Recalculate steady part of the system matrixes if required
        if ( 
                this->_steady_mat_type != 0 
            )
        {
            MPI_TIME_EXEC( this->_build_steady_matrixes<FreqRegimeE::REGULAR>( ); , this->exec_time_build_steady )
            this->_steady_mat_type  = 0;
        }
    }

    // Update frequency index for far-field integrators
    if constexpr( qtf_type == QTFTypeE::QTF_DIFF_CODE )
    {
        this->_qtf_diff_code_integrator.set_frequency_indices( freq_i, freq_j );
    }
    else
    {
        this->_qtf_sum_code_integrator.set_frequency_indices( freq_i, freq_j );
    }

    // Fold database for current frequency and water depth
    cusfloat H = pow2s( w ) * this->_input->water_depth / this->_input->grav_acc;

    fold_database( H );

    // Update integration functor to ne new frequency
    this->_gwfcns_interf_qb.set_ang_freq( w );
    this->_gwfcns_interf_qf.set_ang_freq( w );
    this->_gwfcns_interf_qfa.set_ang_freq( w );
    this->_gwfcns_interf_wl.set_ang_freq( w );
    
    // Re-calculate wave dependent system matrix term
    MPI_TIME_EXEC(  this->_build_wave_matrixes_2<FreqRegimeE::REGULAR>( w );, this->exec_time_build_wave )
                    this->_build_rhs_so( 
                                            qtf_type,
                                            freq_i,
                                            freq_j,
                                            raos_hist,
                                            qtf_body_fields,
                                            qtf_body_fields_wl,
                                            qtf_fs_fields,
                                            qtf_fs_fields_wl,
                                            qtf_annuli_fields
                                        );

    // Calculate system matrixes condition number if required
    if ( this->_is_condition_number )
    {
        cusfloat pf_cond = 0.0;
    
        this->_solver->Cond( this->_pf_gp->sysmat,    pf_cond     );
    
        if ( this->_mpi_config->is_root( ) )
        {
            std::cout << "Condition Number -> PF: " << pf_cond << std::endl;
        }
    }

    // Calculate potential through the potential formulation
    MPI_TIME_EXEC( 
                    ( this->_solver->Solve( this->_pf_gp->sysmat, this->_pf_gp->field_values ) );, 
                    this->exec_time_solve_pf 
                ) 
    MPI_Bcast(
                    this->_pf_gp->field_values,
                    this->_pf_gp->field_values_np,
                    mpi_cuscomplex,
                    this->_mpi_config->proc_root,
                    MPI_COMM_WORLD                
                );
    
}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::_process_far_field(
                                                                QTFTypeE        qtf_type,
                                                                std::size_t     freq_i,
                                                                std::size_t     freq_j,
                                                                cuscomplex*     qb_rhs
                                                            )
{
    // Calculate far field contributions for the required QTF type
    if ( qtf_type == QTFTypeE::QTF_DIFF_CODE )
    {
        this->_qtf_diff_code_integrator.compute_far_field( freq_i, freq_j );
    }
    else if ( qtf_type == QTFTypeE::QTF_SUM_CODE )
    {
        this->_qtf_sum_code_integrator.compute_far_field( freq_i, freq_j );
    }

    // Add far field contributions to the RHS vector
    const std::size_t heads_np = static_cast<std::size_t>( this->_input->heads_np );
    for ( std::size_t ih0=0; ih0<static_cast<std::size_t>(this->_input->heads_np); ih0++ )
    {
        for ( std::size_t ih1=0; ih1<static_cast<std::size_t>(heads_np); ih1++ )
        {
            std::size_t local_index     = ih0 * heads_np + ih1;

            for ( std::size_t i=0; i<static_cast<std::size_t>(this->_mesh_gp->panels_tnp); i++ )
            {
                std::size_t global_index    = ( 
                                                    ih0 * heads_np * this->_mesh_gp->panels_tnp
                                                    +
                                                    ih1 * this->_mesh_gp->panels_tnp
                                                );
                
                
                
                if ( qtf_type == QTFTypeE::QTF_DIFF_CODE )
                {
                    qb_rhs[global_index]        += ( 
                                                        this->_qtf_diff_code_integrator.qc[local_index] 
                                                        + 
                                                        cuscomplex( 0.0, 1.0 ) * this->_qtf_diff_code_integrator.qs[local_index]
                                                    );
                }
                else if ( qtf_type == QTFTypeE::QTF_SUM_CODE )
                {
                    qb_rhs[global_index]        += ( 
                                                        this->_qtf_sum_code_integrator.qc[local_index] 
                                                        + 
                                                        cuscomplex( 0.0, 1.0 ) * this->_qtf_sum_code_integrator.qs[local_index]
                                                    );
                }

            }
        }
    }
}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::set_qtf_annuli_fields(
                                                                    const std::vector<RDDQC*>* annuli_fields
                                                                )
{
    this->_qtf_annuli_fields = annuli_fields;
}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::append_memory_report(
                                                                    std::vector<MemoryReportEntry>& entries,
                                                                    const std::string& prefix
                                                                ) const
{
    if ( this->_sf_gp != nullptr )
    {
        this->_sf_gp->append_memory_report( entries, memory_report_path( prefix, "sf_gp" ) );
    }
    if ( this->_pot_gp != nullptr )
    {
        this->_pot_gp->append_memory_report( entries, memory_report_path( prefix, "pot_gp" ) );
    }
    if ( this->_pf_gp != nullptr )
    {
        this->_pf_gp->append_memory_report( entries, memory_report_path( prefix, "pf_gp" ) );
    }
    if ( this->_ppf_rhs != nullptr && this->_pot_gp != nullptr && this->_input != nullptr )
    {
        std::size_t rhs_len = static_cast<std::size_t>( this->_pot_gp->sysmat_nrows )
                                * static_cast<std::size_t>( this->_input->dofs_np + this->_input->heads_np );
        add_memory_entry(
                            entries,
                            memory_report_path( prefix, "ppf_rhs" ),
                            rhs_len * sizeof( cuscomplex )
                        );
    }
    if ( this->_solver != nullptr )
    {
        this->_solver->append_memory_report( entries, memory_report_path( prefix, "scalapack_solver" ) );
    }
}


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::update_results( SimulationData* sim_data )
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
                                    sim_data->fields_np, 
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


template<std::size_t N, int mode_pf, RecalcSteadyE recalc_steady>
template<QTFTypeE qtf_type>
void FormulationKernelBackend<N, mode_pf, recalc_steady>::update_results_so( SimulationData* sim_data )
{
    // Update data to simulation results
    copy_vector( 
                    sim_data->fields_so_np, 
                    this->_pf_gp->field_values,
                    sim_data->panels_potential_so
                );
    
}