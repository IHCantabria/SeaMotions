
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
#include "gf_intensities.hpp"
#include "../../math/integration.hpp"
#include "../../green/source.hpp"
#include "../../waves/wave_dispersion_base_fo.hpp"
#include "../../waves/waves_common.hpp"


void    calculate_gf_intensity_sysmat(
                                                    Input*          input,
                                                    SclCmpx*        scl,
                                                    MeshGroup*      mesh_gp,
                                                    GWFDnInterface* gwf_interf,
                                                    cusfloat        w,
                                                    cuscomplex*     sysmat_steady,
                                                    cuscomplex*     sysmat,
                                                    cuscomplex*     sources_int
                                   )
{
    /***************************************/
    /******** Fill system matrix  **********/
    /***************************************/
    // Create auxiliar lambda function
    auto wave_fcn   =   [gwf_interf]
                        ( 
                            cusfloat xi, 
                            cusfloat eta,
                            cusfloat x, 
                            cusfloat y, 
                            cusfloat z 
                        )
                        {
                            return (*gwf_interf)( xi, eta, x, y, z );
                        };
    
    // Loop over panels to integrate value
    int         col_count   = 0;
    int         index       = 0;
    cuscomplex  int_value( 0.0, 0.0 );
    cuscomplex  wave_value( 0.0, 0.0 );
    PanelGeom*  panel_j     = nullptr;
    SourceNode* source_i    = nullptr;
    int         row_count   = 0;

    for ( int i=scl->start_col_0; i<scl->end_col_0; i++ )
    {
        // Get memory address of the ith panel
        source_i = mesh_gp->source_nodes[i];
        gwf_interf->set_source_i( source_i, 1.0 );

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
        {
            // Get memory address of the panel jth
            panel_j = mesh_gp->source_nodes[j]->panel;
            gwf_interf->set_source_j( mesh_gp->source_nodes[j] );

            // Integrate green function normal derivative along the current panel
            if ( i == j )
            {
                if ( panel_j->type == DIFFRAC_PANEL_CODE )
                {
                    int_value   =   -cuscomplex( 0.5, 0.0 );
                }
                else if ( panel_j->type == LID_PANEL_CODE )
                {
                    int_value   =   -cuscomplex( 4.0 * PI, 0.0 );
                }
            }
            else
            {
                wave_value      =   adaptive_quadrature_panel(
                                                                source_i->panel,
                                                                wave_fcn,
                                                                input->gfdn_abs_err,
                                                                input->gfdn_rel_err,
                                                                input->is_block_adaption,
                                                                true,
                                                                input->gauss_order
                                                            );
                
                int_value       =   wave_value / 4.0 / PI;
                if ( 
                        panel_j->type == LID_PANEL_CODE
                        &&
                        source_i->panel->type == LID_PANEL_CODE
                    )
                {
                    int_value   = - int_value;
                }
            }
            index           = col_count*scl->num_rows_local+row_count;
            sysmat[index]   = sysmat_steady[index] + int_value;

            // Advance row count
            row_count++;
        }

        // Advance column count
        col_count++;
    }

    /***************************************/
    /***** Fill Hydromechanics RHS  ********/
    /***************************************/

    // Declare local variables to be used
    int         count           = 0;

    // Fill RHS vector
    count       = 0;
    for ( int i=0; i<input->dofs_np; i++ )
    {
        for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
        {
            panel_j = mesh_gp->source_nodes[j]->panel;
            if ( panel_j->type == DIFFRAC_PANEL_CODE )
            {
                sources_int[count] = ( 
                                            mesh_gp->source_nodes[j]->normal_vec[i]
                                            *
                                            mesh_gp->source_nodes[j]->panel->is_move_f
                                        );
            }
            else if ( panel_j->type == LID_PANEL_CODE )
            {
                sources_int[count] = cuscomplex( 0.0, 0.0 );
            }
            count++;
        }
    }

    /***************************************/
    /****** Fill Wave Exciting RHS  ********/
    /***************************************/
    // Define local variables to manage array indexes
                count       = input->dofs_np * mesh_gp->source_nodes_tnp;
    cusfloat    k           = w2k( w, input->water_depth, input->grav_acc );
    cuscomplex  wave_dx     = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_dy     = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_dz     = cuscomplex( 0.0, 0.0 );

    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
        {
            panel_j = mesh_gp->source_nodes[j]->panel;
            if ( panel_j->type == DIFFRAC_PANEL_CODE )
            {
                // Get wave potential derivatives for the panel
                wave_dx             =   wave_potential_fo_space_dx(
                                                                        1.0,
                                                                        w,
                                                                        k,
                                                                        input->water_depth,
                                                                        input->grav_acc,
                                                                        mesh_gp->source_nodes[j]->panel->center[0],
                                                                        mesh_gp->source_nodes[j]->panel->center[1],
                                                                        mesh_gp->source_nodes[j]->panel->center[2],
                                                                        input->heads[i]
                                                                    );

                wave_dy             =   wave_potential_fo_space_dy(
                                                                        1.0,
                                                                        w,
                                                                        k,
                                                                        input->water_depth,
                                                                        input->grav_acc,
                                                                        mesh_gp->source_nodes[j]->panel->center[0],
                                                                        mesh_gp->source_nodes[j]->panel->center[1],
                                                                        mesh_gp->source_nodes[j]->panel->center[2],
                                                                        input->heads[i]
                                                                    );

                wave_dz             =   wave_potential_fo_space_dz(
                                                                        1.0,
                                                                        w,
                                                                        k,
                                                                        input->water_depth,
                                                                        input->grav_acc,
                                                                        mesh_gp->source_nodes[j]->panel->center[0],
                                                                        mesh_gp->source_nodes[j]->panel->center[1],
                                                                        mesh_gp->source_nodes[j]->panel->center[2],
                                                                        input->heads[i]
                                                                    );
                
                // Calculate normal derivative of the wave flow velocities for the jth panel
                sources_int[count]  = -(
                                            wave_dx * mesh_gp->source_nodes[j]->normal_vec[0]
                                            +
                                            wave_dy * mesh_gp->source_nodes[j]->normal_vec[1]
                                            +
                                            wave_dz * mesh_gp->source_nodes[j]->normal_vec[2]
                                        );
            }
            else if ( panel_j->type == LID_PANEL_CODE )
            {
                sources_int[count]  = cuscomplex( 0.0, 0.0 );
            }
            
            count++;
        }
    }

    // Solve system of equations
    scl->Solve( sysmat, sources_int );
}


void    calculate_gf_intensity_steady_sysmat_lin(
                                                    Input*          input,
                                                    SclCmpx*        scl,
                                                    MeshGroup*      mesh_gp,
                                                    cuscomplex*     sysmat
                                                )
{
    /***************************************/
    /******** Fill system matrix  **********/
    /***************************************/
    // Loop over panels to integrate value
    int         col_count   = 0;
    cuscomplex  int_value( 0.0, 0.0 );
    int         row_count   = 0;

    // Define local variables to work with the fast solver
    const int   ndim                    = 3;
    cusfloat    field_point_i[ndim];    clear_vector( ndim, field_point_i );
    PanelGeom*  panel_i                 = nullptr;
    cusfloat    vel_0[ndim];            clear_vector( ndim, vel_0 );
    cusfloat    vel_1[ndim];            clear_vector( ndim, vel_1 );
    cusfloat    vel_2[ndim];            clear_vector( ndim, vel_2 );
    cusfloat    vel_3[ndim];            clear_vector( ndim, vel_3 );
    cusfloat    vel_4[ndim];            clear_vector( ndim, vel_4 );
    cusfloat    vel_5[ndim];            clear_vector( ndim, vel_5 );
    cusfloat    vel_total[ndim];        clear_vector( ndim, vel_total );

    // Define field points to calculate the source influence matrix
    int         field_points_np   = mesh_gp->panels_tnp;
    cusfloat*   field_points      = generate_empty_vector<cusfloat>( 3 * field_points_np );

    for ( int i=0; i<mesh_gp->panels_tnp; i++ )
    {
        copy_vector( 3, mesh_gp->panels[i]->center, &(field_points[3*i]) );
    }

    // Loop over panels and field points to create the steady source matrix
    for ( int i=scl->start_col_0; i<scl->end_col_0; i++ )
    {
        // Get pointer to ith panel
        panel_i = mesh_gp->panels[i];

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=0; j<field_points_np; j++ )
        {
            // Integrate green function normal derivative along the current panel
            if ( i == j )
            {
                int_value       = cuscomplex( 0.0, 0.0 );
            }
            else
            {
                // Reset velocity values
                clear_vector( ndim, vel_0 );
                clear_vector( ndim, vel_1 );
                clear_vector( ndim, vel_2 );
                clear_vector( ndim, vel_3 );
                clear_vector( ndim, vel_4 );
                clear_vector( ndim, vel_5 );
                clear_vector( ndim, vel_total );
                
                // Calcualte velocity corresponding to the r0 source
                calculate_source_velocity_newman(
                                                    mesh_gp->panels[i],
                                                    &(field_points[3*j]), 
                                                    0,
                                                    0, 
                                                    vel_0
                                                );

                // Calculate velocity corresponding to the r1 source
                field_point_i[0]    =   field_points[3*j];
                field_point_i[1]    =   field_points[3*j+1];
                field_point_i[2]    =   field_points[3*j+2] + 2 * input->water_depth;
                calculate_source_velocity_newman(
                                                    mesh_gp->panels_mirror[i],
                                                    field_point_i, 
                                                    0,
                                                    0, 
                                                    vel_1
                                                );
                
                // Calculate velocity corresponding to the r2 source
                field_point_i[0]    =   field_points[3*j];
                field_point_i[1]    =   field_points[3*j+1];
                field_point_i[2]    =   field_points[3*j+2];
                calculate_source_velocity_newman(
                                                    mesh_gp->panels_mirror[i],
                                                    field_point_i, 
                                                    0,
                                                    0, 
                                                    vel_2
                                                );

                // Calculate velocity corresponding to the r3 source
                field_point_i[0]    =   field_points[3*j];
                field_point_i[1]    =   field_points[3*j+1];
                field_point_i[2]    =   field_points[3*j+2] + 2.0 * input->water_depth;
                calculate_source_velocity_newman(
                                                    mesh_gp->panels[i],
                                                    field_point_i, 
                                                    0,
                                                    0, 
                                                    vel_3
                                                );

                // Calculate velocity corresponding to the r4 source
                field_point_i[0]    =   field_points[3*j];
                field_point_i[1]    =   field_points[3*j+1];
                field_point_i[2]    =   -field_points[3*j+2] + 2.0 * input->water_depth;
                calculate_source_velocity_newman(
                                                    mesh_gp->panels_mirror[i],
                                                    field_point_i, 
                                                    0,
                                                    0, 
                                                    vel_4
                                                );

                // Calculate velocity corresponding to the r5 source
                field_point_i[0]    =   field_points[3*j];
                field_point_i[1]    =   field_points[3*j+1];
                field_point_i[2]    =   field_points[3*j+2] + 4.0 * input->water_depth;
                calculate_source_velocity_newman(
                                                    mesh_gp->panels_mirror[i],
                                                    field_point_i, 
                                                    0,
                                                    0, 
                                                    vel_5
                                                );
                
                // Compose total velocity vector
                vel_total[0]    = vel_0[0] + vel_1[0] + vel_2[0] + vel_3[0] + vel_4[0] + vel_5[0];
                vel_total[1]    = vel_0[1] + vel_1[1] + vel_2[1] + vel_3[1] + vel_4[1] + vel_5[1];
                vel_total[2]    = vel_0[2] + vel_1[2] + vel_2[2] + vel_3[2] + vel_4[2] + vel_5[2];
                                        
                int_value = (
                                mesh_gp->source_nodes[j]->normal_vec[0] * vel_total[0]
                                +
                                mesh_gp->source_nodes[j]->normal_vec[1] * vel_total[1]
                                +
                                mesh_gp->source_nodes[j]->normal_vec[2] * vel_total[2]
                            ) / 4.0 / PI;
            }

            if ( 
                    mesh_gp->source_nodes[j]->panel->type == LID_PANEL_CODE
                    &&
                    panel_i->type == LID_PANEL_CODE
                )
            {
                int_value       = - int_value;
            }

            sysmat[col_count*scl->num_rows_local+row_count] = int_value;

            // Advance row count
            row_count++;
        }

        // Advance column count
        col_count++;
    }

    // Delete heap memory associated to this block of code
    mkl_free( field_points );
}


void    calculate_gf_intensity_steady_sysmat_nlin(
                                                    Input*          input,
                                                    SclCmpx*        scl,
                                                    MeshGroup*      mesh_gp,
                                                    GRFDnInterface* grf_interf,
                                                    cuscomplex*     sysmat
                                                )
{
    /***************************************/
    /******** Fill system matrix  **********/
    /***************************************/
    // Loop over panels to integrate value
    int         col_count   = 0;
    cuscomplex  int_value( 0.0, 0.0 );
    int         row_count   = 0;

    
    // Create auxiliar lambda function
    auto steady_fcn =   [grf_interf]
                        ( 
                            cusfloat xi, 
                            cusfloat eta,
                            cusfloat x, 
                            cusfloat y, 
                            cusfloat z 
                        )
                        {
                            return (*grf_interf)( xi, eta, x, y, z );
                        };

    // Declare local variables to work with the adaptive
    // integration scheme
    SourceNode*     source_i    = nullptr;
    
    for ( int i=scl->start_col_0; i<scl->end_col_0; i++ )
    {
        // Get memory address of the ith panel
        source_i = mesh_gp->source_nodes[i];
        grf_interf->set_source_i( source_i );

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
        {
            // Get memory address of the panel jth
            grf_interf->set_source_j( mesh_gp->source_nodes[j] );

            // Integrate green function normal derivative along the current panel
            if ( i == j )
            {
                int_value       = cuscomplex( 0.0, 0.0 );
            }
            else
            {
                int_value       =   adaptive_quadrature_panel(
                                                                source_i->panel,
                                                                steady_fcn,
                                                                input->gfdn_abs_err,
                                                                input->gfdn_rel_err,
                                                                input->is_block_adaption,
                                                                false,
                                                                input->gauss_order
                                                            );
                int_value       =   int_value / 4.0 / PI;
            }

            if ( 
                    mesh_gp->source_nodes[j]->panel->type == LID_PANEL_CODE
                    &&
                    source_i->panel->type == LID_PANEL_CODE
                )
            {
                int_value       = - int_value;
            }

            sysmat[col_count*scl->num_rows_local+row_count] = int_value;

            // Advance row count
            row_count++;
        }

        // Advance column count
        col_count++;
    }
}