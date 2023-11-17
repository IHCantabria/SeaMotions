
// Include local modules
#include "velocities.hpp"

#include "../../green/source.hpp"
#include "../../math/integration.hpp"


void    calculate_raddif_velocity_mat_steady(
                                                Input*      input,
                                                MeshGroup*  mesh_gp,
                                                MLGCmpx*    vel_x_mat,
                                                MLGCmpx*    vel_y_mat,
                                                MLGCmpx*    vel_z_mat
                                            )
{
    /***************************************/
    /******** Fill system matrix  **********/
    /***************************************/
    // Loop over panels to integrate value
    int         col_count       = 0;
    cusfloat*   field_points    = vel_x_mat->field_points;
    int         row_count       = 0;
    int         rows_np         = vel_x_mat->sysmat_nrows;

    // Define local variables to work with the fast solver
    const int   ndim                    = 3;
    cusfloat    field_point_i[ndim];    clear_vector( ndim, field_point_i );
    PanelGeom*  panel_j                 = nullptr;
    PanelGeom*  panel_mirror_j          = nullptr;
    cusfloat    vel_0[ndim];            clear_vector( ndim, vel_0 );
    cusfloat    vel_1[ndim];            clear_vector( ndim, vel_1 );
    cusfloat    vel_2[ndim];            clear_vector( ndim, vel_2 );
    cusfloat    vel_3[ndim];            clear_vector( ndim, vel_3 );
    cusfloat    vel_4[ndim];            clear_vector( ndim, vel_4 );
    cusfloat    vel_5[ndim];            clear_vector( ndim, vel_5 );
    cusfloat    vel_total[ndim];        clear_vector( ndim, vel_total );

    // Loop over panels and field points to create the steady source matrix
    for ( int i=vel_x_mat->start_row; i<vel_x_mat->end_row; i++ )
    {

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=vel_x_mat->start_col; j<vel_x_mat->end_col; j++ )
        {
            // Get pointer to ith panel
            panel_j         = mesh_gp->panels[j];
            panel_mirror_j  = mesh_gp->panels_mirror[j];

            if ( 
                    mesh_gp->panels[i]->type == DIFFRAC_PANEL_CODE
                    &&
                    mesh_gp->panels[j]->type == DIFFRAC_PANEL_CODE
                )
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
                                                    panel_j,
                                                    &(field_points[3*i]), 
                                                    0,
                                                    0, 
                                                    vel_0
                                                );

                // Calculate velocity corresponding to the r1 source
                field_point_i[0]    =   field_points[3*i];
                field_point_i[1]    =   field_points[3*i+1];
                field_point_i[2]    =   field_points[3*i+2] + 2 * input->water_depth;
                calculate_source_velocity_newman(
                                                    panel_mirror_j,
                                                    field_point_i, 
                                                    0,
                                                    0, 
                                                    vel_1
                                                );
                
                // Calculate velocity corresponding to the r2 source
                field_point_i[0]    =   field_points[3*i];
                field_point_i[1]    =   field_points[3*i+1];
                field_point_i[2]    =   field_points[3*i+2];
                calculate_source_velocity_newman(
                                                    panel_mirror_j,
                                                    field_point_i, 
                                                    0,
                                                    0, 
                                                    vel_2
                                                );

                // Calculate velocity corresponding to the r3 source
                field_point_i[0]    =   field_points[3*i];
                field_point_i[1]    =   field_points[3*i+1];
                field_point_i[2]    =   field_points[3*i+2] + 2.0 * input->water_depth;
                calculate_source_velocity_newman(
                                                    panel_j,
                                                    field_point_i, 
                                                    0,
                                                    0, 
                                                    vel_3
                                                );

                // Calculate velocity corresponding to the r4 source
                field_point_i[0]    =   field_points[3*i];
                field_point_i[1]    =   field_points[3*i+1];
                field_point_i[2]    =   -field_points[3*i+2] + 2.0 * input->water_depth;
                calculate_source_velocity_newman(
                                                    panel_mirror_j,
                                                    field_point_i, 
                                                    0,
                                                    0, 
                                                    vel_4
                                                );

                // Calculate velocity corresponding to the r5 source
                field_point_i[0]    =   field_points[3*i];
                field_point_i[1]    =   field_points[3*i+1];
                field_point_i[2]    =   field_points[3*i+2] + 4.0 * input->water_depth;
                calculate_source_velocity_newman(
                                                    panel_mirror_j,
                                                    field_point_i, 
                                                    0,
                                                    0, 
                                                    vel_5
                                                );
                
                // Compose total velocity vector
                vel_total[0]    = vel_0[0] + vel_1[0] + vel_2[0] + vel_3[0] + vel_4[0] + vel_5[0];
                vel_total[1]    = vel_0[1] + vel_1[1] + vel_2[1] + vel_3[1] + vel_4[1] + vel_5[1];
                vel_total[2]    = vel_0[2] + vel_1[2] + vel_2[2] + vel_3[2] + vel_4[2] + vel_5[2];
                                        
                
                vel_x_mat->sysmat_steady[row_count*rows_np+col_count] = vel_total[0];
                vel_y_mat->sysmat_steady[row_count*rows_np+col_count] = vel_total[1];
                vel_z_mat->sysmat_steady[row_count*rows_np+col_count] = vel_total[2];
            }

            // Advance column count
            col_count++;
        }

        // Advance row count
        row_count++;
    }

}


void    calculate_raddif_velocity_mat(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                GWFDxInterface* gwfdx_interf,
                                                GWFDyInterface* gwfdy_interf,
                                                GWFDzInterface* gwfdz_interf,
                                                MLGCmpx*        vel_x_mat,
                                                MLGCmpx*        vel_y_mat,
                                                MLGCmpx*        vel_z_mat
                                        )
{
    /***************************************/
    /******** Fill system matrix  **********/
    /***************************************/

    // Define local variables to work with the fast solver
    int         col_count       = 0;
    int         index                   = 0;
    const int   ndim                    = 3;
    int         row_count       = 0;
    int         rows_np         = vel_x_mat->sysmat_nrows;
    cuscomplex  vel_total[ndim];

    // Define lambda function for the integrations interface
    auto gwfdx_fcn =   [gwfdx_interf]
                        ( 
                            cusfloat xi, 
                            cusfloat eta,
                            cusfloat x, 
                            cusfloat y, 
                            cusfloat z 
                        )
                        {
                            return (*gwfdx_interf)( xi, eta, x, y, z );
                        };
    
    auto gwfdy_fcn =   [gwfdy_interf]
                        ( 
                            cusfloat xi, 
                            cusfloat eta,
                            cusfloat x, 
                            cusfloat y, 
                            cusfloat z 
                        )
                        {
                            return (*gwfdy_interf)( xi, eta, x, y, z );
                        };

    auto gwfdz_fcn =   [gwfdz_interf]
                        ( 
                            cusfloat xi, 
                            cusfloat eta,
                            cusfloat x, 
                            cusfloat y, 
                            cusfloat z 
                        )
                        {
                            return (*gwfdz_interf)( xi, eta, x, y, z );
                        };

    // Loop over panels and field points to create the steady source matrix
    for ( int i=vel_x_mat->start_row; i<vel_x_mat->end_row; i++ )
    {
        // Get memory address of the ith panel
        gwfdx_interf->set_source_i( mesh_gp->source_nodes[i] );
        gwfdy_interf->set_source_i( mesh_gp->source_nodes[i] );
        gwfdz_interf->set_source_i( mesh_gp->source_nodes[i] );

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=vel_x_mat->start_col; j<vel_x_mat->end_col; j++ )
        {
            // Get memory address of the panel jth
            gwfdx_interf->set_source_j( mesh_gp->source_nodes[j] );
            gwfdy_interf->set_source_j( mesh_gp->source_nodes[j] );
            gwfdz_interf->set_source_j( mesh_gp->source_nodes[j] );

            if ( 
                    mesh_gp->panels[i]->type == DIFFRAC_PANEL_CODE
                    &&
                    mesh_gp->panels[j]->type == DIFFRAC_PANEL_CODE
                )
            {
                

                // Integrate green function normal derivative along the current panel
                vel_total[0]    =   adaptive_quadrature_panel(
                                                                mesh_gp->panels[j],
                                                                gwfdx_fcn,
                                                                input->gfdn_abs_err,
                                                                input->gfdn_rel_err,
                                                                input->is_block_adaption,
                                                                false,
                                                                input->gauss_order
                                                            );

                vel_total[1]    =   adaptive_quadrature_panel(
                                                                mesh_gp->panels[j],
                                                                gwfdy_fcn,
                                                                input->gfdn_abs_err,
                                                                input->gfdn_rel_err,
                                                                input->is_block_adaption,
                                                                false,
                                                                input->gauss_order
                                                            );

                vel_total[2]    =   adaptive_quadrature_panel(
                                                                mesh_gp->panels[j],
                                                                gwfdz_fcn,
                                                                input->gfdn_abs_err,
                                                                input->gfdn_rel_err,
                                                                input->is_block_adaption,
                                                                false,
                                                                input->gauss_order
                                                            );

                index   = row_count*rows_np+col_count;
                vel_x_mat->sysmat[index] = vel_x_mat->sysmat_steady[index] + vel_total[0];
                vel_y_mat->sysmat[index] = vel_y_mat->sysmat_steady[index] + vel_total[1];
                vel_z_mat->sysmat[index] = vel_z_mat->sysmat_steady[index] + vel_total[2];
            }

            // Advance column count
            col_count++;
        }

        // Advance row count
        row_count++;
    }

}