
// Include local modules
#include "velocities.hpp"

#include "../../containers/simulation_data.hpp"
#include "../../green/source.hpp"
#include "../../math/integration.hpp"
#include "tools.hpp"


void    calculate_raddif_velocity_mat(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                GWFDxInterface* gwfdx_interf,
                                                GWFDyInterface* gwfdy_interf,
                                                GWFDzInterface* gwfdz_interf,
                                                MLGCmpx*        vel_x_gp,
                                                MLGCmpx*        vel_y_gp,
                                                MLGCmpx*        vel_z_gp
                                        )
{
    /***************************************/
    /******** Fill system matrix  **********/
    /***************************************/

    // Define local variables to work with the fast solver
    int         col_count       = 0;
    int         index           = 0;
    const int   ndim            = 3;
    int         row_count       = 0;
    int         rows_np         = vel_x_gp->sysmat_nrows;
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
    for ( int i=vel_x_gp->start_row; i<vel_x_gp->end_row; i++ )
    {
        // Get memory address of the ith panel
        gwfdx_interf->set_source_i( mesh_gp->source_nodes[i] );
        gwfdy_interf->set_source_i( mesh_gp->source_nodes[i] );
        gwfdz_interf->set_source_i( mesh_gp->source_nodes[i] );

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=vel_x_gp->start_col; j<vel_x_gp->end_col; j++ )
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
                vel_x_gp->sysmat[index] = vel_x_gp->sysmat_steady[index] + vel_total[0];
                vel_y_gp->sysmat[index] = vel_y_gp->sysmat_steady[index] + vel_total[1];
                vel_z_gp->sysmat[index] = vel_z_gp->sysmat_steady[index] + vel_total[2];
            }

            // Advance column count
            col_count++;
        }

        // Advance row count
        row_count++;
    }

}


void    calculate_raddif_velocity_mat_steady(
                                                Input*      input,
                                                MeshGroup*  mesh_gp,
                                                MLGCmpx*    vel_x_gp,
                                                MLGCmpx*    vel_y_gp,
                                                MLGCmpx*    vel_z_gp
                                            )
{
    /***************************************/
    /******** Fill system matrix  **********/
    /***************************************/
    // Loop over panels to integrate value
    int         col_count       = 0;
    cusfloat*   field_points    = vel_x_gp->field_points;
    int         row_count       = 0;
    int         rows_np         = vel_x_gp->sysmat_nrows;

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
    for ( int i=vel_x_gp->start_row; i<vel_x_gp->end_row; i++ )
    {

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=vel_x_gp->start_col; j<vel_x_gp->end_col; j++ )
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
                                        
                
                vel_x_gp->sysmat_steady[row_count*rows_np+col_count] = vel_total[0];
                vel_y_gp->sysmat_steady[row_count*rows_np+col_count] = vel_total[1];
                vel_z_gp->sysmat_steady[row_count*rows_np+col_count] = vel_total[2];
            }

            // Advance column count
            col_count++;
        }

        // Advance row count
        row_count++;
    }

}


void    calculate_velocities_total(
                                        Input*          input,
                                        MpiConfig*      mpi_config,
                                        cusfloat        ang_freq,
                                        cuscomplex*     intensities,
                                        cuscomplex*     raos,
                                        MLGCmpx*        vel_x_gp,
                                        MLGCmpx*        vel_y_gp,
                                        MLGCmpx*        vel_z_gp,
                                        cuscomplex*     vel_x_total,
                                        cuscomplex*     vel_y_total,
                                        cuscomplex*     vel_z_total
                                    )
{
    // Declare auxiliar variables to use in the function
    int index       =   0;
    int index_2     =   0;
    int index_3     =   0;

    /*******************************************************/
    /***  Calculate radiation and diffraction velocities ***/
    /*******************************************************/
    calculate_fields_raddif_lin(
                                    input,
                                    intensities,
                                    vel_x_gp
                                );

    calculate_fields_raddif_lin(
                                    input,
                                    intensities,
                                    vel_y_gp
                                );
    
    calculate_fields_raddif_lin(
                                    input,
                                    intensities,
                                    vel_z_gp
                                );


    /***************************************************************/
    /******** Sum panel velocities from all processes **************/
    /***************************************************************/
    cuscomplex* vel_x_raddif_p0  = nullptr;
    cuscomplex* vel_y_raddif_p0  = nullptr;
    cuscomplex* vel_z_raddif_p0  = nullptr;

    if ( mpi_config->is_root( ) )
    {
        vel_x_raddif_p0   = generate_empty_vector<cuscomplex>( vel_x_gp->fields_np * vel_x_gp->sysmat_nrows );
        vel_y_raddif_p0   = generate_empty_vector<cuscomplex>( vel_x_gp->fields_np * vel_x_gp->sysmat_nrows );
        vel_z_raddif_p0   = generate_empty_vector<cuscomplex>( vel_x_gp->fields_np * vel_x_gp->sysmat_nrows );
    }

    MPI_Reduce(
                    vel_x_gp->field_values,
                    vel_x_raddif_p0,
                    vel_x_gp->fields_np * vel_x_gp->sysmat_nrows,
                    mpi_cuscomplex,
                    MPI_SUM,
                    mpi_config->proc_root,
                    MPI_COMM_WORLD
                );

    MPI_Reduce(
                    vel_y_gp->field_values,
                    vel_y_raddif_p0,
                    vel_y_gp->fields_np * vel_y_gp->sysmat_nrows,
                    mpi_cuscomplex,
                    MPI_SUM,
                    mpi_config->proc_root,
                    MPI_COMM_WORLD
                );

    MPI_Reduce(
                    vel_z_gp->field_values,
                    vel_z_raddif_p0,
                    vel_z_gp->fields_np * vel_z_gp->sysmat_nrows,
                    mpi_cuscomplex,
                    MPI_SUM,
                    mpi_config->proc_root,
                    MPI_COMM_WORLD
                );
    
    /*******************************************************/
    /*******  Add incident wave potential velocities *******/
    /*******************************************************/

    if ( mpi_config->is_root( ) )
    {
        // Calculate incident wave potential
        cusfloat    k   = w2k( 
                                    ang_freq,
                                    input->water_depth,
                                    input->grav_acc
                                );
        
        for ( int i=0; i<input->heads_np; i++ )
        {
            for ( int j=0; j<vel_x_gp->field_points_np; j++ )
            {
                index                   =   vel_x_gp->field_points_np * i + j;
                vel_x_total[index]      =   wave_potential_airy_space_dx(
                                                                            input->wave_amplitude,
                                                                            ang_freq,
                                                                            k,
                                                                            input->water_depth,
                                                                            input->grav_acc,
                                                                            vel_x_gp->field_points[3*j],
                                                                            vel_x_gp->field_points[3*j+1],
                                                                            vel_x_gp->field_points[3*j+2],
                                                                            input->heads[i]
                                                                        );

                vel_y_total[index]      =   wave_potential_airy_space_dy(
                                                                            input->wave_amplitude,
                                                                            ang_freq,
                                                                            k,
                                                                            input->water_depth,
                                                                            input->grav_acc,
                                                                            vel_x_gp->field_points[3*j],
                                                                            vel_x_gp->field_points[3*j+1],
                                                                            vel_x_gp->field_points[3*j+2],
                                                                            input->heads[i]
                                                                        );

                vel_z_total[index]      =   wave_potential_airy_space_dz(
                                                                            input->wave_amplitude,
                                                                            ang_freq,
                                                                            k,
                                                                            input->water_depth,
                                                                            input->grav_acc,
                                                                            vel_x_gp->field_points[3*j],
                                                                            vel_x_gp->field_points[3*j+1],
                                                                            vel_x_gp->field_points[3*j+2],
                                                                            input->heads[i]
                                                                        );
            }
        }
    }

    /*******************************************************/
    /************  Add  diffraction velocities *************/
    /*******************************************************/

    if ( mpi_config->is_root( ) )
    {
        for ( int i=0; i<input->heads_np; i++ )
        {
            for ( int j=0; j<vel_x_gp->field_points_np; j++ )
            {
                index                   = vel_x_gp->field_points_np * i + j;
                index_2                 = ( input->dofs_np + i ) * vel_x_gp->field_points_np + j;
                vel_x_total[index]      += vel_x_raddif_p0[index_2];
                vel_y_total[index]      += vel_y_raddif_p0[index_2];
                vel_z_total[index]      += vel_z_raddif_p0[index_2];
            }
        }
    }

    /***************************************************************/
    /*************** Add Radiation Wave Velocities *****************/
    /***************************************************************/

    if ( mpi_config->is_root( ) )
    {
        for ( int i=0; i<input->heads_np; i++ )
        {
            for ( int j=0; j<vel_x_gp->field_points_nb; j++ )
            {
                for ( int k=0; k<input->dofs_np; k++ )
                {
                    for ( int r=vel_x_gp->field_points_cnp[j]; r<vel_x_gp->field_points_cnp[j+1]; r++ )
                    {
                        index                   = i * vel_x_gp->field_points_np + r;
                        index_2                 = k * vel_x_gp->field_points_np + r;
                        index_3                 = i * ( input->dofs_np * vel_x_gp->field_points_nb ) + j * input->dofs_np + k;
                        vel_x_total[index]      += raos[index_3] * vel_x_raddif_p0[index_2];
                        vel_y_total[index]      += raos[index_3] * vel_y_raddif_p0[index_2];
                        vel_z_total[index]      += raos[index_3] * vel_z_raddif_p0[index_2];
                    }
                }
            }
        }
    }

    /*******************************************************/
    /**************  Deallocate heap memory ****************/
    /*******************************************************/ 
    if ( mpi_config->is_root( ) )
    {
        mkl_free( vel_x_raddif_p0 );
        mkl_free( vel_y_raddif_p0 );
        mkl_free( vel_z_raddif_p0 );
    }

}