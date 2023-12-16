
// Include general usage libraries
#include <fstream>

// Include local modules
#include "velocities.hpp"

#include "../../containers/simulation_data.hpp"
#include "../../green/source.hpp"
#include "../../math/integration.hpp"
#include "../../interfaces/grfdx_interface.hpp"
#include "../../interfaces/grfdy_interface.hpp"
#include "../../interfaces/grfdz_interface.hpp"
#include "../../interfaces/gwfdx_interface.hpp"
#include "../../interfaces/gwfdy_interface.hpp"
#include "../../interfaces/gwfdz_interface.hpp"
#include "tools.hpp"


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
    int         cols_np         = vel_x_gp->sysmat_ncols;
    cusfloat*   field_points    = vel_x_gp->field_points;
    int         row_count       = 0;

    // Define local variables to work with the fast solver
    const int   ndim                    = 3;
    cusfloat    dist_fp_cp              = 0.0;
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

    // Get correction factor for the number of field points if the
    // non adaptive quadrature is being used
    int fc = input->gauss_np_factor_2d( );

    // Loop over panels and field points to create the steady source matrix
    for ( int i=vel_x_gp->start_row; i<vel_x_gp->end_row+1; i++ )
    {

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        col_count = 0;
        for ( int j=vel_x_gp->start_col; j<vel_x_gp->end_col; j++ )
        {
            // Get pointer to ith panel
            panel_j         = mesh_gp->panels[j];
            panel_mirror_j  = mesh_gp->panels_mirror[j];

            if ( 
                    mesh_gp->panels[i/fc]->type == DIFFRAC_PANEL_CODE
                    &&
                    mesh_gp->panels[j]->type == DIFFRAC_PANEL_CODE
                )
            {
                dist_fp_cp = eucledian_dist(
                                                3,
                                                panel_j->center,
                                                &(field_points[3*i])
                                            );

                dist_fp_cp = 5.0;
                
                if ( dist_fp_cp > 1e-6 )
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
                                            
                    
                    vel_x_gp->sysmat_steady[row_count*cols_np+col_count] = vel_total[0] / 4.0 / PI;
                    vel_y_gp->sysmat_steady[row_count*cols_np+col_count] = vel_total[1] / 4.0 / PI;
                    vel_z_gp->sysmat_steady[row_count*cols_np+col_count] = vel_total[2] / 4.0 / PI;
                }
                else
                {
                    vel_x_gp->sysmat_steady[row_count*cols_np+col_count] = - 0.5 * panel_j->normal_vec[0];
                    vel_y_gp->sysmat_steady[row_count*cols_np+col_count] = - 0.5 * panel_j->normal_vec[1];
                    vel_z_gp->sysmat_steady[row_count*cols_np+col_count] = - 0.5 * panel_j->normal_vec[2];
                }
            }

            // Advance column count
            col_count++;
        }

        // Advance row count
        row_count++;
    }

}


void    calculate_raddif_velocity_mat_steady_nlin(
                                                    Input*          input,
                                                    MeshGroup*      mesh_gp,
                                                    MLGCmpx*        vel_x_gp,
                                                    MLGCmpx*        vel_y_gp,
                                                    MLGCmpx*        vel_z_gp
                                                )
{
    // Define potential funcions objects interface
    GRFDxInterface* green_rank_dx   = new   GRFDxInterface( 
                                                                mesh_gp->source_nodes[0],
                                                                mesh_gp->source_nodes[0],
                                                                input->water_depth
                                                            );
    
    GRFDyInterface* green_rank_dy   = new   GRFDyInterface( 
                                                                mesh_gp->source_nodes[0],
                                                                mesh_gp->source_nodes[0],
                                                                input->water_depth
                                                            );

    GRFDzInterface* green_rank_dz   = new   GRFDzInterface( 
                                                                mesh_gp->source_nodes[0],
                                                                mesh_gp->source_nodes[0],
                                                                input->water_depth
                                                            );

    auto            rank_dx_fcn     =   [green_rank_dx]
                                        (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z) -> cuscomplex
                                        {
                                            return (*green_rank_dx)( xi, eta, x, y, z );
                                        };
    
    auto            rank_dy_fcn     =   [green_rank_dy]
                                        (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z) -> cuscomplex
                                        {
                                            return (*green_rank_dy)( xi, eta, x, y, z );
                                        };

    auto            rank_dz_fcn     =   [green_rank_dz]
                                        (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z) -> cuscomplex
                                        {
                                            return (*green_rank_dz)( xi, eta, x, y, z );
                                        };

    // Generate potential matrix
    int         count = 0;
    cuscomplex  vel_x_rank_term( 0.0, 0.0 );
    cuscomplex  vel_y_rank_term( 0.0, 0.0 );
    cuscomplex  vel_z_rank_term( 0.0, 0.0 );
    for ( int i=0; i<vel_x_gp->field_points_np; i++ )
    {
        // Change field point
        green_rank_dx->set_field_point_j(
                                                &(vel_x_gp->field_points[3*i])
                                            );
        green_rank_dy->set_field_point_j(
                                                &(vel_y_gp->field_points[3*i])
                                            );
        green_rank_dz->set_field_point_j(
                                                &(vel_z_gp->field_points[3*i])
                                            );

        for ( int j=vel_x_gp->start_col; j<vel_x_gp->end_col; j++ )
        {
            // Compute steady and wave terms over the panel
            if ( 
                    mesh_gp->panels[i]->type == DIFFRAC_PANEL_CODE
                    &&
                    mesh_gp->panels[j]->type == DIFFRAC_PANEL_CODE
                )
            {
                vel_x_rank_term         = adaptive_quadrature_panel(
                                                                        mesh_gp->panels[j],
                                                                        rank_dx_fcn,
                                                                        input->pot_abs_err,
                                                                        input->pot_rel_err,
                                                                        false,
                                                                        true,
                                                                        input->gauss_order
                                                                    );
                vel_y_rank_term         = adaptive_quadrature_panel(
                                                                        mesh_gp->panels[j],
                                                                        rank_dy_fcn,
                                                                        input->pot_abs_err,
                                                                        input->pot_rel_err,
                                                                        false,
                                                                        true,
                                                                        input->gauss_order
                                                                    );
                vel_z_rank_term         = adaptive_quadrature_panel(
                                                                        mesh_gp->panels[j],
                                                                        rank_dz_fcn,
                                                                        input->pot_abs_err,
                                                                        input->pot_rel_err,
                                                                        false,
                                                                        true,
                                                                        input->gauss_order
                                                                    );
                vel_x_gp->sysmat_steady[count]  = vel_x_rank_term / 4.0 / PI;
                vel_y_gp->sysmat_steady[count]  = vel_y_rank_term / 4.0 / PI;
                vel_z_gp->sysmat_steady[count]  = vel_z_rank_term / 4.0 / PI;

            }
            else
            {
                vel_x_gp->sysmat[count]   = cuscomplex( 0.0, 0.0 );
                vel_y_gp->sysmat[count]   = cuscomplex( 0.0, 0.0 );
                vel_z_gp->sysmat[count]   = cuscomplex( 0.0, 0.0 );
            }
            
            count++;
        }
    }

    // Delete heap memory allocated in the current function
    delete green_rank_dx;
    delete green_rank_dy;
    delete green_rank_dz;

}


void    calculate_raddif_velocity_mat_wave(
                                                    Input*          input,
                                                    MeshGroup*      mesh_gp,
                                                    cusfloat        ang_freq,
                                                    MLGCmpx*        vel_x_gp,
                                                    MLGCmpx*        vel_y_gp,
                                                    MLGCmpx*        vel_z_gp
                                            )
{
    // Define potential funcions objects interface
    GWFDxInterface* green_wave_dx   = new   GWFDxInterface( 
                                                                mesh_gp->source_nodes[0],
                                                                mesh_gp->source_nodes[0],
                                                                ang_freq,
                                                                input->water_depth,
                                                                input->grav_acc
                                                            );
    
    GWFDyInterface* green_wave_dy   = new   GWFDyInterface( 
                                                                mesh_gp->source_nodes[0],
                                                                mesh_gp->source_nodes[0],
                                                                ang_freq,
                                                                input->water_depth,
                                                                input->grav_acc
                                                            );

    GWFDzInterface* green_wave_dz   = new   GWFDzInterface( 
                                                                mesh_gp->source_nodes[0],
                                                                mesh_gp->source_nodes[0],
                                                                ang_freq,
                                                                input->water_depth,
                                                                input->grav_acc
                                                            );

    auto            wave_dx_fcn     =   [green_wave_dx]
                                        (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z) -> cuscomplex
                                        {
                                            return (*green_wave_dx)( xi, eta, x, y, z );
                                        };
    
    auto            wave_dy_fcn     =   [green_wave_dy]
                                        (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z) -> cuscomplex
                                        {
                                            return (*green_wave_dy)( xi, eta, x, y, z );
                                        };

    auto            wave_dz_fcn     =   [green_wave_dz]
                                        (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z) -> cuscomplex
                                        {
                                            return (*green_wave_dz)( xi, eta, x, y, z );
                                        };

    // Generate potential matrix
    int         count       = 0;
    cusfloat    dist_fp_cp  = 0.0;
    cuscomplex  vel_x_wave_term( 0.0, 0.0 );
    cuscomplex  vel_y_wave_term( 0.0, 0.0 );
    cuscomplex  vel_z_wave_term( 0.0, 0.0 );

    // Get correction factor for the number of field points if the
    // non adaptive quadrature is being used
    int fc = input->gauss_np_factor_2d( );

    for ( int i=vel_x_gp->start_row; i<vel_x_gp->end_row+1; i++ )
    {
        // Change field point
        green_wave_dx->set_field_point_j(
                                                &(vel_x_gp->field_points[3*i])
                                            );
        green_wave_dy->set_field_point_j(
                                                &(vel_y_gp->field_points[3*i])
                                            );
        green_wave_dz->set_field_point_j(
                                                &(vel_z_gp->field_points[3*i])
                                            );

        for ( int j=vel_x_gp->start_col; j<vel_x_gp->end_col; j++ )
        {
            // Compute steady and wave terms over the panel
            if ( 
                    mesh_gp->panels[i/fc]->type == DIFFRAC_PANEL_CODE
                    &&
                    mesh_gp->panels[j]->type == DIFFRAC_PANEL_CODE
                )
            {
                dist_fp_cp = eucledian_dist(
                                                3,
                                                mesh_gp->panels[j]->center,
                                                &(vel_x_gp->field_points[3*i])
                                            );
                dist_fp_cp = 5.0;
                if ( dist_fp_cp > 1e-6 )
                {
                    vel_x_wave_term         = adaptive_quadrature_panel(
                                                                            mesh_gp->panels[j],
                                                                            wave_dx_fcn,
                                                                            input->gfdn_abs_err,
                                                                            input->gfdn_rel_err,
                                                                            input->is_block_adaption,
                                                                            false,
                                                                            input->gauss_order
                                                                        );
                    vel_y_wave_term         = adaptive_quadrature_panel(
                                                                            mesh_gp->panels[j],
                                                                            wave_dy_fcn,
                                                                            input->gfdn_abs_err,
                                                                            input->gfdn_rel_err,
                                                                            input->is_block_adaption,
                                                                            false,
                                                                            input->gauss_order
                                                                        );
                    vel_z_wave_term         = adaptive_quadrature_panel(
                                                                            mesh_gp->panels[j],
                                                                            wave_dz_fcn,
                                                                            input->gfdn_abs_err,
                                                                            input->gfdn_rel_err,
                                                                            input->is_block_adaption,
                                                                            false,
                                                                            input->gauss_order
                                                                        );
                    vel_x_gp->sysmat[count] = vel_x_gp->sysmat_steady[count] + vel_x_wave_term / 4.0 / PI;
                    vel_y_gp->sysmat[count] = vel_y_gp->sysmat_steady[count] + vel_y_wave_term / 4.0 / PI;
                    vel_z_gp->sysmat[count] = vel_z_gp->sysmat_steady[count] + vel_z_wave_term / 4.0 / PI;
                }
                else
                {
                    vel_x_gp->sysmat[count] = vel_x_gp->sysmat_steady[count];
                    vel_y_gp->sysmat[count] = vel_y_gp->sysmat_steady[count];
                    vel_z_gp->sysmat[count] = vel_z_gp->sysmat_steady[count];
                }
            }
            else
            {
                vel_x_gp->sysmat[count]   = cuscomplex( 0.0, 0.0 );
                vel_y_gp->sysmat[count]   = cuscomplex( 0.0, 0.0 );
                vel_z_gp->sysmat[count]   = cuscomplex( 0.0, 0.0 );
            }
            
            count++;
        }
    }

    // Delete heap memory allocated in the current function
    delete green_wave_dx;
    delete green_wave_dy;
    delete green_wave_dz;

}


void    calculate_velocities_total(
                                                    Input*          input,
                                                    MpiConfig*      mpi_config,
                                                    MeshGroup*      mesh_gp,
                                                    cusfloat        ang_freq,
                                                    int             ang_freq_num,
                                                    cuscomplex*     intensities,
                                                    cuscomplex*     raos,
                                                    MLGCmpx*        vel_x_gp,
                                                    MLGCmpx*        vel_y_gp,
                                                    MLGCmpx*        vel_z_gp,
                                                    cuscomplex*     vel_x_total,
                                                    cuscomplex*     vel_y_total,
                                                    cuscomplex*     vel_z_total,
                                                    SimulationData* sim_data
                                    )
{
    // Declare auxiliar variables to use in the function
    int index       =   0;
    int index_2     =   0;
    int index_3     =   0;
    
    /*******************************************************/
    /*********  Calculate total velocity matrixes **********/
    /*******************************************************/
    calculate_raddif_velocity_mat_wave(
                                            input,
                                            mesh_gp,
                                            ang_freq,
                                            vel_x_gp,
                                            vel_y_gp,
                                            vel_z_gp
                                        );

    cuscomplex value( 0.0, 0.0 );
    int ngpf    = input->gauss_np_factor_2d( );
    int _idx0   = 0;
    int _idx1   = 0;
    // for ( int i=0; i<vel_x_gp->field_points_np/ngpf; i++ )
    // {
    //     for ( int j=0; j<ngpf; j++ )
    //     {
    //         for ( int k=0; k<vel_x_gp->sysmat_ncols; k++ )
    //         {
    //             _idx0 = i * ngpf * vel_x_gp->sysmat_ncols + k;
    //             _idx1 = i * ngpf * vel_x_gp->sysmat_ncols + j * vel_x_gp->sysmat_ncols + k;
    //             value = vel_x_gp->sysmat[_idx0] - vel_x_gp->sysmat[_idx1];

    //             if ( i == 6 && j == 3 )
    //             {
    //                 double ccc = 0.0;
    //             }

    //             if ( !assert_complex_equality( value, cuscomplex( 0.0, 0.0 ), 1e-3 ) )
    //             {
    //                 std::cout << "vel_x_gp.sysmat not equal at - i: " << i << " - j: " << j << " - k: " << k << " - idx0: " << _idx0 << " - idx1: " << _idx1 << " - Value: " << value << std::endl;
    //                 throw std::runtime_error( "" );
    //             }

    //             value = vel_y_gp->sysmat[_idx0] - vel_y_gp->sysmat[_idx1];
    //             if ( !assert_complex_equality( value, cuscomplex( 0.0, 0.0 ), 1e-3 ) )
    //             {
    //                 std::cout << "vel_y_gp.sysmat not equal at - i: " << i << " - j: " << j << " - k: " << k << " - idx0: " << _idx0 << " - idx1: " << _idx1 << " - Value: " << value << std::endl;
    //                 throw std::runtime_error( "" );
    //             }

    //             value = vel_z_gp->sysmat[_idx0] - vel_z_gp->sysmat[_idx1];
    //             if ( !assert_complex_equality( value, cuscomplex( 0.0, 0.0 ), 1e-3 ) )
    //             {
    //                 std::cout << "vel_z_gp.sysmat not equal at - i: " << i << " - j: " << j << " - k: " << k << " - idx0: " << _idx0 << " - idx1: " << _idx1 << " - Value: " << value << std::endl;
    //                 throw std::runtime_error( "" );
    //             }
    //         }
    //     }
    // }

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
                if ( index == 49 || index == 196 )
                {
                    double abc = 0;
                }

                // if ( index == 49 )
                // {
                //     std::cout << std::endl;
                //     std::cout << "vel_x_total[" << index << "]: " << vel_x_total[index] << " - " << vel_x_total << std::endl;
                //     std::cout << "vel_y_total[" << index << "]: " << vel_y_total[index] << " - " << vel_y_total << std::endl;
                //     std::cout << "vel_z_total[" << index << "]: " << vel_z_total[index] << " - " << vel_z_total << std::endl;
                // }
            }
        }
    }

    std::string     base_path( "E:/sergio/0050_OASIS_SM/_check_potentials/sm_freqs/velocity_total/" );
    std::stringstream ss0; ss0 << "vel_x_ang_freq_" << ang_freq_num << ".dat";
    std::stringstream ss1; ss1 << "vel_y_ang_freq_" << ang_freq_num << ".dat";
    std::stringstream ss2; ss2 << "vel_z_ang_freq_" << ang_freq_num << ".dat";
    std::string     vel_x_fipath = base_path + ss0.str( );
    std::string     vel_y_fipath = base_path + ss1.str( );
    std::string     vel_z_fipath = base_path + ss2.str( );
    std::ofstream   vel_x_outf( vel_x_fipath );
    std::ofstream   vel_y_outf( vel_y_fipath );
    std::ofstream   vel_z_outf( vel_z_fipath );
    CHECK_FILE_UNIT_STATUS( vel_x_outf, vel_x_fipath );
    CHECK_FILE_UNIT_STATUS( vel_y_outf, vel_y_fipath );
    CHECK_FILE_UNIT_STATUS( vel_z_outf, vel_z_fipath );

    vel_x_outf << mesh_gp->panels_tnp << "  " << input->dofs_np + input->heads_np + 2 << std::endl;
    vel_y_outf << mesh_gp->panels_tnp << "  " << input->dofs_np + input->heads_np + 2 << std::endl;
    vel_z_outf << mesh_gp->panels_tnp << "  " << input->dofs_np + input->heads_np + 2 << std::endl;

    int        idx      = 0;
    PanelGeom* panel_i  = nullptr;
    for ( int i=0; i<input->dofs_np+input->heads_np; i++ )
    {
        for ( int j=0; j<vel_x_gp->field_points_np; j++ )
        {
            idx = i * vel_x_gp->field_points_np + j;
            vel_x_outf << vel_x_raddif_p0[idx].real( ) << "  ";
            vel_x_outf << vel_x_raddif_p0[idx].imag( ) << std::endl;

            vel_y_outf << vel_y_raddif_p0[idx].real( ) << "  ";
            vel_y_outf << vel_y_raddif_p0[idx].imag( ) << std::endl;

            vel_z_outf << vel_z_raddif_p0[idx].real( ) << "  ";
            vel_z_outf << vel_z_raddif_p0[idx].imag( ) << std::endl;
        }
    }

    for ( int i=0; i<vel_x_gp->field_points_np; i++ )
    {
        vel_x_outf << vel_x_total[i].real( ) << "  ";
        vel_x_outf << vel_x_total[i].imag( ) << std::endl;

        vel_y_outf << vel_y_total[i].real( ) << "  ";
        vel_y_outf << vel_y_total[i].imag( ) << std::endl;

        vel_z_outf << vel_z_total[i].real( ) << "  ";
        vel_z_outf << vel_z_total[i].imag( ) << std::endl;
    }

    

    /*******************************************************/
    /************  Add  diffraction velocities *************/
    /*******************************************************/

    if ( mpi_config->is_root( ) )
    {
        int idtest =0 ;
        cuscomplex a, b, c;
        for ( int i=0; i<input->heads_np; i++ )
        {
            for ( int j=0; j<vel_x_gp->field_points_np; j++ )
            {
                // index                   = vel_x_gp->sysmat_nrows * i + j;
                // index_2                 = ( input->dofs_np + i ) * vel_x_gp->field_points_np + j;
                // idtest                  = index_2;
                // a                       = vel_x_raddif_p0[index_2];
                // b                       = vel_y_raddif_p0[index_2];
                // c                       = vel_z_raddif_p0[index_2];
                // vel_x_total[index]      += a;
                // vel_y_total[index]      += b;
                // vel_z_total[index]      += b;
                index                   = vel_x_gp->sysmat_nrows * i + j;
                index_2                 = ( input->dofs_np + i ) * vel_x_gp->field_points_np + j;
                idtest                  = index_2;

                if ( index == 49 )
                {
                    a                       = vel_x_raddif_p0[index_2];
                    b                       = vel_y_raddif_p0[index_2];
                    c                       = vel_z_raddif_p0[index_2];
                    vel_x_total[index]      += a;
                    vel_y_total[index]      += b;
                    vel_z_total[index]      += c;

                    // std::cout << std::endl;
                    // std::cout << "vel_x_raddif_p0[" << idtest << "]: " << a << std::endl;
                    // std::cout << "vel_y_raddif_p0[" << idtest << "]: " << b << std::endl;
                    // std::cout << "vel_z_raddif_p0[" << idtest << "]: " << c << std::endl;
                    // std::cout << "vel_x_total[" << index << "]: " << vel_x_total[index] << " - " << vel_x_total << std::endl;
                    // std::cout << "vel_y_total[" << index << "]: " << vel_y_total[index] << " - " << vel_y_total << std::endl;
                    // std::cout << "vel_z_total[" << index << "]: " << vel_z_total[index] << " - " << vel_z_total << std::endl;
                }
                else
                {
                    a                       = vel_x_raddif_p0[index_2];
                    b                       = vel_y_raddif_p0[index_2];
                    c                       = vel_z_raddif_p0[index_2];
                    vel_x_total[index]      += a;
                    vel_y_total[index]      += b;
                    vel_z_total[index]      += c;
                }

                if ( index == 49 || index == 196 )
                {
                    double abc = 0;
                }
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
                for ( int r=vel_x_gp->field_points_cnp[j]; r<vel_x_gp->field_points_cnp[j+1]; r++ )
                {
                    for ( int k=0; k<input->dofs_np; k++ )
                    {
                        index                   = i * vel_x_gp->sysmat_nrows + r;
                        index_2                 = k * vel_x_gp->sysmat_nrows + r;
                        index_3                 = i * ( input->dofs_np * vel_x_gp->field_points_nb ) + j * input->dofs_np + k;
                        vel_x_total[index]      += raos[index_3] * vel_x_raddif_p0[index_2];
                        vel_y_total[index]      += raos[index_3] * vel_y_raddif_p0[index_2];
                        vel_z_total[index]      += raos[index_3] * vel_z_raddif_p0[index_2];

                        
                        // if ( index == 49 )
                        // {
                        //     std::cout << std::endl;
                        //     std::cout << "raos[" << index_3 << "]: " << raos[index_3] << std::endl;
                        //     std::cout << "vel_x_raddif_p0[" << index_2 << "]: " << vel_x_raddif_p0[index_2] << std::endl;
                        //     std::cout << "vel_y_raddif_p0[" << index_2 << "]: " << vel_y_raddif_p0[index_2] << std::endl;
                        //     std::cout << "vel_z_raddif_p0[" << index_2 << "]: " << vel_z_raddif_p0[index_2] << std::endl;
                        //     std::cout << "vel_x_total[" << index << "]: " << vel_x_total[index] << std::endl;
                        //     std::cout << "vel_y_total[" << index << "]: " << vel_y_total[index] << std::endl;
                        //     std::cout << "vel_z_total[" << index << "]: " << vel_z_total[index] << std::endl;
                        // }
                    }

                    if ( index == 49 || index == 196 )
                    {
                        double abc = 0;
                    }
                }
            }
        }

        for ( int i=0; i<vel_x_gp->field_points_np; i++ )
        {
            vel_x_outf << vel_x_total[i].real( ) << "  ";
            vel_x_outf << vel_x_total[i].imag( ) << std::endl;

            vel_y_outf << vel_y_total[i].real( ) << "  ";
            vel_y_outf << vel_y_total[i].imag( ) << std::endl;

            vel_z_outf << vel_z_total[i].real( ) << "  ";
            vel_z_outf << vel_z_total[i].imag( ) << std::endl;
        }

        vel_x_outf.close( );
        vel_y_outf.close( );
        vel_z_outf.close( );

        // for ( int i=0; i<input->heads_np; i++ )
        // {
        //     for ( int j=0; j<vel_x_gp->field_points_nb; j++ )
        //     {
        //         for ( int r=vel_x_gp->field_points_cnp[j]; r<vel_x_gp->field_points_cnp[j+1]; r++ )
        //         {
        //             index                   = i * vel_x_gp->sysmat_nrows + r;
        //             std::cout << endl;
        //             std::cout << "Field Point: " << index << " - " << vel_x_gp->field_points[3*r] << " " << vel_x_gp->field_points[3*r+1] << " " << vel_x_gp->field_points[3*r+2] << std::endl;
        //             std::cout << "Vel_x: " << vel_x_total[index];
        //             std::cout << " - Mag: " << std::abs( vel_x_total[index] );
        //             std::cout << " - Pha: " << 57.3 * std::atan2( vel_x_total[index].imag(), vel_x_total[index].real( ) );
        //             std::cout << std::endl;
        //             std::cout << "Vel_y: " << vel_y_total[index];
        //             std::cout << " - Mag: " << std::abs( vel_y_total[index] );
        //             std::cout << " - Pha: " << 57.3 * std::atan2( vel_y_total[index].imag(), vel_y_total[index].real( ) );
        //             std::cout << std::endl;
        //             std::cout << "Vel_z: " << vel_z_total[index];
        //             std::cout << " - Mag: " << std::abs( vel_z_total[index] );
        //             std::cout << " - Pha: " << 57.3 * std::atan2( vel_z_total[index].imag(), vel_z_total[index].real( ) );
        //             std::cout << std::endl;
        //         }
        //     }
        // }

        // // Write system matrixes to files
        // std::string vel_x_sysmat_fipath( "E:/sergio/0050_OASIS_SM/vel_x_sysmat.dat" );
        // std::string vel_x_sysmat_steady_fipath( "E:/sergio/0050_OASIS_SM/vel_x_sysmat_steady.dat" );

        // std::ofstream of_vel_x_sysmat( vel_x_sysmat_fipath );
        // std::ofstream of_vel_x_sysmat_steady( vel_x_sysmat_steady_fipath );

        // of_vel_x_sysmat << vel_x_gp->sysmat_nrows << " " << vel_x_gp->sysmat_ncols << std::endl;
        // of_vel_x_sysmat_steady << vel_x_gp->sysmat_nrows << " " << vel_x_gp->sysmat_ncols << std::endl;
        // for ( int i=0; i<vel_x_gp->sysmat_nrows; i++ )
        // {
        //     for ( int j=0; j<vel_x_gp->sysmat_ncols; j++ )
        //     {
        //         of_vel_x_sysmat << vel_x_gp->sysmat[i*vel_x_gp->sysmat_ncols+j].real( ) << " " << vel_x_gp->sysmat[i*vel_x_gp->sysmat_ncols+j].imag( ) << std::endl;
        //         of_vel_x_sysmat_steady << vel_x_gp->sysmat_steady[i*vel_x_gp->sysmat_ncols+j].real( ) << " " << vel_x_gp->sysmat_steady[i*vel_x_gp->sysmat_ncols+j].imag( ) << std::endl;
        //     }
        // }

        // of_vel_x_sysmat.close( );
        // of_vel_x_sysmat_steady.close( );

        // // Write results to files
        // std::string vel_x_fipath( "E:/sergio/0050_OASIS_SM/vel_x_data.dat" );
        // std::string vel_y_fipath( "E:/sergio/0050_OASIS_SM/vel_y_data.dat" );
        // std::string vel_z_fipath( "E:/sergio/0050_OASIS_SM/vel_z_data.dat" );

        // std::ofstream of_vel_x( vel_x_fipath );
        // std::ofstream of_vel_y( vel_y_fipath );
        // std::ofstream of_vel_z( vel_z_fipath );

        // std::string space4( "    " );
        // for ( int i=0; i<vel_x_gp->sysmat_nrows; i++ )
        // {
        //     of_vel_x << i+1 << space4;
        //     for ( int j=0; j<3; j++ )
        //     {
        //         of_vel_x << vel_x_gp->field_points[3*i+j] << space4;
        //     }
        //     of_vel_x << vel_x_total[i].real( ) << space4 << vel_x_total[i].imag( ) << space4;
        //     of_vel_x << std::abs( vel_x_total[i] ) << space4 << 57.3 * std::atan2( vel_x_total[i].imag( ), vel_x_total[i].real( ) ) << std::endl;

        //     of_vel_y << i+1 << space4;
        //     for ( int j=0; j<3; j++ )
        //     {
        //         of_vel_y << vel_y_gp->field_points[3*i+j] << space4;
        //     }
        //     of_vel_y << vel_y_total[i].real( ) << space4 << vel_y_total[i].imag( ) << space4;
        //     of_vel_y << std::abs( vel_y_total[i] ) << space4 << 57.3 * std::atan2( vel_y_total[i].imag( ), vel_y_total[i].real( ) ) << std::endl;

        //     of_vel_z << i+1 << space4;
        //     for ( int j=0; j<3; j++ )
        //     {
        //         of_vel_z << vel_z_gp->field_points[3*i+j] << space4;
        //     }
        //     of_vel_z << vel_z_total[i].real( ) << space4 << vel_z_total[i].imag( ) << space4;
        //     of_vel_z << std::abs( vel_z_total[i] ) << space4 << 57.3 * std::atan2( vel_z_total[i].imag( ), vel_z_total[i].real( ) ) << std::endl;
        // }

        // of_vel_x.close( );
        // of_vel_y.close( );
        // of_vel_z.close( );
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