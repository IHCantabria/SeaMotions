
// Include general usage libraries
#include <sstream>

// Include local modules
#include "gf_intensities_v2.hpp"
#include "../../math/integration.hpp"
#include "../../green/source.hpp"
#include "../../green/pulsating_fin_depth_v2.hpp"
#include "../../waves/wave_dispersion_base_fo.hpp"
#include "../../waves/waves_common.hpp"


void    calculate_gf_intensity_sysmat(
                                                    Input*                      input,
                                                    MpiConfig*                  mpi_config,
                                                    SclCmpx*                    scl,
                                                    MeshGroup*                  mesh_gp,
                                                    GWFcnsInterfaceT<NUM_GP2>&  gwf_interf,
                                                    cusfloat                    w,
                                                    MLGCmpx*                    sources_mlg,
                                                    MLGCmpx*                    potential_mlg,
                                                    cuscomplex*                 sources_int
                                   )
{
    /***************************************/
    /******** Fill system matrix  **********/
    /***************************************/
    // Create auxiliar lambda function
    // auto wave_fcn   =   [gwf_interf]
    //                     ( 
    //                         cusfloat xi, 
    //                         cusfloat eta,
    //                         cusfloat x, 
    //                         cusfloat y, 
    //                         cusfloat z 
    //                     )
    //                     {
    //                         return (*gwf_interf)( xi, eta, x, y, z );
    //                     };

    // Clean system matrixes
    clear_vector( sources_mlg->sysmat_nrows*sources_mlg->sysmat_ncols, sources_mlg->sysmat );
    clear_vector( potential_mlg->sysmat_nrows*potential_mlg->sysmat_ncols, potential_mlg->sysmat );
    
    // Loop over panels to integrate value
    int         col_count       = 0;
    cusfloat    dist            = 0.0;
    int         index_sources   = 0;
    int         index_pot       = 0;
    cuscomplex  int_value( 0.0, 0.0 );
    cuscomplex  int_dn_value( 0.0, 0.0 );
    cuscomplex  wave_fcn_value( 0.0, 0.0 );
    cuscomplex  wave_fcn_dn_value( 0.0, 0.0 );
    PanelGeom*  panel_j         = nullptr;
    SourceNode* source_i        = nullptr;
    int         row_count       = 0;
    cusfloat    nu              = pow2s( w ) / input->grav_acc;

    // std::stringstream ss_source_steady;
    // std::stringstream ss_source_wave;
    // std::stringstream ss_source_total;

    // std::stringstream ss_pot_steady;
    // std::stringstream ss_pot_wave;
    // std::stringstream ss_pot_total;

    // ss_source_steady << "sysmat_source_steady_" << mpi_config->proc_rank << ".dat";
    // ss_source_wave << "sysmat_source_wave_" << mpi_config->proc_rank << ".dat";
    // ss_source_total << "sysmat_source_total_" << mpi_config->proc_rank << ".dat";

    // ss_pot_steady << "sysmat_pot_steady_" << mpi_config->proc_rank << ".dat";
    // ss_pot_wave << "sysmat_pot_wave_" << mpi_config->proc_rank << ".dat";
    // ss_pot_total << "sysmat_pot_total_" << mpi_config->proc_rank << ".dat";

    // std::ofstream out_steady( ss_source_steady.str( ) );
    // std::ofstream out_wave( ss_source_wave.str( ) );
    // std::ofstream out_total( ss_source_total.str( ) );
    // std::ofstream out_pot_steady( ss_pot_steady.str( ) );
    // std::ofstream out_pot_wave( ss_pot_wave.str( ) );
    // std::ofstream out_pot_total( ss_pot_total.str( ) );

    // out_steady << sources_mlg->sysmat_nrows << " " << sources_mlg->sysmat_ncols << std::endl;
    // out_wave << sources_mlg->sysmat_nrows << " " << sources_mlg->sysmat_ncols << std::endl;
    // out_total << sources_mlg->sysmat_nrows << " " << sources_mlg->sysmat_ncols << std::endl;
    // out_pot_steady << sources_mlg->sysmat_nrows << " " << sources_mlg->sysmat_ncols << std::endl;
    // out_pot_wave << sources_mlg->sysmat_nrows << " " << sources_mlg->sysmat_ncols << std::endl;
    // out_pot_total << sources_mlg->sysmat_nrows << " " << sources_mlg->sysmat_ncols << std::endl;

    for ( int i=scl->start_col_0; i<scl->end_col_0; i++ )
    {
        // Get memory address of the ith panel
        source_i = mesh_gp->source_nodes[i];
        gwf_interf.set_source_i( source_i, 1.0 );

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
        {
            // Get memory address of the panel jth
            panel_j = mesh_gp->source_nodes[j]->panel;
            gwf_interf.set_source_j( mesh_gp->source_nodes[j] );
            
            // Calculate distance in between field point and source
            dist =   l2_norm( 3, source_i->panel->center, mesh_gp->source_nodes[j]->position );
            dist /=  input->water_depth;
            
            int_value = 0.0;
            int_dn_value = 0.0;
            // Integrate green function normal derivative along the current panel
            if ( i==72 && j == 72 )
            {
                double a = 0.0;
            }

            if ( i == j )
            {
                if ( panel_j->type == DIFFRAC_PANEL_CODE )
                {
                    int_dn_value    = -cuscomplex( 0.5, 0.0 );
                    if ( dist > 1 )
                    {
                        quadrature_panel_t<
                                                PanelGeom, 
                                                GWFcnsInterfaceT<NUM_GP2>, 
                                                john_series<NUM_GP2, G_ON, DGDR_OFF, DGDZ_OFF>, 
                                                NUM_GP
                                            >( 
                                                source_i->panel, 
                                                gwf_interf, 
                                                wave_fcn_value,
                                                wave_fcn_dn_value
                                            );
                    }
                    else
                    {
                        quadrature_panel_t<
                                                PanelGeom, 
                                                GWFcnsInterfaceT<NUM_GP2>, 
                                                wave_term_integral<NUM_GP2, G_ON, DGDR_OFF, DGDZ_OFF, FSLID_OFF>, 
                                                NUM_GP
                                            >( 
                                                source_i->panel, 
                                                gwf_interf, 
                                                wave_fcn_value,
                                                wave_fcn_dn_value
                                            );
                    }
                    int_value       =   wave_fcn_value / 4.0 / PI;
                }
                else if ( panel_j->type == LID_PANEL_CODE )
                {
                    int_dn_value    = -cuscomplex( 1.0, 0.0 );
                    if ( dist > 1 )
                    {
                        quadrature_panel_t<
                                                PanelGeom, 
                                                GWFcnsInterfaceT<NUM_GP2>, 
                                                john_series<NUM_GP2, G_ON, DGDR_OFF, DGDZ_OFF>, 
                                                NUM_GP
                                            >( 
                                                source_i->panel, 
                                                gwf_interf, 
                                                wave_fcn_value,
                                                wave_fcn_dn_value
                                            );
                    }
                    else
                    {
                        quadrature_panel_t<
                                                PanelGeom, 
                                                GWFcnsInterfaceT<NUM_GP2>, 
                                                wave_term_integral<NUM_GP2, G_ON, DGDR_OFF, DGDZ_OFF, FSLID_ON>, 
                                                NUM_GP
                                            >( 
                                                source_i->panel, 
                                                gwf_interf, 
                                                wave_fcn_value,
                                                wave_fcn_dn_value
                                            );
                        wave_fcn_value  += cuscomplex( nu * source_i->panel->free_surface_log_int, 0.0 );
                    }
                    int_value = wave_fcn_value / 4.0 / PI;
                }
                
            }
            else
            {
                wave_fcn_value      = 0.0;
                wave_fcn_dn_value   = 0.0;

                if ( dist > 1 )
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
                                            wave_fcn_dn_value
                                        );
                }
                else
                {
                    quadrature_panel_t<
                                            PanelGeom, 
                                            GWFcnsInterfaceT<NUM_GP2>, 
                                            wave_term_integral<NUM_GP2, G_ON, DGDR_ON, DGDZ_ON, FSLID_OFF>, 
                                            NUM_GP
                                        >( 
                                            source_i->panel, 
                                            gwf_interf, 
                                            wave_fcn_value,
                                            wave_fcn_dn_value
                                        );
                }
                
                int_value       =   wave_fcn_value / 4.0 / PI;
                int_dn_value    =   wave_fcn_dn_value / 4.0 / PI;
                if ( 
                        panel_j->type == LID_PANEL_CODE
                        &&
                        source_i->panel->type == LID_PANEL_CODE
                    )
                {
                    int_value       = 0.0;
                    int_dn_value    = - int_dn_value;
                }
            }

            // Apply the integral value accordingly
            index_sources   = col_count*scl->num_rows_local+row_count;
            index_pot       = row_count*scl->num_cols_local+col_count;
            if ( dist > 1 )
            {
                potential_mlg->sysmat[index_pot]    = int_value;
                sources_mlg->sysmat[index_sources]  = int_dn_value;
            }
            else
            {
                potential_mlg->sysmat[index_pot]    = potential_mlg->sysmat_steady[index_pot] + int_value;
                sources_mlg->sysmat[index_sources]  = sources_mlg->sysmat_steady[index_sources] + int_dn_value;
            }

            // out_steady << sources_mlg->sysmat_steady[index_sources].real( ) << " " << sources_mlg->sysmat_steady[index_sources].imag( ) << std::endl;
            // out_wave << int_dn_value.real( ) << " " << int_dn_value.imag( ) << std::endl;
            // out_total << sources_mlg->sysmat[index_sources].real( ) << " " << sources_mlg->sysmat[index_sources].imag( ) << std::endl;

            // out_pot_steady << potential_mlg->sysmat_steady[index_pot].real( ) << " " << potential_mlg->sysmat_steady[index_pot].imag( ) << std::endl;
            // out_pot_wave << int_value.real( ) << " " << int_value.imag( ) << std::endl;
            // out_pot_total << potential_mlg->sysmat[index_pot].real( ) << " " << potential_mlg->sysmat[index_pot].imag( ) << std::endl;

            // Advance row count
            row_count++;
        }

        // Advance column count
        col_count++;
    }
    // out_steady.close( );
    // out_wave.close( );
    // out_total.close( );

    // out_pot_steady.close( );
    // out_pot_wave.close( );
    // out_pot_total.close( );

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
                                            -1i
                                            *
                                            w
                                            *
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

    // std::stringstream ss_rhs;
    // ss_rhs << "rhs_" << mpi_config->proc_rank << ".dat";
    // std::ofstream out_rhs( ss_rhs.str( ) );
    // int count_1 = 0;
    // for ( int i=0; i<input->dofs_np+input->heads_np; i++ )
    // {
    //     for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
    //     {
            
    //         out_rhs << sources_int[count_1].real( ) << " " << sources_int[count_1].imag( ) << "\n";
    //         count_1++;
    //     }
    // }

    // Solve system of equations
    scl->Solve( sources_mlg->sysmat, sources_int );

    // std::stringstream ss_sources;
    // ss_sources << "sources_" << mpi_config->proc_rank << ".dat";
    // std::ofstream out_sources( ss_sources.str( ) );
    // int count_2 = 0;
    // for ( int i=0; i<input->dofs_np+input->heads_np; i++ )
    // {
    //     for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
    //     {
            
    //         out_sources << sources_int[count_2].real( ) << " " << sources_int[count_2].imag( ) << "\n";
    //         count_2++;
    //     }
    // }

    // if ( mpi_config->is_root( ) )
    // {
    //     int count_3 = 0;
    //     for ( int i=0; i<1; i++ )
    //     {
    //         for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
    //         {
                
    //             std::cout << std::abs( sources_int[count_3] ) << " " << sources_int[count_3].real( ) << " " << sources_int[count_3].imag( ) << "\n";
    //             count_3++;
    //         }
    //     }
    // }
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


// void    calculate_gf_intensity_steady_sysmat_nlin(
//                                                     Input*          input,
//                                                     SclCmpx*        scl,
//                                                     MeshGroup*      mesh_gp,
//                                                     GRFDnInterface* grf_interf,
//                                                     cuscomplex*     sysmat
//                                                 )
// {
//     /***************************************/
//     /******** Fill system matrix  **********/
//     /***************************************/
//     // Loop over panels to integrate value
//     int         col_count   = 0;
//     cuscomplex  int_value( 0.0, 0.0 );
//     int         row_count   = 0;

    
//     // Create auxiliar lambda function
//     auto steady_fcn =   [grf_interf]
//                         ( 
//                             cusfloat xi, 
//                             cusfloat eta,
//                             cusfloat x, 
//                             cusfloat y, 
//                             cusfloat z 
//                         )
//                         {
//                             return (*grf_interf)( xi, eta, x, y, z );
//                         };

//     // Declare local variables to work with the adaptive
//     // integration scheme
//     SourceNode*     source_i    = nullptr;
    
//     for ( int i=scl->start_col_0; i<scl->end_col_0; i++ )
//     {
//         // Get memory address of the ith panel
//         source_i = mesh_gp->source_nodes[i];
//         grf_interf->set_source_i( source_i );

//         // Loop over rows to calcualte the influence of the panel
//         // over each collocation point
//         row_count = 0;
//         for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
//         {
//             // Get memory address of the panel jth
//             grf_interf->set_source_j( mesh_gp->source_nodes[j] );

//             // Integrate green function normal derivative along the current panel
//             if ( i == j )
//             {
//                 int_value       = cuscomplex( 0.0, 0.0 );
//             }
//             else
//             {
//                 int_value       =   adaptive_quadrature_panel(
//                                                                 source_i->panel,
//                                                                 steady_fcn,
//                                                                 input->gfdn_abs_err,
//                                                                 input->gfdn_rel_err,
//                                                                 input->is_block_adaption,
//                                                                 false,
//                                                                 input->gauss_order
//                                                             );
//                 int_value       =   int_value / 4.0 / PI;
//             }

//             if ( 
//                     mesh_gp->source_nodes[j]->panel->type == LID_PANEL_CODE
//                     &&
//                     source_i->panel->type == LID_PANEL_CODE
//                 )
//             {
//                 int_value       = - int_value;
//             }

//             sysmat[col_count*scl->num_rows_local+row_count] = int_value;

//             // Advance row count
//             row_count++;
//         }

//         // Advance column count
//         col_count++;
//     }
// }