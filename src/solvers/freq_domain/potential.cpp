
// Include local modules
#include "potential.hpp"
#include "../../green/source.hpp"
#include "../../interfaces/hmf_interface.hpp"
#include "../../math/integration.hpp"
#include "../../math/math_interface.hpp"
#include "../../waves.hpp"


void    calculate_influence_potmat_steady(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cusfloat*       field_points,
                                                int             field_points_np,
                                                cuscomplex*     inf_pot_mat
                                            )
{
    // Configure MPI
    int source_start_pos    = 0;
    int source_end_pos      = 0;
    mpi_config->get_1d_bounds( 
                                mesh_gp->source_nodes_tnp, 
                                source_start_pos, 
                                source_end_pos 
                            );
    
    // Generate potential matrix
    int         count       = 0;

    if ( input->is_log_sin_ana )
    {
        // Define local auxiliar variables
        const int   ndim                = 3;
        cusfloat    field_point_0[3];   clear_vector( ndim, field_point_0 );
        cusfloat    field_point_1[3];   clear_vector( ndim, field_point_1 );
        cusfloat    field_point_2[3];   clear_vector( ndim, field_point_2 );
        cusfloat    field_point_3[3];   clear_vector( ndim, field_point_3 );
        cusfloat    field_point_4[3];   clear_vector( ndim, field_point_4 );
        cusfloat    field_point_5[3];   clear_vector( ndim, field_point_5 );
        cusfloat    pot_0               = 0.0;
        cusfloat    pot_1               = 0.0;
        cusfloat    pot_2               = 0.0;
        cusfloat    pot_3               = 0.0;
        cusfloat    pot_4               = 0.0;
        cusfloat    pot_5               = 0.0;
        cuscomplex  pot_term            = 0.0;

        for ( int i=0; i<field_points_np; i++ )
        {
            // Calculate field points for the different radius
            copy_vector( ndim, &(field_points[3*i]), field_point_0 );

            field_point_1[0]    = field_points[3*i];
            field_point_1[1]    = field_points[3*i+1];
            field_point_1[2]    = field_points[3*i+2] + 2.0 * input->water_depth;

            copy_vector( ndim, &(field_points[3*i]), field_point_2 );

            field_point_3[0]    = field_points[3*i];
            field_point_3[1]    = field_points[3*i+1];
            field_point_3[2]    = field_points[3*i+2] + 2.0 * input->water_depth;

            field_point_4[0]    = field_points[3*i];
            field_point_4[1]    = field_points[3*i+1];
            field_point_4[2]    = -field_points[3*i+2] + 2.0 * input->water_depth;

            field_point_5[0]    = field_points[3*i];
            field_point_5[1]    = field_points[3*i+1];
            field_point_5[2]    = field_points[3*i+2] + 4.0 * input->water_depth;

            for ( int j=source_start_pos; j<source_end_pos; j++ )
            {
                // Compute steady and wave terms over the panel
                if ( 
                        mesh_gp->source_nodes[i]->panel->type == DIFFRAC_PANEL_CODE
                        &&
                        mesh_gp->source_nodes[j]->panel->type == DIFFRAC_PANEL_CODE
                    )
                {
                    // Calculate potential contribution for r0
                    calculate_source_potential_newman(
                                                            mesh_gp->panels[j],
                                                            field_point_0,
                                                            0, 
                                                            0,
                                                            pot_0
                                                        );

                    // Calculate potential contribution for r1
                    calculate_source_potential_newman(
                                                            mesh_gp->panels_mirror[j], 
                                                            field_point_1,
                                                            0, 
                                                            0, 
                                                            pot_1
                                                        );

                    // Calculate potential contribution for r2
                    calculate_source_potential_newman(
                                                            mesh_gp->panels_mirror[j], 
                                                            field_point_2,
                                                            0, 
                                                            0, 
                                                            pot_2
                                                        );

                    // Calculate potential contribution for r3
                    calculate_source_potential_newman(
                                                            mesh_gp->panels[j], 
                                                            field_point_3,
                                                            0, 
                                                            0, 
                                                            pot_3
                                                        );

                    // Calculate potential contribution for r4
                    calculate_source_potential_newman(
                                                            mesh_gp->panels_mirror[j], 
                                                            field_point_4,
                                                            0, 
                                                            0, 
                                                            pot_4
                                                        );

                    // Calculate potential contribution for r5
                    calculate_source_potential_newman(
                                                            mesh_gp->panels_mirror[j], 
                                                            field_point_5,
                                                            0, 
                                                            0, 
                                                            pot_5
                                                        );

                    // Calculate total potential contribution
                    pot_term = pot_0 + pot_1 + pot_2 + pot_3 + pot_4 + pot_5;

                    inf_pot_mat[count]  =  -pot_term / 4.0 / PI;
                }
                count++;
            }
        }
    }
    else
    {
        // Define potential funcions objects interface
        GRFInterface*   green_interf_steady = new   GRFInterface(
                                                                    mesh_gp->source_nodes[0],
                                                                    0.0,
                                                                    mesh_gp->source_nodes[0]->panel->center,
                                                                    input->water_depth
                                                                );

        auto            steady_fcn          =   [green_interf_steady]
                                                (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z) -> cuscomplex
                                                {
                                                    return (*green_interf_steady)( xi, eta, x, y, z );
                                                };

        // Define local auxiliar variables
        cuscomplex  pot_steady_term( 0.0, 0.0 );

        // Loop over field points and panels to generate the potential steady mat
        for ( int i=0; i<mesh_gp->source_nodes_tnp; i++ )
        {
            // Change field point
            green_interf_steady->set_field_point( 
                                                    mesh_gp->source_nodes[i]->panel->center
                                                );

            for ( int j=source_start_pos; j<source_end_pos; j++ )
            {
                // Change source point
                green_interf_steady->set_source(
                                                    mesh_gp->source_nodes[j],
                                                    1.0
                                                );
                
                // Compute steady and wave terms over the panel
                if ( 
                        mesh_gp->source_nodes[i]->panel->type == DIFFRAC_PANEL_CODE
                        &&
                        mesh_gp->source_nodes[j]->panel->type == DIFFRAC_PANEL_CODE
                    )
                {
                    pot_steady_term = adaptive_quadrature_panel(
                                                                    mesh_gp->source_nodes[j]->panel,
                                                                    steady_fcn,
                                                                    input->pot_abs_err,
                                                                    input->pot_rel_err,
                                                                    input->is_block_adaption,
                                                                    false,
                                                                    input->gauss_order
                                                                );
                    
                    inf_pot_mat[count] = pot_steady_term / 4.0 / PI;

                }
                count++;
            }
        }
    }
}


void    calculate_influence_potmat(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cusfloat        ang_freq,
                                                cuscomplex*     inf_pot_steady,
                                                cusfloat*       field_points,
                                                int             field_points_np,
                                                cuscomplex*     inf_pot_total
                                    )
{
    // Define potential funcions objects interface
    GWFInterface*   green_interf_wave   = new   GWFInterface(
                                                                mesh_gp->source_nodes[0],
                                                                0.0,
                                                                mesh_gp->source_nodes[0]->panel->center,
                                                                ang_freq,
                                                                input->water_depth,
                                                                input->grav_acc
                                                            );

    auto            wave_fcn            =   [green_interf_wave]
                                            (cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z) -> cuscomplex
                                            {
                                                return (*green_interf_wave)( xi, eta, x, y, z );
                                            };
    
    // Configure MPI
    int source_start_pos    = 0;
    int source_end_pos      = 0;
    mpi_config->get_1d_bounds( 
                                mesh_gp->source_nodes_tnp, 
                                source_start_pos, 
                                source_end_pos 
                            );
    
    // Generate potential matrix
    int         count = 0;
    cuscomplex  pot_wave_term( 0.0, 0.0 );
    for ( int i=0; i<field_points_np; i++ )
    {
        // Change field point
        green_interf_wave->set_field_point(
                                                &(field_points[3*i])
                                            );

        for ( int j=source_start_pos; j<source_end_pos; j++ )
        {
            // Change source point
            green_interf_wave->set_source(
                                                mesh_gp->source_nodes[j],
                                                1.0
                                        );
            
            // Compute steady and wave terms over the panel
            if ( 
                    mesh_gp->panels[i]->type == DIFFRAC_PANEL_CODE
                    &&
                    mesh_gp->panels[j]->type == DIFFRAC_PANEL_CODE
                )
            {
                pot_wave_term           = adaptive_quadrature_panel(
                                                                        mesh_gp->panels[j],
                                                                        wave_fcn,
                                                                        input->pot_abs_err,
                                                                        input->pot_rel_err,
                                                                        input->is_block_adaption,
                                                                        false,
                                                                        input->gauss_order
                                                                    );
                inf_pot_total[count]    = inf_pot_steady[count] + pot_wave_term / 4.0 / PI;

            }
            else
            {
                inf_pot_total[count]    = cuscomplex( 0.0, 0.0 );
            }
            
            count++;
        }
    }

}


void    calculate_potpanel_raddif_lin(
                                                Input*          input,
                                                cuscomplex*     inf_pot_mat,
                                                int             rows_np,
                                                int             cols_np,
                                                int             start_col,
                                                cuscomplex*     sources,
                                                cuscomplex*     panel_pot
                                )
{
    // Loop over RHS to compute all the panel potentials the panels potentials
    cuscomplex  alpha( 1.0, 0.0 );
    cuscomplex  beta( 0.0, 0.0 );
    int         icnx = 1;
    int         icny = 1;
    for ( int i=0; i<( input->dofs_np + input->heads_np ); i++ )
    {
        cblas_gemv<cuscomplex>( 
                                    CblasRowMajor,
                                    CblasNoTrans,
                                    rows_np,
                                    cols_np,
                                    &alpha,
                                    inf_pot_mat,
                                    cols_np,
                                    &(sources[i*rows_np+start_col]),
                                    icnx,
                                    &beta,
                                    &(panel_pot[i*rows_np]),
                                    icny
                                );
    }
}


void    calculate_potpanel_raddif_nlin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     sources,
                                                cusfloat        ang_freq
                                        )
{
    // Calculate MPI data chunks
    int elem_end_pos     = 0;
    int elem_start_pos   = 0;
    mpi_config->get_1d_bounds( 
                                    mesh_gp->panels_tnp, 
                                    elem_start_pos, 
                                    elem_end_pos 
                                );
    
    // Create Function to integrate potential value
    GRFInterface*   green_interf_steady = new   GRFInterface(
                                                                mesh_gp->source_nodes[0],
                                                                sources[0],
                                                                mesh_gp->panels[0]->center,
                                                                input->water_depth
                                                            );
    GWFInterface*   green_interf_wave   = new   GWFInterface(
                                                                mesh_gp->source_nodes[0],
                                                                sources[0],
                                                                mesh_gp->panels[0]->center,
                                                                ang_freq,
                                                                input->water_depth,
                                                                input->grav_acc
                                                            );

    auto steady_fcn =   [green_interf_steady]
                        ( 
                            cusfloat    xi,
                            cusfloat    eta,
                            cusfloat    X,
                            cusfloat    Y,
                            cusfloat    Z
                        ) -> cuscomplex
                        {
                            return (*green_interf_steady)( xi, eta, X, Y, Z );
                        };

    auto wave_fcn   =   [green_interf_wave]
                        ( 
                            cusfloat    xi,
                            cusfloat    eta,
                            cusfloat    X,
                            cusfloat    Y,
                            cusfloat    Z
                        ) -> cuscomplex
                        {
                            return (*green_interf_wave)( xi, eta, X, Y, Z );
                        };

    // Loop over panel to get the radiation potential over them
    cuscomplex panel_potential  = cuscomplex( 0.0, 0.0 );
    cuscomplex pot_i_steady     = cuscomplex( 0.0, 0.0 );
    cuscomplex pot_i_wave       = cuscomplex( 0.0, 0.0 );
    for ( int i=0; i<mesh_gp->panels_tnp; i++ )
    {
        // Set new field point for the calculation of the potential
        green_interf_steady->set_field_point( mesh_gp->panels[i]->center );
        green_interf_wave->set_field_point( mesh_gp->panels[i]->center );

        // Clean panel potential
        panel_potential = cuscomplex( 0.0, 0.0 );

        // Loop over source points to calculate potential
        // over ith panel center
        for ( int j=0; j<mesh_gp->source_nodes_tnp; j++ )
        {
            if ( 
                    mesh_gp->panels[i]->type == DIFFRAC_PANEL_CODE
                    &&
                    mesh_gp->source_nodes[j]->panel->type == DIFFRAC_PANEL_CODE
                )
            {
                // Set current source value
                green_interf_steady->set_source(
                                                mesh_gp->source_nodes[j],
                                                sources[j]
                                            );
                green_interf_wave->set_source(
                                                mesh_gp->source_nodes[j],
                                                sources[j]
                                            );
                
                pot_i_steady    = adaptive_quadrature_panel(
                                                                mesh_gp->source_nodes[j]->panel,
                                                                steady_fcn,
                                                                input->pot_abs_err,
                                                                input->pot_rel_err,
                                                                input->is_block_adaption,
                                                                false,
                                                                input->gauss_order
                                                            );
                pot_i_wave      = adaptive_quadrature_panel(
                                                                mesh_gp->source_nodes[j]->panel,
                                                                wave_fcn,
                                                                input->pot_abs_err,
                                                                input->pot_rel_err,
                                                                input->is_block_adaption,
                                                                false,
                                                                input->gauss_order
                                                            );
                panel_potential +=  ( pot_i_steady + pot_i_wave ) /4.0 / PI;

            }
        }
    }

    // Delete heap allocated memory
    delete green_interf_steady;
    delete green_interf_wave;
}


void    calculate_potpanel_total_lin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cusfloat        ang_freq,
                                                cuscomplex*     intensities,
                                                cuscomplex*     pot_steady_sysmat,
                                                cuscomplex*     pot_sysmat,
                                                cusfloat*       pot_fp,
                                                int*            pot_fp_cnp,
                                                int             pot_fp_nb,
                                                cuscomplex*     raos,
                                                cuscomplex*     potpanel_total
                                    )
{
    // Declare auxiliar variables to use in the function
    int index       =   0;
    int index_2     =   0;
    int index_3     =   0;
    int pot_fp_np   =   pot_fp_cnp[pot_fp_nb-1];

    /***************************************************************/
    /************* Radiation-Diffraction Potential *****************/
    /***************************************************************/

    // Calculate potential influence coeffcients matrix
    calculate_influence_potmat(
                                    input,
                                    mpi_config,
                                    mesh_gp,
                                    ang_freq,
                                    pot_steady_sysmat,
                                    pot_fp,
                                    pot_fp_np,
                                    pot_sysmat
                                );

    // Calculate radiation-diffraction panels potential
    cuscomplex* potpanel_raddif = generate_empty_vector<cuscomplex>( ( input->dofs_np + input->heads_np ) * pot_fp_np );

    int col_start_pos    = 0;
    int col_end_pos      = 0;
    mpi_config->get_1d_bounds( 
                                mesh_gp->source_nodes_tnp, 
                                col_start_pos, 
                                col_end_pos 
                            );
    
    calculate_potpanel_raddif_lin(
                                    input,
                                    pot_sysmat,
                                    pot_fp_np,
                                    mesh_gp->source_nodes_tnp,
                                    col_start_pos,
                                    intensities,
                                    potpanel_raddif
                                );

    /***************************************************************/
    /*************** Add Inicident Wave Potential ******************/
    /***************************************************************/

    // Calculate incident wave potential
    cusfloat    k   = w2k( 
                                ang_freq,
                                input->water_depth,
                                input->grav_acc
                            );
    
    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<pot_fp_np; j++ )
        {
            index                   =   pot_fp_np * i + j;
            potpanel_total[index]   =   wave_potential_airy_space(
                                                                    input->wave_amplitude,
                                                                    ang_freq,
                                                                    k,
                                                                    input->water_depth,
                                                                    input->grav_acc,
                                                                    pot_fp[3*j],
                                                                    pot_fp[3*j+1],
                                                                    pot_fp[3*j+2],
                                                                    input->heads[i]
                                                                );
        }
    }

    /***************************************************************/
    /************** Add Diffraction Wave Potential *****************/
    /***************************************************************/

    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<pot_fp_np; j++ )
        {
            index                   = pot_fp_np * i + j;
            index_2                 = ( input->dofs_np + i ) * pot_fp_np + j;
            potpanel_total[index]   = potpanel_raddif[index_2];
        }
    }

    /***************************************************************/
    /*************** Add Radiation Wave Potential *****************/
    /***************************************************************/

    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<pot_fp_nb; j++ )
        {
            for ( int k=0; k<input->dofs_np; k++ )
            {
                for ( int r=pot_fp_cnp[j]; r<pot_fp_cnp[j+1]; r++ )
                {
                    index                   = i * pot_fp_np + r;
                    index_2                 = k * pot_fp_np + r;
                    index_3                 = i * ( input->dofs_np * pot_fp_nb ) + j * input->dofs_np + k;
                    potpanel_total[index]   += raos[index_3] * potpanel_raddif[index_2];
                }
            }
        }
    }

    // Delete heap memory allocated for the current function
    mkl_free( potpanel_raddif );
}