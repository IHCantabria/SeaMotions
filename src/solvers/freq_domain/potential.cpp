
// Include local modules
#include "../../containers/matlin_group.hpp"
#include "../../green/source.hpp"
#include "../../interfaces/hmf_interface.hpp"
#include "../../math/integration.hpp"
#include "../../math/math_interface.hpp"
#include "potential.hpp"
#include "../../waves.hpp"


void    calculate_influence_potmat_steady(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                MLGCmpx*        pot_gp
                                            )
{
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

        for ( int i=0; i<pot_gp->field_points_np; i++ )
        {
            // Calculate field points for the different radius
            copy_vector( ndim, &(pot_gp->field_points[3*i]), field_point_0 );

            field_point_1[0]    = pot_gp->field_points[3*i];
            field_point_1[1]    = pot_gp->field_points[3*i+1];
            field_point_1[2]    = pot_gp->field_points[3*i+2] + 2.0 * input->water_depth;

            copy_vector( ndim, &(pot_gp->field_points[3*i]), field_point_2 );

            field_point_3[0]    = pot_gp->field_points[3*i];
            field_point_3[1]    = pot_gp->field_points[3*i+1];
            field_point_3[2]    = pot_gp->field_points[3*i+2] + 2.0 * input->water_depth;

            field_point_4[0]    = pot_gp->field_points[3*i];
            field_point_4[1]    = pot_gp->field_points[3*i+1];
            field_point_4[2]    = -pot_gp->field_points[3*i+2] + 2.0 * input->water_depth;

            field_point_5[0]    = pot_gp->field_points[3*i];
            field_point_5[1]    = pot_gp->field_points[3*i+1];
            field_point_5[2]    = pot_gp->field_points[3*i+2] + 4.0 * input->water_depth;

            for ( int j=pot_gp->start_col; j<pot_gp->end_col; j++ )
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

                    pot_gp->sysmat_steady[count]  =  -pot_term / 4.0 / PI;
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
                                                    mesh_gp->panels[i]->center
                                                );

            for ( int j=pot_gp->start_col; j<pot_gp->end_col; j++ )
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
                                                                    mesh_gp->panels[j],
                                                                    steady_fcn,
                                                                    input->pot_abs_err,
                                                                    input->pot_rel_err,
                                                                    input->is_block_adaption,
                                                                    false,
                                                                    input->gauss_order
                                                                );
                    
                    pot_gp->sysmat_steady[count] = pot_steady_term / 4.0 / PI;

                }
                count++;
            }
        }
    }
}


void    calculate_influence_potmat(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                cusfloat        ang_freq,
                                                MLGCmpx*        pot_gp
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

    // Generate potential matrix
    int         count = 0;
    cuscomplex  pot_wave_term( 0.0, 0.0 );
    for ( int i=0; i<pot_gp->field_points_np; i++ )
    {
        // Change field point
        green_interf_wave->set_field_point(
                                                &(pot_gp->field_points[3*i])
                                            );

        for ( int j=pot_gp->start_col; j<pot_gp->end_col; j++ )
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
                pot_gp->sysmat[count]   = pot_gp->sysmat_steady[count] + pot_wave_term / 4.0 / PI;

            }
            else
            {
                pot_gp->sysmat[count]   = cuscomplex( 0.0, 0.0 );
            }
            
            count++;
        }
    }

}


void    calculate_potpanel_raddif_nlin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     intensities,
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
                                                                intensities[0],
                                                                mesh_gp->panels[0]->center,
                                                                input->water_depth
                                                            );
    GWFInterface*   green_interf_wave   = new   GWFInterface(
                                                                mesh_gp->source_nodes[0],
                                                                intensities[0],
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
                                                intensities[j]
                                            );
                green_interf_wave->set_source(
                                                mesh_gp->source_nodes[j],
                                                intensities[j]
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
                                                cuscomplex*     raos,
                                                MLGCmpx*        pot_gp,
                                                cuscomplex*     potpanel_total
                                    )
{
    // Declare auxiliar variables to use in the function
    int index       =   0;
    int index_2     =   0;
    int index_3     =   0;

    /***************************************************************/
    /************* Radiation-Diffraction Potential *****************/
    /***************************************************************/

    // Calculate potential influence coeffcients matrix
    calculate_influence_potmat(
                                    input,
                                    mesh_gp,
                                    ang_freq,
                                    pot_gp
                                );

    // Calculate radiation-diffraction panels potential
    int col_start_pos    = 0;
    int col_end_pos      = 0;
    mpi_config->get_1d_bounds( 
                                mesh_gp->source_nodes_tnp, 
                                col_start_pos, 
                                col_end_pos 
                            );
    
    calculate_potpanel_raddif_lin(
                                    input,
                                    intensities,
                                    pot_gp
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
        for ( int j=0; j<pot_gp->field_points_np; j++ )
        {
            index                   =   pot_gp->field_points_np * i + j;
            potpanel_total[index]   =   wave_potential_airy_space(
                                                                    input->wave_amplitude,
                                                                    ang_freq,
                                                                    k,
                                                                    input->water_depth,
                                                                    input->grav_acc,
                                                                    pot_gp->field_points[3*j],
                                                                    pot_gp->field_points[3*j+1],
                                                                    pot_gp->field_points[3*j+2],
                                                                    input->heads[i]
                                                                );
        }
    }

    /***************************************************************/
    /************** Add Diffraction Wave Potential *****************/
    /***************************************************************/

    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<pot_gp->field_points_np; j++ )
        {
            index                   = pot_gp->field_points_np * i + j;
            index_2                 = ( input->dofs_np + i ) * pot_gp->field_points_np + j;
            potpanel_total[index]   = pot_gp->field_values[index_2];
        }
    }

    /***************************************************************/
    /*************** Add Radiation Wave Potential *****************/
    /***************************************************************/

    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<pot_gp->field_points_nb; j++ )
        {
            for ( int k=0; k<input->dofs_np; k++ )
            {
                for ( int r=pot_gp->field_points_cnp[j]; r<pot_gp->field_points_cnp[j+1]; r++ )
                {
                    index                   = i * pot_gp->field_points_np + r;
                    index_2                 = k * pot_gp->field_points_np + r;
                    index_3                 = i * ( input->dofs_np * pot_gp->field_points_nb ) + j * input->dofs_np + k;
                    potpanel_total[index]   += raos[index_3] * pot_gp->field_values[index_2];
                }
            }
        }
    }

}