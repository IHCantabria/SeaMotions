
// Include local modules
#include "potential.hpp"

#include "../../containers/matlin_group.hpp"
#include "../../green/source.hpp"
#include "../../interfaces/hmf_interface.hpp"
#include "../../math/integration.hpp"
#include "../../math/math_interface.hpp"
#include "tools.hpp"
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

    // Delete allocated heap memory data for the
    // current function
    delete green_interf_wave;

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
    calculate_fields_raddif_lin(
                                    input,
                                    intensities,
                                    pot_gp
                                );

    cuscomplex value( 0.0, 0.0 );
    int ngpf    = input->gauss_np_factor_1d( );
    int _idx0   = 0;
    int _idx1   = 0;
    // for ( int i=0; i<pot_gp->field_points_np/ngpf; i++ )
    // {
    //     for ( int j=0; j<ngpf; j++ )
    //     {
    //         for ( int k=0; k<pot_gp->sysmat_ncols; k++ )
    //         {
    //             _idx0 = i * ngpf * pot_gp->sysmat_ncols + k;
    //             _idx1 = i * ngpf * pot_gp->sysmat_ncols + j * pot_gp->sysmat_ncols + k;
    //             value = pot_gp->sysmat[_idx0] - pot_gp->sysmat[_idx1];

    //             if ( i == 6 && j == 3 )
    //             {
    //                 double ccc = 0.0;
    //             }

    //             if ( !assert_complex_equality( value, cuscomplex( 0.0, 0.0 ), 1e-3 ) )
    //             {
    //                 std::cout << "pot_gp.sysmat not equal at - i: " << i << " - j: " << j << " - k: " << k << " - idx0: " << _idx0 << " - idx1: " << _idx1 << " - Value: " << value << std::endl;
    //                 throw std::runtime_error( "" );
    //             }
    //         }
    //     }
    // }

    /***************************************************************/
    /******** Sum panel potentials from all processes **************/
    /***************************************************************/
    cuscomplex* pot_raddif_p0  = nullptr;

    if ( mpi_config->is_root( ) )
    {
        pot_raddif_p0   = generate_empty_vector<cuscomplex>( pot_gp->fields_np * pot_gp->sysmat_nrows );
    }

    MPI_Reduce(
                    pot_gp->field_values,
                    pot_raddif_p0,
                    pot_gp->fields_np * pot_gp->sysmat_nrows,
                    mpi_cuscomplex,
                    MPI_SUM,
                    mpi_config->proc_root,
                    MPI_COMM_WORLD
                );

    /***************************************************************/
    /*************** Add Inicident Wave Potential ******************/
    /***************************************************************/

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
                
                if ( index == 19 )
                {
                    double ab = 0.0;
                }
            }
        }
    }

    /***************************************************************/
    /************** Add Diffraction Wave Potential *****************/
    /***************************************************************/

    if ( mpi_config->is_root( ) )
    {
        for ( int i=0; i<input->heads_np; i++ )
        {
            for ( int j=0; j<pot_gp->field_points_np; j++ )
            {
                index                   = pot_gp->field_points_np * i + j;
                index_2                 = ( input->dofs_np + i ) * pot_gp->field_points_np + j;
                potpanel_total[index]   += pot_raddif_p0[index_2];

                if ( index == 19 )
                {
                    double a = 0.0;
                }
            }
        }
    }

    /***************************************************************/
    /*************** Add Radiation Wave Potential *****************/
    /***************************************************************/

    if ( mpi_config->is_root( ) )
    {
        int idx = 0;
        // for ( int i=0; i<input->dofs_np; i++ )
        // {
        //     for ( int j=0; j<pot_gp->field_points_np; j++ )
        //     {
        //         idx = i*pot_gp->sysmat_nrows+j;
        //         std::cout << "Raddiation Dof: " << i;
        //         std::cout << " - Panel Centre: " << pot_gp->field_points[3*j] << " " << pot_gp->field_points[3*j+1] << " " << pot_gp->field_points[3*j+2];
        //         std::cout << " - Panel Potential: " << pot_raddif_p0[idx];
        //         std::cout << " - Mag: " << std::abs( pot_raddif_p0[idx] );
        //         std::cout << " - Arg: " << 57.3 * std::atan2( pot_raddif_p0[idx].imag( ), pot_raddif_p0[idx].real( ) );
        //         std::cout << " - Index: " << idx;
        //         std::cout << std::endl;
        //     }
        // }

        // for ( int i=0; i<input->heads_np; i++ )
        // {
        //     for ( int j=0; j<pot_gp->field_points_np; j++ )
        //     {
        //         idx = (input->dofs_np+i)*pot_gp->sysmat_nrows+j;
        //         std::cout << "Diffraction Heading: " << i;
        //         std::cout << " - Panel Centre: " << pot_gp->field_points[3*j] << " " << pot_gp->field_points[3*j+1] << " " << pot_gp->field_points[3*j+2];
        //         std::cout << " - Panel Potential: " << pot_raddif_p0[idx];
        //         std::cout << " - Mag: " << std::abs( pot_raddif_p0[idx] );
        //         std::cout << " - Arg: " << 57.3 * std::atan2( pot_raddif_p0[idx].imag( ), pot_raddif_p0[idx].real( ) );
        //         std::cout << " - Index: " << idx;
        //         std::cout << std::endl;
        //     }
        // }

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
                        cuscomplex raos_mult    = raos[index_3] * pot_raddif_p0[index_2];
                        potpanel_total[index]   += raos_mult;
                        std::cout << "index: " << index << " - index_2: " << index_2 << " - index_3: " << index_3;
                        std::cout << " - Raos " << raos[index_3] << " - PotRadiff: " << pot_raddif_p0[index_2];
                        std::cout << " - RaosMultp: " << raos_mult << std::endl;
                    }

                    if ( index == 19 )
                    {
                        double ab = 0.0;
                    }
                }
            }
        }
        // for ( int i=0; i<input->heads_np; i++ )
        // {
        //     for ( int j=0; j<pot_gp->field_points_np; j++ )
        //     {
        //         idx = i*pot_gp->sysmat_nrows+j;
        //         std::cout << "Potpanel total: " << i;
        //         std::cout << " - Panel Centre: " << pot_gp->field_points[3*j] << " " << pot_gp->field_points[3*j+1] << " " << pot_gp->field_points[3*j+2];
        //         std::cout << " - Panel Potential: " << potpanel_total[idx];
        //         std::cout << " - Mag: " << std::abs( potpanel_total[idx] );
        //         std::cout << " - Arg: " << 57.3 * std::atan2( potpanel_total[idx].imag( ), potpanel_total[idx].real( ) );
        //         std::cout << " - Index: " << idx;
        //         std::cout << std::endl;
        //     }
        // }

        // Write results to files
        std::string pot_wl_fipath( "E:/sergio/0050_OASIS_SM/pot_wl_data.dat" );

        std::ofstream of_pot_wl( pot_wl_fipath );

        std::string space4( "    " );
        for ( int i=0; i<pot_gp->sysmat_nrows; i++ )
        {
            of_pot_wl << i+1 << space4;
            for ( int j=0; j<3; j++ )
            {
                of_pot_wl << pot_gp->field_points[3*i+j] << space4;
            }
            of_pot_wl << potpanel_total[i].real( ) << space4 << potpanel_total[i].imag( ) << space4;
            of_pot_wl << std::abs( potpanel_total[i] ) << space4 << 57.3 * std::atan2( potpanel_total[i].imag( ), potpanel_total[i].real( ) ) << std::endl;

        }

        of_pot_wl.close( );
    }


    /*******************************************************/
    /**************  Deallocate heap memory ****************/
    /*******************************************************/
    
    if ( mpi_config->is_root( ) )
    {
        mkl_free( pot_raddif_p0 );
    }

}