
// Include general usage libraries
#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision
#include <sstream>

// Include local modules
#include "gf_intensities_v3.hpp"
#include "../../math/integration.hpp"
#include "../../math/random.hpp"
#include "../../green/source.hpp"
#include "../../green/pulsating_fin_depth_v2.hpp"
#include "../../waves/wave_dispersion_base_fo.hpp"
#include "../../waves/waves_common.hpp"


void    calculate_gf_intensity_sysmat_v3(
                                                    Input*                      input,
                                                    MpiConfig*                  mpi_config,
                                                    SclCmpx*                    scl,
                                                    MeshGroup*                  mesh_gp,
                                                    GWFcnsInterfaceT<NUM_GP2>&  gwf_interf,
                                                    cusfloat                    w,
                                                    MLGCmpx*                    source_form_mlg,
                                                    MLGCmpx*                    pot_form_mlg,
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
    clear_vector( source_form_mlg->sysmat_nrows*source_form_mlg->sysmat_ncols, source_form_mlg->sysmat );
    clear_vector( pot_form_mlg->sysmat_nrows*pot_form_mlg->sysmat_ncols, pot_form_mlg->sysmat );
    clear_vector( potential_mlg->sysmat_nrows*potential_mlg->sysmat_ncols, potential_mlg->sysmat );
    
    // Loop over panels to integrate value
    int         col_count       = 0;
    cusfloat    dist            = 0.0;
    cusfloat    distn           = 0.0;
    int         index_sources   = 0;
    int         index_pot       = 0;
    bool        is_john         = false;
    cuscomplex  int_value( 0.0, 0.0 );
    cuscomplex  int_dn_sf_value( 0.0, 0.0 );
    cuscomplex  int_dn_pf_value( 0.0, 0.0 );
    cuscomplex  wave_fcn_value( 0.0, 0.0 );
    cuscomplex  wave_fcn_dn_sf_value( 0.0, 0.0 );
    cuscomplex  wave_fcn_dn_pf_value( 0.0, 0.0 );
    PanelGeom*  panel_j         = nullptr;
    SourceNode* source_i        = nullptr;
    int         row_count       = 0;
    cusfloat    nu              = pow2s( w ) / input->grav_acc;
    cusfloat    log_sing_val    = 0.0;

    cuscomplex* potential_rhs   = generate_empty_vector<cuscomplex>( potential_mlg->sysmat_nrows * ( input->dofs_np + input->heads_np ) );
    int         dofs_offset     = input->dofs_np * potential_mlg->sysmat_nrows;

    // std::stringstream ss_source_form_steady;
    // std::stringstream ss_source_form_wave;
    // std::stringstream ss_source_form_total;

    // std::stringstream ss_pot_form_steady;
    // std::stringstream ss_pot_form_wave;
    // std::stringstream ss_pot_form_total;

    // std::stringstream ss_pot_steady;
    // std::stringstream ss_pot_wave;
    // std::stringstream ss_pot_total;

    // ss_source_form_steady << "sysmat_source_form_steady_" << mpi_config->proc_rank << "_" << w << ".dat";
    // ss_source_form_wave << "sysmat_source_form_wave_" << mpi_config->proc_rank << "_" << w << ".dat";
    // ss_source_form_total << "sysmat_source_form_total_" << mpi_config->proc_rank << "_" << w << ".dat";

    // ss_pot_form_steady << "sysmat_pot_form_steady_" << mpi_config->proc_rank << "_" << w << ".dat";
    // ss_pot_form_wave << "sysmat_pot_form_wave_" << mpi_config->proc_rank << "_" << w << ".dat";
    // ss_pot_form_total << "sysmat_pot_form_total_" << mpi_config->proc_rank << "_" << w << ".dat";

    // ss_pot_steady << "sysmat_pot_steady_" << mpi_config->proc_rank << "_" << w << ".dat";
    // ss_pot_wave << "sysmat_pot_wave_" << mpi_config->proc_rank << "_" << w << ".dat";
    // ss_pot_total << "sysmat_pot_total_" << mpi_config->proc_rank << "_" << w << ".dat";

    // std::ofstream out_source_form_steady( ss_source_form_steady.str( ) );
    // std::ofstream out_source_form_wave( ss_source_form_wave.str( ) );
    // std::ofstream out_source_form_total( ss_source_form_total.str( ) );

    // std::ofstream out_pot_form_steady( ss_pot_form_steady.str( ) );
    // std::ofstream out_pot_form_wave( ss_pot_form_wave.str( ) );
    // std::ofstream out_pot_form_total( ss_pot_form_total.str( ) );

    // std::ofstream out_pot_steady( ss_pot_steady.str( ) );
    // std::ofstream out_pot_wave( ss_pot_wave.str( ) );
    // std::ofstream out_pot_total( ss_pot_total.str( ) );

    // out_source_form_steady << source_form_mlg->sysmat_nrows << " " << source_form_mlg->sysmat_ncols << "\n";
    // out_source_form_wave << source_form_mlg->sysmat_nrows << " " << source_form_mlg->sysmat_ncols << "\n";
    // out_source_form_total << source_form_mlg->sysmat_nrows << " " << source_form_mlg->sysmat_ncols << "\n";

    // out_pot_form_steady << pot_form_mlg->sysmat_nrows << " " << pot_form_mlg->sysmat_ncols << "\n";
    // out_pot_form_wave << pot_form_mlg->sysmat_nrows << " " << pot_form_mlg->sysmat_ncols << "\n";
    // out_pot_form_total << pot_form_mlg->sysmat_nrows << " " << pot_form_mlg->sysmat_ncols << "\n";

    // out_pot_steady << source_form_mlg->sysmat_nrows << " " << source_form_mlg->sysmat_ncols << "\n";
    // out_pot_wave << source_form_mlg->sysmat_nrows << " " << source_form_mlg->sysmat_ncols << "\n";
    // out_pot_total << source_form_mlg->sysmat_nrows << " " << source_form_mlg->sysmat_ncols << "\n";
    
    // std::ofstream out_debug( "out_debug.out" );
    
    potential_mlg->clear_field_values( );

    cuscomplex  wave_dx     = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_dy     = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_dz     = cuscomplex( 0.0, 0.0 );

    cusfloat    k           = w2k( w, input->water_depth, input->grav_acc );


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
            // distn    =  l2_norm( 3, source_i->panel->center, mesh_gp->source_nodes[j]->position );
            distn    =  std::sqrt( 
                                    pow2s( mesh_gp->source_nodes[j]->position[0] - source_i->panel->center[0] )
                                    +
                                    pow2s( mesh_gp->source_nodes[j]->position[1] - source_i->panel->center[1] )
                                );
            dist     =  distn / input->water_depth;
            is_john = dist > 1.0;
            
            int_value       = 0.0;
            int_dn_sf_value = 0.0;
            int_dn_pf_value = 0.0;
            // Integrate green function normal derivative along the current panel
            if ( i==72 && j == 72 )
            {
                double a = 0.0;
            }

            if ( i == j )
            {
                if ( panel_j->type == DIFFRAC_PANEL_CODE )
                {
                    int_dn_sf_value    = cuscomplex( 0.5, 0.0 );
                    int_dn_pf_value    = cuscomplex( 0.5, 0.0 );
                    if ( is_john )
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
                                                wave_fcn_dn_sf_value,
                                                wave_fcn_dn_pf_value
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
                                                wave_fcn_dn_sf_value,
                                                wave_fcn_dn_pf_value
                                            );
                    }
                    int_value       =   wave_fcn_value / 4.0 / PI;
                }
                else if ( panel_j->type == LID_PANEL_CODE )
                {
                    int_dn_sf_value    = -cuscomplex( 1.0, 0.0 );
                    int_dn_pf_value    = -cuscomplex( 1.0, 0.0 );
                    if ( is_john )
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
                                                wave_fcn_dn_sf_value,
                                                wave_fcn_dn_pf_value
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
                                                wave_fcn_dn_sf_value,
                                                wave_fcn_dn_pf_value
                                            );
                        log_sing_val    = 2.0 * ( LOG2_GAMMA - std::log( nu ) - source_i->panel->free_surface_log_int ) * source_i->panel->area;
                        // std::cout << "log_value: " << source_i->panel->free_surface_log_int << " - " << log_sing_val << std::endl;
                        wave_fcn_value  += cuscomplex( nu * log_sing_val, 0.0 );
                    }
                    int_value = wave_fcn_value / 4.0 / PI;
                    // int_value = 0.0;
                }
                
            }
            else
            {
                wave_fcn_value          = 0.0;
                wave_fcn_dn_sf_value    = 0.0;
                wave_fcn_dn_pf_value    = 0.0;

                if ( is_john )
                {
                    // std::cout << "IsJohn... i: " << i << " - j: " << j << " - A: " << dist << std::endl;
                    // std::cout << "A: " << dist << " - dist: " << distn << " - h: " << input->water_depth << std::endl;
                    // std::cout << "Panel I: " << source_i->panel->center[0] << "," << source_i->panel->center[1] << "," << source_i->panel->center[2] << std::endl;
                    // std::cout << "Panel J: " << panel_j->center[0] << "," << panel_j->center[1] << "," << panel_j->center[2] << std::endl;
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
                                            wave_fcn_dn_pf_value
                                        );
                }
                else
                {
                    // bool verbose = ( col_count==983 && row_count==1099 );

                    quadrature_panel_t<
                                            PanelGeom, 
                                            GWFcnsInterfaceT<NUM_GP2>, 
                                            wave_term_integral<NUM_GP2, G_ON, DGDR_ON, DGDZ_ON, FSLID_OFF>, 
                                            NUM_GP
                                        >( 
                                            source_i->panel, 
                                            gwf_interf, 
                                            wave_fcn_value,
                                            wave_fcn_dn_sf_value,
                                            wave_fcn_dn_pf_value
                                        );
                }


                
                int_value       =   wave_fcn_value / 4.0 / PI;
                int_dn_sf_value =   wave_fcn_dn_sf_value / 4.0 / PI;
                int_dn_pf_value =   wave_fcn_dn_pf_value / 4.0 / PI;

                // if ( col_count==983 && row_count==1099 )
                // {
                //     constexpr int NN = 1;
                //     cusfloat R[NN];
                //     cusfloat z[NN];
                //     cusfloat zeta[NN];
                //     cuscomplex G[NN];
                //     cuscomplex dG_dr[NN];
                //     cuscomplex dG_dz[NN];
                //     cuscomplex dG_dzeta[NN];
                //     cusfloat h = input->water_depth;
                //     cusfloat g = input->grav_acc;
                //     WaveDispersionFONK wave_data( w, h, g );
                //     wave_data.update_full( w, h, g );
                //     BesselFactoryVecUpTo<NN> bessel_factory;

                //     cusfloat dx = panel_j->center[0] - source_i->position[0];
                //     cusfloat dy = panel_j->center[1] - source_i->position[1];
                //     R[0] = std::sqrt( pow2s( dx ) + pow2s( dy ) );
                //     z[0] = panel_j->center[2];
                //     zeta[0] = source_i->position[2];

                //     wave_term_integral<NN, G_ON, DGDR_ON, DGDZ_ON, FSLID_OFF>( R, z, zeta, h, bessel_factory, wave_data, G, dG_dr, dG_dz, dG_dzeta );
                //     std::cout << "OUT (983-1099)" << std::endl;
                //     std::cout << "\t- R:            " << R[0] << std::endl;
                //     std::cout << "\t- z:            " << z[0] << std::endl;
                //     std::cout << "\t- zeta          " << zeta[0] << std::endl;
                //     std::cout << "\t- int_value:    " << wave_fcn_value << std::endl;
                //     std::cout << "\t- G:            " << G[0] << std::endl;
                //     std::cout << "\t- dG_dr:        " << dG_dr[0] << std::endl;
                //     std::cout << "\t- dG_dz:        " << dG_dz[0] << std::endl;
                //     std::cout << "\t- dG_dzeta:     " << dG_dzeta[0] << std::endl;
                //     std::cout << "\t- G*A:          " << G[0] * source_i->panel->area << std::endl;
                //     std::cout << "\t- dG_dr*A:      " << dG_dr[0] * source_i->panel->area << std::endl;
                //     std::cout << "\t- dG_dz*A:      " << dG_dz[0] * source_i->panel->area << std::endl;
                //     std::cout << "\t- dG_dzeta*A:   " << dG_dzeta[0] * source_i->panel->area << std::endl;
                //     std::cout << "\t- wave_fcn:     " << wave_fcn_value << std::endl;
                //     std::cout << "\t- Area:         " << source_i->panel->area << std::endl;
                //     std::cout << std::endl;

                //     std::cout << "PANEL_I: " << std::endl;
                //     std::cout << (*source_i->panel);
                //     std::cout << std::endl;

                //     std::cout << "PANEL_J: " << std::endl;
                //     std::cout << (*panel_j);
                //     std::cout << std::endl;
                // }
                // if ( 
                //         panel_j->type == LID_PANEL_CODE
                //         &&
                //         source_i->panel->type == LID_PANEL_CODE
                //     )
                // {
                //     int_value           = 0.0;
                //     int_dn_sf_value     = + int_dn_sf_value;
                //     int_dn_pf_value     = + int_dn_pf_value;
                // }
            }

            // if ( panel_j->type == LID_PANEL_CODE )
            // {
            //     int_value = cuscomplex( 0.0, 0.0 );
            // }

            // Apply the integral value accordingly
            index_sources   = col_count*scl->num_rows_local+row_count;
            index_pot       = row_count*scl->num_cols_local+col_count;
            if ( is_john )
            {
                potential_mlg->sysmat[index_pot]        = int_value;
                source_form_mlg->sysmat[index_sources]  = int_dn_sf_value;
                pot_form_mlg->sysmat[index_sources]     = int_dn_pf_value;
            }
            else
            {
                potential_mlg->sysmat[index_pot]        = potential_mlg->sysmat_steady[index_pot] + int_value;
                source_form_mlg->sysmat[index_sources]  = source_form_mlg->sysmat_steady[index_sources] + int_dn_sf_value;
                pot_form_mlg->sysmat[index_sources]     = pot_form_mlg->sysmat_steady[index_sources] + int_dn_pf_value;

                // if ( i==1140 && j==1330 )
                // {
                //     std::cout << "int_value_sf: " << source_form_mlg->sysmat_steady[index_sources] << " - " << int_dn_sf_value << std::endl;
                //     std::cout << "int_value_pf: " << pot_form_mlg->sysmat_steady[index_sources] << " - " << int_dn_pf_value << std::endl;
                // }
            }


            // if ( 
            //         ( source_i->panel->type == LID_PANEL_CODE )
            //         &&
            //         ( panel_j->type == LID_PANEL_CODE )
            //     )
            // {
            //     if ( 
            //             std::abs( nu*potential_mlg->sysmat[index_pot] + pot_form_mlg->sysmat[index_sources] ) > 1e-4
            //             &&
            //             !is_john
            //         )
            //     {
            //         if ( col_count==1140 && row_count==1330 )
            //         {
            //             std::cout << "OUT - int_value:    " << potential_mlg->sysmat_steady[index_pot] << " - " << int_value << std::endl;
            //             std::cout << "OUT - int_value_sf: " << source_form_mlg->sysmat_steady[index_sources] << " - " << int_dn_sf_value << std::endl;
            //             std::cout << "OUT - int_value_pf: " << pot_form_mlg->sysmat_steady[index_sources] << " - " << int_dn_pf_value << std::endl;

            //             std::cout << "PANEL_I:" << std::endl;
            //             std::cout << *(source_i->panel);

            //             std::cout << "PANEL_J:" << std::endl;
            //             std::cout << *(panel_j) << std::endl;

            //         }
            //         out_debug << "I: " << col_count << " - J: " << row_count;
            //         out_debug << " - SX: " << source_i->position[0];
            //         out_debug << " - SY: " << source_i->position[1];
            //         out_debug << " - SZ: " << source_i->position[2];
            //         out_debug << " - PX: " << panel_j->center[0];
            //         out_debug << " - PY: " << panel_j->center[1];
            //         out_debug << " - PZ: " << panel_j->center[2];
            //         out_debug << " - Diff: " << nu*potential_mlg->sysmat[index_pot] << " - " << pot_form_mlg->sysmat[index_sources];
            //         out_debug << " - PotSteady: " << potential_mlg->sysmat_steady[index_pot] << " - PotWave: " << int_value;
            //         out_debug << " - PFSteady: " << pot_form_mlg->sysmat_steady[index_sources] << " - PFWave: " << int_dn_pf_value;
            //         out_debug << " - Area: " << source_i->panel->area;
            //         out_debug << " - Nx: " << source_i->normal_vec[0] << " - Ny: " << source_i->normal_vec[1] << " - Nz: " << source_i->normal_vec[2] << std::endl;
            //     }
            //     // pot_form_mlg->sysmat[index_sources] = +nu*potential_mlg->sysmat[index_pot];
            // }

            if ( source_i->panel->type == DIFFRAC_PANEL_CODE )
            {
                for ( int id=0; id<input->dofs_np; id++ )
                {
                    potential_rhs[id*potential_mlg->sysmat_nrows+j] += (
                                                                            mesh_gp->source_nodes[i]->normal_vec[id]
                                                                            *
                                                                            mesh_gp->source_nodes[i]->panel->is_move_f
                                                                        ) * potential_mlg->sysmat[index_pot];
                }
                for ( int id=0; id<input->heads_np; id++ )
                {
                    // panel_j = mesh_gp->source_nodes[j]->panel;
                    // Get wave potential derivatives for the panel
                    wave_dx             =   wave_potential_fo_space_dx(
                                                                            1.0,
                                                                            w,
                                                                            k,
                                                                            input->water_depth,
                                                                            input->grav_acc,
                                                                            mesh_gp->source_nodes[i]->panel->center[0],
                                                                            mesh_gp->source_nodes[i]->panel->center[1],
                                                                            mesh_gp->source_nodes[i]->panel->center[2],
                                                                            input->heads[id]
                                                                        );

                    wave_dy             =   wave_potential_fo_space_dy(
                                                                            1.0,
                                                                            w,
                                                                            k,
                                                                            input->water_depth,
                                                                            input->grav_acc,
                                                                            mesh_gp->source_nodes[i]->panel->center[0],
                                                                            mesh_gp->source_nodes[i]->panel->center[1],
                                                                            mesh_gp->source_nodes[i]->panel->center[2],
                                                                            input->heads[id]
                                                                        );

                    wave_dz             =   wave_potential_fo_space_dz(
                                                                            1.0,
                                                                            w,
                                                                            k,
                                                                            input->water_depth,
                                                                            input->grav_acc,
                                                                            mesh_gp->source_nodes[i]->panel->center[0],
                                                                            mesh_gp->source_nodes[i]->panel->center[1],
                                                                            mesh_gp->source_nodes[i]->panel->center[2],
                                                                            input->heads[id]
                                                                        );
                    
                    // Calculate normal derivative of the wave flow velocities for the jth panel
                    potential_rhs[dofs_offset+id*potential_mlg->sysmat_nrows+j] += -(
                                                                                        wave_dx * mesh_gp->source_nodes[i]->normal_vec[0]
                                                                                        +
                                                                                        wave_dy * mesh_gp->source_nodes[i]->normal_vec[1]
                                                                                        +
                                                                                        wave_dz * mesh_gp->source_nodes[i]->normal_vec[2]
                                                                                    ) * potential_mlg->sysmat[index_pot];
                }
            }

            // out_source_form_steady << std::fixed;
            // out_source_form_steady << std::setprecision( 16 ) << source_form_mlg->sysmat_steady[index_sources].real( ) << " " << source_form_mlg->sysmat_steady[index_sources].imag( ) << "\n";
            // out_source_form_wave << std::fixed;
            // out_source_form_wave << std::setprecision( 16 ) << int_dn_sf_value.real( ) << " " << int_dn_sf_value.imag( ) << "\n";
            // out_source_form_total << std::fixed;
            // out_source_form_total << std::setprecision( 16 ) << source_form_mlg->sysmat[index_sources].real( ) << " " << source_form_mlg->sysmat[index_sources].imag( ) << "\n";

            // out_pot_form_steady << std::fixed;
            // out_pot_form_steady << std::setprecision( 16 ) << pot_form_mlg->sysmat_steady[index_sources].real( ) << " " << pot_form_mlg->sysmat_steady[index_sources].imag( ) << "\n";
            // out_pot_form_wave << std::fixed;
            // out_pot_form_wave << std::setprecision( 16 ) << int_dn_pf_value.real( ) << " " << int_dn_pf_value.imag( ) << "\n";
            // out_pot_form_total << std::fixed;
            // out_pot_form_total << std::setprecision( 16 ) << pot_form_mlg->sysmat[index_sources].real( ) << " " << pot_form_mlg->sysmat[index_sources].imag( ) << "\n";

            // out_pot_steady << std::fixed;
            // out_pot_steady << std::setprecision( 16 ) << potential_mlg->sysmat_steady[index_pot].real( ) << " " << potential_mlg->sysmat_steady[index_pot].imag( ) << "\n";
            // out_pot_wave << std::fixed;
            // out_pot_wave << std::setprecision( 16 ) << int_value.real( ) << " " << int_value.imag( ) << "\n";
            // out_pot_total << std::fixed;
            // out_pot_total << std::setprecision( 16 ) << potential_mlg->sysmat[index_pot].real( ) << " " << potential_mlg->sysmat[index_pot].imag( ) << "\n";

            // Advance row count
            row_count++;
        }

        // Advance column count
        col_count++;
    }
    // out_source_form_steady.close( );
    // out_source_form_wave.close( );
    // out_source_form_total.close( );

    // out_pot_form_steady.close( );
    // out_pot_form_wave.close( );
    // out_pot_form_total.close( );

    // out_pot_steady.close( );
    // out_pot_wave.close( );
    // out_pot_total.close( );

    // out_debug.close( );

    // Potential rhs
    // cuscomplex* rhs_pot = generate_empty_vector<cuscomplex>( potential_mlg->sysmat_nrows * input->dofs_np );

    // for ( int i=0; i<input->dofs_np; i++ )
    // {
    //     for ( int j=0; j<potential_mlg->sysmat_nrows; j++ )
    //     {
    //         rhs_pot[i*potential_mlg->sysmat_nrows+j] = mesh_gp->panels[j]->normal_vec[i];
    //     }
    // }

    // clear_vector( potential_mlg->sysmat_nrows * ( input->dofs_np + input->heads_np ), potential_rhs );
    // cuscomplex  alpha( 1.0, 0.0 );
    // cuscomplex  beta( 0.0, 0.0 );
    // int         icnx = 1;
    // int         icny = 1;
    // for ( int i=0; i<input->dofs_np; i++ )
    // {
    //     cblas_gemv<cuscomplex>( 
    //                                 CblasRowMajor,
    //                                 CblasNoTrans,
    //                                 potential_mlg->sysmat_nrows,
    //                                 potential_mlg->sysmat_ncols,
    //                                 &alpha,
    //                                 potential_mlg->sysmat,
    //                                 potential_mlg->sysmat_ncols,
    //                                 &(rhs_pot[i*potential_mlg->sysmat_nrows+potential_mlg->start_col]),
    //                                 icnx,
    //                                 &beta,
    //                                 &(potential_rhs[i*potential_mlg->sysmat_nrows]),
    //                                 icny
    //                             );
    // }

    /***************************************/
    /***** Fill Hydromechanics RHS  ********/
    /***************************************/

    MPI_Allreduce(
                        potential_rhs,
                        pot_form_mlg->field_values,
                        mesh_gp->panels_tnp * ( input->dofs_np + input->heads_np ),
                        mpi_cuscomplex,
                        MPI_SUM,
                        MPI_COMM_WORLD
                    );

    // std::stringstream ss_rhs_pf;
    // ss_rhs_pf << "rhs_pf_" << mpi_config->proc_rank << "_" << w << ".dat";
    // std::ofstream out_rhs_pf( ss_rhs_pf.str( ) );
    // int count_1 = 0;
    // for ( int i=0; i<input->dofs_np+input->heads_np; i++ )
    // {
    //     for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
    //     {
            
    //         out_rhs_pf << pot_form_mlg->field_values[count_1].real( ) << " " << pot_form_mlg->field_values[count_1].imag( ) << "\n";
    //         count_1++;
    //     }
    // }
    // out_rhs_pf.close( );

    cusfloat potential_cond   = 0.0;
    cusfloat pot_form_cond    = 0.0;
    cusfloat source_form_cond = 0.0;

    scl->Cond( potential_mlg->sysmat,   potential_cond      );
    scl->Cond( pot_form_mlg->sysmat,    pot_form_cond       );
    scl->Cond( source_form_mlg->sysmat, source_form_cond    );

    if ( mpi_config->is_root( ) )
    {
        std::cout << "Condition Number -> PF: " << pot_form_cond << " - SF: " << source_form_cond << " - POT: " << potential_cond << std::endl;
    }
                    
    scl->Solve( pot_form_mlg->sysmat, pot_form_mlg->field_values );

    // int     nrhs    = input->dofs_np+input->heads_np;
    // int     info    = 0;
    // int*    ipiv    = generate_empty_vector<int>( pot_form_mlg->sysmat_nrows );
    // gesv<cuscomplex>( 
    //                     &(pot_form_mlg->sysmat_nrows),
    //                     &(nrhs),
    //                     pot_form_mlg->sysmat,
    //                     &(pot_form_mlg->sysmat_nrows),
    //                     ipiv,
    //                     pot_form_mlg->field_values,
    //                     &(pot_form_mlg->sysmat_nrows),
    //                     &info
    //                 );


    // if ( info != 0 )
    // {
    //     std::cerr << "ERROR - Calculating RAO" << std::endl;
    //     std::cerr << "Error solving system of equations - Info: " << info << std::endl;
    //     throw std::runtime_error( "" );
    // }

    // std::stringstream ss_pot_pf;
    // ss_pot_pf << "pot_pf_" << mpi_config->proc_rank << "_" << w << ".dat";
    // std::ofstream out_pot_pf( ss_pot_pf.str( ) );
    // count_1 = 0;
    // for ( int i=0; i<input->dofs_np+input->heads_np; i++ )
    // {
    //     for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
    //     {
            
    //         out_pot_pf << pot_form_mlg->field_values[count_1].real( ) << " " << pot_form_mlg->field_values[count_1].imag( ) << "\n";
    //         count_1++;
    //     }
    // }
    // out_pot_pf.close( );
    
    // Declare local variables to be used
    int count = 0;

    // Fill RHS vector
    count       = 0;
    for ( int i=0; i<input->dofs_np; i++ )
    {
        for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
        {
            panel_j = mesh_gp->source_nodes[j]->panel;
            // if ( i==2 )
            // {
            //     std::cout << "J: " << j << " - PanelType: " << panel_j->type << std::endl;
            // }
            if ( panel_j->type == DIFFRAC_PANEL_CODE )
            {
                sources_int[count] = ( 
                                            mesh_gp->source_nodes[j]->normal_vec[i]
                                            *
                                            mesh_gp->source_nodes[j]->panel->is_move_f
                                        );

                // if ( i == 2 )
                // {
                //     std::cout << "Count: " << count;
                //     std::cout << " - J: " << j;
                //     std::cout << " - X: " << mesh_gp->source_nodes[j]->position[0];
                //     std::cout << " - Y: " << mesh_gp->source_nodes[j]->position[1];
                //     std::cout << " - Z: " << mesh_gp->source_nodes[j]->position[2];
                //     std::cout << " - X: " << mesh_gp->source_nodes[j]->normal_vec[0];
                //     std::cout << " - Y: " << mesh_gp->source_nodes[j]->normal_vec[1];
                //     std::cout << " - Z: " << mesh_gp->source_nodes[j]->normal_vec[2];
                //     std::cout << " - RHS: " << sources_int[count];
                //     std::cout << std::endl;
                // }
            }
            else if ( panel_j->type == LID_PANEL_CODE )
            {
                // if ( i == 2 )
                // {
                //     std::cout << "J: " << j << " - Nz: " << mesh_gp->source_nodes[j]->normal_vec[i] << std::endl;
                // }
                sources_int[count] = -nu*potential_mlg->field_values[i*potential_mlg->sysmat_nrows+j];
            }
            count++;
        }
    }

    /***************************************/
    /****** Fill Wave Exciting RHS  ********/
    /***************************************/
    // Define local variables to manage array indexes
    count       = input->dofs_np * mesh_gp->source_nodes_tnp;
    // cuscomplex  wave_dx     = cuscomplex( 0.0, 0.0 );
    // cuscomplex  wave_dy     = cuscomplex( 0.0, 0.0 );
    // cuscomplex  wave_dz     = cuscomplex( 0.0, 0.0 );

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
                sources_int[count]  = -nu*potential_mlg->field_values[(input->dofs_np+i)*potential_mlg->sysmat_nrows+j];
            }
            
            count++;
        }
    }

    // std::cout << "StartRow: " << scl->start_row_0 << " - EndRow: " << scl->end_row_0 << std::endl;
    // std::stringstream ss_rhs;
    // ss_rhs << "rhs_sf_" << mpi_config->proc_rank << "_" << w << ".dat";
    // std::ofstream out_rhs( ss_rhs.str( ) );
    // count_1 = 0;
    // for ( int i=0; i<input->dofs_np+input->heads_np; i++ )
    // {
    //     for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
    //     {
            
    //         out_rhs << sources_int[count_1].real( ) << " " << sources_int[count_1].imag( ) << "\n";
    //         count_1++;
    //     }
    // }
    // out_rhs.close( );

    // Solve system of equations
    scl->Solve( source_form_mlg->sysmat, sources_int );

    // clear_vector( source_form_mlg->sysmat_nrows, ipiv );

    // gesv<cuscomplex>( 
    //                     &(source_form_mlg->sysmat_nrows),
    //                     &(nrhs),
    //                     source_form_mlg->sysmat,
    //                     &(source_form_mlg->sysmat_nrows),
    //                     ipiv,
    //                     sources_int,
    //                     &(source_form_mlg->sysmat_nrows),
    //                     &info
    //                 );

    // if ( info != 0 )
    // {
    //     std::cerr << "ERROR - Calculating RAO" << std::endl;
    //     std::cerr << "Error solving system of equations - Info: " << info << std::endl;
    //     throw std::runtime_error( "" );
    // }

    // mkl_free( ipiv );


    // std::stringstream ss_sources;
    // ss_sources << "sources_" << mpi_config->proc_rank << "_" << w << ".dat";
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
    // out_sources.close( );

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

    mkl_free( potential_rhs );
}


void    calculate_gf_intensity_steady_sysmat_lin_v3(
                                                    Input*          input,
                                                    SclCmpx*        scl,
                                                    MeshGroup*      mesh_gp,
                                                    cuscomplex*     source_form_sysmat,
                                                    cuscomplex*     pot_form_sysmat
                                                )
{
    /***************************************/
    /******** Fill system matrix  **********/
    /***************************************/
    // Loop over panels to integrate value
    int         col_count   = 0;
    cuscomplex  int_value_sf( 0.0, 0.0 );
    cuscomplex  int_value_pf( 0.0, 0.0 );
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
    cusfloat    vel_total_pf[ndim];     clear_vector( ndim, vel_total_pf );
    cusfloat    vel_total_sf[ndim];     clear_vector( ndim, vel_total_sf );

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
                int_value_sf = cuscomplex( 0.0, 0.0 );
                int_value_pf = cuscomplex( 0.0, 0.0 );
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
                clear_vector( ndim, vel_total_pf );
                clear_vector( ndim, vel_total_sf );
                
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
                vel_total_sf[0] = - ( vel_0[0] + vel_1[0] + vel_2[0] + vel_3[0] + vel_4[0] + vel_5[0] );
                vel_total_sf[1] = - ( vel_0[1] + vel_1[1] + vel_2[1] + vel_3[1] + vel_4[1] + vel_5[1] );
                vel_total_sf[2] = - ( vel_0[2] + vel_1[2] + vel_2[2] + vel_3[2] - vel_4[2] + vel_5[2] );

                vel_total_pf[0] = - vel_total_sf[0];
                vel_total_pf[1] = - vel_total_sf[1];
                vel_total_pf[2] = - ( - vel_0[2] + vel_1[2] + vel_2[2] - vel_3[2] + vel_4[2] + vel_5[2] );
                                        
                int_value_sf    = (
                                        mesh_gp->source_nodes[j]->normal_vec[0] * vel_total_sf[0]
                                        +
                                        mesh_gp->source_nodes[j]->normal_vec[1] * vel_total_sf[1]
                                        +
                                        mesh_gp->source_nodes[j]->normal_vec[2] * vel_total_sf[2]
                                    ) / 4.0 / PI;

                int_value_pf    = (
                                        mesh_gp->source_nodes[i]->normal_vec[0] * vel_total_pf[0]
                                        +
                                        mesh_gp->source_nodes[i]->normal_vec[1] * vel_total_pf[1]
                                        +
                                        mesh_gp->source_nodes[i]->normal_vec[2] * vel_total_pf[2]
                                    ) / 4.0 / PI;

                // if ( i==1140 && j==1330 )
                // {
                //     std::cout << "int_value_sf: " << int_value_sf << std::endl;
                //     std::cout << "int_value_pf: " << int_value_pf << std::endl;
                // }
            }

            // if ( 
            //         mesh_gp->source_nodes[j]->panel->type == LID_PANEL_CODE
            //         &&
            //         panel_i->type == LID_PANEL_CODE
            //     )
            // {
            //     int_value       = + int_value;
            // }

            source_form_sysmat[col_count*scl->num_rows_local+row_count] = int_value_sf;
            pot_form_sysmat[col_count*scl->num_rows_local+row_count]    = int_value_pf;

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