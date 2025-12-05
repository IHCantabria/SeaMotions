
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
#include "hydromechanics.hpp"
#include "../../math/integration.hpp"


void    calculate_hydromechanic_coeffs_lin( 
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     panels_pot,
                                                cusfloat        ang_freq,
                                                cusfloat*       added_mass,
                                                cusfloat*       damping_rad,
                                                cuscomplex*     panels_press
                                            )
{
    // Define local variables
    int         dofs_np         = input->dofs_np;
    int         elem_end_pos    = 0;
    int         elem_start_pos  = 0;
    int         index           = 0;
    int         index_1         = 0;
    cuscomplex  press_i         = cuscomplex( 0.0, 0.0 );
    cusfloat    rho_w           = input->water_density;

    // Get hydromechanic coefficients for all the degrees of freedom
    for ( int ib=0; ib<mesh_gp->meshes_np; ib++ )
    {
         // Calculate MPI data chunks
        elem_end_pos    = 0;
        elem_start_pos  = 0;
        mpi_config->get_1d_bounds( 
                                        mesh_gp->panels_np[ib], 
                                        elem_start_pos, 
                                        elem_end_pos 
                                    );
        
        for ( int id=0; id<input->dofs_np; id++ )
        {
            for ( int jd=0; jd<input->dofs_np; jd++ )
            {
                // Get current matrix index
                for ( int jb=0; jb<mesh_gp->meshes_np; jb++ )
                {
                    index = (
                                ib * ( dofs_np * dofs_np * mesh_gp->meshes_np )
                                +
                                id * ( dofs_np * mesh_gp->meshes_np )
                                + 
                                jb * dofs_np 
                                + 
                                jd
                            );
                    added_mass[index]   = 0.0;
                    damping_rad[index]  = 0.0;
                    for ( int ie=elem_start_pos; ie<elem_end_pos; ie++ )
                    {
                        if ( mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->type == DIFFRAC_PANEL_CODE )
                        {
                            // Integrate pressure over panel
                            index_1                 = mesh_gp->panels_cnp[jb] + ie + mesh_gp->panels_tnp * id;
                            press_i                 = ( 
                                                            panels_pot[index_1] 
                                                            * 
                                                            mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->normal_vec[jd] 
                                                            * 
                                                            mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->area
                                                        );
                            added_mass[index]       +=  rho_w * press_i.real( );
                            damping_rad[index]      +=  rho_w * press_i.imag( ) * ang_freq;
                            panels_press[index_1]   = - cuscomplex( 0.0, 1.0 ) * rho_w * ang_freq * panels_pot[index_1] ;
                        }
                    }

                }
            }
        }
    }

    // std::ofstream               out_file( "D:/sergio/0058_Isobara/2_Simulations/seamotions_files/dp0/AQW-3/AQW/AQ/Analysis/sm_potentials.dat" );
    // cusfloat rhow               = 1025.0;
    // cusfloat added_mass_surge   = 0.0;
    // cusfloat damping_surge      = 0.0;
    // for ( int i=0; i<mesh_gp->panels_tnp; i++ )
    // {
    //     press_i             = rhow * panels_pot[i] * mesh_gp->panels[i]->normal_vec[0] * mesh_gp->panels[i]->area;
    //     added_mass_surge    += press_i.imag( ) / ang_freq;
    //     damping_surge       -= press_i.real( );

    //     out_file << panels_pot[i].real( ) << "\t";
    //     out_file << panels_pot[i].imag( ) << "\t";
    //     out_file << press_i.real( ) << "\t";
    //     out_file << press_i.imag( ) << "\t";
    //     out_file << mesh_gp->panels[i]->center[0] << "\t";
    //     out_file << mesh_gp->panels[i]->center[1] << "\t";
    //     out_file << mesh_gp->panels[i]->center[2] << "\t";
    //     out_file << mesh_gp->panels[i]->normal_vec[0] << "\t";
    //     out_file << mesh_gp->panels[i]->normal_vec[1] << "\t";
    //     out_file << mesh_gp->panels[i]->normal_vec[2] << "\t";
    //     out_file << mesh_gp->panels[i]->area << "\t";
    //     out_file << "\n";
    // }

    // out_file.close( );

    // std::cout << "Added Mass:   " << added_mass_surge << " " << added_mass[0] << std::endl;
    // std::cout << "Damping:      " << damping_surge << " " << damping_rad[0] << std::endl;

    // std::cout << "ADDED MASS MATRIX: " << std::endl;
    // int ndofs = 6;
    // for ( int i=0; i<ndofs; i++ )
    // {
    //     std::cout << added_mass[i*ndofs];
    //     for ( int j=1; j<ndofs; j++ )
    //     {
    //         std::cout << " " << added_mass[i*ndofs+j];
    //     }
    //     std::cout << "\n";
    // }
}


// void    calculate_hydromechanic_coeffs_nlin(
//                                                 Input*          input,
//                                                 MpiConfig*      mpi_config,
//                                                 MeshGroup*      mesh_gp,
//                                                 HMFInterface*   hmf_interf,
//                                                 cusfloat        ang_freq,
//                                                 cusfloat*       added_mass,
//                                                 cusfloat*       damping_rad
//                                             )
// {
//     // Define local variables
//     int         dofs_np = input->dofs_np;
//     cusfloat    rho_w   = input->water_density;

//     // Generate lambda function for the integration
//     auto target_fcn = [hmf_interf]
//                     (
//                         cusfloat xi, 
//                         cusfloat eta, 
//                         cusfloat x,
//                         cusfloat y,
//                         cusfloat z
//                     )
//                     {
//                         return (*hmf_interf)( xi, eta, x, y, z );
//                     };

//     // Allocate space for pressure vector
//     int max_panels  = 0;
//     for ( int i=0; i<mesh_gp->meshes_np; i++ )
//     {
//         // Get ith mesh panels
//         if ( mesh_gp->panels_np[i] > max_panels )
//         {
//             max_panels = mesh_gp->panels_np[i];
//         }
//     }
//     cuscomplex* pressure    = generate_empty_vector<cuscomplex>( mesh_gp->meshes_np * dofs_np * max_panels );

//     // Loop over first dimension of degrees of freedrom
//     int         elem_end_pos    = 0;
//     int         elem_start_pos  = 0;
//     int         index           = 0;
//     int         index_1         = 0;
//     cuscomplex  press_i         = 0.0;

//     for ( int ib=0; ib<mesh_gp->meshes_np; ib++ )
//     {
//         // Calculate MPI data chunks
//         elem_end_pos    = 0;
//         elem_start_pos  = 0;
//         mpi_config->get_1d_bounds( 
//                                         mesh_gp->panels_np[ib], 
//                                         elem_start_pos, 
//                                         elem_end_pos 
//                                     );

//         for ( int id=0; id<dofs_np; id++ )
//         {
//             for ( int jb=0; jb<mesh_gp->meshes_np; jb++ )
//             {
//                 // Set ith dof
//                 hmf_interf->set_start_index_i( 
//                                                     mesh_gp->source_nodes_tnp*id,
//                                                     mesh_gp->source_nodes_cnp[jb],
//                                                     mesh_gp->source_nodes_cnp[jb+1]
//                                                 );

//                 // Loop over panels to integrate the wave radiation
//                 // pressure along the floating object external shape
//                 for ( int ie=elem_start_pos; ie<elem_end_pos; ie++ )
//                 {
//                     // Set new panel
//                     hmf_interf->set_panel( mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie] );

//                     // Integrate pressure over panel
//                     index           = ( max_panels * mesh_gp->meshes_np ) * id + max_panels * jb + ie;
//                     if ( mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->type == DIFFRAC_PANEL_CODE )
//                     {
//                         pressure[index] = adaptive_quadrature_panel(
//                                                                         mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie],
//                                                                         target_fcn,
//                                                                         input->press_abs_err,
//                                                                         input->press_rel_err,
//                                                                         input->is_block_adaption,
//                                                                         true,
//                                                                         input->gauss_order
//                                                                     );
//                     }
//                 }
//             }
//         }

//         // Get hydromechanic coefficients for all the degrees of freedom
//         for ( int id=0; id<dofs_np; id++ )
//         {
//             for ( int jd=0; jd<dofs_np; jd++ )
//             {
//                 // Get current matrix index
//                 for ( int jb=0; jb<mesh_gp->meshes_np; jb++ )
//                 {
//                     index = (
//                                 ib * ( dofs_np * dofs_np * mesh_gp->meshes_np )
//                                 +
//                                 id * ( dofs_np * mesh_gp->meshes_np )
//                                 + 
//                                 jb * dofs_np 
//                                 + 
//                                 jd
//                             );
//                     added_mass[index]   = 0.0;
//                     damping_rad[index]  = 0.0;
//                     for ( int ie=elem_start_pos; ie<elem_end_pos; ie++ )
//                     {
//                         // Integrate pressure over panel
//                         index_1             = ( max_panels * mesh_gp->meshes_np ) * id + max_panels * jb + ie;
//                         press_i             = pressure[index_1] * mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->normal_vec[jd];
//                         added_mass[index]   -=  rho_w * press_i.real( );
//                         damping_rad[index]  -=  rho_w * ang_freq * press_i.imag( );
//                     }

//                 }
//             }
//         }
//     }

//     // Deallocate local allocated heap memory
//     mkl_free( pressure );
// }
