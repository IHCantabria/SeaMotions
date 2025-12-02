
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
#include "../../cli_header_banner.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../containers/mpi_timer.hpp"
#include "../../interfaces/stab_interface_t.hpp"
#include "../../mesh/stability_mesh.hpp"
#include "../../mesh/mesh.hpp"
#include "stab_input.hpp"
#include "../../tools.hpp"


template<typename T, int NGP>
void        quadrature_panel_hydrostat(
                                            T*                  panel,
                                            cusfloat*           press_force,
                                            bool                verbose
                                        )
{
    // Calculate function value
    // target_fcn.template operator()<Kernel>( 
    //                                             panel->xl, 
    //                                             panel->yl, 
    //                                             panel->gauss_points_global_x, 
    //                                             panel->gauss_points_global_y, 
    //                                             panel->gauss_points_global_z,
    //                                             verbose
    //                                         );

    // result_G        = 0.0;
    // result_G_dn_sf  = 0.0;
    // result_G_dn_pf  = 0.0;
    // for ( int i=0; i<NGP*NGP; i++ )
    // {
    //     GAUSS_2D_LOOP( result_G,        target_fcn.G            );
    //     GAUSS_2D_LOOP( result_G_dn_sf,  target_fcn.dG_dn_sf     );
    //     GAUSS_2D_LOOP( result_G_dn_pf,  target_fcn.dG_dn_pf     );

    //     if ( verbose )
    //     {
    //         std::cout << "Wx: " << GaussPointsT<2,NGP>::weights_x[i];
    //         std::cout << " - Wy: " << GaussPointsT<2,NGP>::weights_y[i];
    //         std::cout << " - fcn: " << target_fcn.G[i];
    //         std::cout << " - Jac: " << panel->jac_det_gauss_points[i] << std::endl;
    //     }

    // }

}


// void calculate_hydrostatic_forces( 
//                                         Input*          input, 
//                                         StabilityMesh*  mesh, 
//                                         cusfloat*       gforce, 
//                                         cusfloat*       gmoment 
//                                 )
// {
    // // Define local variables
    // cusfloat lforce[3];
    // cusfloat lmoment[3];

    // Define calculation interface
    // StabInterfaceT<NUM_GP2> stab_interf(  );

    // // Loop over panels to calculate global forces and moments
    // for ( int i=0; i<mesh->elems_np; i++ )
    // {
    //     if ( mesh->[i] < 1 )
    //     {

    //     }
    // }
// }


int main( int argc, char* argv[ ] )
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 1))
    {
        return 1;
    }

    std::string case_fopath(argv[1]);

    /*****************************************/
    /****** Initialize MPI environment *******/
    /*****************************************/
    MPI_Init( NULL, NULL );

    // Get total number of processors
    int procs_total = 0;
    MPI_Comm_size(
                    MPI_COMM_WORLD,
                    &procs_total
                );

    // Get current process rank
    int proc_rank = 0;
    MPI_Comm_rank(
                    MPI_COMM_WORLD,
                    &proc_rank
                );

    // Create MPI Configuration system
    MpiConfig mpi_config( proc_rank, procs_total, MPI_ROOT_PROC_ID, MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );

    /*****************************************/
    /******** Print Header Section ***********/
    /*****************************************/
    cli_header_banner<true>( case_fopath, "Stability" );


    /*****************************************/
    /************ Read Input data ************/
    /*****************************************/
    StabInput input( case_fopath );

    
    /*****************************************/
    /*********** Stability kernel ************/
    /*****************************************/
    MpiTimer case_timer;

    // // Read input mesh
    // std::string mesh_fipath = "d:/sergio/0050_OASIS_SM/SeaMotionsStabValidation_files/user_files/0_Box/box_stability_tri.poly";
    // std::string body_name   = "box";
    // cusfloat    cog[3]      = { 0.0, 0.0, 5.0 };
    cusfloat    draft       = 4.1;

    Mesh mesh_ref( input.mesh_fipath, input.body_name, input.cog, false , 0 );
    StabilityMesh mesh_mov( input.mesh_fipath, input.body_name, input.cog, false , 0, draft );

    mesh_mov.move( 0.0, 0.0, -draft, 0.0, 0.0, 0.0 );
    mesh_mov.check_underwater_panels( );

    mesh_mov.write_underwater_panels( 
                                        case_fopath,
                                        input.mesh_finame + "_tri"
                                    );

    // // Calculate hydrostatic forces
    // cusfloat force[3];
    // cusfloat moment[3];
    // calculate_hydrostatic_forces( &input, &mesh_mov, force, moment );

    /*****************************************/
    /**** Close program and final actions ****/
    /*****************************************/

    // Print Elapsed time
    if ( mpi_config.is_root( ) )
    {
        std::cout << std::endl << std::endl;
        std::cout << " -> Seamotions (Frequency) finished!" << std::endl;
        std::cout << " ---> Elapsed wall time for calculation [s]: " << case_timer << std::endl;
    }

    return 0;
}
