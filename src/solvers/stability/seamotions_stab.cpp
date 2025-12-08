
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
#include "../../containers/initial_stability.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../containers/mpi_timer.hpp"
#include "../../interfaces/hydrostatic_force_interface.hpp"
#include "../../interfaces/stab_interface_t.hpp"
#include "../../math/integration.hpp"
#include "../../math/nonlinear_solvers_math.hpp"
#include "../../mesh/rigid_body_mesh.hpp"
#include "../../mesh/mesh.hpp"
#include "stab_input.hpp"
#include "../../tools.hpp"
#include "../../math/newmark.hpp"


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
//                                         RigidBodyMesh*  mesh, 
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

constexpr int _NDOF = 6;
template<typename T>
class HydrostaticForces
{
private:
    /* Define class private attributes */ 
    RigidBodyMesh*  _mesh           = nullptr;  // Storage pointer of the mesh object used to represent the floater external surface
    T               _force_interf   = nullptr;  // Storage pointer of the functor interface used to calculate the force over the panel
    cusfloat        _weight         = 0.0;      // Weight of the floater to be accounted on external force vector

public:
    /* Define class constructor */
    HydrostaticForces( 
                            RigidBodyMesh*  mesh_in,
                            T               force_interf_in,
                            cusfloat        weight_in
                        )
    {
        // Storage input arguments
        this->_force_interf = force_interf_in;
        this->_mesh         = mesh_in;
        this->_weight       = weight_in;
    }

    /* Define class overloaded operators */

    // Overload = operator to have a functor behaviour
    // that matches the time solver interface
    void    operator()   ( 
                            cusfloat    time,
                            cusfloat    time_step,
                            cusfloat*   pos,
                            cusfloat*   ,
                            cusfloat*   ,
                            cusfloat*   rhs
                        ) const
    {
        /*  1.  Move and rotate mesh around centre of gravity to 
                be the state predicted by the stepper */
        this->_mesh->move( pos[0], pos[1], pos[2], pos[3], pos[4], pos[5] );

        /*  2.  Calculate hydrostatic forces and moments around the centre of 
                gravity*/
        
        //  2.1 Clean input RHS vector in order to acount for old spurious data
        clear_vector( _NDOF, rhs );

        //  2.2 Calculate hydrostatic forces
        int count_normal = 0;
        cusfloat rprev = 0.0;
        for ( int i=0; i<this->_mesh->get_elems_np( ); i++ )
        {
            ( *this->_force_interf )( this->_mesh->get_panel( i ), rhs );
            if ( this->_mesh->get_panel( i )->normal_vec[2] > 0.1 )
            {
                count_normal++;
                // std::cout << " --> Normal pointing upwards" << std::endl;
            }
            // std::cout << "rhs[2]: " << rhs[2] << " - Prev: " << rprev << " - Diff: " << rhs[2] - rprev << " - Normal: " << this->_mesh->get_panel( i )->normal_vec[0] << ", " << this->_mesh->get_panel( i )->normal_vec[1] << ", " << this->_mesh->get_panel( i )->normal_vec[2] << std::endl;
            rprev = rhs[2];
        }

        // std::cout << "Num Normals: " << count_normal << std::endl;
        // std::cout << "Hydro Force: " << rhs[2];
        /*  3.  Add body weight to vertical force */
        rhs[2] -= this->_weight;
        // std::cout << " - Weight: " << this->_weight << " - Diff:" << rhs[2] << std::endl;

    }
};


template<typename T>
class HydrostaticForcesNLin
{
private:
    /* Define class private attributes */ 
    RigidBodyMesh*  _mesh           = nullptr;  // Storage pointer of the mesh object used to represent the floater external surface
    T               _force_interf   = nullptr;  // Storage pointer of the functor interface used to calculate the force over the panel
    cusfloat        _pos_init[6]    ;           // Storage initial position of the mesh. It will be used to refresh mesh state during execution except for heave
    cusfloat        _weight         = 0.0;      // Weight of the floater to be accounted on external force vector

public:
    /* Define class constructor */
    HydrostaticForcesNLin( 
                            RigidBodyMesh*  mesh_in,
                            cusfloat*       pos_in,
                            T               force_interf_in,
                            cusfloat        weight_in
                        )
    {
        // Storage input arguments
        this->_force_interf = force_interf_in;
        this->_mesh         = mesh_in;
        this->_weight       = weight_in;

        copy_vector( 6, pos_in, this->_pos_init );
    }

    /* Define class overloaded operators */

    // Overload = operator to have a functor behaviour
    // that matches the time solver interface
    cusfloat operator()   ( 
                            cusfloat heave
                        ) const
    {
        /*  1.  Move and rotate mesh around centre of gravity to 
                be the state predicted by the stepper */
        this->_mesh->move( 
                            this->_pos_init[0], 
                            this->_pos_init[1], 
                            heave, 
                            this->_pos_init[3], 
                            this->_pos_init[4], 
                            this->_pos_init[5] 
                        );

        /*  2.  Cut mesh along the free surface */
        this->_mesh->check_underwater_panels( );

        /*  3.  Calculate hydrostatic forces and moments around the centre of 
                gravity*/
        
        //  3.1 Clean input RHS vector in order to acount for old spurious data
        // clear_vector( _NDOF, rhs );
        cusfloat rhs[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

        //  3.2 Calculate hydrostatic forces
        PanelGeom* paneli = nullptr;
        int count_normal = 0;
        cusfloat rprev = 0.0;
        for ( int i=0; i<this->_mesh->get_elems_np( ); i++ )
        {
            paneli = this->_mesh->get_panel( i );
            if ( paneli->location_zone == -1 )
            {
                ( *this->_force_interf )( this->_mesh->get_panel( i ), rhs );
            }
            // if ( this->_mesh->get_panel( i )->normal_vec[2] > 0.1 )
            // {
                // count_normal++;
                // std::cout << " --> Normal pointing upwards" << std::endl;
            // }
            // std::cout << "rhs[2]: " << rhs[2] << " - Prev: " << rprev << " - Diff: " << rhs[2] - rprev << " - Normal: " << this->_mesh->get_panel( i )->normal_vec[0] << ", " << this->_mesh->get_panel( i )->normal_vec[1] << ", " << this->_mesh->get_panel( i )->normal_vec[2] << std::endl;
            // rprev = rhs[2];
        }

        // std::cout << "Num Normals: " << count_normal << std::endl;
        // std::cout << "Hydro Force: " << rhs[2];
        /*  3.  Add body weight to vertical force */
        rhs[2] -= this->_weight;
        // std::cout << " - Weight: " << this->_weight << " - Diff:" << rhs[2] << std::endl;

        return rhs[2];
    }
};


void calculate_hydrostatic_properties( 
                                        StabInput*      input, 
                                        RigidBodyMesh*  mesh
                                    )
{
    for ( std::size_t i=0; i<input->draft_hs.size( ); i++ )
    {
        // Adjust mesh to the draft
        mesh->move( 0.0, 0.0, -input->draft_hs[i], 0.0, 0.0, 0.0 );

        // Cut along the free surface in order to properly 
        // calculate initial stability parameters
        mesh->check_underwater_panels( );

        // Calcualte initial estability parameters
    }
}


void test_newmark( StabInput* input, RigidBodyMesh* mesh )
{
    // Use mesh bounding box to have an estimation of the
    // dynamical properties of object to set up dynamical
    // simulation to find the equilibrim
    cusfloat    lx      = mesh->x_max - mesh->x_min;
    cusfloat    ly      = mesh->y_max - mesh->y_min;
    cusfloat    lz      = mesh->z_max - mesh->z_min;
    cusfloat    area    = lx * ly;
    cusfloat    volume  = area * lz;
    cusfloat    mass    = volume * input->water_density;
    cusfloat    ixx     = mass * ( pow2s( ly ) + pow2s( lz ) ) / 12.0;
    cusfloat    iyy     = mass * ( pow2s( lx ) + pow2s( lz ) ) / 12.0;
    cusfloat    izz     = std::sqrt( pow2s( ixx ) + pow2s( iyy ) );
    cusfloat    k22     = input->water_density * input->grav_acc * area;
    cusfloat    k33     = input->water_density * input->grav_acc * volume * 10.0;
    cusfloat    k44     = input->water_density * input->grav_acc * volume * 10.0;
    cusfloat    d22     = 2.0 * std::sqrt( mass * k22 );
    cusfloat    d33     = 2.0 * std::sqrt( mass * k33 );
    cusfloat    d44     = 2.0 * std::sqrt( mass * k44 );
    cusfloat    t22     = 2.0 * PI * std::sqrt( mass / k22 );
    cusfloat    t33     = 2.0 * PI * std::sqrt( ixx  / k33 );
    cusfloat    t44     = 2.0 * PI * std::sqrt( iyy  / k44 );

    // Create matrixes in dense form
    cusfloat    mass_d[36]  = { 
                                    mass, 0.0, 0.0, 0.0, 0.0, 0.0,
                                    0.0, mass, 0.0, 0.0, 0.0, 0.0,
                                    0.0, 0.0, mass, 0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0, ixx, 0.0, 0.0,
                                    0.0, 0.0, 0.0, 0.0, iyy, 0.0,
                                    0.0, 0.0, 0.0, 0.0, 0.0, izz
                                };
    
    cusfloat    stiff_d[36] = { 
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                                };

    cusfloat    damp_d[36]  = { 
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                    0.0, 0.0, d22, 0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0, d33, 0.0, 0.0,
                                    0.0, 0.0, 0.0, 0.0, d44, 0.0,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                                };

    // Create system matrixes in CSR format
    CSRMatrix*  mass_mat        = new CSRMatrix( 6, mass_d     );
    CSRMatrix*  damp_mat        = new CSRMatrix( 6, damp_d     );
    CSRMatrix*  stiff_mat       = new CSRMatrix( 6, stiff_d    );

    // Create restrictions vector
    int         restrictions[6] = { 1, 1, 0, 1, 1, 1 };
    cusfloat    y0_pos[6]       = { 0.0, 0.0, -1.72, 0.0, 0.0, 0.0 };
    cusfloat    y0_vel[6]       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    cusfloat    y0_acc[6]       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    // Calculate maximum time of simulation based on natural periods
    // estimated
    cusfloat    c2          = restrictions[ 2 ];
    cusfloat    c3          = restrictions[ 3 ];
    cusfloat    c4          = restrictions[ 4 ];
    cusfloat    max_time    = ( c2 * t22 + c3 * t33 + c4 * t44 ) / ( c2 + c3 + c4 );

    // Create functor to calculate hydrostatic forces
    HydrostaticForceInterface<NUM_GP>   hydrostat_force_interf( 
                                                                    input->water_density, 
                                                                    input->grav_acc 
                                                                );

    HydrostaticForces                   hydrostatic_force(          
                                                                    mesh,
                                                                    &hydrostat_force_interf, 
                                                                    input->mass * input->grav_acc 
                                                                );

    // Create Newmark-Beta instance
    NewmarkBeta nwk(
                        &hydrostatic_force,
                        mass_mat,
                        stiff_mat,
                        damp_mat,
                        0.01,
                        0.0,
                        y0_pos,
                        y0_vel,
                        y0_acc,
                        restrictions
                    );

    // Loop over time until reach equilibrium
    while ( true )
    {
        // Advance one step in time
        nwk.step();

        // std::cout << "Time: " << nwk.time << " - rhs: " << nwk.rhs[2] << std::endl;
        if ( std::abs( nwk.rhs[2] ) < 0.1 )
        {
            break;
        }
        // std::cout << " - X: " << nwk.y_pos[0];
        // std::cout << " - Y: " << nwk.y_pos[1];
        // std::cout << " - Z: " << nwk.y_pos[2];
        // std::cout << " - RX: " << nwk.y_pos[3];
        // std::cout << " - RY: " << nwk.y_pos[4];
        // std::cout << " - RZ: " << nwk.y_pos[5] << std::endl;
        // std::cout << " - VX: " << nwk.y_vel[0];
        // std::cout << " - VY: " << nwk.y_vel[1];
        // std::cout << " - VZ: " << nwk.y_vel[2];
        // std::cout << " - VRX: " << nwk.y_vel[3];
        // std::cout << " - VRY: " << nwk.y_vel[4];
        // std::cout << " - VRZ: " << nwk.y_vel[5] << std::endl;
        // std::cout << " - AX: " << nwk.y_acc[0];
        // std::cout << " - AY: " << nwk.y_acc[1];
        // std::cout << " - AZ: " << nwk.y_acc[2];
        // std::cout << " - ARX: " << nwk.y_acc[3];
        // std::cout << " - ARY: " << nwk.y_acc[4];
        // std::cout << " - ARZ: " << nwk.y_acc[5] << std::endl;
        // std::cout << std::endl;

        // Check for time limit
        if ( nwk.time >=  max_time*10 )
        {
            break;
        }
    }

    std::cout << "Time: " << nwk.time << " - Z: " << nwk.y_pos[2] << " - rhs: " << nwk.rhs[2] << std::endl;

    // Delete heap allocations
    delete mass_mat;
    delete stiff_mat;
    delete damp_mat;

}


void test_bisection( StabInput* input, RigidBodyMesh* mesh )
{
    // Create restrictions vector
    cusfloat    y0_pos[6]       = { 0.0, 0.0, 0.0, 0.0, 5.0/57.3, 0.0 };

    // Create functor to calculate hydrostatic forces
    HydrostaticForceInterface<NUM_GP>   hydrostat_force_interf( 
                                                                    input->water_density, 
                                                                    input->grav_acc 
                                                                );

    HydrostaticForcesNLin               hydrostatic_force(          
                                                                    mesh,
                                                                    y0_pos,
                                                                    &hydrostat_force_interf, 
                                                                    input->mass * input->grav_acc 
                                                                );

    // Use bisection method to find equilibrium position
    cusfloat    lz  = mesh->z_max - mesh->z_min;
    cusfloat    sol = 0.0;
    int         info = 0;
    bisection( hydrostatic_force, -lz, lz, 0.01, 1e-6, 100, true, sol, info );

    InitialStability<NUM_GP, RigidBodyMesh> init_stab( 
                                                            input->water_density,
                                                            input->grav_acc,
                                                            sol,
                                                            input->cog,
                                                            input->rad_gyr,
                                                            mesh
                                                        );

    init_stab.print( );

}


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
    std::string auto_flush_fopath( "S:/seamotions_validation/0_seamotions/1_H50/0_Box/1_results/_track_mesh_changes" );
    cusfloat    draft       = 5.0;

    Mesh mesh_ref( input.mesh_fipath, input.body_name, input.cog, false , 0, false );
    RigidBodyMesh mesh_mov( input.mesh_fipath, input.body_name, input.cog, false , 0, draft, auto_flush_fopath );

    
    mesh_mov.write( 
                        case_fopath
                );

    mesh_mov.move( 0.0, 0.0, -draft, 30.0/57.3, 0.0, 0.0 );
    mesh_mov.check_underwater_panels( );

    mesh_mov.write_underwater_panels( 
                                        case_fopath,
                                        input.mesh_finame + "_tri"
                                    );

    InitialStability<NUM_GP, RigidBodyMesh> init_stab( 
                                                            input.water_density,
                                                            input.grav_acc,
                                                            draft,
                                                            input.cog,
                                                            input.rad_gyr,
                                                            &mesh_mov
                                                        );
    init_stab.print( );

    // test_newmark( &input, &mesh_mov );
    test_bisection( &input, &mesh_mov );
    
    // calculate_hydrostatic_properties( &input, &mesh_mov );

    // InitialStabilityParams isp;
    // calculate_initial_stability_paraters( draft, &mesh_mov, &isp );

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
        std::cout << " -> Seamotions (Stability) finished!" << std::endl;
        std::cout << " ---> Elapsed wall time for calculation [s]: " << case_timer << std::endl;
    }

    return 0;
}
