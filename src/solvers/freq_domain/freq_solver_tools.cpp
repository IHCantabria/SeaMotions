
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

// Include general usage scientific libraries
#include <complex>
#include <fstream>
#include "mkl_cblas.h"

// Include local modules
#include "freq_solver_tools.hpp"

#include "../../containers/matlin_group.hpp"
#include "../../containers/simulation_data.hpp"
#include "../../green/kochin.hpp"
// #include "../../interfaces/grf_interface.hpp"
// #include "../../interfaces/gwfdx_interface.hpp"
// #include "../../interfaces/gwfdy_interface.hpp"
// #include "../../interfaces/gwfdz_interface.hpp"
#include "../../interfaces/gwfcns_interface_t.hpp"
#include "qtf.hpp"
#include "qtf_indirect_method.hpp"
#include "../../green/integrals_database.hpp"
#include "../../solvers/freq_domain/diffraction.hpp"
#include "../../solvers/freq_domain/froude_krylov.hpp"
#include "../../solvers/freq_domain/hydromechanics.hpp"
#include "../../solvers/freq_domain/formulation_kernel_backend.hpp"
// #include "../../solvers/freq_domain/potential.hpp"
#include "../../solvers/freq_domain/panel_fields.hpp"
#include "../../solvers/freq_domain/raos.hpp"
#include "../../solvers/freq_domain/tools.hpp"
// #include "../../solvers/freq_domain/velocities.hpp"
// #include "../../solvers/freq_domain/wave_elevation.hpp"
#include "../../waves/wave_dispersion_fo.hpp"
#include "../../waves/waves_common.hpp"

#include "../../green/source.hpp"


void    calculate_influence_potmat_steady(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                MLGCmpx*        pot_gp
                                            )
{
    // Generate potential matrix
    int         count       = 0;

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

    // std::cout << "pot_gp->field_points_np: " << pot_gp->field_points_np << "\n";
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

        // std::cout << "pot_gp->start_col: " << pot_gp->start_col << " - pot_gp->end_col: " << pot_gp->end_col << "\n";
        for ( int j=pot_gp->start_col; j<pot_gp->end_col+1; j++ )
        {
            // // Compute steady and wave terms over the panel
            // if ( 
            //         mesh_gp->source_nodes[i]->panel->type == DIFFRAC_PANEL_CODE
            //         ||
            //         mesh_gp->source_nodes[j]->panel->type == LID_PANEL_CODE
            //     )
            // {
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
            pot_gp->sysmat_steady[count]  =  pot_term / 4.0 / PI;
            // std::cout << count << pot_gp->sysmat_steady[count] << std::endl;
            // }

            // if ( j==1140 && i==1330 )
            // {
            //     std::cout << "POT OUT: " << pot_gp->sysmat_steady[count] << std::endl;
            //     std::cout << "Pot_0: " << pot_0 << std::endl;
            //     std::cout << "Pot_1: " << pot_1 << std::endl;
            //     std::cout << "Pot_2: " << pot_2 << std::endl;
            //     std::cout << "Pot_3: " << pot_3 << std::endl;
            //     std::cout << "Pot_4: " << pot_4 << std::endl;
            //     std::cout << "Pot_5: " << pot_5 << std::endl;
            //     std::cout << "FP0: " << field_point_0[0] << ", " << field_point_0[1] << ", " << field_point_0[2] << std::endl;
            //     std::cout << "FP1: " << field_point_1[0] << ", " << field_point_1[1] << ", " << field_point_1[2] << std::endl;
            //     std::cout << "FP2: " << field_point_2[0] << ", " << field_point_2[1] << ", " << field_point_2[2] << std::endl;
            //     std::cout << "FP3: " << field_point_3[0] << ", " << field_point_3[1] << ", " << field_point_3[2] << std::endl;
            //     std::cout << "FP4: " << field_point_4[0] << ", " << field_point_4[1] << ", " << field_point_4[2] << std::endl;
            //     std::cout << "FP5: " << field_point_5[0] << ", " << field_point_5[1] << ", " << field_point_5[2] << std::endl;
            //     std::cout << std::endl;
            //     std::cout << "PANEL_I: " << std::endl;
            //     std::cout << *(mesh_gp->panels[j]) << std::endl;

            //     std::cout << "PANEL_J: " << std::endl;
            //     std::cout << *(mesh_gp->panels[i]) << std::endl;
            // }

            count++;
        }
    }
    
}


void    calculate_freq_domain_coeffs(
                                                MpiConfig*      mpi_config,
                                                Input*          input,
                                                Hydrostatics**  hydrostatics,
                                                Output*         output
                                    )
{
    /****************************************************/
    /*** Create total mesh for the interacting bodies ***/
    /****************************************************/

    // Group all meshes in a vector
    Mesh** all_meshes = new Mesh*[input->bodies_np];
    for ( int i=0; i<input->bodies_np; i++ )
    {
        all_meshes[i] = input->bodies[i]->mesh;
    }

    // Create new mesh from the meshes of all objects
    std::cout << "Mesh group!" << std::endl;
    MeshGroup*  mesh_gp         = new MeshGroup( 
                                                all_meshes,
                                                input->bodies_np,
                                                input->is_wl_points
                                            );

    Mesh**      fs_mesh         = nullptr;
    MeshGroup*  mesh_fs_qtf_gp  = nullptr;
    std::cout << "Mesh group done!" << std::endl;
    if ( input->is_fs_qtf )
    {
        // Create vector with the QTF FS mesh to feed the
        // MeshGroup object constructor
        fs_mesh         = new Mesh*[1];
        fs_mesh[0]      = input->bodies[0]->mesh_fs_qtf;

        // Create partition circle water line if required
        if ( input->out_qtf_so_model == 2 )
        {
            fs_mesh[0]->detect_pc_points( input->wl_det_prec );
        }

        // Create mesh group
        mesh_fs_qtf_gp  = new MeshGroup(
                                            fs_mesh,
                                            1,
                                            false
                                        );
    }

    if ( input->is_log_sin_ana )
    {
        mesh_gp->define_mirror_panels( );
    }

    /****************************************************/
    /********* Create Scalapack solver instance *********/
    /****************************************************/
    SclCmpx scl( 
                    mesh_gp->source_nodes_tnp,
                    input->dofs_np + input->heads_np,
                    mpi_config->procs_total,
                    mpi_config->proc_rank,
                    mpi_config->proc_root,
                    MPI_COMM_WORLD
                );

    /****************************************************/
    /***************** Launch solver ********************/
    /****************************************************/
    // if ( input->is_fast_solver )
    // {
    freq_domain_linear_solver(
                                    input,
                                    mpi_config,
                                    mesh_gp,
                                    mesh_fs_qtf_gp,
                                    &scl,
                                    hydrostatics,
                                    output
                            );
    // }
    // else
    // {
        // freq_domain_nonlinear_solver(
        //                                 input,
        //                                 mpi_config,
        //                                 mesh_gp,
        //                                 &scl,
        //                                 hydrostatics,
        //                                 output
        //                             );
    // }
    

    delete []   all_meshes;
    delete      mesh_fs_qtf_gp;
    delete      mesh_gp;
    delete []   fs_mesh;
}




void    calculate_global_hydstiffness(
                                                Input*          input,
                                                Hydrostatics**  hydrostatics,
                                                cusfloat*       hydstiffness
                                    )
{
    int index   = 0;
    for ( int i=0; i<input->bodies_np; i++ )
    {
        for ( int j=0; j<input->dofs_np; j++ )
        {
            for ( int k=0; k<input->dofs_np; k++ )
            {
                index               =   (
                                            i * input->bodies_np *  pow2s( input->dofs_np ) 
                                            + 
                                            i * input->dofs_np 
                                            + 
                                            j * input->dofs_np * input->bodies_np
                                            +
                                            k
                                        );
                hydstiffness[index] = hydrostatics[i]->hydstiffmat[j*input->dofs_np+k];
            }
        }
    }
}


void    calculate_global_structural_mass(
                                                Input*          input,
                                                cusfloat*       structural_mass_p0
                                        )
{
    // Allocate space for the body ith structural mass
    // matrix
    cusfloat*   body_mass   = generate_empty_vector<cusfloat>( pow2s( input->dofs_np ) );
    int index = 0;
    for ( int i=0; i<input->bodies_np; i++ )
    {
        // Clear body matrix to not get spurious data from 
        // the previous body definition
        clear_vector( pow2s( input->dofs_np ), body_mass );

        // Define body mass matrix
        body_mass[0] = input->bodies[i]->mass;  // Surge
        body_mass[7] = input->bodies[i]->mass;  // Sway
        body_mass[14] = input->bodies[i]->mass; // Heave

        if ( input->bodies[i]->interia_by_rad )
        {
            body_mass[21] = input->bodies[i]->mass * pow2s( input->bodies[i]->rad_inertia[0] ); // Roll
            body_mass[28] = input->bodies[i]->mass * pow2s( input->bodies[i]->rad_inertia[1] ); // Pitch
            body_mass[35] = input->bodies[i]->mass * pow2s( input->bodies[i]->rad_inertia[2] ); // Yaw
        }
        else
        {
            body_mass[21] = input->bodies[i]->inertia[0]; // Roll
            body_mass[22] = input->bodies[i]->inertia[1]; // Roll - Pitch
            body_mass[23] = input->bodies[i]->inertia[2]; // Roll - Yaw
            body_mass[27] = input->bodies[i]->inertia[1]; // Pitch - Roll
            body_mass[28] = input->bodies[i]->inertia[3]; // Pitch
            body_mass[29] = input->bodies[i]->inertia[4]; // Pitch - Yaw
            body_mass[33] = input->bodies[i]->inertia[2]; // Yaw - Roll
            body_mass[34] = input->bodies[i]->inertia[4]; // Yaw - Pitch-
            body_mass[35] = input->bodies[i]->inertia[5]; // Yaw
        }

        // Copy ith body structural mass matrix to the global matrix assembly
        for ( int j=0; j<input->dofs_np; j++ )
        {
            for ( int k=0; k<input->dofs_np; k++ )
            {
                index                       =   (
                                                    i * input->bodies_np *  pow2s( input->dofs_np ) 
                                                    + 
                                                    i * input->dofs_np 
                                                    + 
                                                    j * input->dofs_np * input->bodies_np
                                                    +
                                                    k
                                                );
                structural_mass_p0[index]   = body_mass[j*input->dofs_np+k];
            }
        }
    }

    // Deallocate heap memory space for variables in the
    // current function
    mkl_free( body_mass );
}


void    freq_domain_linear_solver(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                MeshGroup*      mesh_fs_qtf_gp,
                                                SclCmpx*        scl,
                                                Hydrostatics**  hydrostatics,
                                                Output*         output
                                )
{
    /****************************************************************/
    /************ Allocate space for the simulation data ************/
    /****************************************************************/
    if ( mpi_config->is_root( ) ) std::cout << "Starting Linear solver..." << std::endl;
    // Allocate space for the intensities, hydromechanics and wave exciting forces
    SimulationData* sim_data    = new SimulationData(
                                                        input,
                                                        mpi_config,
                                                        mesh_gp->meshes_np,
                                                        input->dofs_np,
                                                        input->heads_np,
                                                        scl->num_rows
                                                    );

    // Allocate kernel object
    FormulationKernelBackend<NUM_GP, PF_OFF> formulation_kernel( input, mpi_config, mesh_gp );

    // Allocate space for the system matrixes of the required potentials
    int         ipm_cols_np     = 0;
    int         ipm_sc          = 0;
    int         ipm_ed          = 0;

    ipm_sc      = scl->start_col - 1;
    ipm_ed      = ipm_sc + scl->num_cols_local - 1;
    ipm_cols_np = scl->num_cols_local;
    if ( mpi_config->is_root( ) ) std::cout << "Starting MatLinGroup..." << std::endl;
    MatLinGroup<cuscomplex>*    source_form_gp  = new MatLinGroup<cuscomplex>(
                                                                                mesh_gp->panels_tnp,
                                                                                ipm_cols_np,
                                                                                mesh_gp->meshes_np,
                                                                                ( input->dofs_np + input->heads_np ),
                                                                                0,
                                                                                mesh_gp->panels_tnp-1,
                                                                                ipm_sc,
                                                                                ipm_ed,
                                                                                true
                                                                            );

    MatLinGroup<cuscomplex>*    pot_form_gp     = new MatLinGroup<cuscomplex>(
                                                                                mesh_gp->panels_tnp,
                                                                                ipm_cols_np,
                                                                                mesh_gp->meshes_np,
                                                                                ( input->dofs_np + input->heads_np ),
                                                                                0,
                                                                                mesh_gp->panels_tnp-1,
                                                                                ipm_sc,
                                                                                ipm_ed,
                                                                                true
                                                                            );

    MatLinGroup<cuscomplex>*    potpanel_lin_gp = new MatLinGroup<cuscomplex>(
                                                                                mesh_gp->panels_tnp,
                                                                                ipm_cols_np,
                                                                                mesh_gp->meshes_np,
                                                                                ( input->dofs_np + input->heads_np ),
                                                                                0,
                                                                                mesh_gp->panels_tnp-1,
                                                                                ipm_sc,
                                                                                ipm_ed,
                                                                                true
                                                                            );
    
    MatLinGroup<cuscomplex>*    qtf_wl_we_gp        = nullptr;
    MatLinGroup<cuscomplex>*    vel_x_body_gp       = nullptr;
    MatLinGroup<cuscomplex>*    vel_y_body_gp       = nullptr;
    MatLinGroup<cuscomplex>*    vel_z_body_gp       = nullptr;
    if ( input->is_calc_mdrift )
    {
        sim_data->add_mean_drift_data( 
                                            mesh_gp->panels_raddif_tnp,
                                            mesh_gp->panels_wl_tnp,
                                            input->gauss_np_factor_2d( ),
                                            input->gauss_np_factor_1d( )
                                        );
        qtf_wl_we_gp            = new MatLinGroup<cuscomplex>(
                                                                mesh_gp->panels_wl_tnp * input->gauss_np_factor_1d( ),
                                                                ipm_cols_np,
                                                                mesh_gp->meshes_np,
                                                                ( input->dofs_np + input->heads_np ),
                                                                0,
                                                                ( mesh_gp->panels_wl_tnp * input->gauss_np_factor_1d( ) ) - 1,
                                                                ipm_sc,
                                                                ipm_ed,
                                                                false
                                                            );
        vel_x_body_gp           = new MatLinGroup<cuscomplex>(
                                                                mesh_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ),
                                                                ipm_cols_np,
                                                                mesh_gp->meshes_np,
                                                                ( input->dofs_np + input->heads_np ),
                                                                0,
                                                                ( mesh_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ) ) - 1,
                                                                ipm_sc,
                                                                ipm_ed,
                                                                false
                                                            );
        vel_y_body_gp           = new MatLinGroup<cuscomplex>(
                                                                mesh_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ),
                                                                ipm_cols_np,
                                                                mesh_gp->meshes_np,
                                                                ( input->dofs_np + input->heads_np ),
                                                                0,
                                                                ( mesh_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ) ) - 1,
                                                                ipm_sc,
                                                                ipm_ed,
                                                                false
                                                            );
        vel_z_body_gp           = new MatLinGroup<cuscomplex>(
                                                                mesh_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ),
                                                                ipm_cols_np,
                                                                mesh_gp->meshes_np,
                                                                ( input->dofs_np + input->heads_np ),
                                                                0,
                                                                ( mesh_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ) ) - 1,
                                                                ipm_sc,
                                                                ipm_ed,
                                                                false
                                                            );
    }

    MatLinGroup<cuscomplex>*    qtf_body_pot_gp     = nullptr;
    MatLinGroup<cuscomplex>*    qtf_fs_pot_gp       = nullptr;
    MatLinGroup<cuscomplex>*    qtf_fs_vel_x_gp     = nullptr;
    MatLinGroup<cuscomplex>*    qtf_fs_vel_y_gp     = nullptr;
    MatLinGroup<cuscomplex>*    qtf_fs_vel_z_gp     = nullptr;
    MatLinGroup<cuscomplex>*    qtf_pc_pot_gp       = nullptr;
    MatLinGroup<cuscomplex>*    qtf_pc_vel_x_gp     = nullptr;
    MatLinGroup<cuscomplex>*    qtf_pc_vel_y_gp     = nullptr;
    MatLinGroup<cuscomplex>*    qtf_pc_vel_z_gp     = nullptr;
    MatLinGroup<cuscomplex>*    qtf_wl_vel_x_gp     = nullptr;
    MatLinGroup<cuscomplex>*    qtf_wl_vel_y_gp     = nullptr;
    MatLinGroup<cuscomplex>*    qtf_wl_vel_z_gp     = nullptr;
    if ( input->out_qtf )
    {
        // Add data variables to storage the QTF terms values
        sim_data->add_qtf_data(
                                    input->angfreqs_np
                                );
        
        // Add data variables storage for the calculation of the QTF terms
        sim_data->add_qtf_base_data(
                                        mesh_gp->panels_raddif_tnp,
                                        input->gauss_np_factor_2d( ),
                                        mesh_gp->panels_wl_tnp,
                                        input->gauss_np_factor_1d( ),
                                        input->angfreqs_np
                                    );

        // Add data variables storage for the calculation of the second order
        // potential using the indirect method
        if ( input->out_qtf_so_model == 1 )
        {
            sim_data->add_qtf_indirect_data(
                                                mesh_gp->panels_raddif_tnp,
                                                input->gauss_np_factor_2d( ),
                                                mesh_fs_qtf_gp->panels_raddif_tnp,
                                                input->gauss_np_factor_2d( ),
                                                mesh_gp->panels_wl_tnp,
                                                input->gauss_np_factor_1d( ),
                                                input->angfreqs_np
                                            );

            qtf_body_pot_gp         = new MatLinGroup<cuscomplex>(
                                                                    mesh_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ),
                                                                    ipm_cols_np,
                                                                    mesh_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );

            qtf_fs_pot_gp           = new MatLinGroup<cuscomplex>(
                                                                    mesh_fs_qtf_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ),
                                                                    ipm_cols_np,
                                                                    mesh_fs_qtf_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_fs_qtf_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );

            qtf_fs_vel_x_gp         = new MatLinGroup<cuscomplex>(
                                                                    mesh_fs_qtf_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ),
                                                                    ipm_cols_np,
                                                                    mesh_fs_qtf_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_fs_qtf_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );

            qtf_fs_vel_y_gp         = new MatLinGroup<cuscomplex>(
                                                                    mesh_fs_qtf_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ),
                                                                    ipm_cols_np,
                                                                    mesh_fs_qtf_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_fs_qtf_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );

            qtf_fs_vel_z_gp         = new MatLinGroup<cuscomplex>(
                                                                    mesh_fs_qtf_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ),
                                                                    ipm_cols_np,
                                                                    mesh_fs_qtf_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_fs_qtf_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );

            qtf_wl_vel_x_gp         = new MatLinGroup<cuscomplex>(
                                                                    mesh_gp->panels_wl_tnp * input->gauss_np_factor_1d( ),
                                                                    ipm_cols_np,
                                                                    mesh_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_gp->panels_wl_tnp * input->gauss_np_factor_1d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );

            qtf_wl_vel_y_gp         = new MatLinGroup<cuscomplex>(
                                                                    mesh_gp->panels_wl_tnp * input->gauss_np_factor_1d( ),
                                                                    ipm_cols_np,
                                                                    mesh_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_gp->panels_wl_tnp * input->gauss_np_factor_1d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );

            qtf_wl_vel_z_gp         = new MatLinGroup<cuscomplex>(
                                                                    mesh_gp->panels_wl_tnp * input->gauss_np_factor_1d( ),
                                                                    ipm_cols_np,
                                                                    mesh_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_gp->panels_wl_tnp * input->gauss_np_factor_1d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );
        }

        // Define new fields to storage data for the direct QTF model
        if ( input->out_qtf_so_model == 2 )
        {
            sim_data->add_qtf_direct_data(
                                                mesh_gp->panels_raddif_tnp,
                                                input->gauss_np_factor_2d( ),
                                                mesh_fs_qtf_gp->panels_raddif_tnp,
                                                input->gauss_np_factor_2d( ),
                                                mesh_gp->panels_wl_tnp,
                                                input->gauss_np_factor_1d( ),
                                                mesh_fs_qtf_gp->panels_wl_tnp,
                                                input->gauss_np_factor_1d( ),
                                                input->angfreqs_np
                                            );

            qtf_fs_pot_gp           = new MatLinGroup<cuscomplex>(
                                                                    mesh_fs_qtf_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ),
                                                                    ipm_cols_np,
                                                                    mesh_fs_qtf_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_fs_qtf_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );

            qtf_fs_vel_x_gp         = new MatLinGroup<cuscomplex>(
                                                                    mesh_fs_qtf_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ),
                                                                    ipm_cols_np,
                                                                    mesh_fs_qtf_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_fs_qtf_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );

            qtf_fs_vel_y_gp         = new MatLinGroup<cuscomplex>(
                                                                    mesh_fs_qtf_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ),
                                                                    ipm_cols_np,
                                                                    mesh_fs_qtf_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_fs_qtf_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );

            qtf_fs_vel_z_gp         = new MatLinGroup<cuscomplex>(
                                                                    mesh_fs_qtf_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ),
                                                                    ipm_cols_np,
                                                                    mesh_fs_qtf_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_fs_qtf_gp->panels_raddif_tnp * input->gauss_np_factor_2d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );

            qtf_pc_pot_gp           = new MatLinGroup<cuscomplex>(
                                                                    mesh_fs_qtf_gp->panels_wl_tnp * input->gauss_np_factor_1d( ),
                                                                    ipm_cols_np,
                                                                    mesh_fs_qtf_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_fs_qtf_gp->panels_wl_tnp * input->gauss_np_factor_1d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );

            qtf_pc_vel_x_gp         = new MatLinGroup<cuscomplex>(
                                                                    mesh_fs_qtf_gp->panels_wl_tnp * input->gauss_np_factor_1d( ),
                                                                    ipm_cols_np,
                                                                    mesh_fs_qtf_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_fs_qtf_gp->panels_wl_tnp * input->gauss_np_factor_1d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );

            qtf_pc_vel_y_gp         = new MatLinGroup<cuscomplex>(
                                                                    mesh_fs_qtf_gp->panels_wl_tnp * input->gauss_np_factor_1d( ),
                                                                    ipm_cols_np,
                                                                    mesh_fs_qtf_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_fs_qtf_gp->panels_wl_tnp * input->gauss_np_factor_1d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );

            qtf_pc_vel_z_gp         = new MatLinGroup<cuscomplex>(
                                                                    mesh_fs_qtf_gp->panels_wl_tnp * input->gauss_np_factor_1d( ),
                                                                    ipm_cols_np,
                                                                    mesh_fs_qtf_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_fs_qtf_gp->panels_wl_tnp * input->gauss_np_factor_1d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );
            
            qtf_wl_vel_x_gp         = new MatLinGroup<cuscomplex>(
                                                                    mesh_gp->panels_wl_tnp * input->gauss_np_factor_1d( ),
                                                                    ipm_cols_np,
                                                                    mesh_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_gp->panels_wl_tnp * input->gauss_np_factor_1d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );

            qtf_wl_vel_y_gp         = new MatLinGroup<cuscomplex>(
                                                                    mesh_gp->panels_wl_tnp * input->gauss_np_factor_1d( ),
                                                                    ipm_cols_np,
                                                                    mesh_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_gp->panels_wl_tnp * input->gauss_np_factor_1d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );

            qtf_wl_vel_z_gp         = new MatLinGroup<cuscomplex>(
                                                                    mesh_gp->panels_wl_tnp * input->gauss_np_factor_1d( ),
                                                                    ipm_cols_np,
                                                                    mesh_gp->meshes_np,
                                                                    ( input->dofs_np + input->heads_np ),
                                                                    0,
                                                                    ( mesh_gp->panels_wl_tnp * input->gauss_np_factor_1d( ) ) - 1,
                                                                    ipm_sc,
                                                                    ipm_ed,
                                                                    false
                                                                );
        }
    }

    /****************************************************************/
    /************** Define field points for evaluation **************/
    /****************************************************************/

    // Define field points to calculate potential influence matrix
    int _count_pot_np   = 0;
    for ( int i=0; i<mesh_gp->panels_tnp; i++ )
    {
        if ( mesh_gp->panels[i]->type == DIFFRAC_PANEL_CODE || mesh_gp->panels[i]->type == LID_PANEL_CODE )
        {
            copy_vector( 3, mesh_gp->panels[i]->center, &(potpanel_lin_gp->field_points[3*_count_pot_np]) );
            _count_pot_np++;
        }
    }
    copy_vector( mesh_gp->meshes_np+1, mesh_gp->panels_raddif_cnp, potpanel_lin_gp->field_points_cnp );

    // Define field points to calculate the potential on the WL to evaluate
    // second order forces
    if ( input->is_calc_mdrift )
    {
        define_gauss_points_wl(
                                    input,
                                    mesh_gp,
                                    qtf_wl_we_gp
                                );

        // Define field points for the evaluation of the velocities field
        define_gauss_points_diffrac_panels( input, mesh_gp, vel_x_body_gp );
        define_gauss_points_diffrac_panels( input, mesh_gp, vel_y_body_gp );
        define_gauss_points_diffrac_panels( input, mesh_gp, vel_z_body_gp );

    }

    if ( input->out_qtf_so_model == 1 )
    {
        define_gauss_points_diffrac_panels( 
                                                input, 
                                                mesh_gp, 
                                                qtf_body_pot_gp 
                                            );

        define_gauss_points_diffrac_panels( 
                                                input, 
                                                mesh_fs_qtf_gp, 
                                                qtf_fs_pot_gp 
                                            );

        define_gauss_points_diffrac_panels( 
                                                input, 
                                                mesh_fs_qtf_gp, 
                                                qtf_fs_vel_x_gp 
                                            );

        define_gauss_points_diffrac_panels( 
                                                input, 
                                                mesh_fs_qtf_gp, 
                                                qtf_fs_vel_y_gp 
                                            );

        define_gauss_points_diffrac_panels( 
                                                input, 
                                                mesh_fs_qtf_gp, 
                                                qtf_fs_vel_z_gp 
                                            );

        define_gauss_points_wl(
                                    input,
                                    mesh_gp,
                                    qtf_wl_vel_x_gp
                                );

        define_gauss_points_wl(
                                    input,
                                    mesh_gp,
                                    qtf_wl_vel_y_gp
                                );

        define_gauss_points_wl(
                                    input,
                                    mesh_gp,
                                    qtf_wl_vel_z_gp
                                );
    }

    if ( input->out_qtf_so_model == 2)
    {
        // Define gauss points over the FS
        define_gauss_points_diffrac_panels( 
                                                input, 
                                                mesh_fs_qtf_gp, 
                                                qtf_fs_pot_gp 
                                            );

        define_gauss_points_diffrac_panels( 
                                                input, 
                                                mesh_fs_qtf_gp, 
                                                qtf_fs_vel_x_gp 
                                            );

        define_gauss_points_diffrac_panels( 
                                                input, 
                                                mesh_fs_qtf_gp, 
                                                qtf_fs_vel_y_gp 
                                            );

        define_gauss_points_diffrac_panels( 
                                                input, 
                                                mesh_fs_qtf_gp, 
                                                qtf_fs_vel_z_gp 
                                            );
        
        // Define WL points over the PC
        define_gauss_points_wl(
                                    input,
                                    mesh_fs_qtf_gp,
                                    qtf_pc_pot_gp
                                );

        define_gauss_points_wl(
                                    input,
                                    mesh_fs_qtf_gp,
                                    qtf_pc_vel_x_gp
                                );

        define_gauss_points_wl(
                                    input,
                                    mesh_fs_qtf_gp,
                                    qtf_pc_vel_y_gp
                                );

        define_gauss_points_wl(
                                    input,
                                    mesh_fs_qtf_gp,
                                    qtf_pc_vel_z_gp
                                );

        // Define WL points over the body
        define_gauss_points_wl(
                                    input,
                                    mesh_gp,
                                    qtf_wl_vel_x_gp
                                );

        define_gauss_points_wl(
                                    input,
                                    mesh_gp,
                                    qtf_wl_vel_y_gp
                                );

        define_gauss_points_wl(
                                    input,
                                    mesh_gp,
                                    qtf_wl_vel_z_gp
                                );
    }

    /****************************************************************/
    /*************** Create Green function interface ****************/
    /****************************************************************/
    GWFcnsInterfaceT<NUM_GP*NUM_GP> gwfcns_interf(
                                                    mesh_gp->source_nodes[0],
                                                    mesh_gp->source_nodes[0],
                                                    input->angfreqs[0],
                                                    input->water_depth,
                                                    input->grav_acc
                                                );

    /****************************************************************/
    /************ Calculate global structural mass matrix ***********/
    /****************************************************************/
    if ( mpi_config->is_root( ) )
    {
        calculate_global_structural_mass(
                                            input,
                                            sim_data->structural_mass_p0
                                        );
    }


    /****************************************************************/
    /************* Calculate global hydrostatic matrix **************/
    /****************************************************************/
    if ( mpi_config->is_root( ) )
    {
        calculate_global_hydstiffness(
                                            input,
                                            hydrostatics,
                                            sim_data->hydrostiff_p0
                                        );
    }

    /****************************************************************/
    /***** Calculate steady contributions to system matrixes ********/
    /****************************************************************/

    // // Calculate steady part of the sources influence matrix
    // calculate_gf_intensity_steady_sysmat_lin(
    //                                             input,
    //                                             scl,
    //                                             mesh_gp,
    //                                             sysmat_steady
    //                                         );
    
    // // Calculate steady part of the potential influence matrix matrix
    // calculate_influence_potmat_steady(
    //                                         input,
    //                                         mesh_gp,
    //                                         potpanel_lin_gp
    //                                 );
    
    // // Calculate steady part of the potential influence matrix to calculate
    // // the mean drift
    // if ( input->is_calc_mdrift )
    // {
    //     calculate_influence_potmat_steady(
    //                                                 input,
    //                                                 mesh_gp,
    //                                                 qtf_wl_we_gp
    //                                         );

    //     calculate_raddif_velocity_mat_steady(
    //                                                 input,
    //                                                 mesh_gp,
    //                                                 vel_x_body_gp,
    //                                                 vel_y_body_gp,
    //                                                 vel_z_body_gp
    //                                         );
    // }

    // if ( input->out_qtf_so_model == 1 )
    // {
    //     // Calculate potential over the body
    //     calculate_influence_potmat_steady(
    //                                             input,
    //                                             mesh_gp,
    //                                             qtf_body_pot_gp
    //                                         );

    //     // Calculate potential over the FS
    //     calculate_influence_potmat_steady(
    //                                             input,
    //                                             mesh_fs_qtf_gp,
    //                                             qtf_fs_pot_gp
    //                                         );

    //     // Calculate velocities over the FS
    //     calculate_raddif_velocity_mat_steady(
    //                                             input,
    //                                             mesh_fs_qtf_gp,
    //                                             qtf_fs_vel_x_gp,
    //                                             qtf_fs_vel_y_gp,
    //                                             qtf_fs_vel_z_gp
    //                                         );

    //     // Calculate velocities over the WL
    //     calculate_raddif_velocity_mat_steady(
    //                                             input,
    //                                             mesh_gp,
    //                                             qtf_wl_vel_x_gp,
    //                                             qtf_wl_vel_y_gp,
    //                                             qtf_wl_vel_z_gp
    //                                         );
    // }

    // if ( input->out_qtf_so_model == 2 )
    // {
    //     // Calcualte velocities over the FS
    //     calculate_influence_potmat_steady(
    //                                             input,
    //                                             mesh_fs_qtf_gp,
    //                                             qtf_fs_pot_gp
    //                                         );

    //     calculate_raddif_velocity_mat_steady(
    //                                             input,
    //                                             mesh_fs_qtf_gp,
    //                                             qtf_fs_vel_x_gp,
    //                                             qtf_fs_vel_y_gp,
    //                                             qtf_fs_vel_z_gp
    //                                         );

    //     // Calculate potential ove the PC
    //     calculate_influence_potmat_steady(
    //                                             input,
    //                                             mesh_fs_qtf_gp,
    //                                             qtf_pc_pot_gp
    //                                         );

    //     // Calculate velocities over the PC
    //     calculate_raddif_velocity_mat_steady(
    //                                             input,
    //                                             mesh_gp,
    //                                             qtf_pc_vel_x_gp,
    //                                             qtf_pc_vel_y_gp,
    //                                             qtf_pc_vel_z_gp
    //                                         );
        
    //     // Calculate velocities over the Body WL
    //     calculate_raddif_velocity_mat_steady(
    //                                             input,
    //                                             mesh_gp,
    //                                             qtf_wl_vel_x_gp,
    //                                             qtf_wl_vel_y_gp,
    //                                             qtf_wl_vel_z_gp
    //                                         );
    // }

    /****************************************************************/
    /******* Calculate wave contributions to system matrixes ********/
    /********** to obtain the hydrodynamic coefficients *************/
    /****************************************************************/

    // Loop over frequencies
    for ( int i=0; i<input->angfreqs_np; i++ )
    {
        MPI_Barrier( MPI_COMM_WORLD );
        double freq_tstart = MPI_Wtime( );


        if ( mpi_config->is_root( ) ) std::cout << "Calculating angular frequency: " << input->angfreqs[i] << std::endl;

        formulation_kernel.solve( input->angfreqs[i] );
        formulation_kernel.update_results( sim_data );

        // Recalculate wave properties for the current 
        // angular frequency
        // gwf_interf->set_ang_freq( input->angfreqs[i] );
        // gwf_dn_interf->set_ang_freq( input->angfreqs[i] );
        
        // gwfcns_interf.set_ang_freq( input->angfreqs[i] );

        // // Fold for current frequency and water depth
        // cusfloat H = pow2s( input->angfreqs[i] ) * input->water_depth / input->grav_acc;
        // fold_database( H );

        // // Calculate sources intensity
        // MPI_Barrier( MPI_COMM_WORLD );
        // double mat_tstart = MPI_Wtime( );
        // calculate_gf_intensity_sysmat_v3(
        //                                     input,
        //                                     mpi_config,
        //                                     scl,
        //                                     mesh_gp,
        //                                     gwfcns_interf,
        //                                     input->angfreqs[i],
        //                                     source_form_gp,
        //                                     pot_form_gp,
        //                                     potpanel_lin_gp,
        //                                     sim_data->intensities
        //                                 );
        // double mat_tend = MPI_Wtime( );
        // // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - Calculating system matrixes... -> Done!" << std::endl;
        
        // // Gather source values from each processor
        // // MPI_Barrier( MPI_COMM_WORLD );
        // // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - BroadCast 1" << std::endl;
        // MPI_Bcast(
        //             sim_data->intensities,
        //             scl->num_rows * ( input->dofs_np + input->heads_np ),
        //             mpi_cuscomplex,
        //             mpi_config->proc_root,
        //             MPI_COMM_WORLD                
        //         );

        // MPI_Bcast(
        //             pot_form_gp->field_values,
        //             scl->num_rows * ( input->dofs_np + input->heads_np ),
        //             mpi_cuscomplex,
        //             mpi_config->proc_root,
        //             MPI_COMM_WORLD                
        //         );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - BroadCast 1 -> Done!" << std::endl;

        // Calculate potential influence coeffcients matrix
        // calculate_influence_field_mat(
        //                                 input,
        //                                 mesh_gp,
        //                                 gwf_interf,
        //                                 potpanel_lin_gp
        //                             );

        // Calculate panels potential
        // MPI_Barrier( MPI_COMM_WORLD );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - calculate_fields_raddif_lin" << std::endl;
        // calculate_fields_raddif_lin(
        //                                 input,
        //                                 sim_data->intensities,
        //                                 potpanel_lin_gp
        //                             );

        // // std::stringstream ss_pot_sf;
        // // ss_pot_sf << "pot_sf_" << mpi_config->proc_rank << "_" << input->angfreqs[i] << ".dat";
        // // std::ofstream out_pot_sf( ss_pot_sf.str( ) );
        // // int count_1 = 0;
        // // for ( int i=0; i<input->dofs_np+input->heads_np; i++ )
        // // {
        // //     for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
        // //     {
                
        // //         out_pot_sf << potpanel_lin_gp->field_values[count_1].real( ) << " " << potpanel_lin_gp->field_values[count_1].imag( ) << "\n";
        // //         count_1++;
        // //     }
        // // }
        // // out_pot_sf.close( );
        // // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - calculate_fields_raddif_lin -> Done!" << std::endl;

        // clear_vector( 
        //                 mesh_gp->panels_tnp * ( input->dofs_np + input->heads_np ),
        //                 sim_data->panels_potential
        //             );

        // // MPI_Barrier( MPI_COMM_WORLD );
        // // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - MPI_Allreduce" << std::endl;
        // MPI_Allreduce(
        //                 potpanel_lin_gp->field_values,
        //                 sim_data->panels_potential,
        //                 mesh_gp->panels_tnp * ( input->dofs_np + input->heads_np ),
        //                 mpi_cuscomplex,
        //                 MPI_SUM,
        //                 MPI_COMM_WORLD
        //             );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - MPI_Allreduce -> Done!" << std::endl;

        // for ( int k=0; k<(input->dofs_np+input->heads_np); k++ )
        // {
        //     for ( int i=0; i<mesh_gp->panels_tnp; i++ )
        //     {
        //         if ( mesh_gp->panels[i]->type == LID_PANEL_CODE )
        //         {
        //             sim_data->panels_potential[k*mesh_gp->panels_tnp+i] = cuscomplex( 0.0, 0.0 );
        //         }
        //     }
        // }

        // Calculate added mass and damping
        // MPI_Barrier( MPI_COMM_WORLD );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - calculate_hydromechanic_coeffs_lin" << std::endl;
        calculate_hydromechanic_coeffs_lin( 
                                                input,
                                                mpi_config,
                                                mesh_gp,
                                                sim_data->panels_potential,
                                                input->angfreqs[i],
                                                sim_data->added_mass,
                                                sim_data->damping_rad,
                                                sim_data->panels_pressure
                                            );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - calculate_hydromechanic_coeffs_lin -> Done!" << std::endl;

        // Calculate diffraction forces
        // MPI_Barrier( MPI_COMM_WORLD );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - calculate_diffraction_forces_lin" << std::endl;
        calculate_diffraction_forces_lin(
                                                input,
                                                mpi_config,
                                                mesh_gp,
                                                sim_data->panels_potential,
                                                input->angfreqs[i],
                                                sim_data->wave_diffrac,
                                                sim_data->panels_pressure
                                        );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - calculate_diffraction_forces_lin -> Done!" << std::endl;
        
        
        // Calculate Froude-Krylov forces
        // MPI_Barrier( MPI_COMM_WORLD );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - calculate_froude_krylov_fo" << std::endl;
        calculate_froude_krylov_fo(
                                    input,
                                    mpi_config,
                                    mesh_gp,
                                    input->angfreqs[i],
                                    sim_data->froude_krylov
                                );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - calculate_froude_krylov_fo -> Done!" << std::endl;

        // Join data from all processors
        // MPI_Barrier( MPI_COMM_WORLD );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - reduce_added_mass" << std::endl;
        MPI_Reduce(
                            sim_data->added_mass,
                            sim_data->added_mass_p0,
                            sim_data->hydmech_np,
                            mpi_cusfloat,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                        );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - reduce_added_mass -> Done!" << std::endl;

        // MPI_Barrier( MPI_COMM_WORLD );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - reduce_damping_rad" << std::endl;
        MPI_Reduce(
                            sim_data->damping_rad,
                            sim_data->damping_rad_p0,
                            sim_data->hydmech_np,
                            mpi_cusfloat,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                    );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - reduce_damping_rad -> Done!" << std::endl;

        // MPI_Barrier( MPI_COMM_WORLD );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - reduce_wave_exc" << std::endl;
        MPI_Reduce(
                            sim_data->froude_krylov,
                            sim_data->froude_krylov_p0,
                            sim_data->wave_exc_np,
                            mpi_cuscomplex,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                        );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - reduce_wave_exc -> Done!" << std::endl;

        // MPI_Barrier( MPI_COMM_WORLD );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - reduce_wave_diffrac" << std::endl;
        MPI_Reduce(
                            sim_data->wave_diffrac,
                            sim_data->wave_diffrac_p0,
                            sim_data->wave_exc_np,
                            mpi_cuscomplex,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                        );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - reduce_wave_diffrac -> Done!" << std::endl;

        // Calculate total wave exciting forces
        // MPI_Barrier( MPI_COMM_WORLD );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - Calculate total wave exciting forces" << std::endl;
        if ( mpi_config->is_root( ) )
        {
            sv_add(
                        sim_data->wave_exc_np,
                        sim_data->wave_diffrac_p0,
                        sim_data->froude_krylov_p0,
                        sim_data->wave_exc_p0
                    );
        }
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - Calculate total wave exciting forces -> Done!" << std::endl;

        // Calculate raos
        // MPI_Barrier( MPI_COMM_WORLD );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - calculate_raos" << std::endl;
        if ( mpi_config->is_root( ) )
        {
            calculate_raos(
                                input,
                                sim_data->structural_mass_p0,
                                sim_data->added_mass_p0,
                                sim_data->damping_rad_p0,
                                sim_data->hydrostiff_p0,
                                sim_data->wave_diffrac_p0,
                                sim_data->froude_krylov_p0,
                                input->angfreqs[i],
                                sim_data->raos
                            );
        }
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - calculate_raos -> Done!" << std::endl;

        // Calculate mean drift
        // if ( input->is_calc_mdrift )
        // {
        //     // Update integration interfaces status to the current angular frequency
        //     gwf_dx_interf->set_ang_freq( input->angfreqs[i] );
        //     gwf_dy_interf->set_ang_freq( input->angfreqs[i] );
        //     gwf_dz_interf->set_ang_freq( input->angfreqs[i] );

        //     // Broadcast RAOs values to be available in all the processes
        //     MPI_Bcast(
        //                 sim_data->raos,
        //                 sim_data->wave_exc_np,
        //                 mpi_cuscomplex,
        //                 mpi_config->proc_root,
        //                 MPI_COMM_WORLD    
        //             );
            
        //     // Calculate total potential at the target WL points
        //     // calculate_potpanel_total_lin(
        //     //                                 input,
        //     //                                 mpi_config,
        //     //                                 mesh_gp,
        //     //                                 input->angfreqs[i],
        //     //                                 sim_data->intensities,
        //     //                                 sim_data->raos,
        //     //                                 qtf_wl_we_gp,
        //     //                                 sim_data->qtf_wl_we_pot_total
        //     //                             );

            // calculate_fields_lin(
            //                         input,
            //                         mpi_config,
            //                         mesh_gp,
            //                         gwf_interf,
            //                         wave_potential_fo_space,
            //                         input->angfreqs[i],
            //                         sim_data->intensities,
            //                         sim_data->raos,
            //                         qtf_wl_we_gp,
            //                         sim_data->mdrift_wl_pot_fk,
            //                         sim_data->mdrift_wl_pot_raddif,
            //                         sim_data->mdrift_wl_pot_total
            //                     );
            
        //     // Calculate velocities over panels
            // calculate_velocities_total(
            //                             input,
            //                             mpi_config,
            //                             mesh_gp,
            //                             input->angfreqs[i],
            //                             sim_data->intensities,
            //                             sim_data->raos,
            //                             vel_x_body_gp,
            //                             vel_y_body_gp,
            //                             vel_z_body_gp,
            //                             sim_data->mdrift_press_vel_x,
            //                             sim_data->mdrift_press_vel_y,
            //                             sim_data->mdrift_press_vel_z
            //                         );
            
        //     calculate_fields_lin(
        //                             input,
        //                             mpi_config,
        //                             mesh_gp,
        //                             gwf_dx_interf,
        //                             wave_potential_fo_space_dx,
        //                             input->angfreqs[i],
        //                             sim_data->intensities,
        //                             sim_data->raos,
        //                             vel_x_body_gp,
        //                             sim_data->mdrift_body_vel_x_fk,
        //                             sim_data->mdrift_body_vel_x_raddif,
        //                             sim_data->mdrift_body_vel_x_total
        //                         );
            
        //     calculate_fields_lin(
        //                             input,
        //                             mpi_config,
        //                             mesh_gp,
        //                             gwf_dy_interf,
        //                             wave_potential_fo_space_dy,
        //                             input->angfreqs[i],
        //                             sim_data->intensities,
        //                             sim_data->raos,
        //                             vel_y_body_gp,
        //                             sim_data->mdrift_body_vel_y_fk,
        //                             sim_data->mdrift_body_vel_y_raddif,
        //                             sim_data->mdrift_body_vel_y_total
        //                         );

        //     calculate_fields_lin(
        //                             input,
        //                             mpi_config,
        //                             mesh_gp,
        //                             gwf_dz_interf,
        //                             wave_potential_fo_space_dz,
        //                             input->angfreqs[i],
        //                             sim_data->intensities,
        //                             sim_data->raos,
        //                             vel_z_body_gp,
        //                             sim_data->mdrift_body_vel_z_fk,
        //                             sim_data->mdrift_body_vel_z_raddif,
        //                             sim_data->mdrift_body_vel_z_total
        //                         );

        //     if ( mpi_config->is_root( ) )
        //     {
        //         // Calculate relative wave elevation
        //         calculate_relative_wave_elevation_lin(
        //                                                     input,
        //                                                     qtf_wl_we_gp,
        //                                                     sim_data->mdrift_wl_pot_total,
        //                                                     input->angfreqs[i],
        //                                                     sim_data->raos,
        //                                                     sim_data->mdrift_wl_rel_we
        //                                                 );

        //         // Calculate mean drift forces
        //         calculate_qtf_terms_force(
        //                                         input,
        //                                         mesh_gp,
        //                                         QTF_DIFF_CODE,
        //                                         sim_data->mdrift_wl_rel_we,
        //                                         sim_data->mdrift_wl_rel_we,
        //                                         sim_data->raos,
        //                                         sim_data->raos,
        //                                         sim_data->mdrift_body_vel_x_total,
        //                                         sim_data->mdrift_body_vel_y_total,
        //                                         sim_data->mdrift_body_vel_z_total,
        //                                         sim_data->mdrift_body_vel_x_total,
        //                                         sim_data->mdrift_body_vel_y_total,
        //                                         sim_data->mdrift_body_vel_z_total,
        //                                         input->angfreqs[i],
        //                                         input->angfreqs[i],
        //                                         sim_data->mdrift,
        //                                         sim_data->mdrift_wl,
        //                                         sim_data->mdrift_bern,
        //                                         sim_data->mdrift_acc,
        //                                         sim_data->mdrift_mom,
        //                                         qtf_wl_we_gp,
        //                                         vel_x_body_gp,
        //                                         false
        //                                     );

        //         // Check if QTF calcuation is required to storage panel information
        //         if ( input->out_qtf )
        //         {
        //             sim_data->storage_qtf_base_freq(
        //                                                 i,
        //                                                 sim_data->mdrift_body_vel_x_total,
        //                                                 sim_data->mdrift_body_vel_y_total,
        //                                                 sim_data->mdrift_body_vel_z_total,
        //                                                 sim_data->raos,
        //                                                 sim_data->mdrift_wl_rel_we
        //                                             );

        //             if ( input->out_qtf_so_model == 1 )
        //             {
        //                 // Calculate potential field over the body
        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_gp,
        //                                         gwf_interf,
        //                                         wave_potential_fo_space,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_body_pot_gp,
        //                                         sim_data->mdrift_body_pot_fk,
        //                                         sim_data->mdrift_body_pot_raddif,
        //                                         sim_data->mdrift_body_pot_total
        //                                     );

        //                 // Calculate potential field over the fs
        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_fs_qtf_gp,
        //                                         gwf_interf,
        //                                         wave_potential_fo_space,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_fs_pot_gp,
        //                                         sim_data->mdrift_fs_pot_fk,
        //                                         sim_data->mdrift_fs_pot_raddif,
        //                                         sim_data->mdrift_fs_pot_total
        //                                     );

        //                 // Calculate velocities field over the free surface
        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_fs_qtf_gp,
        //                                         gwf_dx_interf,
        //                                         wave_potential_fo_space_dx,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_fs_vel_x_gp,
        //                                         sim_data->mdrift_fs_vel_x_fk,
        //                                         sim_data->mdrift_fs_vel_x_raddif,
        //                                         sim_data->mdrift_fs_vel_x_total
        //                                     );

        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_fs_qtf_gp,
        //                                         gwf_dy_interf,
        //                                         wave_potential_fo_space_dy,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_fs_vel_y_gp,
        //                                         sim_data->mdrift_fs_vel_y_fk,
        //                                         sim_data->mdrift_fs_vel_y_raddif,
        //                                         sim_data->mdrift_fs_vel_y_total
        //                                     );

        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_fs_qtf_gp,
        //                                         gwf_dz_interf,
        //                                         wave_potential_fo_space_dz,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_fs_vel_z_gp,
        //                                         sim_data->mdrift_fs_vel_z_fk,
        //                                         sim_data->mdrift_fs_vel_z_raddif,
        //                                         sim_data->mdrift_fs_vel_z_total
        //                                     );

        //                 // Calculate velocities field over the WL
        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_gp,
        //                                         gwf_dx_interf,
        //                                         wave_potential_fo_space_dx,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_wl_vel_x_gp,
        //                                         sim_data->mdrift_wl_vel_x_fk,
        //                                         sim_data->mdrift_wl_vel_x_raddif,
        //                                         sim_data->mdrift_wl_vel_x_total
        //                                     );
                        
        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_gp,
        //                                         gwf_dy_interf,
        //                                         wave_potential_fo_space_dy,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_wl_vel_y_gp,
        //                                         sim_data->mdrift_wl_vel_y_fk,
        //                                         sim_data->mdrift_wl_vel_y_raddif,
        //                                         sim_data->mdrift_wl_vel_y_total
        //                                     );

        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_gp,
        //                                         gwf_dz_interf,
        //                                         wave_potential_fo_space_dz,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_wl_vel_z_gp,
        //                                         sim_data->mdrift_wl_vel_z_fk,
        //                                         sim_data->mdrift_wl_vel_z_raddif,
        //                                         sim_data->mdrift_wl_vel_z_total
        //                                     );

        //                 // Calculate Kochin coefficients for the far field integral
        //                 calculate_kochin_pert_coeffs(
        //                                                 input,
        //                                                 mesh_gp,
        //                                                 i,
        //                                                 sim_data
        //                                             );

        //                 calculate_kochin_rad_coeffs(
        //                                                 input,
        //                                                 mesh_gp,
        //                                                 i,
        //                                                 sim_data
        //                                             );
                        
        //                 // Storage parameters necessary for the second order potential calculation
        //                 // using the indirect method
        //                 sim_data->storage_qtf_indirect_freq(
        //                                                         i,
        //                                                         sim_data->mdrift_body_pot_raddif,
        //                                                         sim_data->mdrift_body_vel_x_raddif,
        //                                                         sim_data->mdrift_body_vel_y_raddif,
        //                                                         sim_data->mdrift_body_vel_z_raddif,
        //                                                         sim_data->mdrift_fs_pot_fk,
        //                                                         sim_data->mdrift_fs_pot_raddif,
        //                                                         sim_data->mdrift_fs_pot_total,
        //                                                         sim_data->mdrift_fs_vel_x_fk,
        //                                                         sim_data->mdrift_fs_vel_y_fk,
        //                                                         sim_data->mdrift_fs_vel_z_fk,
        //                                                         sim_data->mdrift_fs_vel_x_total,
        //                                                         sim_data->mdrift_fs_vel_y_total,
        //                                                         sim_data->mdrift_fs_vel_z_total,
        //                                                         sim_data->mdrift_wl_pot_raddif,
        //                                                         sim_data->mdrift_wl_vel_x_total,
        //                                                         sim_data->mdrift_wl_vel_y_total,
        //                                                         sim_data->mdrift_wl_vel_z_total,
        //                                                         sim_data->mdrift_kochin_pert_cos,
        //                                                         sim_data->mdrift_kochin_pert_sin,
        //                                                         sim_data->mdrift_kochin_rad_cos,
        //                                                         sim_data->mdrift_kochin_rad_sin
        //                                                     );
        //             }

        //             if ( input->out_qtf_so_model == 2 )
        //             {
        //                 // Calculate potential field over the fs
        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_fs_qtf_gp,
        //                                         gwf_interf,
        //                                         wave_potential_fo_space,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_fs_pot_gp,
        //                                         sim_data->mdrift_fs_pot_fk,
        //                                         sim_data->mdrift_fs_pot_raddif,
        //                                         sim_data->mdrift_fs_pot_total
        //                                     );

        //                 // Calculate velocities field over the free surface
        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_fs_qtf_gp,
        //                                         gwf_dx_interf,
        //                                         wave_potential_fo_space_dx,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_fs_vel_x_gp,
        //                                         sim_data->mdrift_fs_vel_x_fk,
        //                                         sim_data->mdrift_fs_vel_x_raddif,
        //                                         sim_data->mdrift_fs_vel_x_total
        //                                     );

        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_fs_qtf_gp,
        //                                         gwf_dy_interf,
        //                                         wave_potential_fo_space_dy,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_fs_vel_y_gp,
        //                                         sim_data->mdrift_fs_vel_y_fk,
        //                                         sim_data->mdrift_fs_vel_y_raddif,
        //                                         sim_data->mdrift_fs_vel_y_total
        //                                     );

        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_fs_qtf_gp,
        //                                         gwf_dz_interf,
        //                                         wave_potential_fo_space_dz,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_fs_vel_z_gp,
        //                                         sim_data->mdrift_fs_vel_z_fk,
        //                                         sim_data->mdrift_fs_vel_z_raddif,
        //                                         sim_data->mdrift_fs_vel_z_total
        //                                     );
                        
        //                 // Calculate potential field over the PC
        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_gp,
        //                                         gwf_dx_interf,
        //                                         wave_potential_fo_space,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_pc_pot_gp,
        //                                         sim_data->mdrift_pc_pot_fk,
        //                                         sim_data->mdrift_pc_pot_raddif,
        //                                         sim_data->mdrift_pc_pot_total
        //                                     );
                        
        //                 // Calculate velocities field over the PC
        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_gp,
        //                                         gwf_dx_interf,
        //                                         wave_potential_fo_space_dx,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_pc_vel_x_gp,
        //                                         sim_data->mdrift_pc_vel_x_fk,
        //                                         sim_data->mdrift_pc_vel_x_raddif,
        //                                         sim_data->mdrift_pc_vel_x_total
        //                                     );
                        
        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_gp,
        //                                         gwf_dy_interf,
        //                                         wave_potential_fo_space_dy,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_pc_vel_y_gp,
        //                                         sim_data->mdrift_pc_vel_y_fk,
        //                                         sim_data->mdrift_pc_vel_y_raddif,
        //                                         sim_data->mdrift_pc_vel_y_total
        //                                     );

        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_gp,
        //                                         gwf_dz_interf,
        //                                         wave_potential_fo_space_dz,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_pc_vel_z_gp,
        //                                         sim_data->mdrift_pc_vel_z_fk,
        //                                         sim_data->mdrift_pc_vel_z_raddif,
        //                                         sim_data->mdrift_pc_vel_z_total
        //                                     );
                        
        //                 // Calculate velocities field over the WL
        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_gp,
        //                                         gwf_dx_interf,
        //                                         wave_potential_fo_space_dx,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_wl_vel_x_gp,
        //                                         sim_data->mdrift_wl_vel_x_fk,
        //                                         sim_data->mdrift_wl_vel_x_raddif,
        //                                         sim_data->mdrift_wl_vel_x_total
        //                                     );
                        
        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_gp,
        //                                         gwf_dy_interf,
        //                                         wave_potential_fo_space_dy,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_wl_vel_y_gp,
        //                                         sim_data->mdrift_wl_vel_y_fk,
        //                                         sim_data->mdrift_wl_vel_y_raddif,
        //                                         sim_data->mdrift_wl_vel_y_total
        //                                     );

        //                 calculate_fields_lin(
        //                                         input,
        //                                         mpi_config,
        //                                         mesh_gp,
        //                                         gwf_dz_interf,
        //                                         wave_potential_fo_space_dz,
        //                                         input->angfreqs[i],
        //                                         sim_data->intensities,
        //                                         sim_data->raos,
        //                                         qtf_wl_vel_z_gp,
        //                                         sim_data->mdrift_wl_vel_z_fk,
        //                                         sim_data->mdrift_wl_vel_z_raddif,
        //                                         sim_data->mdrift_wl_vel_z_total
        //                                     );

        //                 // Storage results
        //                 sim_data->storage_qtf_direct_freq( i );
        //             }
        //         }
        //     }
        // }

        // Gather sources
        // MPI_Barrier( MPI_COMM_WORLD );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - Gather_intensities" << std::endl;
        double save_tstart = MPI_Wtime( );
        // if ( input->out_sources )
        // {
        //     MPI_Gather( 
        //                     sim_data->intensities,
        //                     scl->num_rows * ( input->dofs_np + input->heads_np ),
        //                     mpi_cuscomplex,
        //                     sim_data->intensities_p0,
        //                     scl->num_rows * ( input->dofs_np + input->heads_np ),
        //                     mpi_cuscomplex,
        //                     mpi_config->proc_root,
        //                     MPI_COMM_WORLD
        //                 );
        // }
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - Gather_intensities -> Done!" << std::endl;

        // Gather panels pressure if any
        // MPI_Barrier( MPI_COMM_WORLD );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - Reduce_panels_pressure" << std::endl;
        if ( input->out_pressure )
        {
            MPI_Reduce(
                            sim_data->panels_pressure,
                            sim_data->panels_pressure_p0,
                            scl->num_rows * ( input->dofs_np + input->heads_np ),
                            mpi_cuscomplex,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                        );
        }
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - Reduce_panels_pressure -> Done!" << std::endl;

        // Output values to disk
        // MPI_Barrier( MPI_COMM_WORLD );
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - Save_to_disk" << std::endl;
        if ( mpi_config->is_root( ) )
        {
            if ( input->out_hydmech )
            {
                output->save_hydromechanics_format(
                                                        i,
                                                        _DN_ADDED_MASS,
                                                        sim_data->added_mass_p0
                                                    );

                output->save_hydromechanics_format(
                                                        i,
                                                        _DN_DAMPING_RAD,
                                                        sim_data->damping_rad_p0
                                                    );
            }

            if ( input->out_diffrac )
            {
                output->save_wave_exciting_format(
                                                    i,
                                                    _DN_DIFFRAC,
                                                    sim_data->wave_diffrac_p0
                                                );
            }

            if ( input->out_fk )
            {
                output->save_wave_exciting_format(
                                                    i,
                                                    _DN_FK,
                                                    sim_data->froude_krylov_p0
                                                );
            }

            if ( input->out_wex )
            {
                output->save_wave_exciting_format(
                                                    i,
                                                    _DN_WEX,
                                                    sim_data->wave_exc_p0
                                                );
            }

            if ( input->out_raos )
            {
                output->save_wave_exciting_format(
                                                    i,
                                                    _DN_RAO,
                                                    sim_data->raos
                                                );
            }

            // Storage sources
            if ( input->out_sources )
            {
                output->save_fields_data( 
                                            i,
                                            _DN_SRC_INT,
                                            sim_data->intensities
                                        );
            }

            // Storage panels potential
            if ( input->out_potential )
            {
                output->save_fields_data( 
                                            i,
                                            _DN_POT_INT,
                                            sim_data->panels_potential
                                        );
            }

            // Storage panels pressure
            if ( input->out_pressure )
            {
                output->save_fields_data( 
                                            i,
                                            _DN_PRESS_INT,
                                            sim_data->panels_pressure_p0
                                        );
            }

            if ( input->is_calc_mdrift )
            {
                output->save_wave_exciting_format(
                                                    i,
                                                    _DN_MDRIFT,
                                                    sim_data->mdrift
                                                );

                if ( input->out_qtf_comp )
                {
                    output->save_wave_exciting_format(
                                                        i,
                                                        _DN_MDRIFT_WL,
                                                        sim_data->mdrift_wl
                                                    );
                    
                    output->save_wave_exciting_format(
                                                        i,
                                                        _DN_MDRIFT_BERN,
                                                        sim_data->mdrift_bern
                                                    );

                    output->save_wave_exciting_format(
                                                        i,
                                                        _DN_MDRIFT_ACC,
                                                        sim_data->mdrift_acc
                                                    );
                    
                    output->save_wave_exciting_format(
                                                        i,
                                                        _DN_MDRIFT_MOM,
                                                        sim_data->mdrift_mom
                                                    );
                }
            }

        }
        // std::cout << "Proc. Rank: " << mpi_config->proc_rank << " - Save_to_disk -> Done!" << std::endl;
        double save_tend = MPI_Wtime( );

        MPI_Barrier( MPI_COMM_WORLD );
        double freq_tend = MPI_Wtime( );

        if ( mpi_config->is_root( ) )
        {
            std::cout << "Execution time [s]: " << ( freq_tend - freq_tstart ) << std::endl;
            // std::cout << "Execution Mat time [s]: " << ( mat_tend - mat_tstart ) << std::endl;
            std::cout << "Execution Save time [s]: " << ( save_tend - save_tstart ) << std::endl;
        }
    }

    /****************************************************************/
    /*************** Calculate second order forces ******************/
    /****************************************************************/
    // if ( input->out_qtf )
    // {
    //     if ( mpi_config->is_root( ) )
    //     {
    //         // Loop over frequencies to get the second order force matrix matrix
    //         for ( int i=0; i<input->angfreqs_np; i++ )
    //         {
    //             for ( int j=0; j<input->angfreqs_np; j++ )
    //             {
    //                 if ( i != j )
    //                 {
    //                     std::cout << "i: " << i << " - j: " << j << " - ang_freq_i: " << input->angfreqs[i];
    //                     std::cout << " - ang_freq_j: " << input->angfreqs[j] << " - Proc: " << mpi_config->proc_rank << "\n";
    //                     if ( input->out_qtf_so_model == 0 )
    //                     {
    //                         // Calculate second order force using Pinkster model
    //                         calculate_pinkster(
    //                                                 input,
    //                                                 mpi_config,
    //                                                 mesh_gp,
    //                                                 input->angfreqs[i],
    //                                                 input->angfreqs[j],
    //                                                 sim_data->qtf_diff_secord_force
    //                                             );
    //                     }
    //                     else if ( input->out_qtf_so_model == 1 )
    //                     {
    //                         calculate_secord_force_indirect(
    //                                                             input,
    //                                                             mesh_gp,
    //                                                             input->angfreqs[i],
    //                                                             input->angfreqs[j],
    //                                                             QTF_DIFF_CODE,
    //                                                             qtf_body_pot_gp,
    //                                                             qtf_fs_pot_gp,
    //                                                             qtf_wl_vel_x_gp,
    //                                                             sim_data
    //                                                         );
    //                     }

    //                     // Distribute Pinkster force along the second order force
    //                     // matrix format
    //                     qtf_distribute_matrix_data(
    //                                                     input,
    //                                                     i,
    //                                                     j,
    //                                                     sim_data->qtf_diff_secord_force,
    //                                                     sim_data->qtf_diff_freqs,
    //                                                     1,
    //                                                     0
    //                                                 );
    //                     if ( input->out_qtf_comp )
    //                     {
    //                         qtf_distribute_matrix_data(
    //                                                     input,
    //                                                     i,
    //                                                     j,
    //                                                     sim_data->qtf_diff_secord_force,
    //                                                     sim_data->qtf_diff_secord_force_freqs,
    //                                                     1,
    //                                                     0
    //                                                 );
    //                     }
    //                 }

    //                 // Calculate QTF diff force terms
    //                 calculate_qtf_terms_force(
    //                                                 input,
    //                                                 mesh_gp,
    //                                                 QTF_DIFF_CODE,
    //                                                 &(sim_data->qtf_wl_rel_we_total_freq[i*sim_data->qtf_wl_heads_np]),
    //                                                 &(sim_data->qtf_wl_rel_we_total_freq[j*sim_data->qtf_wl_heads_np]),
    //                                                 &(sim_data->qtf_raos_freq[i*sim_data->wave_exc_np]),
    //                                                 &(sim_data->qtf_raos_freq[j*sim_data->wave_exc_np]),
    //                                                 &(sim_data->qtf_body_vel_x_total_freq[i*sim_data->qtf_body_heads_np]),
    //                                                 &(sim_data->qtf_body_vel_y_total_freq[i*sim_data->qtf_body_heads_np]),
    //                                                 &(sim_data->qtf_body_vel_z_total_freq[i*sim_data->qtf_body_heads_np]),
    //                                                 &(sim_data->qtf_body_vel_x_total_freq[j*sim_data->qtf_body_heads_np]),
    //                                                 &(sim_data->qtf_body_vel_y_total_freq[j*sim_data->qtf_body_heads_np]),
    //                                                 &(sim_data->qtf_body_vel_z_total_freq[j*sim_data->qtf_body_heads_np]),
    //                                                 input->angfreqs[i],
    //                                                 input->angfreqs[j],
    //                                                 sim_data->qtf,
    //                                                 sim_data->qtf_diff_wl,
    //                                                 sim_data->qtf_diff_bern,
    //                                                 sim_data->qtf_diff_acc,
    //                                                 sim_data->qtf_diff_mom,
    //                                                 qtf_wl_we_gp,
    //                                                 vel_x_body_gp,
    //                                                 true
    //                                             );

    //                 qtf_distribute_matrix_data(
    //                                                 input,
    //                                                 i,
    //                                                 j,
    //                                                 sim_data->qtf,
    //                                                 sim_data->qtf_diff_freqs,
    //                                                 0,
    //                                                 1
    //                                             );

    //                 if ( input->out_qtf_comp )
    //                 {
    //                     qtf_distribute_matrix_data(
    //                                                 input,
    //                                                 i,
    //                                                 j,
    //                                                 sim_data->qtf_diff_wl,
    //                                                 sim_data->qtf_diff_wl_freqs,
    //                                                 0,
    //                                                 0
    //                                             );
                        
    //                     qtf_distribute_matrix_data(
    //                                                 input,
    //                                                 i,
    //                                                 j,
    //                                                 sim_data->qtf_diff_bern,
    //                                                 sim_data->qtf_diff_bern_freqs,
    //                                                 0,
    //                                                 0
    //                                             );

    //                     qtf_distribute_matrix_data(
    //                                                 input,
    //                                                 i,
    //                                                 j,
    //                                                 sim_data->qtf_diff_acc,
    //                                                 sim_data->qtf_diff_acc_freqs,
    //                                                 0,
    //                                                 0
    //                                             );

    //                     qtf_distribute_matrix_data(
    //                                                 input,
    //                                                 i,
    //                                                 j,
    //                                                 sim_data->qtf_diff_mom,
    //                                                 sim_data->qtf_diff_mom_freqs,
    //                                                 0,
    //                                                 0
    //                                             );
    //                 }

    //                 // Calculate QTF sum force terms
    //                 calculate_qtf_terms_force(
    //                                                 input,
    //                                                 mesh_gp,
    //                                                 QTF_SUM_CODE,
    //                                                 &(sim_data->qtf_wl_rel_we_total_freq[i*sim_data->qtf_wl_heads_np]),
    //                                                 &(sim_data->qtf_wl_rel_we_total_freq[j*sim_data->qtf_wl_heads_np]),
    //                                                 &(sim_data->qtf_raos_freq[i*sim_data->wave_exc_np]),
    //                                                 &(sim_data->qtf_raos_freq[j*sim_data->wave_exc_np]),
    //                                                 &(sim_data->qtf_body_vel_x_total_freq[i*sim_data->qtf_body_heads_np]),
    //                                                 &(sim_data->qtf_body_vel_y_total_freq[i*sim_data->qtf_body_heads_np]),
    //                                                 &(sim_data->qtf_body_vel_z_total_freq[i*sim_data->qtf_body_heads_np]),
    //                                                 &(sim_data->qtf_body_vel_x_total_freq[j*sim_data->qtf_body_heads_np]),
    //                                                 &(sim_data->qtf_body_vel_y_total_freq[j*sim_data->qtf_body_heads_np]),
    //                                                 &(sim_data->qtf_body_vel_z_total_freq[j*sim_data->qtf_body_heads_np]),
    //                                                 input->angfreqs[i],
    //                                                 input->angfreqs[j],
    //                                                 sim_data->qtf,
    //                                                 sim_data->qtf_sum_wl,
    //                                                 sim_data->qtf_sum_bern,
    //                                                 sim_data->qtf_sum_acc,
    //                                                 sim_data->qtf_sum_mom,
    //                                                 qtf_wl_we_gp,
    //                                                 vel_x_body_gp,
    //                                                 true
    //                                             );

    //                 qtf_distribute_matrix_data(
    //                                                 input,
    //                                                 i,
    //                                                 j,
    //                                                 sim_data->qtf,
    //                                                 sim_data->qtf_sum_freqs,
    //                                                 0,
    //                                                 1
    //                                             );

    //                 if ( input->out_qtf_comp )
    //                 {
    //                     qtf_distribute_matrix_data(
    //                                                 input,
    //                                                 i,
    //                                                 j,
    //                                                 sim_data->qtf_sum_wl,
    //                                                 sim_data->qtf_sum_wl_freqs,
    //                                                 0,
    //                                                 0
    //                                             );
                        
    //                     qtf_distribute_matrix_data(
    //                                                 input,
    //                                                 i,
    //                                                 j,
    //                                                 sim_data->qtf_sum_bern,
    //                                                 sim_data->qtf_sum_bern_freqs,
    //                                                 0,
    //                                                 0
    //                                             );

    //                     qtf_distribute_matrix_data(
    //                                                 input,
    //                                                 i,
    //                                                 j,
    //                                                 sim_data->qtf_sum_acc,
    //                                                 sim_data->qtf_sum_acc_freqs,
    //                                                 0,
    //                                                 0
    //                                             );

    //                     qtf_distribute_matrix_data(
    //                                                 input,
    //                                                 i,
    //                                                 j,
    //                                                 sim_data->qtf_sum_mom,
    //                                                 sim_data->qtf_sum_mom_freqs,
    //                                                 0,
    //                                                 0
    //                                             );
    //                 }

    //             }
    //         }

    //         // Save data into the disk file
    //         output->save_qtf_format( 
    //                                     "qtf_diff",
    //                                     sim_data->qtf_diff_freqs
    //                                 );

    //         output->save_qtf_format(
    //                                     "qtf_sum",
    //                                     sim_data->qtf_sum_freqs
    //                                 );

    //         if ( input->out_qtf_comp )
    //         {
    //             // Storage QTF acceleration term
    //             output->save_qtf_format(
    //                                         "qtf_diff_acc",
    //                                         sim_data->qtf_diff_acc_freqs
    //                                     );
                
    //             output->save_qtf_format(
    //                                         "qtf_sum_acc",
    //                                         sim_data->qtf_sum_acc_freqs
    //                                     );

    //             // Storage QTF bernoulli term
    //             output->save_qtf_format(
    //                                         "qtf_diff_bern",
    //                                         sim_data->qtf_diff_bern_freqs
    //                                     );
                
    //             output->save_qtf_format(
    //                                         "qtf_sum_bern",
    //                                         sim_data->qtf_sum_bern_freqs
    //                                     );

    //             // Storage QTF momentum term
    //             output->save_qtf_format(
    //                                         "qtf_diff_mom",
    //                                         sim_data->qtf_diff_mom_freqs
    //                                     );
                
    //             output->save_qtf_format(
    //                                         "qtf_sum_mom",
    //                                         sim_data->qtf_sum_mom_freqs
    //                                     );

    //             // Storage QTF second order potential term
    //             output->save_qtf_format(
    //                                         "qtf_diff_sop",
    //                                         sim_data->qtf_diff_secord_force_freqs
    //                                     );
                
    //             output->save_qtf_format(
    //                                         "qtf_sum_sop",
    //                                         sim_data->qtf_sum_secord_force_freqs
    //                                     );

    //             // Storage QTF wl term
    //             output->save_qtf_format(
    //                                         "qtf_diff_wl",
    //                                         sim_data->qtf_diff_wl_freqs
    //                                     );
                
    //             output->save_qtf_format(
    //                                         "qtf_sum_wl",
    //                                         sim_data->qtf_sum_wl_freqs
    //                                     );
                
    //         }
    //     }
    // }

    /****************************************************************/
    /************** Delete heap memory allocated data ***************/
    /****************************************************************/

    // Delete integration interfaces
    MPI_Barrier( MPI_COMM_WORLD );
    // delete grf_dn_interf;
    // delete gwf_dn_interf;

    // Delete simulation data
    delete sim_data;
    
    // Delete sysmtem matrixes
    // mkl_free( sysmat );
    // mkl_free( sysmat_steady );
    
    delete source_form_gp;
    delete pot_form_gp;
    delete potpanel_lin_gp;

    if ( input->is_calc_mdrift )
    {
        // delete  gwf_interf;
        // delete  gwf_dx_interf;
        // delete  gwf_dy_interf;
        // delete  gwf_dz_interf;
        delete  qtf_wl_we_gp;
        delete  vel_x_body_gp;
        delete  vel_y_body_gp;
        delete  vel_z_body_gp;
    }

    if ( input->out_qtf_so_model == 1 )
    {
        delete qtf_body_pot_gp;
        delete qtf_fs_pot_gp;
        delete qtf_fs_vel_x_gp;
        delete qtf_fs_vel_y_gp;
        delete qtf_fs_vel_z_gp;
        delete qtf_wl_vel_x_gp;
        delete qtf_wl_vel_y_gp;
        delete qtf_wl_vel_z_gp;
    }

    if ( input->out_qtf_so_model == 1 )
    {
        delete qtf_fs_pot_gp;
        delete qtf_fs_vel_x_gp;
        delete qtf_fs_vel_y_gp;
        delete qtf_fs_vel_z_gp;
        delete qtf_pc_pot_gp;
        delete qtf_pc_vel_x_gp;
        delete qtf_pc_vel_y_gp;
        delete qtf_pc_vel_z_gp;
        delete qtf_wl_vel_x_gp;
        delete qtf_wl_vel_y_gp;
        delete qtf_wl_vel_z_gp;
    }
}


// void    freq_domain_nonlinear_solver(
//                                                 Input*          input,
//                                                 MpiConfig*      mpi_config,
//                                                 MeshGroup*      mesh_gp,
//                                                 SclCmpx*        scl,
//                                                 Hydrostatics**  hydrostatics,
//                                                 Output*         output
//                                     )
// {
//     /****************************************************/
//     /****** Allocate space for the simulation data ******/
//     /****************************************************/

//     // Allocate space for the intensities, hydromechanics and wave exciting forces
//     SimulationData* sim_data    = new SimulationData(
//                                                         input,
//                                                         mpi_config,
//                                                         mesh_gp->meshes_np,
//                                                         input->dofs_np,
//                                                         input->heads_np,
//                                                         scl->num_rows_local,
//                                                         scl->num_rows
//                                                     );

//     // Allocate space for the intensities influence system matrix
//     cuscomplex* sysmat          = generate_empty_vector<cuscomplex>( scl->num_rows_local * scl->num_cols_local );
//     cuscomplex* sysmat_steady   = generate_empty_vector<cuscomplex>( scl->num_rows_local * scl->num_cols_local );

//     /****************************************************/
//     /********* Create Green function interface *********/
//     /****************************************************/
//     GRFDnInterface* grf_dn_interf   = new   GRFDnInterface(
//                                                                 mesh_gp->source_nodes[0],
//                                                                 mesh_gp->source_nodes[0],
//                                                                 input->water_depth
//                                                             );

//     GWFDnInterface* gwf_dn_interf   = new   GWFDnInterface(
//                                                                 mesh_gp->source_nodes[0],
//                                                                 mesh_gp->source_nodes[0],
//                                                                 input->angfreqs[0],
//                                                                 input->water_depth,
//                                                                 input->grav_acc
//                                                             );
    
//     HMFInterface*   hmf_interf      = new   HMFInterface(
//                                                                 mesh_gp->source_nodes,
//                                                                 sim_data->intensities,
//                                                                 mesh_gp->panels[0],
//                                                                 0,
//                                                                 0,
//                                                                 0,
//                                                                 0,
//                                                                 input->angfreqs[0],
//                                                                 input
//                                                         );

//     /****************************************************/
//     /***** Calculate global structural mass matrix ******/
//     /****************************************************/
//     if ( mpi_config->is_root( ) )
//     {
//         calculate_global_structural_mass(
//                                             input,
//                                             sim_data->structural_mass_p0
//                                         );
//     }


//     /****************************************************/
//     /******* Calculate global hydrostatic matrix ********/
//     /****************************************************/
//     if ( mpi_config->is_root( ) )
//     {
//         calculate_global_hydstiffness(
//                                             input,
//                                             hydrostatics,
//                                             sim_data->hydrostiff_p0
//                                         );
//     }

//     /****************************************************/
//     /******* Calculate hydrodynamic coefficients ********/
//     /****************************************************/

//     // Calculate intensities system matrix steady contribution
//     calculate_gf_intensity_steady_sysmat_nlin(
//                                                     input,
//                                                     scl,
//                                                     mesh_gp,
//                                                     grf_dn_interf,
//                                                     sysmat
//                                                 );

//     // Loop over frequencies
//     for ( int i=0; i<input->angfreqs_np; i++ )
//     {
//         MPI_Barrier( MPI_COMM_WORLD );
//         double freq_tstart = MPI_Wtime( );

//         std::cout << "Calculating angular frequency: " << input->angfreqs[i] << std::endl;

//         // Recalculate wave properties for the current 
//         // angular frequency
//         hmf_interf->set_ang_freq( input->angfreqs[i] );
//         gwf_dn_interf->set_ang_freq( input->angfreqs[i] );

//         // Calculate sources intensity
//         MPI_Barrier( MPI_COMM_WORLD );
//         calculate_gf_intensity_sysmat(
//                                             input,
//                                             scl,
//                                             mesh_gp,
//                                             gwf_dn_interf,
//                                             input->angfreqs[i],
//                                             sysmat_steady,
//                                             sysmat,
//                                             sim_data->intensities
//                                         );
        
//         // Gather source values from each processor
//         MPI_Bcast(
//                     sim_data->intensities,
//                     scl->num_rows * ( input->dofs_np + input->heads_np ),
//                     mpi_cuscomplex,
//                     mpi_config->proc_root,
//                     MPI_COMM_WORLD                
//                 );

//         // Update sources values for the integration objects
//         hmf_interf->set_source_values( sim_data->intensities );

//         // Calculate panel potentials
//         calculate_potpanel_raddif_nlin(
//                                                 input,
//                                                 mpi_config,
//                                                 mesh_gp,
//                                                 sim_data->intensities,
//                                                 input->angfreqs[i]
//                                         );

//         // Calculate added mass and damping coefficients
//         calculate_hydromechanic_coeffs_nlin( 
//                                                 input,
//                                                 mpi_config,
//                                                 mesh_gp,
//                                                 hmf_interf,
//                                                 input->angfreqs[i],
//                                                 sim_data->added_mass,
//                                                 sim_data->damping_rad
//                                             );
        
//         // Calculate diffraction forces
//         calculate_diffraction_forces_nlin(
//                                                 input,
//                                                 mpi_config,
//                                                 mesh_gp,
//                                                 hmf_interf,
//                                                 input->angfreqs[i],
//                                                 sim_data->wave_diffrac
//                                         );

//         // Calculate Froude-Krylov forces
//         calculate_froude_krylov_fo(
//                                     input,
//                                     mpi_config,
//                                     mesh_gp,
//                                     input->angfreqs[i],
//                                     sim_data->froude_krylov
//                                 );

//         // Join data from all processors
//         MPI_Reduce(
//                             sim_data->added_mass,
//                             sim_data->added_mass_p0,
//                             sim_data->hydmech_np,
//                             mpi_cusfloat,
//                             MPI_SUM,
//                             mpi_config->proc_root,
//                             MPI_COMM_WORLD
//                         );

//         MPI_Reduce(
//                             sim_data->damping_rad,
//                             sim_data->damping_rad_p0,
//                             sim_data->hydmech_np,
//                             mpi_cusfloat,
//                             MPI_SUM,
//                             mpi_config->proc_root,
//                             MPI_COMM_WORLD
//                     );

//         MPI_Reduce(
//                             sim_data->froude_krylov,
//                             sim_data->froude_krylov_p0,
//                             sim_data->wave_exc_np,
//                             mpi_cuscomplex,
//                             MPI_SUM,
//                             mpi_config->proc_root,
//                             MPI_COMM_WORLD
//                         );

//         MPI_Reduce(
//                             sim_data->wave_diffrac,
//                             sim_data->wave_diffrac_p0,
//                             sim_data->wave_exc_np,
//                             mpi_cuscomplex,
//                             MPI_SUM,
//                             mpi_config->proc_root,
//                             MPI_COMM_WORLD
//                         );

//         // Calculate total wave exciting forces
//         if ( mpi_config->is_root( ) )
//         {
//             sv_add(
//                         sim_data->wave_exc_np,
//                         sim_data->wave_diffrac_p0,
//                         sim_data->froude_krylov_p0,
//                         sim_data->wave_exc_p0
//                     );
//         }

//         // Calculate raos
//         if ( mpi_config->is_root( ) )
//         {
//             calculate_raos(
//                                 input,
//                                 sim_data->structural_mass_p0,
//                                 sim_data->added_mass_p0,
//                                 sim_data->damping_rad_p0,
//                                 sim_data->hydrostiff_p0,
//                                 sim_data->wave_diffrac_p0,
//                                 sim_data->froude_krylov_p0,
//                                 input->angfreqs[i],
//                                 sim_data->raos
//                             );
//         }

//         // Output values to disk
//         if ( mpi_config->is_root( ) )
//         {
//             if ( input->out_hydmech )
//             {
//                 output->save_hydromechanics_format(
//                                                         i,
//                                                         _DN_ADDED_MASS,
//                                                         sim_data->added_mass_p0
//                                                     );

//                 output->save_hydromechanics_format(
//                                                         i,
//                                                         _DN_DAMPING_RAD,
//                                                         sim_data->damping_rad_p0
//                                                     );
//             }

//             if ( input->out_diffrac )
//             {
//                 output->save_wave_exciting_format(
//                                                     i,
//                                                     _DN_DIFFRAC,
//                                                     sim_data->wave_diffrac_p0
//                                                 );
//             }

//             if ( input->out_fk )
//             {
//                 output->save_wave_exciting_format(
//                                                     i,
//                                                     _DN_FK,
//                                                    sim_data->froude_krylov_p0
//                                                 );
//             }

//             if ( input->out_wex )
//             {
//                 output->save_wave_exciting_format(
//                                                     i,
//                                                     _DN_WEX,
//                                                     sim_data->wave_exc_p0
//                                                 );
//             }

//             if ( input->out_raos )
//             {
//                 output->save_wave_exciting_format(
//                                                     i,
//                                                     _DN_RAO,
//                                                     sim_data->raos
//                                                 );
//             }

//         }

//         MPI_Barrier( MPI_COMM_WORLD );
//         double freq_tend = MPI_Wtime( );

//         if ( mpi_config->is_root( ) )
//         {
//             std::cout << "Execution time [s]: " << ( freq_tend - freq_tstart ) << std::endl;
//         }
//     }

//     /****************************************************/
//     /******** Delete heap memory allocated data *********/
//     /****************************************************/
//     MPI_Barrier( MPI_COMM_WORLD );

//     // Delete integration interfaces
//     delete grf_dn_interf;
//     delete gwf_dn_interf;
//     delete hmf_interf;

//     // Delete system matrix
//     mkl_free( sysmat );
//     mkl_free( sysmat_steady );
    
//     // Delete mesh group object
//     delete mesh_gp;
// }
