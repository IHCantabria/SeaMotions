
// Include general usage scientific libraries
#include <complex>
#include <fstream>
#include "mkl_cblas.h"

// Include local modules
#include "freq_solver_tools.hpp"
#include "../../containers/matlin_group.hpp"
#include "../../containers/simulation_data.hpp"
#include "../../interfaces/grf_interface.hpp"
#include "../../solvers/freq_domain/diffraction.hpp"
#include "../../solvers/freq_domain/gf_intensities.hpp"
#include "../../solvers/freq_domain/hydromechanics.hpp"
#include "../../solvers/freq_domain/potential.hpp"
#include "../../solvers/freq_domain/raos.hpp"
#include "../../solvers/freq_domain/froude_krylov.hpp"
#include "../../solvers/freq_domain/wave_elevation.hpp"
#include "../../waves.hpp"


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
    MeshGroup*  mesh_gp =   new   MeshGroup( 
                                                all_meshes,
                                                input->bodies_np,
                                                input->is_wl_points
                                            );

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
    /****** Allocate space for the simulation data ******/
    /****************************************************/

    // Allocate space for the sources, hydromechanics and wave exciting forces
    int         hydmech_np          = pow2s( input->dofs_np * mesh_gp->meshes_np );
    int         wave_exc_np         = input->heads_np * mesh_gp->meshes_np * input->dofs_np;

    cusfloat*   added_mass          = generate_empty_vector<cusfloat>( hydmech_np );
    cusfloat*   damping_rad         = generate_empty_vector<cusfloat>( hydmech_np );
    cuscomplex* froude_krylov       = generate_empty_vector<cuscomplex>( wave_exc_np );
    cuscomplex* raos                = generate_empty_vector<cuscomplex>( wave_exc_np );
    cuscomplex* sources             = generate_empty_vector<cuscomplex>( ( input->dofs_np + input->heads_np ) * scl.num_rows_local );
    cuscomplex* sysmat              = generate_empty_vector<cuscomplex>( scl.num_rows_local * scl.num_cols_local );
    cuscomplex* sysmat_steady       = generate_empty_vector<cuscomplex>( scl.num_rows_local * scl.num_cols_local );
    cuscomplex* wave_diffrac        = generate_empty_vector<cuscomplex>( wave_exc_np );

    cusfloat*   added_mass_p0       = generate_empty_vector<cusfloat>( hydmech_np );
    cusfloat*   damping_rad_p0      = generate_empty_vector<cusfloat>( hydmech_np );
    cuscomplex* froude_krylov_p0    = generate_empty_vector<cuscomplex>( wave_exc_np );
    cusfloat*   hydrostiff_p0       = generate_empty_vector<cusfloat>( hydmech_np );
    cusfloat*   structural_mass_p0  = generate_empty_vector<cusfloat>( hydmech_np );
    cuscomplex* wave_diffrac_p0     = generate_empty_vector<cuscomplex>( wave_exc_np );
    cuscomplex* wave_exc_p0         = generate_empty_vector<cuscomplex>( wave_exc_np );

    // Define memory allocations for the constant source case ( Fast Mode )
    int         ipm_cols_np             = 0;
    int         ipm_sc                  = 0;
    int         ipm_ed                  = 0;
    cuscomplex* mdwe                    = nullptr;
    cuscomplex* mdwe_panel_pot          = nullptr;
    cusfloat*   mdwe_pot_fp             = nullptr;
    cusfloat*   mdwe_cog_to_fp          = nullptr;
    int         mdwe_pot_fp_np          = 0;
    cuscomplex* mdwe_potpanel_total     = nullptr;
    cuscomplex* panel_pot               = nullptr;
    cuscomplex* panel_pot_p0            = nullptr;
    cusfloat*   pot_fp                  = nullptr;
    int         pot_fp_np               = 0;
    cuscomplex* pot_smat                = nullptr;
    cuscomplex* pot_steady_smat         = nullptr;
    cuscomplex* pot_steady_mdwe_smat    = nullptr;
    cuscomplex* pot_mdwe_smat           = nullptr;

    MatLinGroup<cuscomplex>*    mdrift_we       = nullptr;
    MatLinGroup<cuscomplex>*    potpanel_lin    = nullptr;
    MatLinGroup<cuscomplex>*    vel_body        = nullptr;

    if ( input->is_fast_solver )
    {
        // Define column ranges in function of the number of processors
        mpi_config->get_1d_bounds( 
                                            mesh_gp->source_nodes_tnp, 
                                            ipm_sc, 
                                            ipm_ed 
                                        );
        ipm_cols_np             = ipm_ed - ipm_sc;

        // Allocate space for the influence potential matrixes
        pot_smat            = generate_empty_vector<cuscomplex>( mesh_gp->panels_tnp * ipm_cols_np );
        panel_pot           = generate_empty_vector<cuscomplex>( ( input->dofs_np + input->heads_np ) * mesh_gp->panels_tnp );
        panel_pot_p0        = generate_empty_vector<cuscomplex>( ( input->dofs_np + input->heads_np ) * mesh_gp->panels_tnp );
        pot_steady_smat     = generate_empty_vector<cuscomplex>( mesh_gp->panels_tnp * ipm_cols_np );

        if ( input->out_mdrift )
        {
            mdwe_potpanel_total     = generate_empty_vector<cuscomplex>( mesh_gp->panels_wl_tnp );
            pot_steady_mdwe_smat    = generate_empty_vector<cuscomplex>( mesh_gp->panels_wl_tnp * ipm_cols_np );
            pot_mdwe_smat           = generate_empty_vector<cuscomplex>( mesh_gp->panels_wl_tnp * ipm_cols_np );

            vel_body = new MatLinGroup<cuscomplex>(
                                                        mesh_gp->panels_tnp,
                                                        mesh_gp->panels_tnp,
                                                        3
                                                    );
        }
        
    }

    /****************************************************/
    /********* Create Green function interface *********/
    /****************************************************/
    GRFDnInterface* grf_dn_interf   = new   GRFDnInterface(
                                                                mesh_gp->source_nodes[0],
                                                                mesh_gp->source_nodes[0],
                                                                input->water_depth
                                                            );

    GWFDnInterface* gwf_dn_interf   = new   GWFDnInterface(
                                                                mesh_gp->source_nodes[0],
                                                                mesh_gp->source_nodes[0],
                                                                input->angfreqs[0],
                                                                input->water_depth,
                                                                input->grav_acc
                                                            );
    
    HMFInterface*   hmf_interf      = new   HMFInterface(
                                                                mesh_gp->source_nodes,
                                                                sources,
                                                                mesh_gp->panels[0],
                                                                0,
                                                                0,
                                                                0,
                                                                0,
                                                                input->angfreqs[0],
                                                                input
                                                        );

    /****************************************************/
    /***** Calculate global structural mass matrix ******/
    /****************************************************/
    if ( mpi_config->is_root( ) )
    {
        calculate_global_structural_mass(
                                            input,
                                            structural_mass_p0
                                        );
    }


    /****************************************************/
    /******* Calculate global hydrostatic matrix ********/
    /****************************************************/
    if ( mpi_config->is_root( ) )
    {
        calculate_global_hydstiffness(
                                            input,
                                            hydrostatics,
                                            hydrostiff_p0
                                        );
    }

    /****************************************************/
    /******* Calculate hydrodynamic coefficients ********/
    /****************************************************/

    // Calculate steady contribution to the sources intensity
    // for the fast solver
    if ( input->is_fast_solver )
    {
        // Calculate steady part of the sources influence matrix
        double source_steady_t0 = MPI_Wtime( );
        calculate_gf_intensity_steady_sysmat(
                                                input,
                                                &scl,
                                                mesh_gp,
                                                grf_dn_interf,
                                                sysmat_steady
                                            );
        double source_steady_t1 = MPI_Wtime( );
        std::cout << "Time integration sources steady: " << source_steady_t1 - source_steady_t0 << std::endl;

        // Define field points to calculate potential influence matrix
        pot_fp_np     = mesh_gp->diffrac_panels_np;
        pot_fp        = generate_empty_vector<cusfloat>( 3 * pot_fp_np );

        int _count_pot_np   = 0;
        for ( int i=0; i<mesh_gp->panels_tnp; i++ )
        {
            if ( mesh_gp->panels[i]->type == DIFFRAC_PANEL_CODE )
            {
                copy_vector( 3, mesh_gp->panels[i]->center, &(pot_fp[3*_count_pot_np]) );
                _count_pot_np++;
            }
        }

        if ( 
                input->out_mdrift
                ||
                input->out_qtf
            )
        {
            // Allocate memory for the variables associated to the mean drift
            // and QTFs for the fast solver
            mdwe                = generate_empty_vector<cuscomplex>( mdwe_pot_fp_np );
            mdwe_pot_fp_np      = mesh_gp->panels_wl_tnp;
            mdwe_pot_fp         = generate_empty_vector<cusfloat>( 3 * mdwe_pot_fp_np );
            mdwe_cog_to_fp      = generate_empty_vector<cusfloat>( 3 * mdwe_pot_fp_np );

            // Get all WL line centers to be more accesible through a vector
            for ( int i=0; i<mesh_gp->panels_wl_tnp; i++ )
            {
                copy_vector( 3, mesh_gp->panels_wl[i]->center_wl, &(mdwe_pot_fp[3*i]) );
            }

            // Get the radius from the WL line center to the body COG
            for ( int i=0; i<mesh_gp->meshes_np; i++ )
            {
                for ( int j=mesh_gp->panels_wl_cnp[i]; j<mesh_gp->panels_wl_cnp[i+1]; j++ )
                {
                    sv_sub(
                                3,
                                &(mdwe_pot_fp[3*j]),
                                input->bodies[i]->cog,
                                &(mdwe_cog_to_fp[3*j])
                            );
                }
            }

        }

        // Calculate steady part of the potential influence matrix matrix
        double pot_steady_t0 = MPI_Wtime( );
        calculate_influence_potmat_steady(
                                                input,
                                                mpi_config,
                                                mesh_gp,
                                                pot_fp,
                                                pot_fp_np,
                                                pot_steady_smat
                                        );
        double pot_steady_t1 = MPI_Wtime( );
        std::cout << "Time integration potential steady: " << pot_steady_t1 - pot_steady_t0 << std::endl;

        // Calculate steady parto of the potential influence matrix to calculate
        // the mean drift
        calculate_influence_potmat_steady(
                                                input,
                                                mpi_config,
                                                mesh_gp,
                                                mdwe_pot_fp,
                                                mdwe_pot_fp_np,
                                                pot_steady_mdwe_smat
                                            );

        // Calculate 

    }

    // Loop over frequencies
    for ( int i=0; i<input->angfreqs_np; i++ )
    {
        MPI_Barrier( MPI_COMM_WORLD );
        double freq_tstart = MPI_Wtime( );

        std::cout << "Calculating angular frequency: " << input->angfreqs[i] << std::endl;
        // Recalculate wave properties for the current 
        // angular frequency
        hmf_interf->set_ang_freq( input->angfreqs[i] );
        gwf_dn_interf->set_ang_freq( input->angfreqs[i] );

        // Calculate sources intensity
        MPI_Barrier( MPI_COMM_WORLD );
        calculate_gf_intensity_sysmat(
                                            input,
                                            &scl,
                                            mesh_gp,
                                            gwf_dn_interf,
                                            input->angfreqs[i],
                                            sysmat_steady,
                                            sysmat,
                                            sources
                                        );
        
        // Gather source values from each processor
        MPI_Bcast(
                    sources,
                    scl.num_rows * ( input->dofs_np + input->heads_np ),
                    mpi_cuscomplex,
                    mpi_config->proc_root,
                    MPI_COMM_WORLD                
                );

        // int pshift = 10;
        // if ( mpi_config->is_root( ) )
        // {
        //     std::cout << "Sources Intensity: " << std::endl;
        //     for ( int i=0; i<5; i++ )
        //     {
        //         std::cout << "Source[" <<  i << "]: " << sources[pshift+i] << " - " << std::abs( sources[pshift+i] ) << " - " << std::arg( sources[pshift+i] ) << std::endl;
        //     }
        // }

        // Update sources values for the integration objects
        hmf_interf->set_source_values( sources );

        // Calculate potential matrix if any
        if ( input->is_fast_solver )
        {
            // Calculate potential influence coeffcients matrix
            calculate_influence_potmat(
                                            input,
                                            mpi_config,
                                            mesh_gp,
                                            input->angfreqs[i],
                                            pot_steady_smat,
                                            pot_fp,
                                            pot_fp_np,
                                            pot_smat
                                        );

            // Calculate panels potential
            calculate_potpanel_raddif_lin(
                                                    input,
                                                    pot_smat,
                                                    mesh_gp->panels_tnp,
                                                    ipm_cols_np,
                                                    ipm_sc,
                                                    sources,
                                                    panel_pot
                                            );

            MPI_Allreduce(
                            panel_pot,
                            panel_pot_p0,
                            mesh_gp->panels_tnp * ( input->dofs_np + input->heads_np ),
                            mpi_cuscomplex,
                            MPI_SUM,
                            MPI_COMM_WORLD
                        );

            // std::cout << "Panel Potential: " << std::endl;
            // for ( int i=0; i<5; i++ )
            // {
            //     std::cout << "Potential[" <<  i << "]: " << panel_pot[pshift+i] << " - " << std::abs( panel_pot[pshift+i] ) << " - " << std::arg( panel_pot[pshift+i] ) << std::endl;
            // }

            // Calculate added mass and damping
            calculate_hydromechanic_coeffs_lin( 
                                                    input,
                                                    mpi_config,
                                                    mesh_gp,
                                                    panel_pot_p0,
                                                    input->angfreqs[i],
                                                    added_mass,
                                                    damping_rad
                                                );

            // Calculate diffraction forces
            calculate_diffraction_forces_lin(
                                                    input,
                                                    mpi_config,
                                                    mesh_gp,
                                                    panel_pot_p0,
                                                    input->angfreqs[i],
                                                    wave_diffrac
                                            );
        }
        else
        {
            // Calculate panel potentials
            calculate_potpanel_raddif_nlin(
                                                    input,
                                                    mpi_config,
                                                    mesh_gp,
                                                    sources,
                                                    input->angfreqs[i]
                                            );

            // Calculate added mass and damping coefficients
            calculate_hydromechanic_coeffs_nlin( 
                                                    input,
                                                    mpi_config,
                                                    mesh_gp,
                                                    hmf_interf,
                                                    input->angfreqs[i],
                                                    added_mass,
                                                    damping_rad
                                                );
            
            // Calculate diffraction forces
            calculate_diffraction_forces_nlin(
                                                    input,
                                                    mpi_config,
                                                    mesh_gp,
                                                    hmf_interf,
                                                    input->angfreqs[i],
                                                    wave_diffrac
                                            );

        }

        // Calculate Froude-Krylov forces
        calculate_froude_krylov(
                                    input,
                                    mpi_config,
                                    mesh_gp,
                                    input->angfreqs[i],
                                    froude_krylov
                                );

        // Join data from all processors
        MPI_Reduce(
                            added_mass,
                            added_mass_p0,
                            hydmech_np,
                            mpi_cusfloat,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                        );

        MPI_Reduce(
                            damping_rad,
                            damping_rad_p0,
                            hydmech_np,
                            mpi_cusfloat,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                    );

        MPI_Reduce(
                            froude_krylov,
                            froude_krylov_p0,
                            wave_exc_np,
                            mpi_cuscomplex,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                        );

        MPI_Reduce(
                            wave_diffrac,
                            wave_diffrac_p0,
                            wave_exc_np,
                            mpi_cuscomplex,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                        );

        // Calculate total wave exciting forces
        if ( mpi_config->is_root( ) )
        {
            sv_add(
                        wave_exc_np,
                        wave_diffrac_p0,
                        froude_krylov_p0,
                        wave_exc_p0
                    );
        }

        // Calculate raos
        if ( mpi_config->is_root( ) )
        {
            calculate_raos(
                                input,
                                structural_mass_p0,
                                added_mass_p0,
                                damping_rad_p0,
                                hydrostiff_p0,
                                wave_diffrac_p0,
                                froude_krylov_p0,
                                input->angfreqs[i],
                                raos
                            );
        }

        if ( input->out_mdrift )
        {
            // Broadcast RAOs values to be available in all the processes
            MPI_Bcast(
                        raos,
                        wave_exc_np,
                        mpi_cuscomplex,
                        mpi_config->proc_root,
                        MPI_COMM_WORLD    
                    );
            
            if ( input->is_fast_solver )
            {

                // Calculate total potential at the target WL points
                calculate_potpanel_total_lin(
                                                input,
                                                mpi_config,
                                                mesh_gp,
                                                input->angfreqs[i],
                                                sources,
                                                pot_steady_mdwe_smat,
                                                pot_mdwe_smat,
                                                mdwe_pot_fp,
                                                mesh_gp->panels_wl_cnp,
                                                mesh_gp->meshes_np,
                                                raos,
                                                mdwe_potpanel_total
                                            );
                
                // Calculate relative wave elevation
                calculate_relative_wave_elevation_lin(
                                                            input,
                                                            mpi_config,
                                                            mdwe_cog_to_fp,
                                                            mesh_gp->panels_wl_cnp,
                                                            mesh_gp->meshes_np,
                                                            mdwe_potpanel_total,
                                                            input->angfreqs[i],
                                                            raos,
                                                            mdwe
                                                        );

                
            }

        }

        // Output values to disk
        if ( mpi_config->is_root( ) )
        {
            if ( input->out_hydmech )
            {
                output->save_hydromechanics_format(
                                                        i,
                                                        _DN_ADDED_MASS,
                                                        added_mass_p0
                                                    );

                output->save_hydromechanics_format(
                                                        i,
                                                        _DN_DAMPING_RAD,
                                                        damping_rad_p0
                                                    );
            }

            if ( input->out_diffrac )
            {
                output->save_wave_exciting_format(
                                                    i,
                                                    _DN_DIFFRAC,
                                                    wave_diffrac_p0
                                                );
            }

            if ( input->out_fk )
            {
                output->save_wave_exciting_format(
                                                    i,
                                                    _DN_FK,
                                                    froude_krylov_p0
                                                );
            }

            if ( input->out_wex )
            {
                output->save_wave_exciting_format(
                                                    i,
                                                    _DN_WEX,
                                                    wave_exc_p0
                                                );
            }

            if ( input->out_raos )
            {
                output->save_wave_exciting_format(
                                                    i,
                                                    _DN_RAO,
                                                    raos
                                                );
            }

        }

        MPI_Barrier( MPI_COMM_WORLD );
        double freq_tend = MPI_Wtime( );

        if ( mpi_config->is_root( ) )
        {
            std::cout << "Execution time [s]: " << ( freq_tend - freq_tstart ) << std::endl;
        }
    }

    // Delete heap memory allocated data
    MPI_Barrier( MPI_COMM_WORLD );
    delete grf_dn_interf;
    delete gwf_dn_interf;
    delete hmf_interf;

    mkl_free( added_mass_p0 );
    mkl_free( damping_rad_p0 );
    mkl_free( froude_krylov_p0 );
    mkl_free( hydrostiff_p0 );
    mkl_free( structural_mass_p0 );
    mkl_free( wave_diffrac_p0 );
    mkl_free( wave_exc_p0 );
    
    mkl_free( added_mass );
    mkl_free( damping_rad );
    mkl_free( froude_krylov );
    mkl_free( raos );
    mkl_free( sources );
    mkl_free( sysmat );
    mkl_free( sysmat_steady );
    mkl_free( wave_diffrac );

    if ( input->is_fast_solver )
    {
        mkl_free( pot_smat );
        mkl_free( pot_steady_smat );
        mkl_free( panel_pot_p0 );
        mkl_free( panel_pot );
        mkl_free( pot_fp );

        if ( input->out_mdrift )
        {
            mkl_free( mdwe );
            mkl_free( mdwe_cog_to_fp );
            mkl_free( mdwe_pot_fp );
            mkl_free( mdwe_potpanel_total );
            mkl_free( pot_mdwe_smat );
            mkl_free( pot_steady_mdwe_smat );

            delete vel_body_gp;
        }
    }
    
    delete mesh_gp;
    delete [] all_meshes;
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


void    linear_solver(
                            Input*          input,
                            MpiConfig*      mpi_config,
                            MeshGroup*      mesh_gp,
                            SclCmpx*        scl,
                            Hydrostatics**  hydrostatics,
                            Output*         output
                    )
{
    /****************************************************************/
    /************ Allocate space for the simulation data ************/
    /****************************************************************/

    // Allocate space for the intensities, hydromechanics and wave exciting forces
    SimulationData* sim_data    = new SimulationData(
                                                        mesh_gp->meshes_np,
                                                        input->dofs_np,
                                                        input->heads_np,
                                                        scl->num_rows_local,
                                                        scl->num_cols_local,
                                                        scl->num_rows,
                                                        mpi_config
                                                    );

    // Allocate space for the intensities influence system matrix
    cuscomplex* sysmat          = generate_empty_vector<cuscomplex>( scl->num_rows_local * scl->num_cols_local );
    cuscomplex* sysmat_steady   = generate_empty_vector<cuscomplex>( scl->num_rows_local * scl->num_cols_local );
    
    // Allocate space for the system matrixes of the required potentials
    int         ipm_cols_np     = 0;
    int         ipm_sc          = 0;
    int         ipm_ed          = 0;

    mpi_config->get_1d_bounds( 
                                mesh_gp->panels_tnp, 
                                ipm_sc, 
                                ipm_ed 
                            );
    ipm_cols_np = ipm_ed - ipm_sc;

    MatLinGroup<cuscomplex>*    potpanel_lin_gp = new MatLinGroup<cuscomplex>(
                                                                                mesh_gp->panels_raddif_tnp,
                                                                                ipm_cols_np,
                                                                                mesh_gp->meshes_np,
                                                                                ( input->dofs_np + input->heads_np ),
                                                                                0,
                                                                                mesh_gp->panels_raddif_tnp-1,
                                                                                ipm_sc,
                                                                                ipm_ed
                                                                            );
    MatLinGroup<cuscomplex>*    mdrift_we_gp        = nullptr;
    MatLinGroup<cuscomplex>*    vel_body_gp         = nullptr;
    if ( input->out_mdrift )
    {
        sim_data->add_mean_drift_data( mesh_gp->panels_wl_tnp );

        mdrift_we_gp            = new MatLinGroup<cuscomplex>(
                                                                mesh_gp->panels_wl_tnp,
                                                                ipm_cols_np,
                                                                mesh_gp->meshes_np,
                                                                ( input->dofs_np + input->heads_np ),
                                                                0,
                                                                mesh_gp->panels_wl_tnp,
                                                                ipm_sc,
                                                                ipm_ed
                                                            );

        vel_body_gp             = new MatLinGroup<cuscomplex>(
                                                                mesh_gp->panels_tnp,
                                                                ipm_cols_np,
                                                                mesh_gp->meshes_np,
                                                                ( input->dofs_np + input->heads_np ),
                                                                0,
                                                                mesh_gp->panels_tnp-1,
                                                                ipm_sc,
                                                                ipm_ed
                                                            );
    }

    /****************************************************************/
    /************** Define field points for evaluation **************/
    /****************************************************************/

    // Define field points to calculate potential influence matrix
    int _count_pot_np   = 0;
    for ( int i=0; i<mesh_gp->panels_tnp; i++ )
    {
        if ( mesh_gp->panels[i]->type == DIFFRAC_PANEL_CODE )
        {
            copy_vector( 3, mesh_gp->panels[i]->center, &(potpanel_lin_gp->field_points[3*_count_pot_np]) );
            _count_pot_np++;
        }
    }
    copy_vector( mesh_gp->meshes_np, mesh_gp->panels_raddif_cnp, potpanel_lin_gp->field_points_cnp );

    // Define field points to calculate the potential on the WL to evaluate
    // second order forces
    if ( 
            input->out_mdrift
            ||
            input->out_qtf
        )
    {
        // Get all WL line centers to be more accesible through a vector
        for ( int i=0; i<mesh_gp->panels_wl_tnp; i++ )
        {
            copy_vector( 3, mesh_gp->panels_wl[i]->center_wl, &(mdrift_we_gp->field_points[3*i]) );
        }
        copy_vector( mesh_gp->meshes_np, mesh_gp->panels_wl_cnp, mdrift_we_gp->field_points_cnp );

        // Get the radius from the WL line center to the body COG
        for ( int i=0; i<mesh_gp->meshes_np; i++ )
        {
            for ( int j=mesh_gp->panels_wl_cnp[i]; j<mesh_gp->panels_wl_cnp[i+1]; j++ )
            {
                sv_sub(
                            3,
                            &(mdrift_we_gp->field_points[3*j]),
                            input->bodies[i]->cog,
                            &(mdrift_we_gp->cog_to_field_points[3*j])
                        );
            }
        }

    }

    /****************************************************************/
    /*************** Create Green function interface ****************/
    /****************************************************************/
    GRFDnInterface* grf_dn_interf   = new   GRFDnInterface(
                                                                mesh_gp->source_nodes[0],
                                                                mesh_gp->source_nodes[0],
                                                                input->water_depth
                                                            );

    GWFDnInterface* gwf_dn_interf   = new   GWFDnInterface(
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

    // Calculate steady part of the sources influence matrix
    double source_steady_t0 = MPI_Wtime( );
    calculate_gf_intensity_steady_sysmat(
                                            input,
                                            scl,
                                            mesh_gp,
                                            grf_dn_interf,
                                            sysmat_steady
                                        );
    double source_steady_t1 = MPI_Wtime( );
    std::cout << "Time integration sources steady: " << source_steady_t1 - source_steady_t0 << std::endl;

    

    // Calculate steady part of the potential influence matrix matrix
    double pot_steady_t0 = MPI_Wtime( );
    calculate_influence_potmat_steady(
                                            input,
                                            mpi_config,
                                            mesh_gp,
                                            potpanel_lin_gp
                                    );
    double pot_steady_t1 = MPI_Wtime( );
    std::cout << "Time integration potential steady: " << pot_steady_t1 - pot_steady_t0 << std::endl;

    // Calculate steady parto of the potential influence matrix to calculate
    // the mean drift
    calculate_influence_potmat_steady(
                                            input,
                                            mpi_config,
                                            mesh_gp,
                                            mdrift_we_gp
                                        );

    /****************************************************************/
    /******* Calculate wave contributions to system matrixes ********/
    /********** to obtain the hydrodynamic coefficients *************/
    /****************************************************************/

    // Loop over frequencies
    for ( int i=0; i<input->angfreqs_np; i++ )
    {
        MPI_Barrier( MPI_COMM_WORLD );
        double freq_tstart = MPI_Wtime( );

        std::cout << "Calculating angular frequency: " << input->angfreqs[i] << std::endl;

        // Recalculate wave properties for the current 
        // angular frequency
        gwf_dn_interf->set_ang_freq( input->angfreqs[i] );

        // Calculate sources intensity
        MPI_Barrier( MPI_COMM_WORLD );
        calculate_gf_intensity_sysmat(
                                            input,
                                            scl,
                                            mesh_gp,
                                            gwf_dn_interf,
                                            input->angfreqs[i],
                                            sysmat_steady,
                                            sysmat,
                                            sim_data->intensities
                                        );
        
        // Gather source values from each processor
        MPI_Bcast(
                    sim_data->intensities,
                    scl->num_rows * ( input->dofs_np + input->heads_np ),
                    mpi_cuscomplex,
                    mpi_config->proc_root,
                    MPI_COMM_WORLD                
                );


        // Calculate potential influence coeffcients matrix
        calculate_influence_potmat(
                                        input,
                                        mpi_config,
                                        mesh_gp,
                                        input->angfreqs[i],
                                        potpanel_lin_gp
                                    );

        // Calculate panels potential
        calculate_potpanel_raddif_lin(
                                                input,
                                                sim_data->intensities,
                                                potpanel_lin_gp
                                        );

        MPI_Allreduce(
                        potpanel_lin_gp->field_values,
                        sim_data->panels_potential,
                        mesh_gp->panels_tnp * ( input->dofs_np + input->heads_np ),
                        mpi_cuscomplex,
                        MPI_SUM,
                        MPI_COMM_WORLD
                    );

        // Calculate added mass and damping
        calculate_hydromechanic_coeffs_lin( 
                                                input,
                                                mpi_config,
                                                mesh_gp,
                                                sim_data->panels_potential,
                                                input->angfreqs[i],
                                                sim_data->added_mass,
                                                sim_data->damping_rad
                                            );

        // Calculate diffraction forces
        calculate_diffraction_forces_lin(
                                                input,
                                                mpi_config,
                                                mesh_gp,
                                                sim_data->panels_potential,
                                                input->angfreqs[i],
                                                sim_data->wave_diffrac
                                        );
        
        // Calculate Froude-Krylov forces
        calculate_froude_krylov(
                                    input,
                                    mpi_config,
                                    mesh_gp,
                                    input->angfreqs[i],
                                    sim_data->froude_krylov
                                );

        // Join data from all processors
        MPI_Reduce(
                            sim_data->added_mass,
                            sim_data->added_mass_p0,
                            sim_data->hydmech_np,
                            mpi_cusfloat,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                        );

        MPI_Reduce(
                            sim_data->damping_rad,
                            sim_data->damping_rad_p0,
                            sim_data->hydmech_np,
                            mpi_cusfloat,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                    );

        MPI_Reduce(
                            sim_data->froude_krylov,
                            sim_data->froude_krylov_p0,
                            sim_data->wave_exc_np,
                            mpi_cuscomplex,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                        );

        MPI_Reduce(
                            sim_data->wave_diffrac,
                            sim_data->wave_diffrac_p0,
                            sim_data->wave_exc_np,
                            mpi_cuscomplex,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                        );

        // Calculate total wave exciting forces
        if ( mpi_config->is_root( ) )
        {
            sv_add(
                        sim_data->wave_exc_np,
                        sim_data->wave_diffrac_p0,
                        sim_data->froude_krylov_p0,
                        sim_data->wave_exc_p0
                    );
        }

        // Calculate raos
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

        if ( input->out_mdrift )
        {
            // Broadcast RAOs values to be available in all the processes
            MPI_Bcast(
                        sim_data->raos,
                        sim_data->wave_exc_np,
                        mpi_cuscomplex,
                        mpi_config->proc_root,
                        MPI_COMM_WORLD    
                    );
            
            // Calculate total potential at the target WL points
            calculate_potpanel_total_lin(
                                            input,
                                            mpi_config,
                                            mesh_gp,
                                            input->angfreqs[i],
                                            sim_data->intensities,
                                            sim_data->raos,
                                            mdrift_we_gp,
                                            sim_data->mdrift_we_pot_total
                                        );
            
            // Calculate relative wave elevation
            calculate_relative_wave_elevation_lin(
                                                        input,
                                                        mpi_config,
                                                        potpanel_lin_gp,
                                                        sim_data->mdrift_we_pot_total,
                                                        input->angfreqs[i],
                                                        sim_data->raos,
                                                        sim_data->mdrift_rel_we
                                                    );
        }

        // Output values to disk
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

        }

        MPI_Barrier( MPI_COMM_WORLD );
        double freq_tend = MPI_Wtime( );

        if ( mpi_config->is_root( ) )
        {
            std::cout << "Execution time [s]: " << ( freq_tend - freq_tstart ) << std::endl;
        }
    }

    /****************************************************************/
    /************** Delete heap memory allocated data ***************/
    /****************************************************************/

    // Delete integration interfaces
    MPI_Barrier( MPI_COMM_WORLD );
    delete grf_dn_interf;
    delete gwf_dn_interf;

    // Delete simulation data
    delete sim_data;
    
    // Delete sysmtem matrixes
    mkl_free( sysmat );
    mkl_free( sysmat_steady );
    
    delete potpanel_lin_gp;

    if ( input->out_mdrift )
    {
        delete      mdrift_we_gp;
        delete      vel_body_gp;
    }
    
    // Delete mesh group data
    delete mesh_gp;
}


void    nonlinear_solver(
                            Input*          input,
                            MpiConfig*      mpi_config,
                            MeshGroup*      mesh_gp,
                            Hydrostatics**  hydrostatics,
                            Output*         output
                        )
{
    /****************************************************/
    /****** Allocate space for the simulation data ******/
    /****************************************************/

    // Allocate space for the sources, hydromechanics and wave exciting forces
    int         hydmech_np          = pow2s( input->dofs_np * mesh_gp->meshes_np );
    int         wave_exc_np         = input->heads_np * mesh_gp->meshes_np * input->dofs_np;

    cusfloat*   added_mass          = generate_empty_vector<cusfloat>( hydmech_np );
    cusfloat*   damping_rad         = generate_empty_vector<cusfloat>( hydmech_np );
    cuscomplex* froude_krylov       = generate_empty_vector<cuscomplex>( wave_exc_np );
    cuscomplex* raos                = generate_empty_vector<cuscomplex>( wave_exc_np );
    cuscomplex* sources             = generate_empty_vector<cuscomplex>( ( input->dofs_np + input->heads_np ) * scl.num_rows_local );
    cuscomplex* sysmat              = generate_empty_vector<cuscomplex>( scl.num_rows_local * scl.num_cols_local );
    cuscomplex* sysmat_steady       = generate_empty_vector<cuscomplex>( scl.num_rows_local * scl.num_cols_local );
    cuscomplex* wave_diffrac        = generate_empty_vector<cuscomplex>( wave_exc_np );

    cusfloat*   added_mass_p0       = generate_empty_vector<cusfloat>( hydmech_np );
    cusfloat*   damping_rad_p0      = generate_empty_vector<cusfloat>( hydmech_np );
    cuscomplex* froude_krylov_p0    = generate_empty_vector<cuscomplex>( wave_exc_np );
    cusfloat*   hydrostiff_p0       = generate_empty_vector<cusfloat>( hydmech_np );
    cusfloat*   structural_mass_p0  = generate_empty_vector<cusfloat>( hydmech_np );
    cuscomplex* wave_diffrac_p0     = generate_empty_vector<cuscomplex>( wave_exc_np );
    cuscomplex* wave_exc_p0         = generate_empty_vector<cuscomplex>( wave_exc_np );

    // Define memory allocations for the constant source case ( Fast Mode )
    int         ipm_cols_np             = 0;
    int         ipm_sc                  = 0;
    int         ipm_ed                  = 0;
    cuscomplex* mdwe                    = nullptr;
    cuscomplex* mdwe_panel_pot          = nullptr;
    cusfloat*   mdwe_pot_fp             = nullptr;
    cusfloat*   mdwe_cog_to_fp          = nullptr;
    int         mdwe_pot_fp_np          = 0;
    cuscomplex* mdwe_potpanel_total     = nullptr;
    cuscomplex* panel_pot               = nullptr;
    cuscomplex* panel_pot_p0            = nullptr;
    cusfloat*   pot_fp                  = nullptr;
    int         pot_fp_np               = 0;
    cuscomplex* pot_smat                = nullptr;
    cuscomplex* pot_steady_smat         = nullptr;
    cuscomplex* pot_steady_mdwe_smat    = nullptr;
    cuscomplex* pot_mdwe_smat           = nullptr;

    MatLinGroup<cuscomplex>*    mdrift_we       = nullptr;
    MatLinGroup<cuscomplex>*    potpanel_lin    = nullptr;
    MatLinGroup<cuscomplex>*    vel_body        = nullptr;

    if ( input->is_fast_solver )
    {
        // Define column ranges in function of the number of processors
        mpi_config->get_1d_bounds( 
                                            mesh_gp->source_nodes_tnp, 
                                            ipm_sc, 
                                            ipm_ed 
                                        );
        ipm_cols_np             = ipm_ed - ipm_sc;

        // Allocate space for the influence potential matrixes
        pot_smat            = generate_empty_vector<cuscomplex>( mesh_gp->panels_tnp * ipm_cols_np );
        panel_pot           = generate_empty_vector<cuscomplex>( ( input->dofs_np + input->heads_np ) * mesh_gp->panels_tnp );
        panel_pot_p0        = generate_empty_vector<cuscomplex>( ( input->dofs_np + input->heads_np ) * mesh_gp->panels_tnp );
        pot_steady_smat     = generate_empty_vector<cuscomplex>( mesh_gp->panels_tnp * ipm_cols_np );

        if ( input->out_mdrift )
        {
            mdwe_potpanel_total     = generate_empty_vector<cuscomplex>( mesh_gp->panels_wl_tnp );
            pot_steady_mdwe_smat    = generate_empty_vector<cuscomplex>( mesh_gp->panels_wl_tnp * ipm_cols_np );
            pot_mdwe_smat           = generate_empty_vector<cuscomplex>( mesh_gp->panels_wl_tnp * ipm_cols_np );

            vel_body = new MatLinGroup<cuscomplex>(
                                                        mesh_gp->panels_tnp,
                                                        mesh_gp->panels_tnp,
                                                        3
                                                    );
        }
        
    }

    /****************************************************/
    /********* Create Green function interface *********/
    /****************************************************/
    GRFDnInterface* grf_dn_interf   = new   GRFDnInterface(
                                                                mesh_gp->source_nodes[0],
                                                                mesh_gp->source_nodes[0],
                                                                input->water_depth
                                                            );

    GWFDnInterface* gwf_dn_interf   = new   GWFDnInterface(
                                                                mesh_gp->source_nodes[0],
                                                                mesh_gp->source_nodes[0],
                                                                input->angfreqs[0],
                                                                input->water_depth,
                                                                input->grav_acc
                                                            );
    
    HMFInterface*   hmf_interf      = new   HMFInterface(
                                                                mesh_gp->source_nodes,
                                                                sources,
                                                                mesh_gp->panels[0],
                                                                0,
                                                                0,
                                                                0,
                                                                0,
                                                                input->angfreqs[0],
                                                                input
                                                        );

    /****************************************************/
    /***** Calculate global structural mass matrix ******/
    /****************************************************/
    if ( mpi_config->is_root( ) )
    {
        calculate_global_structural_mass(
                                            input,
                                            structural_mass_p0
                                        );
    }


    /****************************************************/
    /******* Calculate global hydrostatic matrix ********/
    /****************************************************/
    if ( mpi_config->is_root( ) )
    {
        calculate_global_hydstiffness(
                                            input,
                                            hydrostatics,
                                            hydrostiff_p0
                                        );
    }

    /****************************************************/
    /******* Calculate hydrodynamic coefficients ********/
    /****************************************************/

    // Calculate steady contribution to the sources intensity
    // for the fast solver
    if ( input->is_fast_solver )
    {
        // Calculate steady part of the sources influence matrix
        double source_steady_t0 = MPI_Wtime( );
        calculate_gf_intensity_steady_sysmat(
                                                input,
                                                &scl,
                                                mesh_gp,
                                                grf_dn_interf,
                                                sysmat_steady
                                            );
        double source_steady_t1 = MPI_Wtime( );
        std::cout << "Time integration sources steady: " << source_steady_t1 - source_steady_t0 << std::endl;

        // Define field points to calculate potential influence matrix
        pot_fp_np     = mesh_gp->diffrac_panels_np;
        pot_fp        = generate_empty_vector<cusfloat>( 3 * pot_fp_np );

        int _count_pot_np   = 0;
        for ( int i=0; i<mesh_gp->panels_tnp; i++ )
        {
            if ( mesh_gp->panels[i]->type == DIFFRAC_PANEL_CODE )
            {
                copy_vector( 3, mesh_gp->panels[i]->center, &(pot_fp[3*_count_pot_np]) );
                _count_pot_np++;
            }
        }

        if ( 
                input->out_mdrift
                ||
                input->out_qtf
            )
        {
            // Allocate memory for the variables associated to the mean drift
            // and QTFs for the fast solver
            mdwe                = generate_empty_vector<cuscomplex>( mdwe_pot_fp_np );
            mdwe_pot_fp_np      = mesh_gp->panels_wl_tnp;
            mdwe_pot_fp         = generate_empty_vector<cusfloat>( 3 * mdwe_pot_fp_np );
            mdwe_cog_to_fp      = generate_empty_vector<cusfloat>( 3 * mdwe_pot_fp_np );

            // Get all WL line centers to be more accesible through a vector
            for ( int i=0; i<mesh_gp->panels_wl_tnp; i++ )
            {
                copy_vector( 3, mesh_gp->panels_wl[i]->center_wl, &(mdwe_pot_fp[3*i]) );
            }

            // Get the radius from the WL line center to the body COG
            for ( int i=0; i<mesh_gp->meshes_np; i++ )
            {
                for ( int j=mesh_gp->panels_wl_cnp[i]; j<mesh_gp->panels_wl_cnp[i+1]; j++ )
                {
                    sv_sub(
                                3,
                                &(mdwe_pot_fp[3*j]),
                                input->bodies[i]->cog,
                                &(mdwe_cog_to_fp[3*j])
                            );
                }
            }

        }

        // Calculate steady part of the potential influence matrix matrix
        double pot_steady_t0 = MPI_Wtime( );
        calculate_influence_potmat_steady(
                                                input,
                                                mpi_config,
                                                mesh_gp,
                                                pot_fp,
                                                pot_fp_np,
                                                pot_steady_smat
                                        );
        double pot_steady_t1 = MPI_Wtime( );
        std::cout << "Time integration potential steady: " << pot_steady_t1 - pot_steady_t0 << std::endl;

        // Calculate steady parto of the potential influence matrix to calculate
        // the mean drift
        calculate_influence_potmat_steady(
                                                input,
                                                mpi_config,
                                                mesh_gp,
                                                mdwe_pot_fp,
                                                mdwe_pot_fp_np,
                                                pot_steady_mdwe_smat
                                            );

        // Calculate 

    }

    // Loop over frequencies
    for ( int i=0; i<input->angfreqs_np; i++ )
    {
        MPI_Barrier( MPI_COMM_WORLD );
        double freq_tstart = MPI_Wtime( );

        std::cout << "Calculating angular frequency: " << input->angfreqs[i] << std::endl;
        // Recalculate wave properties for the current 
        // angular frequency
        hmf_interf->set_ang_freq( input->angfreqs[i] );
        gwf_dn_interf->set_ang_freq( input->angfreqs[i] );

        // Calculate sources intensity
        MPI_Barrier( MPI_COMM_WORLD );
        calculate_gf_intensity_sysmat(
                                            input,
                                            &scl,
                                            mesh_gp,
                                            gwf_dn_interf,
                                            input->angfreqs[i],
                                            sysmat_steady,
                                            sysmat,
                                            sources
                                        );
        
        // Gather source values from each processor
        MPI_Bcast(
                    sources,
                    scl.num_rows * ( input->dofs_np + input->heads_np ),
                    mpi_cuscomplex,
                    mpi_config->proc_root,
                    MPI_COMM_WORLD                
                );

        // int pshift = 10;
        // if ( mpi_config->is_root( ) )
        // {
        //     std::cout << "Sources Intensity: " << std::endl;
        //     for ( int i=0; i<5; i++ )
        //     {
        //         std::cout << "Source[" <<  i << "]: " << sources[pshift+i] << " - " << std::abs( sources[pshift+i] ) << " - " << std::arg( sources[pshift+i] ) << std::endl;
        //     }
        // }

        // Update sources values for the integration objects
        hmf_interf->set_source_values( sources );

        // Calculate potential matrix if any
        if ( input->is_fast_solver )
        {
            // Calculate potential influence coeffcients matrix
            calculate_influence_potmat(
                                            input,
                                            mpi_config,
                                            mesh_gp,
                                            input->angfreqs[i],
                                            pot_steady_smat,
                                            pot_fp,
                                            pot_fp_np,
                                            pot_smat
                                        );

            // Calculate panels potential
            calculate_potpanel_raddif_lin(
                                                    input,
                                                    pot_smat,
                                                    mesh_gp->panels_tnp,
                                                    ipm_cols_np,
                                                    ipm_sc,
                                                    sources,
                                                    panel_pot
                                            );

            MPI_Allreduce(
                            panel_pot,
                            panel_pot_p0,
                            mesh_gp->panels_tnp * ( input->dofs_np + input->heads_np ),
                            mpi_cuscomplex,
                            MPI_SUM,
                            MPI_COMM_WORLD
                        );

            // std::cout << "Panel Potential: " << std::endl;
            // for ( int i=0; i<5; i++ )
            // {
            //     std::cout << "Potential[" <<  i << "]: " << panel_pot[pshift+i] << " - " << std::abs( panel_pot[pshift+i] ) << " - " << std::arg( panel_pot[pshift+i] ) << std::endl;
            // }

            // Calculate added mass and damping
            calculate_hydromechanic_coeffs_lin( 
                                                    input,
                                                    mpi_config,
                                                    mesh_gp,
                                                    panel_pot_p0,
                                                    input->angfreqs[i],
                                                    added_mass,
                                                    damping_rad
                                                );

            // Calculate diffraction forces
            calculate_diffraction_forces_lin(
                                                    input,
                                                    mpi_config,
                                                    mesh_gp,
                                                    panel_pot_p0,
                                                    input->angfreqs[i],
                                                    wave_diffrac
                                            );
        }
        else
        {
            // Calculate panel potentials
            calculate_potpanel_raddif_nlin(
                                                    input,
                                                    mpi_config,
                                                    mesh_gp,
                                                    sources,
                                                    input->angfreqs[i]
                                            );

            // Calculate added mass and damping coefficients
            calculate_hydromechanic_coeffs_nlin( 
                                                    input,
                                                    mpi_config,
                                                    mesh_gp,
                                                    hmf_interf,
                                                    input->angfreqs[i],
                                                    added_mass,
                                                    damping_rad
                                                );
            
            // Calculate diffraction forces
            calculate_diffraction_forces_nlin(
                                                    input,
                                                    mpi_config,
                                                    mesh_gp,
                                                    hmf_interf,
                                                    input->angfreqs[i],
                                                    wave_diffrac
                                            );

        }

        // Calculate Froude-Krylov forces
        calculate_froude_krylov(
                                    input,
                                    mpi_config,
                                    mesh_gp,
                                    input->angfreqs[i],
                                    froude_krylov
                                );

        // Join data from all processors
        MPI_Reduce(
                            added_mass,
                            added_mass_p0,
                            hydmech_np,
                            mpi_cusfloat,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                        );

        MPI_Reduce(
                            damping_rad,
                            damping_rad_p0,
                            hydmech_np,
                            mpi_cusfloat,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                    );

        MPI_Reduce(
                            froude_krylov,
                            froude_krylov_p0,
                            wave_exc_np,
                            mpi_cuscomplex,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                        );

        MPI_Reduce(
                            wave_diffrac,
                            wave_diffrac_p0,
                            wave_exc_np,
                            mpi_cuscomplex,
                            MPI_SUM,
                            mpi_config->proc_root,
                            MPI_COMM_WORLD
                        );

        // Calculate total wave exciting forces
        if ( mpi_config->is_root( ) )
        {
            sv_add(
                        wave_exc_np,
                        wave_diffrac_p0,
                        froude_krylov_p0,
                        wave_exc_p0
                    );
        }

        // Calculate raos
        if ( mpi_config->is_root( ) )
        {
            calculate_raos(
                                input,
                                structural_mass_p0,
                                added_mass_p0,
                                damping_rad_p0,
                                hydrostiff_p0,
                                wave_diffrac_p0,
                                froude_krylov_p0,
                                input->angfreqs[i],
                                raos
                            );
        }

        if ( input->out_mdrift )
        {
            // Broadcast RAOs values to be available in all the processes
            MPI_Bcast(
                        raos,
                        wave_exc_np,
                        mpi_cuscomplex,
                        mpi_config->proc_root,
                        MPI_COMM_WORLD    
                    );
            
            if ( input->is_fast_solver )
            {

                // Calculate total potential at the target WL points
                calculate_potpanel_total_lin(
                                                input,
                                                mpi_config,
                                                mesh_gp,
                                                input->angfreqs[i],
                                                sources,
                                                pot_steady_mdwe_smat,
                                                pot_mdwe_smat,
                                                mdwe_pot_fp,
                                                mesh_gp->panels_wl_cnp,
                                                mesh_gp->meshes_np,
                                                raos,
                                                mdwe_potpanel_total
                                            );
                
                // Calculate relative wave elevation
                calculate_relative_wave_elevation_lin(
                                                            input,
                                                            mpi_config,
                                                            mdwe_cog_to_fp,
                                                            mesh_gp->panels_wl_cnp,
                                                            mesh_gp->meshes_np,
                                                            mdwe_potpanel_total,
                                                            input->angfreqs[i],
                                                            raos,
                                                            mdwe
                                                        );

                
            }

        }

        // Output values to disk
        if ( mpi_config->is_root( ) )
        {
            if ( input->out_hydmech )
            {
                output->save_hydromechanics_format(
                                                        i,
                                                        _DN_ADDED_MASS,
                                                        added_mass_p0
                                                    );

                output->save_hydromechanics_format(
                                                        i,
                                                        _DN_DAMPING_RAD,
                                                        damping_rad_p0
                                                    );
            }

            if ( input->out_diffrac )
            {
                output->save_wave_exciting_format(
                                                    i,
                                                    _DN_DIFFRAC,
                                                    wave_diffrac_p0
                                                );
            }

            if ( input->out_fk )
            {
                output->save_wave_exciting_format(
                                                    i,
                                                    _DN_FK,
                                                    froude_krylov_p0
                                                );
            }

            if ( input->out_wex )
            {
                output->save_wave_exciting_format(
                                                    i,
                                                    _DN_WEX,
                                                    wave_exc_p0
                                                );
            }

            if ( input->out_raos )
            {
                output->save_wave_exciting_format(
                                                    i,
                                                    _DN_RAO,
                                                    raos
                                                );
            }

        }

        MPI_Barrier( MPI_COMM_WORLD );
        double freq_tend = MPI_Wtime( );

        if ( mpi_config->is_root( ) )
        {
            std::cout << "Execution time [s]: " << ( freq_tend - freq_tstart ) << std::endl;
        }
    }

    // Delete heap memory allocated data
    MPI_Barrier( MPI_COMM_WORLD );
    delete grf_dn_interf;
    delete gwf_dn_interf;
    delete hmf_interf;

    mkl_free( added_mass_p0 );
    mkl_free( damping_rad_p0 );
    mkl_free( froude_krylov_p0 );
    mkl_free( hydrostiff_p0 );
    mkl_free( structural_mass_p0 );
    mkl_free( wave_diffrac_p0 );
    mkl_free( wave_exc_p0 );
    
    mkl_free( added_mass );
    mkl_free( damping_rad );
    mkl_free( froude_krylov );
    mkl_free( raos );
    mkl_free( sources );
    mkl_free( sysmat );
    mkl_free( sysmat_steady );
    mkl_free( wave_diffrac );
    
    delete mesh_gp;
}
