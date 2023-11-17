
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
    /***************** Launch solver ********************/
    /****************************************************/
    if ( input->is_fast_solver )
    {
        freq_domain_linear_solver(
                                        input,
                                        mpi_config,
                                        mesh_gp,
                                        &scl,
                                        hydrostatics,
                                        output
                                );
    }
    else
    {
        freq_domain_nonlinear_solver(
                                        input,
                                        mpi_config,
                                        mesh_gp,
                                        &scl,
                                        hydrostatics,
                                        output
                                    );
    }
    

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


void    freq_domain_linear_solver(
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
    calculate_gf_intensity_steady_sysmat_lin(
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
                                            mesh_gp,
                                            potpanel_lin_gp
                                    );
    double pot_steady_t1 = MPI_Wtime( );
    std::cout << "Time integration potential steady: " << pot_steady_t1 - pot_steady_t0 << std::endl;

    // Calculate steady parto of the potential influence matrix to calculate
    // the mean drift
    if ( input->out_mdrift )
    {
        calculate_influence_potmat_steady(
                                                input,
                                                mesh_gp,
                                                mdrift_we_gp
                                            );
    }

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


void    freq_domain_nonlinear_solver(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                SclCmpx*        scl,
                                                Hydrostatics**  hydrostatics,
                                                Output*         output
                                    )
{
    /****************************************************/
    /****** Allocate space for the simulation data ******/
    /****************************************************/

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
                                                                sim_data->intensities,
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
                                            sim_data->structural_mass_p0
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
                                            sim_data->hydrostiff_p0
                                        );
    }

    /****************************************************/
    /******* Calculate hydrodynamic coefficients ********/
    /****************************************************/

    // Calculate intensities system matrix steady contribution
    calculate_gf_intensity_steady_sysmat_lin(
                                                input,
                                                scl,
                                                mesh_gp,
                                                grf_dn_interf,
                                                sysmat_steady
                                            );

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

        // Update sources values for the integration objects
        hmf_interf->set_source_values( sim_data->intensities );

        // Calculate panel potentials
        calculate_potpanel_raddif_nlin(
                                                input,
                                                mpi_config,
                                                mesh_gp,
                                                sim_data->intensities,
                                                input->angfreqs[i]
                                        );

        // Calculate added mass and damping coefficients
        calculate_hydromechanic_coeffs_nlin( 
                                                input,
                                                mpi_config,
                                                mesh_gp,
                                                hmf_interf,
                                                input->angfreqs[i],
                                                sim_data->added_mass,
                                                sim_data->damping_rad
                                            );
        
        // Calculate diffraction forces
        calculate_diffraction_forces_nlin(
                                                input,
                                                mpi_config,
                                                mesh_gp,
                                                hmf_interf,
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

    /****************************************************/
    /******** Delete heap memory allocated data *********/
    /****************************************************/
    MPI_Barrier( MPI_COMM_WORLD );

    // Delete integration interfaces
    delete grf_dn_interf;
    delete gwf_dn_interf;
    delete hmf_interf;

    // Delete system matrix
    mkl_free( sysmat );
    mkl_free( sysmat_steady );
    
    // Delete mesh group object
    delete mesh_gp;
}
