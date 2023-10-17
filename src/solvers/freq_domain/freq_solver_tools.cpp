
// Include general usage scientific libraries
#include <complex>
#include <fstream>

// Include local modules
#include "freq_solver_tools.hpp"
#include "../../interfaces/grf_interface.hpp"
#include "../../math/integration.hpp"
#include "../../waves.hpp"


void    calculate_diffraction_forces(
                                            Input*          input,
                                            MpiConfig*      mpi_config,
                                            MeshGroup*      mesh_gp,
                                            HMFInterface*   hmf_interf,
                                            cusfloat        w,
                                            cuscomplex*     wave_diffrac
                                    )
{
    std::cout << "Prock Rank Entry: " << mpi_config->proc_rank << std::endl;
    // Define local variables
    int         dofs_np = input->dofs_np;
    cusfloat    rho_w   = input->water_density;

    // Generate lambda function for the integration
    auto target_fcn = [hmf_interf]
                    (
                        cusfloat xi, 
                        cusfloat eta, 
                        cusfloat x,
                        cusfloat y,
                        cusfloat z
                    )
                    {
                        return (*hmf_interf)( xi, eta, x, y, z );
                    };

    // Allocate space for pressure vector
    int max_panels  = 0;
    for ( int i=0; i<mesh_gp->meshes_np; i++ )
    {
        // Get ith mesh panels
        if ( mesh_gp->panels_np[i] > max_panels )
        {
            max_panels = mesh_gp->panels_np[i];
        }
    }
    cuscomplex* pressure    = generate_empty_vector<cuscomplex>( mesh_gp->meshes_np * max_panels );

    // Loop over first dimension of degrees of freedrom
    GaussPoints gp( 2 );
    int         elem_end_pos    = 0;
    int         elem_start_pos  = 0;
    int         index           = 0;
    int         index_1         = 0;
    cuscomplex  press_i         = 0.0;
    // cuscomplex  pressure = 0.0;
    for ( int ih=0; ih<input->heads_np; ih++ )
    {
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

            // Set start index for the sources evaluation
            hmf_interf->set_start_index_i( 
                                                mesh_gp->source_nodes_tnp * ( dofs_np + ih ),
                                                0,
                                                mesh_gp->source_nodes_tnp
                                            );

            // Loop over panels to integrate the wave radiation
            // pressure along the floating object external shape
            for ( int ie=elem_start_pos; ie<elem_end_pos; ie++ )
            {
                // Set new panel
                hmf_interf->set_panel( mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie] );

                // Integrate pressure over panel
                index           = max_panels * ib + ie;
                pressure[index] = adaptive_quadrature_panel(
                                                                mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie],
                                                                target_fcn,
                                                                1.0,
                                                                &gp,
                                                                true
                                                            );
            }
        }

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

            for ( int id=0; id<dofs_np; id++ )
            {
                index = (
                            ih * ( dofs_np * mesh_gp->meshes_np )
                            +
                            ib * dofs_np
                            + 
                            id
                        );

                wave_diffrac[index]   = 0.0;
                for ( int ie=elem_start_pos; ie<elem_end_pos; ie++ )
                {
                    // Integrate pressure over panel
                    index_1             = max_panels * ib + ie;
                    press_i             = pressure[index_1] * mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->normal_vec[id];
                    wave_diffrac[index] += cuscomplex( 0.0, -w * rho_w ) * press_i;
                }
            }
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
    MeshGroup*  mesh_gp =   new   MeshGroup( 
                                                all_meshes,
                                                input->bodies_np
                                            );

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
    int         hydmech_np      = pow2s( input->dofs_np * mesh_gp->meshes_np );
    int         wave_exc_np     = input->heads_np * mesh_gp->meshes_np * input->dofs_np;

    cusfloat*   added_mass      = generate_empty_vector<cusfloat>( hydmech_np );
    cusfloat*   damping_rad     = generate_empty_vector<cusfloat>( hydmech_np );
    cuscomplex* froude_krylov   = generate_empty_vector<cuscomplex>( wave_exc_np );
    cuscomplex* sources         = generate_empty_vector<cuscomplex>( ( input->dofs_np + input->heads_np ) * scl.num_rows_local );
    cuscomplex* sysmat          = generate_empty_vector<cuscomplex>( scl.num_rows_local * scl.num_cols_local );
    cuscomplex* wave_diffrac    = generate_empty_vector<cuscomplex>( wave_exc_np );

    cusfloat*   added_mass_p0       = nullptr;
    cusfloat*   damping_rad_p0      = nullptr;
    cuscomplex* froude_krylov_p0    = nullptr;
    cusfloat*   hydrostiff_p0       = nullptr;
    cuscomplex* raos_p0             = nullptr;
    cusfloat*   structural_mass_p0  = nullptr;
    cuscomplex* wave_diffrac_p0     = nullptr;
    cuscomplex* wave_exc_p0         = nullptr;

    added_mass_p0       = generate_empty_vector<cusfloat>( hydmech_np );
    damping_rad_p0      = generate_empty_vector<cusfloat>( hydmech_np );
    froude_krylov_p0    = generate_empty_vector<cuscomplex>( wave_exc_np );
    hydrostiff_p0       = generate_empty_vector<cusfloat>( hydmech_np );
    raos_p0             = generate_empty_vector<cuscomplex>( wave_exc_np );
    structural_mass_p0  = generate_empty_vector<cusfloat>( hydmech_np );
    wave_diffrac_p0     = generate_empty_vector<cuscomplex>( wave_exc_np );
    wave_exc_p0         = generate_empty_vector<cuscomplex>( wave_exc_np );

    /****************************************************/
    /********* Create Green function interface *********/
    /****************************************************/
    GWFDnInterface* green_dn_interf = new   GWFDnInterface(
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
                                                                input->water_depth,
                                                                input->grav_acc
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

    // Loop over frequencies
    for ( int i=0; i<input->angfreqs_np; i++ )
    {
        MPI_Barrier( MPI_COMM_WORLD );
        double freq_tstart = MPI_Wtime( );

        std::cout << "Calculating angular frequency: " << input->angfreqs[i] << std::endl;
        // Recalculate wave properties for the current 
        // angular frequency
        hmf_interf->set_ang_freq( input->angfreqs[i] );
        green_dn_interf->set_ang_freq( input->angfreqs[i] );

        // Calculate sources intensity
        MPI_Barrier( MPI_COMM_WORLD );
        calculate_sources_intensity(
                                        input,
                                        &scl,
                                        mesh_gp,
                                        green_dn_interf,
                                        input->angfreqs[i],
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

        // Update sources values for the integration objects
        hmf_interf->set_source_values( sources );

        calculate_panel_potentials(
                                            input,
                                            mpi_config,
                                            mesh_gp,
                                            sources,
                                            input->angfreqs[i]
                                    );

        // Calculate added mass and damping coefficients
        calculate_hydromechanic_coeffs( 
                                            input,
                                            mpi_config,
                                            mesh_gp,
                                            hmf_interf,
                                            input->angfreqs[i],
                                            added_mass,
                                            damping_rad
                                        );
        
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

        // Calculate Froude-Krylov forces
        calculate_froude_krylov(
                                            input,
                                            mpi_config,
                                            mesh_gp,
                                            input->angfreqs[i],
                                            froude_krylov
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

        // Calculate diffraction forces
        MPI_Barrier( MPI_COMM_WORLD );
        calculate_diffraction_forces(
                                            input,
                                            mpi_config,
                                            mesh_gp,
                                            hmf_interf,
                                            input->angfreqs[i],
                                            wave_diffrac
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
                                raos_p0
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
                                                    raos_p0
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
    delete green_dn_interf;
    delete hmf_interf;

    mkl_free( added_mass_p0 );
    mkl_free( damping_rad_p0 );
    mkl_free( froude_krylov_p0 );
    mkl_free( hydrostiff_p0 );
    mkl_free( raos_p0 );
    mkl_free( structural_mass_p0 );
    mkl_free( wave_diffrac_p0 );
    mkl_free( wave_exc_p0 );
    
    mkl_free( added_mass );
    mkl_free( damping_rad );
    mkl_free( froude_krylov );
    mkl_free( sources );
    mkl_free( sysmat );
    mkl_free( wave_diffrac );

    
    delete mesh_gp;
    delete [] all_meshes;
}


void    calculate_froude_krylov(
                                            Input*          input,
                                            MpiConfig*      mpi_config,
                                            MeshGroup*      mesh_gp,
                                            cusfloat        ang_freq,
                                            cuscomplex*     froude_krylov
                                )
{
    // Calculate wave number
    cusfloat k      = w2k( ang_freq, input->water_depth, input->grav_acc );
    cusfloat rho_w  = input->water_density;

    // Generate pressure vector to storage the pressures over the panels
    int max_panels  = 0;
    for ( int i=0; i<mesh_gp->meshes_np; i++ )
    {
        // Get ith mesh panels
        if ( mesh_gp->panels_np[i] > max_panels )
        {
            max_panels = mesh_gp->panels_np[i];
        }
    }
    cuscomplex* pressure = generate_empty_vector<cuscomplex>( max_panels * mesh_gp->meshes_np );

    // Loop around headings to get the Froude-Krylov force
    // for each of them
    GaussPoints gp( 2 );
    int         elem_end_pos    = 0;
    int         elem_start_pos  = 0;
    int         index           = 0;
    int         index_1         = 0;
    cuscomplex  press_i         = cuscomplex( 0.0, 0.0 );

    for ( int ih=0; ih<input->heads_np; ih++ )
    {
        // Definition of the target function to integrate incident
        // wave potential for the desised heading
        auto    target_fcn  =   [ 
                                    input,
                                    ang_freq, 
                                    k,
                                    ih
                                ]
                                (
                                    cusfloat ,
                                    cusfloat ,
                                    cusfloat x,
                                    cusfloat y,
                                    cusfloat z
                                )
                                {
                                    return wave_potential_airy_space(
                                                                        1.0,
                                                                        ang_freq,
                                                                        k,
                                                                        input->water_depth,
                                                                        input->grav_acc,
                                                                        x,
                                                                        y,
                                                                        z,
                                                                        input->heads[ih]
                                                                    );
                                };
        
        // Loop around bodies to calculate the pressure of the incident wave
        // over each panel of the mesh
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

            // Loop over panels integrating incident wave pressure
            for ( int ie=elem_start_pos; ie<elem_end_pos; ie++ )
            {
                // Integrate pressure over panel
                index           = max_panels * ib + ie;
                pressure[index] = adaptive_quadrature_panel(
                                                                mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie],
                                                                target_fcn,
                                                                0.001,
                                                                &gp,
                                                                false
                                                            );
                // std::cout << "Froude-Krylov - Index: " << index << " - Index_1: " << mesh_gp->panels_cnp[ib]+ie;
                // std::cout << " - ie: " << ie << " - ib: " << ib << " - pressure: " << pressure[index] << std::endl << std::flush;
            }
        }

        // Loop over bodies to calculate the projection of the Froude-Krylov force
        // over each DOF
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
                // Generate index for the current Froude-Krylov force and
                // to clean the memory space
                index                   = (
                                                ih * ( mesh_gp->meshes_np * input->dofs_np )
                                                +
                                                ib * input->dofs_np
                                                +
                                                id
                                            );
                
                froude_krylov[index]    = 0.0;
                // Loop over panels to sum all the pressures
                for ( int ie=elem_start_pos; ie<elem_end_pos; ie++ )
                {
                    // Integrate pressure over panel
                    index_1                 = max_panels * ib + ie;
                    press_i                 = pressure[index_1] * mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->normal_vec[id];
                    froude_krylov[index]    += cuscomplex( 0.0, -ang_freq * rho_w ) * press_i;
                    // std::cout << "IB: " << ib << " - Dof: " << id << " - pressure: " << pressure[index_1] << " - normal: " << mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->normal_vec[id] << std::endl;
                }
                // std::cout << "Heading: " << input->heads[ih] << " - IB: " << ib << " - Dof: " << id << " - Froude-Krylov: " << froude_krylov[index] << std::endl;
            }
        }
    }

    // Deallocate local heap memory
    mkl_free( pressure );
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


void    calculate_hydromechanic_coeffs(
                                            Input*          input,
                                            MpiConfig*      mpi_config,
                                            MeshGroup*      mesh_gp,
                                            HMFInterface*   hmf_interf,
                                            cusfloat        ang_freq,
                                            cusfloat*       added_mass,
                                            cusfloat*       damping_rad
                                        )
{
    // Define local variables
    int         dofs_np = input->dofs_np;
    cusfloat    rho_w   = input->water_density;

    // Generate lambda function for the integration
    auto target_fcn = [hmf_interf]
                    (
                        cusfloat xi, 
                        cusfloat eta, 
                        cusfloat x,
                        cusfloat y,
                        cusfloat z
                    )
                    {
                        return (*hmf_interf)( xi, eta, x, y, z );
                    };

    // Allocate space for pressure vector
    int max_panels  = 0;
    for ( int i=0; i<mesh_gp->meshes_np; i++ )
    {
        // Get ith mesh panels
        if ( mesh_gp->panels_np[i] > max_panels )
        {
            max_panels = mesh_gp->panels_np[i];
        }
    }
    cuscomplex* pressure    = generate_empty_vector<cuscomplex>( mesh_gp->meshes_np * dofs_np * max_panels );

    // Loop over first dimension of degrees of freedrom
    GaussPoints gp( 2 );
    int         elem_end_pos    = 0;
    int         elem_start_pos  = 0;
    int         index           = 0;
    int         index_1         = 0;
    cuscomplex  press_i         = 0.0;
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

        for ( int id=0; id<dofs_np; id++ )
        {
            for ( int jb=0; jb<mesh_gp->meshes_np; jb++ )
            {
                // Set ith dof
                hmf_interf->set_start_index_i( 
                                                    mesh_gp->source_nodes_tnp*id,
                                                    mesh_gp->source_nodes_cnp[jb],
                                                    mesh_gp->source_nodes_cnp[jb+1]
                                                );

                // Loop over panels to integrate the wave radiation
                // pressure along the floating object external shape
                for ( int ie=elem_start_pos; ie<elem_end_pos; ie++ )
                {
                    // Set new panel
                    hmf_interf->set_panel( mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie] );

                    // Integrate pressure over panel
                    index           = ( max_panels * mesh_gp->meshes_np ) * id + max_panels * jb + ie;
                    pressure[index] = adaptive_quadrature_panel(
                                                                    mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie],
                                                                    target_fcn,
                                                                    1.0,
                                                                    &gp,
                                                                    true
                                                                );
                }
            }
        }

        // Get hydromechanic coefficients for all the degrees of freedom
        for ( int id=0; id<dofs_np; id++ )
        {
            for ( int jd=0; jd<dofs_np; jd++ )
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
                        // Integrate pressure over panel
                        index_1             = ( max_panels * mesh_gp->meshes_np ) * id + max_panels * jb + ie;
                        press_i             = pressure[index_1] * mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->normal_vec[jd];
                        added_mass[index]   -=  rho_w * press_i.imag( ) / ang_freq;
                        damping_rad[index]  -=  rho_w * press_i.real( );
                    }

                }
            }
        }
    }

    // Deallocate local allocated heap memory
    mkl_free( pressure );
}


void    calculate_panel_potentials(
                                            Input*          input,
                                            MpiConfig*      mpi_config,
                                            MeshGroup*      mesh_gp,
                                            cuscomplex*     sources,
                                            cusfloat        ang_freq
                                    )
{
    GaussPoints gp( 10 );

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
    cuscomplex panel_potential  = complex( 0.0, 0.0 );
    cuscomplex pot_i_steady     = complex( 0.0, 0.0 );
    cuscomplex pot_i_wave       = complex( 0.0, 0.0 );
    for ( int i=0; i<mesh_gp->panels_tnp; i++ )
    {
        // Set new field point for the calculation of the potential
        green_interf_steady->set_field_point( mesh_gp->panels[i]->center );
        green_interf_wave->set_field_point( mesh_gp->panels[i]->center );

        // Clean panel potential
        panel_potential = complex( 0.0, 0.0 );

        // Loop over source points to calculate potential
        // over ith panel center
        for ( int j=0; j<mesh_gp->source_nodes_tnp; j++ )
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
                                                            0.01,
                                                            &gp
                                                        );
            pot_i_wave      = adaptive_quadrature_panel(
                                                            mesh_gp->source_nodes[j]->panel,
                                                            wave_fcn,
                                                            0.01,
                                                            &gp
                                                        );
            panel_potential +=  ( pot_i_steady + pot_i_wave ) /4.0 / PI;
        }
    }

    // Delete heap allocated memory
    delete green_interf_steady;
    delete green_interf_wave;
}


void    calculate_sources_intensity(
                                            Input*          input,
                                            SclCmpx*        scl,
                                            MeshGroup*      mesh_gp,
                                            GWFDnInterface* green_interf,
                                            cusfloat        w,
                                            cuscomplex*     sysmat,
                                            cuscomplex*     sources_int
                                   )
{
    GaussPoints gp = GaussPoints( 4 );

    /***************************************/
    /******** Fill system matrix  **********/
    /***************************************/
    // Create auxiliar lambda function
    auto lmb_fcn = [green_interf]
                    ( 
                        cusfloat xi, 
                        cusfloat eta,
                        cusfloat x, 
                        cusfloat y, 
                        cusfloat z 
                    )
                    {
                        return (*green_interf)( xi, eta, x, y, z );
                    };

    // Loop over panels to integrate value
    int         col_count   = 0;
    cuscomplex  int_value( 0.0, 0.0 );
    SourceNode* source_i    = nullptr;
    int         row_count   = 0;

    for ( int i=scl->start_col_0; i<scl->end_col_0; i++ )
    {
        // Get memory address of the ith panel
        source_i = mesh_gp->source_nodes[i];
        green_interf->set_source_i( source_i );

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
        {
            // Get memory address of the panel jth
            green_interf->set_source_j( mesh_gp->source_nodes[j] );

            // Integrate green function normal derivative along the current panel
            if ( i == j )
            {
                int_value   =   -std::complex( 0.5, 0.0 );
            }
            else
            {
                int_value   =   adaptive_quadrature_panel(
                                                                source_i->panel,
                                                                lmb_fcn,
                                                                1.0,
                                                                &gp,
                                                                false
                                                            );
                int_value   =   int_value / 4.0 / PI;
            }

            sysmat[col_count*scl->num_rows_local+row_count] = int_value;

            // Advance row count
            row_count++;
        }

        // Advance column count
        col_count++;
    }

    /***************************************/
    /***** Fill Hydromechanics RHS  ********/
    /***************************************/

    // Declare local variables to be used
    int         count           = 0;

    // Fill RHS vector
    count       = 0;
    // outfile << "RHS" << std::endl;
    for ( int i=0; i< 6; i++ )
    {
        for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
        {
            sources_int[count] = complex( 0.0, w * mesh_gp->source_nodes[j]->normal_vec[i] );
            count++;
        }
    }

    /***************************************/
    /****** Fill Wave Exciting RHS  ********/
    /***************************************/
    // Define local variables to manage array indexes
                count       = input->dofs_np * mesh_gp->source_nodes_tnp;
    cusfloat    k           = w2k( w, input->water_depth, input->grav_acc );
    cuscomplex  wave_dx     = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_dy     = cuscomplex( 0.0, 0.0 );
    cuscomplex  wave_dz     = cuscomplex( 0.0, 0.0 );

    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
        {
            // Get wave potential derivatives for the panel
            wave_dx             =   wave_potential_airy_space_dx(
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

            wave_dy             =   wave_potential_airy_space_dy(
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

            wave_dz             =   wave_potential_airy_space_dz(
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
            
            count++;
        }
    }

    // Solve system of equations
    scl->Solve( sysmat, sources_int );
}


void    calculate_raos(
                                            Input*          input,
                                            cusfloat*       structural_mass,
                                            cusfloat*       added_mass,
                                            cusfloat*       damping_rad,
                                            cusfloat*       hydstiffness,
                                            cuscomplex*     wave_diffrac,
                                            cuscomplex*     froude_krylov,
                                            cusfloat        ang_freq,
                                            cuscomplex*     rao
                        )
{
    // Allocate space to hold the system matrix
    cuscomplex* sysmat  = generate_empty_vector<cuscomplex>( pow2s( input->bodies_np * input->dofs_np ) );

    // Clear input rao vector to avoid problems with residual data
    clear_vector( input->bodies_np * input->dofs_np * input->heads_np, rao );

    // Fill in matrix system
    int         index   = 0;
    for ( int i=0; i<input->bodies_np; i++ )
    {
        for ( int j=0; j<input->bodies_np; j++ )
        {
            for ( int k=0; k<input->dofs_np; k++ )
            {
                for ( int m=0; m<input->dofs_np; m++ )
                {
                    // Get structural mass
                    index = (
                                i * input->bodies_np * pow2s( input->dofs_np )
                                +
                                j * input->dofs_np
                                +
                                k * input->bodies_np * input->dofs_np
                                +
                                m
                            );
                    sysmat[index] = cuscomplex( 
                                                    -pow2s( ang_freq ) * ( structural_mass[index] + added_mass[index] ) + hydstiffness[index],
                                                    ang_freq * damping_rad[index]
                                                );
                }
            }
        }
    }

    // Fill in right hand side vector
    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<input->bodies_np; j++ )
        {
            for ( int k=0; k<input->dofs_np; k++ )
            {
                index = (
                            i * input->bodies_np * input->dofs_np
                            +
                            j * input->dofs_np
                            +
                            k
                        );
                rao[index] = wave_diffrac[index] + froude_krylov[index];
            }
        }
    }

    // Solve system of equations
    int     rows_np = input->bodies_np * input->dofs_np;
    int     info    = 0;
    int*    ipiv    = generate_empty_vector<int>( rows_np );
    gesv<cuscomplex>( 
                        &rows_np,
                        &(input->heads_np),
                        sysmat,
                        &(rows_np),
                        ipiv,
                        rao,
                        &(rows_np),
                        &info
                    );

    mkl_free( ipiv );

    if ( info != 0 )
    {
        std::cerr << "ERROR - Calculating RAO" << std::endl;
        std::cerr << "Error solving system of equations - Info: " << info << std::endl;
        throw std::runtime_error( "" );
    }

    // Deallocate local heap memory
    mkl_free( sysmat );

}