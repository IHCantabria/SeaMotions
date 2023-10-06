
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
                std::cout << "Potential[" << ib << ", " << ie << "]: " << pressure[index];
                std::cout << " - I1: " << mesh_gp->panels_cnp[ib]+ie << " - index: " << index <<   " - Proc.Rank: " << mpi_config->proc_rank << std::endl << std::flush;
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
                    // std::cout << "Pressure[" << ie << "]: " << w*rho_w*press_i << std::endl;
                    wave_diffrac[index] += cuscomplex( 0.0, -w * rho_w ) * press_i;
                }
                // std::cout << "Heading: " << input->heads[ih] << " - IB: " << ib << " - Dof: " << id << " - Wave Diffrac: " << wave_diffrac[index] << std::endl;
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
                    MPI_COMM_WORLD
                );

    /****************************************************/
    /****** Allocate space for the simulation data ******/
    /****************************************************/
    int         hydmech_np      = pow2s( input->dofs_np * mesh_gp->meshes_np );
    int         wave_exc_np     = input->heads_np * mesh_gp->meshes_np * input->dofs_np;

    cusfloat*   added_mass      = generate_empty_vector<cusfloat>( hydmech_np );
    cuscomplex* all_sources     = generate_empty_vector<cuscomplex>( ( input->dofs_np + input->heads_np ) * mesh_gp->source_nodes_tnp );
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
                                                                all_sources,
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
        std::cout << "Calculating angular frequency: " << input->angfreqs[i] << std::endl;
        // Recalculate wave properties for the current 
        // angular frequency
        hmf_interf->set_ang_freq( input->angfreqs[i] );
        green_dn_interf->set_ang_freq( input->angfreqs[i] );

        // Calculate sources intensity
        MPI_Barrier( MPI_COMM_WORLD );
        calculate_sources_intensity(
                                        input,
                                        mpi_config,
                                        &scl,
                                        mesh_gp,
                                        green_dn_interf,
                                        input->angfreqs[i],
                                        sysmat,
                                        sources
                                    );
        
        // Gather source values from each processor
        MPI_Allgather(
                        sources,
                        scl.num_rows_local * ( input->dofs_np + input->heads_np ),
                        mpi_cuscomplex,
                        all_sources,
                        scl.num_rows * ( input->dofs_np + input->heads_np ),
                        mpi_cuscomplex,
                        MPI_COMM_WORLD
                    );
        // for ( int j=0; j<7; j++ )
        // {
        //     scl.GetGlobalRhs( &(sources[j*scl.num_rows_local]), &(all_sources[j*scl.num_rows_local]) );
        // }

        // Update sources values for the integration objects
        hmf_interf->set_source_values( all_sources );

        std::cout << "Scl.NumRows: " << scl.num_rows << " - Scl.NumRowsLocal: " << scl.num_rows_local << std::endl;
        for ( int j=0; j<scl.num_rows*(input->dofs_np+input->heads_np); j++ )
        {
            // std::cout << "Panel[" << j << "]: " << all_sources[j] << " - " << sources[j] << " - " << std::abs( all_sources[j] ) << " - " << std::arg( all_sources[j] )*180.0/PI << std::endl;
            std::cout << "Panel[" << j << "]: " << all_sources[j] << std::endl << std::flush;
        }

        // if ( mpi_config->procs_total == 1 )
        // {
        //     std::ofstream outfile( "my_coeffs.dat" );

        //     outfile << scl.num_rows*(input->dofs_np+input->heads_np) << std::endl;
        //     for ( int j=0; j<scl.num_rows*(input->dofs_np+input->heads_np); j++ )
        //     {
        //         outfile << all_sources[j].real( ) << "  " << all_sources[j].imag( ) << std::endl;
        //     }

        //     outfile.close( );
        // }

        // if ( mpi_config->procs_total > 1 )
        // {
        //     std::ifstream infile( "my_coeffs.dat" );

        //     cusfloat a, b;
        //     cuscomplex c;
        //     int rows_np = 0;
        //     infile >> rows_np;
        //     for ( int j=0; j<rows_np; j++ )
        //     {
        //         infile >> a >> b;
        //         c = all_sources[j];
        //         all_sources[j] = cuscomplex( a, b );
        //         if ( mpi_config->proc_rank == 0 )
        //         {
        //             std::cout << "Index: " << j << " - Diff: " << c-all_sources[j] << " - AllSources: " << all_sources[j] << std::endl << std::flush;
        //         }
        //     }
        // }

        // for ( int j=0; j<scl.num_rows*(input->dofs_np+input->heads_np); j++ )
        // {
        //     all_sources[j] = 1.0;
        // }

        calculate_panel_potentials(
                                            input,
                                            mpi_config,
                                            mesh_gp,
                                            all_sources,
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
        std::cout << "(P" << mpi_config->proc_rank << ") " << "Reducing Added Mass..." << std::flush;
        MPI_Reduce(
                        added_mass,
                        added_mass_p0,
                        hydmech_np,
                        mpi_cusfloat,
                        MPI_SUM,
                        mpi_config->proc_root,
                        MPI_COMM_WORLD
                    );
        std::cout << "(P" << mpi_config->proc_rank << ") " << "--> Done!" << std::endl << std::flush;

        std::cout << "(P" << mpi_config->proc_rank << ") " << "Reducing Damping..." << std::flush;
        MPI_Reduce(
                        damping_rad,
                        damping_rad_p0,
                        hydmech_np,
                        mpi_cusfloat,
                        MPI_SUM,
                        mpi_config->proc_root,
                        MPI_COMM_WORLD
                    );
        std::cout << "(P" << mpi_config->proc_rank << ") " << "--> Done!" << std::endl;

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

        if ( mpi_config->is_root( ) )
        {
            for ( int i=0; i<input->heads_np; i++ )
            {
                for ( int j=0; j<input->dofs_np; j++ )
                {
                    std::cout << "FK - HD: " << i << " - Dof: " << j << " - Value: " << froude_krylov_p0[input->dofs_np*i+j];
                    std::cout << " |Value|: " << std::abs( froude_krylov_p0[input->dofs_np*i+j] ) << std::endl << std::flush;
                }
            }
        }

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
        
        if ( mpi_config->is_root( ) )
        {
            for ( int i=0; i<input->heads_np; i++ )
            {
                for ( int j=0; j<input->dofs_np; j++ )
                {
                    std::cout << "WD - HD: " << i << " - Dof: " << j << " - Value: " << wave_diffrac_p0[input->dofs_np*i+j];
                    std::cout << " |Value|: " << std::abs( wave_diffrac_p0[input->dofs_np*i+j] ) << std::endl << std::flush;
                }
            }
        }

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
                                mpi_config,
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
                                                        added_mass_p0
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

    }

    // Delete heap memory allocated data
    mkl_free( added_mass_p0 );
    mkl_free( damping_rad_p0 );
    mkl_free( froude_krylov_p0 );
    mkl_free( hydrostiff_p0 );
    mkl_free( raos_p0 );
    mkl_free( structural_mass_p0 );
    mkl_free( wave_diffrac_p0 );
    mkl_free( wave_exc_p0 );
    
    mkl_free( added_mass );
    mkl_free( all_sources );
    mkl_free( damping_rad );
    mkl_free( froude_krylov );
    mkl_free( sources );
    mkl_free( sysmat );
    mkl_free( wave_diffrac );

    delete green_dn_interf;
    delete hmf_interf;
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
    // cuscomplex  pressure = 0.0;
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
                // std::cout << "start index: " << mesh_gp->source_nodes_tnp*id + mesh_gp->source_nodes_cnp[jb] << std::endl; 
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
                    // std::cout << "Pressure [" << index << "]: " << pressure[index] << std::endl;
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
                        // std::cout << "press_i: " << press_i << " - Normal: " << mesh_gp->panels[ie]->normal_vec[jd] << std::endl;
                        added_mass[index]   -=  rho_w * press_i.imag( ) / ang_freq;
                        damping_rad[index]  -=  rho_w * press_i.real( );
                    }

                    // std::cout << "Body: " << ib << " - " << jb << " - Dof i: " << id << " - Dof j: " << jd << " - AddedMass: " << added_mass[index] << " - Damping: " << damping_rad[index] << " - Index: " << index << std::endl;

                }
            }
        }
    }

    // Deallocate local allocated heap memory
    // std::cout << "(P" << mpi_config->proc_rank << ") " << "Deallocating HydMech Pressure..." << std::flush;
    mkl_free( pressure );
    // std::cout << "(P" << mpi_config->proc_rank << ") " << "-->Done!" << std::endl << std::flush;
}


void    calculate_panel_potentials(
                                            Input*          input,
                                            MpiConfig*      mpi_config,
                                            MeshGroup*      mesh_gp,
                                            cuscomplex*     all_sources,
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
                                                                all_sources[0],
                                                                mesh_gp->panels[0]->center,
                                                                input->water_depth
                                                            );
    GWFInterface*   green_interf_wave   = new   GWFInterface(
                                                                mesh_gp->source_nodes[0],
                                                                all_sources[0],
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

        // std::cout << "Field Point: " << mesh_gp->panels[i]->center[0] << " - " << mesh_gp->panels[i]->center[1] << " - " << mesh_gp->panels[i]->center[2] << std::endl;

        // Clean panel potential
        panel_potential = complex( 0.0, 0.0 );

        // Loop over source points to calculate potential
        // over ith panel center
        for ( int j=0; j<mesh_gp->source_nodes_tnp; j++ )
        {
            // Set current source value
            green_interf_steady->set_source(
                                            mesh_gp->source_nodes[j],
                                            all_sources[j]
                                        );
            green_interf_wave->set_source(
                                            mesh_gp->source_nodes[j],
                                            all_sources[j]
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
            // std::cout << "pot_i[" << j << "]: " << pot_i/4.0/PI << " - Source Value: " << all_sources[j] << std::endl;
        }
        // std::cout << "Potential[" << i << "]: " << panel_potential << " - " << std::abs( panel_potential ) << " - " << std::arg( panel_potential )*180.0/PI << std::endl;
        // throw std::runtime_error( "" );
    }

    // Delete heap allocated memory
    delete green_interf_steady;
    delete green_interf_wave;
}


void    calculate_sources_intensity(
                                            Input*          input,
                                            MpiConfig*      mpi_config,
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
    std::cout << "GWFDnInterface: " << green_interf << std::endl;
    std::cout << "Creating system matrix..." << std::endl;
    int         col_count   = 0;
    cuscomplex  int_value( 0.0, 0.0 );
    SourceNode* source_i    = nullptr;
    int         row_count   = 0;

    // std::ofstream outfile( "matrix.dat" );
    // outfile << "Num.Rows: " << scl->num_rows << std::endl;
    // MPI_Barrier( MPI_COMM_WORLD );
    int count_total = 0;
    for ( int i=scl->start_col_0; i<scl->end_col_0; i++ )
    {
        // std::cout << "Source[i]: " << i << " - " << mesh_gp->source_nodes[i] << std::endl;
        // Get memory address of the ith panel
        source_i = mesh_gp->source_nodes[i];
        green_interf->set_source_i( source_i );

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
        {
            // std::cout << "Panel[j]: " << j << " - " << mesh->source_nodes[j] << std::endl;
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
                                                                0.001,
                                                                &gp,
                                                                false
                                                            );
                int_value   =   int_value / 4.0 / PI;
            }

            // if ( std::abs( int_value.real( ) ) < 1e-10 )
            // {
            //     std::cout << "I: " << i << " - J: " << j << " - Zero value" << std::endl;
            // }
            std::cout << "index: " << col_count*scl->num_rows_local+row_count <<  " - int_value: " << int_value << std::endl;
            sysmat[col_count*scl->num_rows_local+row_count] = int_value;
            // outfile << int_value.real( ) << " " << int_value.imag( ) << std::endl;
            // std::cout << "Index Count: " << col_count*scl->num_rows_local+row_count << std::endl;
            count_total += 1;
            // Advance row count
            row_count++;
        }
        std::cout << std::endl;
        // Advance column count
        col_count++;
    }

    // MPI_Barrier( MPI_COMM_WORLD );
    int _rows_np = scl->end_row - scl->start_row;
    int _cols_np = scl->end_col - scl->start_col;
    std::cout << "--> Done!" << std::endl;
    std::cout << "count_total: " << count_total << std::endl;
    std::cout << "Num.Rows: " << _rows_np << std::endl;
    std::cout << "Num.Cols: " << _cols_np << std::endl;

    /***************************************/
    /***** Fill Hydromechanics RHS  ********/
    /***************************************/
    std::cout << "Calculating RHS..." << std::endl;
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
            // std::cout << "Source Node [" << j << "]: " << mesh_gp->source_nodes[j]->normal_vec[i] << std::endl;
            // std::cout << " - Center Panel:  " << mesh_gp->source_nodes[j]->position[0];
            // std::cout << " - " << mesh_gp->source_nodes[j]->position[1];
            // std::cout << " - " << mesh_gp->source_nodes[j]->position[2];
            // std::cout << std::endl;
            // std::cout << " - Normal Vec: " << mesh_gp->source_nodes[j]->normal_vec[0];
            // std::cout << " - " << mesh_gp->source_nodes[j]->normal_vec[1];
            // std::cout << " - " << mesh_gp->source_nodes[j]->normal_vec[2];
            // std::cout << std::endl;
            // outfile << sources_int[start_pos+count].real( ) << " " << sources_int[start_pos+count].imag( ) << std::endl;
            std::cout << "sources_int[" << count <<"]: " << sources_int[count] << std::endl;
            count++;
        }
        std::cout << std::endl;
        // break;
    }
    // outfile.close( );

    std::cout << "--> Done!" << std::endl;

    /***************************************/
    /****** Fill Wave Exciting RHS  ********/
    /***************************************/
    // Define local variables to manage array indexes
                count       = 60;
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
            if ( mpi_config->is_root( ) && count > 59 )
            {
                std::cout << "sources_int[" << count << "]: " << sources_int[count] << std::endl;
            }
            count++;
        }
    }

    // Solve system of equations
    std::cout << "Calculating system of equations..." << std::endl;
    scl->Solve( sysmat, sources_int );
    std::cout << "--> Done!" << std::endl;

    for ( int i=0; i<70; i++ )
    {
        if ( mpi_config->is_root( ) && count > 0 )
        {
            std::cout << "sol[" << i << "]: " << sources_int[i] << std::endl;
        }
    }

    // count       = input->dofs_np * mesh_gp->source_nodes_tnp;
    // for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
    // {
    //     std::cout << "Sources[" << count+j << "]: " << sources_int[count+j] << std::endl;
    // }
    
}


void    calculate_raos(
                                            Input*          input,
                                            MpiConfig*      mpi_config,
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
    cusfloat    mass_i = 0.0;
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
                    if ( k==0 & m==0 )
                    {
                        std::cout << "mass_i: " << mass_i << std::endl;
                        std::cout << "added_mass: " << added_mass[index] << std::endl;
                        std::cout << "daming: " << damping_rad[index] << std::endl;
                        std::cout << "hydstiffness: " << hydstiffness[index] << std::endl;
                        std::cout << "sysmat[index]: " << sysmat[index] << std::endl;
                    }
                }
            }
        }
    }

    for ( int i=0; i<6; i++ )
    {
        std::cout << "sysmat[" << i << ", :]: ";
        for ( int j=0; j<6; j++ )
        {
            std::cout << sysmat[6*i+j] << " ";
        }
        std::cout << std::endl;
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
                if ( k == 0 )
                {
                    std::cout << "Wave Diffrac: " << wave_diffrac[index] << std::endl;
                    std::cout << "Froude Krylov: " << froude_krylov[index] << std::endl;
                    std::cout << "Wave Total: " << rao[index] << std::endl;
                }
            }
        }
    }

    std::cout << "Fw: " << std::endl;
    for ( int i=0; i<6; i++ )
    {
        std::cout << "Fw[" << i << "]: " << rao[i] << std::endl;
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

    for ( int i=0; i<input->heads_np; i++ )
    {
        std::cout << "Heading: " << input->heads[i] << std::endl;
        for ( int j=0; j<input->bodies_np; j++ )
        {
            std::cout << " -> Bodies: " << j << std::endl;
            for ( int k=0; k<input->dofs_np; k++ )
            {
                index = (
                            i * input->bodies_np * input->dofs_np
                            +
                            j * input->dofs_np
                            +
                            k
                        );
                std::cout << " --> Dof[" << k << "]: " << rao[index] << " - Mag: " << std::abs( rao[index] ) << " - Phase: " << std::arg( rao[index] ) << std::endl;
            }
        }
    }

    // Deallocate local heap memory
    mkl_free( sysmat );

}