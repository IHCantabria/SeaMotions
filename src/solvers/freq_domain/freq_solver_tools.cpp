
// Include general usage scientific libraries
#include <complex>
#include <fstream>
#include "mkl_cblas.h"

// Include local modules
#include "freq_solver_tools.hpp"
#include "../../interfaces/grf_interface.hpp"
#include "../../math/integration.hpp"
#include "../../waves.hpp"


void    calculate_diffraction_forces_lin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     panel_pot,
                                                cusfloat        w,
                                                cuscomplex*     wave_diffrac
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
    
    // Loop over headings to define the wave diffraction forces
    for ( int ih=0; ih<input->heads_np; ih++ )
    {
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
                    index_1             = (
                                                input->dofs_np * mesh_gp->panels_tnp 
                                                +
                                                ih * mesh_gp->panels_tnp
                                                +
                                                mesh_gp->panels_cnp[ib]
                                                +
                                                ie
                                            );
                    press_i             = (
                                                panel_pot[index_1]
                                                *
                                                mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->normal_vec[id]
                                                *
                                                mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->area
                                            );
                    wave_diffrac[index] += cuscomplex( 0.0, -w * rho_w ) * press_i;
                }
            }
        }
    }
}


void    calculate_diffraction_forces_nlin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                HMFInterface*   hmf_interf,
                                                cusfloat        w,
                                                cuscomplex*     wave_diffrac
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
                                                                input->press_abs_err,
                                                                input->press_rel_err,
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

    // Allocate space for the sources, hydromechanics and wave exciting forces
    int         hydmech_np          = pow2s( input->dofs_np * mesh_gp->meshes_np );
    int         wave_exc_np         = input->heads_np * mesh_gp->meshes_np * input->dofs_np;

    cusfloat*   added_mass          = generate_empty_vector<cusfloat>( hydmech_np );
    cusfloat*   damping_rad         = generate_empty_vector<cusfloat>( hydmech_np );
    cuscomplex* froude_krylov       = generate_empty_vector<cuscomplex>( wave_exc_np );
    cuscomplex* sources             = generate_empty_vector<cuscomplex>( ( input->dofs_np + input->heads_np ) * scl.num_rows_local );
    cuscomplex* sysmat              = generate_empty_vector<cuscomplex>( scl.num_rows_local * scl.num_cols_local );
    cuscomplex* sysmat_steady       = generate_empty_vector<cuscomplex>( scl.num_rows_local * scl.num_cols_local );
    cuscomplex* wave_diffrac        = generate_empty_vector<cuscomplex>( wave_exc_np );

    cusfloat*   added_mass_p0       = generate_empty_vector<cusfloat>( hydmech_np );
    cusfloat*   damping_rad_p0      = generate_empty_vector<cusfloat>( hydmech_np );
    cuscomplex* froude_krylov_p0    = generate_empty_vector<cuscomplex>( wave_exc_np );
    cusfloat*   hydrostiff_p0       = generate_empty_vector<cusfloat>( hydmech_np );
    cuscomplex* raos_p0             = generate_empty_vector<cuscomplex>( wave_exc_np );
    cusfloat*   structural_mass_p0  = generate_empty_vector<cusfloat>( hydmech_np );
    cuscomplex* wave_diffrac_p0     = generate_empty_vector<cuscomplex>( wave_exc_np );
    cuscomplex* wave_exc_p0         = generate_empty_vector<cuscomplex>( wave_exc_np );

    // Define memory allocations for the constant source case ( Fast Mode )
    cuscomplex* inf_pot_mat         = nullptr;
    cuscomplex* inf_pot_steady_mat  = nullptr;
    int         ipm_cols_np         = 0;
    int         ipm_sc              = 0;
    int         ipm_ed              = 0;
    cuscomplex* panel_pot           = nullptr;
    cuscomplex* panel_pot_p0        = nullptr;
    if ( input->poly_order == 0 )
    {
        // Define column ranges in function of the number of processors
        mpi_config->get_1d_bounds( 
                                            mesh_gp->source_nodes_tnp, 
                                            ipm_sc, 
                                            ipm_ed 
                                        );
        ipm_cols_np         = ipm_ed - ipm_sc;

        // Allocate space for the influence potential matrix
        inf_pot_mat         = generate_empty_vector<cuscomplex>( mesh_gp->panels_tnp * ipm_cols_np );
        panel_pot           = generate_empty_vector<cuscomplex>( ( input->dofs_np + input->heads_np ) * mesh_gp->panels_tnp );
        panel_pot_p0        = generate_empty_vector<cuscomplex>( ( input->dofs_np + input->heads_np ) * mesh_gp->panels_tnp );

        // Calculate steady contribution of the influence potential matrix
        inf_pot_steady_mat  = generate_empty_vector<cuscomplex>( mesh_gp->panels_tnp * ipm_cols_np );
        calculate_influence_potential_steady(
                                                    input,
                                                    mpi_config,
                                                    mesh_gp,
                                                    inf_pot_steady_mat
                                            );
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
                                                                input->water_depth,
                                                                input->grav_acc,
                                                                input->press_abs_err,
                                                                input->press_rel_err
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
    calculate_sources_sysmat_steady(
                                        input,
                                        &scl,
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
        calculate_sources_intensity(
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


        if ( mpi_config->is_root( ) )
        {
            std::cout << "Sources" << std::endl;
            print_vector( 5, sources, 1, 6 );
        }

        // Update sources values for the integration objects
        hmf_interf->set_source_values( sources );

        // Calculate potential matrix if any
        if ( input->poly_order == 0 )
        {
            // Calculate potential influence coeffcients matrix
            calculate_influence_potential_total(
                                                    input,
                                                    mpi_config,
                                                    mesh_gp,
                                                    input->angfreqs[i],
                                                    inf_pot_steady_mat,
                                                    inf_pot_mat
                                                );

            // Calculate panels potential
            calculate_panel_potentials_lin(
                                                    input,
                                                    inf_pot_mat,
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
            calculate_panel_potentials_nlin(
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
    delete grf_dn_interf;
    delete gwf_dn_interf;
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
    mkl_free( sysmat_steady );
    mkl_free( wave_diffrac );

    if ( input->poly_order == 0 )
    {
        mkl_free( inf_pot_mat );
        mkl_free( inf_pot_steady_mat );
        mkl_free( panel_pot_p0 );
        mkl_free( panel_pot );
    }
    
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
                                                                input->press_abs_err,
                                                                input->press_rel_err,
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


void    calculate_hydromechanic_coeffs_lin( 
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     panels_pot,
                                                cusfloat        ang_freq,
                                                cusfloat*       added_mass,
                                                cusfloat*       damping_rad
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
                        // Integrate pressure over panel
                        index_1             = mesh_gp->panels_cnp[jb] + ie + mesh_gp->panels_tnp * id;
                        press_i             = ( 
                                                    panels_pot[index_1] 
                                                    * 
                                                    mesh_gp->panels[mesh_gp->panels_cnp[ib]+ie]->normal_vec[jd] 
                                                    * 
                                                    mesh_gp->panels[mesh_gp->panels_cnp[jb]+ie]->area
                                                );
                        added_mass[index]   -=  rho_w * press_i.imag( ) / ang_freq;
                        damping_rad[index]  -=  rho_w * press_i.real( );
                    }

                }
            }
        }
    }
}


void    calculate_hydromechanic_coeffs_nlin(
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
                                                                    input->press_abs_err,
                                                                    input->press_rel_err,
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


void    calculate_influence_potential_steady(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     inf_pot_mat
                                            )
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
    
    // Configure MPI
    int source_start_pos    = 0;
    int source_end_pos      = 0;
    mpi_config->get_1d_bounds( 
                                mesh_gp->source_nodes_tnp, 
                                source_start_pos, 
                                source_end_pos 
                            );
    
    // Generate potential matrix
    int         count = 0;
    cuscomplex  pot_steady_term( 0.0, 0.0 );
    for ( int i=0; i<mesh_gp->source_nodes_tnp; i++ )
    {
        // Change field point
        green_interf_steady->set_field_point( 
                                                mesh_gp->source_nodes[i]->panel->center
                                            );

        for ( int j=source_start_pos; j<source_end_pos; j++ )
        {
            // Change source point
            green_interf_steady->set_source(
                                                mesh_gp->source_nodes[j],
                                                1.0
                                            );
            
            // Compute steady and wave terms over the panel
            pot_steady_term = adaptive_quadrature_panel(
                                                            mesh_gp->source_nodes[j]->panel,
                                                            steady_fcn,
                                                            input->press_abs_err,
                                                            input->press_rel_err
                                                        );
            
            inf_pot_mat[count] = pot_steady_term / 4.0 / PI;
            count++;
        }
    }

}


void    calculate_influence_potential_total(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cusfloat        ang_freq,
                                                cuscomplex*     inf_pot_steady,
                                                cuscomplex*     inf_pot_total
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
    
    // Configure MPI
    int source_start_pos    = 0;
    int source_end_pos      = 0;
    mpi_config->get_1d_bounds( 
                                mesh_gp->source_nodes_tnp, 
                                source_start_pos, 
                                source_end_pos 
                            );
    
    // Generate potential matrix
    int         count = 0;
    cuscomplex  pot_wave_term( 0.0, 0.0 );
    for ( int i=0; i<mesh_gp->source_nodes_tnp; i++ )
    {
        // Change field point
        green_interf_wave->set_field_point(
                                                mesh_gp->source_nodes[i]->panel->center
                                            );

        for ( int j=source_start_pos; j<source_end_pos; j++ )
        {
            // Change source point
            green_interf_wave->set_source(
                                                mesh_gp->source_nodes[j],
                                                1.0
                                        );
            
            // Compute steady and wave terms over the panel
            pot_wave_term           = adaptive_quadrature_panel(
                                                                    mesh_gp->source_nodes[j]->panel,
                                                                    wave_fcn,
                                                                    input->press_abs_err,
                                                                    input->press_rel_err
                                                                );
            
            inf_pot_total[count]    = inf_pot_steady[count] + pot_wave_term / 4.0 / PI;
            count++;
        }
    }

}


void    calculate_panel_potentials_lin(
                                                Input*          input,
                                                cuscomplex*     inf_pot_mat,
                                                int             rows_np,
                                                int             cols_np,
                                                int             start_col,
                                                cuscomplex*     sources,
                                                cuscomplex*     panel_pot
                                        )
{
    // Loop over RHS to compute all the panel potentials the panels potentials
    cuscomplex  alpha( 1.0, 0.0 );
    cuscomplex  beta( 0.0, 0.0 );
    int         icnx = 1;
    int         icny = 1;
    for ( int i=0; i<( input->dofs_np + input->heads_np ); i++ )
    {
        cblas_gemv<cuscomplex>( 
                                    CblasRowMajor,
                                    CblasNoTrans,
                                    rows_np,
                                    cols_np,
                                    &alpha,
                                    inf_pot_mat,
                                    cols_np,
                                    &(sources[i*rows_np+start_col]),
                                    icnx,
                                    &beta,
                                    &(panel_pot[i*rows_np]),
                                    icny
                                );
    }
}


void    calculate_panel_potentials_nlin(
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
                                                            input->press_abs_err,
                                                            input->press_rel_err
                                                        );
            pot_i_wave      = adaptive_quadrature_panel(
                                                            mesh_gp->source_nodes[j]->panel,
                                                            wave_fcn,
                                                            input->press_abs_err,
                                                            input->press_rel_err
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
                                                GWFDnInterface* gwf_interf,
                                                cusfloat        w,
                                                cuscomplex*     sysmat_steady,
                                                cuscomplex*     sysmat,
                                                cuscomplex*     sources_int
                                   )
{
    /***************************************/
    /******** Fill system matrix  **********/
    /***************************************/
    // Create auxiliar lambda function
    auto wave_fcn   =   [gwf_interf]
                        ( 
                            cusfloat xi, 
                            cusfloat eta,
                            cusfloat x, 
                            cusfloat y, 
                            cusfloat z 
                        )
                        {
                            return (*gwf_interf)( xi, eta, x, y, z );
                        };
    
    // Loop over panels to integrate value
    int         col_count   = 0;
    int         index       = 0;
    cuscomplex  int_value( 0.0, 0.0 );
    cuscomplex  wave_value( 0.0, 0.0 );
    SourceNode* source_i    = nullptr;
    int         row_count   = 0;

    for ( int i=scl->start_col_0; i<scl->end_col_0; i++ )
    {
        // Get memory address of the ith panel
        source_i = mesh_gp->source_nodes[i];
        gwf_interf->set_source_i( source_i );

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
        {
            // Get memory address of the panel jth
            gwf_interf->set_source_j( mesh_gp->source_nodes[j] );

            // Integrate green function normal derivative along the current panel
            if ( i == j )
            {
                int_value   =   -std::complex( 0.5, 0.0 );
            }
            else
            {
                wave_value      =   adaptive_quadrature_panel(
                                                                source_i->panel,
                                                                wave_fcn,
                                                                input->gfdn_abs_err,
                                                                input->gfdn_rel_err,
                                                                false
                                                            );
                int_value       =   wave_value / 4.0 / PI;
            }
            index           = col_count*scl->num_rows_local+row_count;
            sysmat[index]   = sysmat_steady[index] + int_value;

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


void    calculate_sources_sysmat_steady(
                                                Input*          input,
                                                SclCmpx*        scl,
                                                MeshGroup*      mesh_gp,
                                                GRFDnInterface* grf_interf,
                                                cuscomplex*     sysmat
                                        )
{
    /***************************************/
    /******** Fill system matrix  **********/
    /***************************************/
    // Create auxiliar lambda function
    auto steady_fcn =   [grf_interf]
                        ( 
                            cusfloat xi, 
                            cusfloat eta,
                            cusfloat x, 
                            cusfloat y, 
                            cusfloat z 
                        )
                        {
                            return (*grf_interf)( xi, eta, x, y, z );
                        };
    
    // Loop over panels to integrate value
    int         col_count   = 0;
    cuscomplex  int_value( 0.0, 0.0 );
    cuscomplex  steady_value( 0.0, 0.0 );
    SourceNode* source_i    = nullptr;
    int         row_count   = 0;

    for ( int i=scl->start_col_0; i<scl->end_col_0; i++ )
    {
        // Get memory address of the ith panel
        source_i = mesh_gp->source_nodes[i];
        grf_interf->set_source_i( source_i );

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
        {
            // Get memory address of the panel jth
            grf_interf->set_source_j( mesh_gp->source_nodes[j] );

            // Integrate green function normal derivative along the current panel
            if ( i == j )
            {
                int_value       = complex( 0.0, 0.0 );
            }
            else
            {
                steady_value    =   adaptive_quadrature_panel(
                                                                source_i->panel,
                                                                steady_fcn,
                                                                input->gfdn_abs_err,
                                                                input->gfdn_rel_err,
                                                                false
                                                            );
                int_value       =   steady_value / 4.0 / PI;
            }

            sysmat[col_count*scl->num_rows_local+row_count] = int_value;

            // Advance row count
            row_count++;
        }

        // Advance column count
        col_count++;
    }
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