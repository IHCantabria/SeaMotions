
// Include general usage scientific libraries
#include <complex>
#include <fstream>

// Include local modules
#include "freq_solver_tools.hpp"
#include "../../math/integration.hpp"


void calculate_freq_domain_coeffs(
                                    MpiConfig*      mpi_config,
                                    Input*          input,
                                    Hydrostatics*   hydrostatics,
                                    Output*         output
                                )
{
    // /****************************************************/
    // /*** Create total mesh for the interacting bodies ***/
    // /****************************************************/

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
                    input->dofs_np,
                    mpi_config->procs_total,
                    mpi_config->proc_rank
                );

    /****************************************************/
    /****** Allocate space for the simulation data ******/
    /****************************************************/
    cusfloat*   added_mass  = generate_empty_vector<cusfloat>( pow2s( input->dofs_np * mesh_gp->meshes_np ) );
    cuscomplex* all_sources = generate_empty_vector<cuscomplex>( input->dofs_np * mesh_gp->source_nodes_tnp );
    cusfloat*   damping     = generate_empty_vector<cusfloat>( pow2s( input->dofs_np * mesh_gp->meshes_np ) );
    cuscomplex* sources     = generate_empty_vector<cuscomplex>( input->dofs_np * mesh_gp->meshes_np * scl.num_rows_local );
    cuscomplex* sysmat      = generate_empty_vector<cuscomplex>( scl.num_rows_local * scl.num_cols_local );

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
                                                                mesh_gp->source_nodes_tnp,
                                                                mesh_gp->panels[0],
                                                                0,
                                                                0,
                                                                input->angfreqs[0],
                                                                input->water_depth,
                                                                input->grav_acc
                                                        );

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
                        scl.num_rows_local,
                        mpi_cuscomplex,
                        all_sources,
                        scl.num_rows,
                        mpi_cuscomplex,
                        MPI_COMM_WORLD
                    );

        // Update sources values for the integration objects
        hmf_interf->set_source_values( all_sources );

        // Calculate added mass and damping coefficients
        calculate_hydromechanic_coeffs( 
                                            mpi_config,
                                            mesh_gp,
                                            input->dofs_np,
                                            hmf_interf,
                                            added_mass,
                                            damping
                                        );

        // Loop over headings to calculate diffraction and 
        // Froude-Krylov forces
        for ( int j=0; j<input->heads_np; j++ )
        {
            // Calculate wave exciting forces

            // Calculate raos
        }
        MPI_Barrier( MPI_COMM_WORLD );

    }

    // Delete heap memory allocated data
    
    mkl_free( added_mass );
    mkl_free( all_sources );
    mkl_free( damping );
    mkl_free( sources );
    mkl_free( sysmat );

    delete green_dn_interf;
    delete hmf_interf;
    delete [] all_meshes;
}


void    calculate_hydromechanic_coeffs(
                                            MpiConfig*      mpi_config,
                                            MeshGroup*      mesh_gp,
                                            int             dofs_np,
                                            HMFInterface*   hmf_interf,
                                            cusfloat*       added_mass,
                                            cusfloat*       damping
                                        )
{
    // Calculate MPI data chunks
    int elem_end_pos     = 0;
    int elem_start_pos   = 0;
    mpi_config->get_1d_bounds( 
                                    mesh_gp->panels_tnp, 
                                    elem_start_pos, 
                                    elem_end_pos 
                                );

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

    // Loop over first dimension of degrees of freedrom
    GaussPoints gp( 1 );
    cuscomplex  pressure = 0.0;
    for ( int ib=0; ib<mesh_gp->meshes_np; ib++ )
    {
        for ( int jb=0; jb<mesh_gp->meshes_np; jb++ )
        {
            for ( int id=0; id<dofs_np; id++ )
            {
                // Set ith dof
                hmf_interf->set_start_index_i( mesh_gp->source_nodes_cnp[ib]*dofs_np + id*mesh_gp->source_nodes_np[ib] );

                // Loop over second dimension of degrees of freedrom
                for ( int jd=0; jd<dofs_np; jd++ )
                {
                    // Set jth dof
                    hmf_interf->set_dof_j( jd );

                    // Loop over panels to integrate the wave radiation
                    // pressure along the floating object external shape
                    for ( int ie=elem_start_pos; ie<elem_end_pos; ie++ )
                    {
                        // Set new panel
                        hmf_interf->set_panel( mesh_gp->panels[ie] );

                        // Integrate pressure over panel
                        pressure += adaptive_quadrature_panel(
                                                                mesh_gp->panels[ie],
                                                                target_fcn,
                                                                1000.0,
                                                                &gp
                                                            );
                    }
                }
            }
        }
    }
}


void    calculate_sources_intensity(
                                        SclCmpx*        scl,
                                        MeshGroup*      mesh_gp,
                                        GWFDnInterface* green_interf,
                                        cusfloat        w,
                                        cuscomplex*     sysmat,
                                        cuscomplex*     sources_int
                                   )
{
    GaussPoints gp = GaussPoints( 1 );

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
        std::cout << "Source[i]: " << i << " - " << mesh_gp->source_nodes[i] << std::endl;
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
                                                                1000,
                                                                &gp
                                                            );
                int_value   =   int_value / 4.0 / PI;
            }
            
            sysmat[col_count*scl->num_rows_local+row_count] = int_value;
            // outfile << int_value.real( ) << " " << int_value.imag( ) << std::endl;
            // std::cout << "Index Count: " << col_count*scl->num_rows_local+row_count << std::endl;
            count_total += 1;
            // Advance row count
            row_count++;
        }
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
    /************* Fill RHS  ***************/
    /***************************************/
    std::cout << "Calculating RHS..." << std::endl;
    // Declare local variables to be used
    int         count           = 0;
    int         m               = 0; 
    int         start_pos       = 0; 

    // Fill RHS vector
    count       = 0;
    for ( int i=0; i< 6; i++ )
    {
        count = 0;
        for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
        {
            sources_int[start_pos+count] = w * mesh_gp->source_nodes[j]->normal_vec[i];
            // outfile << sources_int[start_pos+count].real( ) << " " << sources_int[start_pos+count].imag( ) << std::endl;
            count++;
        }
    }
    // outfile.close( );

    std::cout << "--> Done!" << std::endl; 

    // Solve system of equations
    std::cout << "Calculating system of equations..." << std::endl;
    scl->Solve( sysmat, sources_int );
    std::cout << "--> Done!" << std::endl;
    
}