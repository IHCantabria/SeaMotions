
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

    // // Group all meshes in a vector
    // Mesh** all_meshes = new Mesh*[input->bodies_np];
    // for ( int i=0; i<input->bodies_np; i++ )
    // {
    //     all_meshes[i] = input->bodies[i]->mesh;
    // }

    // // Create new mesh from the meshes of all objects

    Mesh* all_mesh = input->bodies[0]->mesh;

    all_mesh->define_source_nodes( 
                                    input->poly_order,
                                    input->bodies[0]->cog
                                );
    
    /****************************************************/
    /********* Create Scalapack solver instance *********/
    /****************************************************/
    SclCmpx scl( 
                    all_mesh->source_nodes_np,
                    input->dofs_np,
                    mpi_config->procs_total,
                    mpi_config->proc_rank
                );

    /****************************************************/
    /********* Create Green function interface *********/
    /****************************************************/
    GWFDnInterface* green_interf    = new GWFDnInterface(
                                                            all_mesh->source_nodes[0],
                                                            all_mesh->source_nodes[0],
                                                            input->angfreqs[0],
                                                            input->water_depth,
                                                            input->grav_acc
                                                        );

    /****************************************************/
    /******* Calculate hydrodynamic coefficients ********/
    /****************************************************/

    // Allocate space for intermediate data
    cusfloat*   added_mass  = generate_empty_vector<cusfloat>( input->dofs_np * input->dofs_np );
    cusfloat*   damping     = generate_empty_vector<cusfloat>( input->dofs_np * input->dofs_np );
    cuscomplex* sources     = generate_empty_vector<cuscomplex>( input->dofs_np * all_mesh->elems_np );
    cuscomplex* sysmat      = generate_empty_vector<cuscomplex>( scl.num_rows_local * scl.num_cols_local );

    cuscomplex*  all_sources;
    if ( mpi_config->is_root( ) )
    {
        all_sources = generate_empty_vector<cuscomplex>( all_mesh->elems_np );
    }

    // Loop over frequencies
    for ( int i=0; i<input->angfreqs_np; i++ )
    {
        std::cout << "Calculating angular frequency: " << input->angfreqs[i] << std::endl;
        // Recalculate wave properties for the current 
        // angular frequency
        green_interf->set_ang_freq( input->angfreqs[i] );

        // Calculate sources intensity
        MPI_Barrier( MPI_COMM_WORLD );
        calculate_sources_intensity( 
                                        &scl,
                                        all_mesh,
                                        green_interf,
                                        input->angfreqs[i],
                                        sysmat,
                                        sources
                                    );
        
        // Gather source values from each processor
        MPI_Gather(
                        sources,
                        scl.num_rows_local,
                        mpi_cuscomplex,
                        all_sources,
                        all_mesh->elems_np,
                        mpi_cuscomplex,
                        mpi_config->proc_root,
                        MPI_COMM_WORLD            
                    );

        // Calculate hydromechanic coefficients
        if ( mpi_config->is_root( ) )
        {
            // Calculate added mass and damping coefficients
            calculate_hydromechanic_coeffs( 
                                                all_mesh,
                                                green_interf,
                                                all_sources,
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
        }
        MPI_Barrier( MPI_COMM_WORLD );

    }

    // Delete heap memory allocated data
    delete green_interf;
    mkl_free( added_mass );
    mkl_free( damping );
    mkl_free( sources );
    mkl_free( sysmat );

    if ( mpi_config->is_root( ) )
    {
        mkl_free( all_sources );
    }

}


void    calculate_hydromechanic_coeffs(
                                            Mesh*           mesh,
                                            GWFDnInterface* green_interf,
                                            cuscomplex*     sources,
                                            cusfloat*       added_mass,
                                            cusfloat*       damping
                                        )
{
    // Loop over panels to calculate pressure over them

    for ( int i=0; i<6; i++ )
    {

    }
}


void    calculate_sources_intensity(
                                        SclCmpx*        scl,
                                        Mesh*           mesh,
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
        std::cout << "Source[i]: " << i << " - " << mesh->source_nodes[i] << std::endl;
        // Get memory address of the ith panel
        source_i = mesh->source_nodes[i];
        green_interf->set_source_i( source_i );

        // Loop over rows to calcualte the influence of the panel
        // over each collocation point
        row_count = 0;
        for ( int j=scl->start_row_0; j<scl->end_row_0; j++ )
        {
            // std::cout << "Panel[j]: " << j << " - " << mesh->source_nodes[j] << std::endl;
            // Get memory address of the panel jth
            green_interf->set_source_j( mesh->source_nodes[j] );

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

    // Fill RHS vector for Surge motion
    count       = 0;
    m           = 0;
    start_pos   = m * scl->num_rows_local;
    std::cout << "RHS" << std::endl;
    for ( int i=scl->start_row_0; i<scl->end_row_0; i++ )
    {
        sources_int[start_pos+count] = w * mesh->source_nodes[i]->normal_vec[0];
        // outfile << sources_int[start_pos+count].real( ) << " " << sources_int[start_pos+count].imag( ) << std::endl;
        count++;
    }

    // outfile.close( );

    // Fill RHS vector for Sway motion
    count       = 0;
    m           = 1;
    start_pos   = m * mesh->elems_np;
    for ( int i=scl->start_row_0; i<scl->end_row_0; i++ )
    {
        sources_int[start_pos+count] = w * mesh->source_nodes[i]->normal_vec[1];
        count++;
    }

    // Fill RHS vector for Heave motion
    count       = 0;
    m           = 2;
    start_pos   = m * mesh->elems_np;
    for ( int i=scl->start_row_0; i<scl->end_row_0; i++ )
    {
        sources_int[start_pos+count] = w * mesh->source_nodes[i]->normal_vec[2];
        count++;
    }

    // Fill RHS vector for Roll motion
    count       = 0;
    m           = 3;
    start_pos   = m * mesh->elems_np;
    for ( int i=scl->start_row_0; i<scl->end_row_0; i++ )
    {
        // Get current panel memory location
        source_i = mesh->source_nodes[i];

        // Calculate proyection of the normal velocity
        sources_int[start_pos+count] = w * source_i->normal_vec[3];

        count++;
    }

    // Fill RHS vector for Pitch motion
    count       = 0;
    m           = 4;
    start_pos   = m * mesh->elems_np;
    for ( int i=scl->start_row_0; i<scl->end_row_0; i++ )
    {
        // Get current panel memory location
        source_i = mesh->source_nodes[i];

        // Calculate proyection of the normal velocity
        sources_int[start_pos+count] = w * source_i->normal_vec[4];

        count++;
    }

    // Fill RHS vector for Yaw motion
    count       = 0;
    m           = 5;
    start_pos   = m * mesh->elems_np;
    for ( int i=scl->start_row_0; i<scl->end_row_0; i++ )
    {
        // Get current panel memory location
        source_i = mesh->source_nodes[i];

        // Calculate proyection of the normal velocity
        sources_int[start_pos+count] = w * source_i->normal_vec[5];

        count++;
    }

    std::cout << "--> Done!" << std::endl; 

    // Solve system of equations
    std::cout << "Calculating system of equations..." << std::endl;
    scl->Solve( sysmat, sources_int );
    std::cout << "--> Done!" << std::endl;
    
}