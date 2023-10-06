
// Include general usage libraries
#include "mpi.h"
#include <string>

// Include local modules
#include "../../containers/mpi_config.hpp"
#include "freq_solver_tools.hpp"
#include "../../hydrostatics.hpp"
#include "../../inout/output.hpp"
#include "../../inout/reader.hpp"
#include "../../tools.hpp"



int main( int argc, char* argv[] )
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 1))
    {
        return 1;
    }

    std::string case_fopath(argv[1]);

    /*****************************************/
    /************ Read Input data ************/
    /*****************************************/
    Input* input = read_input_files( case_fopath );

    /*****************************************/
    /****** Initialize MPI environment *******/
    /*****************************************/

    // Init MPI
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

    // Declare container to hold MPI configuration
    MpiConfig mpi_config( 
                            proc_rank,
                            procs_total,
                            0,
                            MPI_COMM_WORLD
                        );


    /*****************************************/
    /**** Launch Hydrostatics Calculation ****/
    /*****************************************/
    MPI_Barrier( MPI_COMM_WORLD );
    double start = MPI_Wtime( );
    Hydrostatics** hydrostatics = new Hydrostatics*[input->bodies_np];
    for ( int i=0; i<input->bodies_np; i++ )
    {
        hydrostatics[i] =   new Hydrostatics( 
                                                input->bodies[i]->mesh,
                                                input->water_density,
                                                input->grav_acc,
                                                input->bodies[i]->mass,
                                                input->bodies[i]->cog,
                                                input->bodies[i]->rad_inertia,
                                                &mpi_config
                                            );
    }
    MPI_Barrier( MPI_COMM_WORLD );
    double end = MPI_Wtime( );

    if ( mpi_config.is_root( ) )
    {
        for ( int i=0 ; i<input->bodies_np; i++ )
        {
            hydrostatics[i]->print( );
        }
        std::cout << "Execution time [s]: " << (end - start) << std::endl;
    }

    /*****************************************/
    /********* Launch Output System **********/
    /*****************************************/
    Output* output = nullptr;
    
    if ( mpi_config.is_root( ) )
    {
        output = new Output( input );
    }

    /*****************************************/
    /****** Storage Initial Parameters *******/
    /*****************************************/

    if ( mpi_config.is_root( ) )
    {
        // Storage frequency set
        output->save_frequencies( input->freqs );

        // Storage headings set
        output->save_headings( input->heads.data( ) );

        // Storage structural mass
        if ( input->out_struct_mass )
        {
            output->save_structural_mass( );
        }

        // Storage hydrostatic stiffness matrix
        if ( input->out_hydstiff )
        {
            output->save_hydstiffness( hydrostatics );
        }

    }

    /*****************************************/
    /***** Calculate Source Distribution *****/
    /*****************************************/
    calculate_freq_domain_coeffs(
                                    &mpi_config,
                                    input, 
                                    hydrostatics,
                                    output 
                                );

    /*****************************************/
    /********* Close MPI environment *********/
    /*****************************************/
    MPI_Finalize( );

    /*****************************************/
    /**** Delete heap memory allocations *****/
    /*****************************************/
    delete input;
    delete output;

    for ( int i=0; i<input->bodies_np; i++ )
    {
        delete hydrostatics[i];
    }
    delete [] hydrostatics;

    std::cout << "Seamotions Freqcuency ended!" << std::endl;

    return 0;
}