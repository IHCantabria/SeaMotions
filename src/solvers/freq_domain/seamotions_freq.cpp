
// Include general usage libraries
#include <iostream>
#include "mpi.h"
#include <string>

// Include local modules
#include "../../containers/mpi_timer.hpp"
#include "frequency_solver.hpp"
#include "../../cli_header_banner.hpp"
#include "../../inout/reader.hpp"
#include "../../tools.hpp"
#include "../../version.hpp"



int main( int argc, char* argv[] )
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
    cli_header_banner<true>( case_fopath, "Frequency" );


    /*****************************************/
    /************ Read Input data ************/
    /*****************************************/
    Input input;
    read_input_files( &input, case_fopath );
    

    /*****************************************/
    /****** Initialize Frequency Solver ******/
    /*****************************************/
    MpiTimer case_timer;

    FrequencySolver<NUM_GP, PF_ON> freq_solver( &input, &mpi_config );
    
    // /*****************************************/
    // /******** First Order Solution ***********/
    // /*****************************************/
    freq_solver.calculate_first_order( );

    // /*****************************************/
    // /********* Close MPI environment *********/
    // /*****************************************/
    case_timer.stop( );
    MPI_Finalize( );

    /*****************************************/
    /**** Close program and final actions ****/
    /*****************************************/

    // Print Elapsed time
    if ( mpi_config.is_root( ) )
    {
        std::cout << std::endl << std::endl;
        std::cout << " -> Seamotions (Frequency) finished!" << std::endl;
        std::cout << " ---> Elapsed wall time for calculation [s]: " << case_timer << std::endl;
    }

    return 0;
}