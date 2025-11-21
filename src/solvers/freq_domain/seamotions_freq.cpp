
// Include general usage libraries
#include "mpi.h"
#include <string>

// Include local modules
#include "../../containers/mpi_timer.hpp"
#include "freq_solver_tools.hpp"
#include "frequency_solver.hpp"
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
    MPI_Init( NULL, NULL );

    /*****************************************/
    /****** Initialize Frequency Solver ******/
    /*****************************************/
    MpiTimer case_timer;

    FrequencySolver<NUM_GP, PF_ON> freq_solver( input );
    
    /*****************************************/
    /******** First Order Solution ***********/
    /*****************************************/
    freq_solver.calculate_first_order( );

    /*****************************************/
    /********* Close MPI environment *********/
    /*****************************************/
    case_timer.stop( );
    MPI_Finalize( );

    /*****************************************/
    /**** Close program and final actions ****/
    /*****************************************/

    // Print Elapsed time
    if ( freq_solver.mpi_config->is_root( ) )
    {
        std::cout << "Elapsed wall time for calculation [s]: " << case_timer << std::endl;
        std::cout << "Seamotions Freqcuency ended!" << std::endl;
    }

    return 0;
}