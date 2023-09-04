
// Include general usage libraries
#include "mpi.h"
#include <string>

// Include local modules
#include "./containers/mpi_config.hpp"
#include "hydrostatics.hpp"
#include "./inout/reader.hpp"
#include "tools.hpp"


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

    // Declare container to hold MPI configuration
    MpiConfig mpi_config;

    // Get total number of processors
    MPI_Comm_size(
                    MPI_COMM_WORLD,
                    &(mpi_config.procs_total)
                );

    // Get current process rank
    MPI_Comm_rank(
                    MPI_COMM_WORLD,
                    &(mpi_config.proc_rank)
                );

    /*****************************************/
    /**** Launch Hydrostatics Calculation ****/
    /*****************************************/
    MPI_Barrier( MPI_COMM_WORLD );
    double start = MPI_Wtime( );
    int body_num = 0;
    Hydrostatics hydrostatics( 
                                    input->bodies[body_num]->mesh,
                                    input->water_density,
                                    input->grav_acc,
                                    input->bodies[body_num]->mass,
                                    input->bodies[body_num]->cog,
                                    input->bodies[body_num]->rad_inertia,
                                    &mpi_config
                                );
    MPI_Barrier( MPI_COMM_WORLD );
    double end = MPI_Wtime( );

    if ( mpi_config.is_root( ) )
    {
        hydrostatics.print( );
        std::cout << "Execution time [s]: " << (end - start) << std::endl;
    }

    /*****************************************/
    /********* Launch Output System **********/
    /*****************************************/


    /*****************************************/
    /***** Calculate Source Distribution *****/
    /*****************************************/

    

    /*****************************************/
    /********* Close MPI environment *********/
    /*****************************************/
    MPI_Finalize( );

    return 0;
}