
// Include general usage libraries
#include <iostream>
#include "mpi.h"


int main( void )
{
    // Initializa MPI environment
    MPI_Init( NULL, NULL );

    // Get number total number of processes
    int world_size = 0;
    MPI_Comm_size( 
                    MPI_COMM_WORLD,
                    &world_size
                );

    // Get current process rank
    int proc_rank = 0;
    MPI_Comm_rank(
                    MPI_COMM_WORLD,
                    &proc_rank
                );

    // Master process
    int number = 0;
    if ( proc_rank == 0 )
    {
        number = -999;
        MPI_Send(
                    &number,
                    1,
                    MPI_INT,
                    1,
                    0,
                    MPI_COMM_WORLD
                );
    }
    else if ( proc_rank == 1 )
    {
        MPI_Recv(
                    &number,
                    1,
                    MPI_INT,
                    0,
                    0,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                );
        std::cout << "Message received at process 1. Number value: " << number << std::endl;
    }

    // Close MPI environment
    MPI_Finalize( );

    return 0;
}