
// Include general usage libraries
#include <iostream>
#include "mpi.h"


int main( void )
{
    // Initialize MPI environment
    MPI_Init( NULL, NULL );

    // Get the total number of processors
    int proc_total = 0;
    MPI_Comm_size(
                    MPI_COMM_WORLD,
                    &proc_total
                );

    // Get the current rank of the processor
    int proc_rank = 0;
    MPI_Comm_rank(
                    MPI_COMM_WORLD,
                    &proc_rank
                );
    
    if ( proc_rank == 0 )
    {
        int data = -999;
        MPI_Bcast(
                    &data,
                    1,
                    MPI_INT,
                    0,
                    MPI_COMM_WORLD
                );
        std::cout << "Data sent from master process" << std::endl;
    }

    if ( proc_rank > 0 )
    {
        // Receive data from master process
        int recv_data = 0;
        // sMPI_Recv(
        //             &recv_data,
        //             1,
        //             MPI_INT,
        //             0,
        //             0,
        //             MPI_COMM_WORLD,
        //             MPI_STATUS_IGNORE
        //         );

        std::cout << "Message received from master process at processor: " << proc_rank << " - data: " << recv_data << std::endl;
    }


    // Close MPI environment
    MPI_Finalize( );

    return 0;
}