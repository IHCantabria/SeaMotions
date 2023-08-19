
// Include general usage libraries
#include <iostream>
#include "mpi.h"

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
#include "../../src/config.hpp"
#include "../../src/math/math_tools.hpp"


int main( void )
{
    // Initialize MPI environment
    MPI_Init( NULL, NULL );

    // Get total number of processors
    int proc_total = 0;
    MPI_Comm_size(
                    MPI_COMM_WORLD,
                    &proc_total
                );

    // Get the rank of the current processor
    int proc_rank = 0;
    MPI_Comm_rank(
                    MPI_COMM_WORLD,
                    &proc_rank
                );
    
    // Send message from processor 0
    if ( proc_rank == 0 )
    {
        // Define information to send
        int         np      = 10;
        cusfloat*   data    = generate_empty_vector<cusfloat>( np );

        for ( int i=0; i<np; i++ )
        {
            data[i] = i;
        }

        // Send information
        MPI_Send(
                    data,
                    np,
                    mpi_cusfloat,
                    1,
                    0,
                    MPI_COMM_WORLD
                );

        // Delete heap memory allocate in this code block
        mkl_free( data );
    }

    // Receive message at processor 1
    if ( proc_rank == 1 )
    {
        // Get information about the size of the matrix
        MPI_Status status;
        MPI_Probe(
                        0,
                        0,
                        MPI_COMM_WORLD,
                        &status
                    );
        
        // Get data size from status object
        int data_size = 0;
        MPI_Get_count(
                        &status,
                        mpi_cusfloat,
                        &data_size
                    );

        // Print the size of the matrix to receive
        std::cout << "Size of the vector to receive: " << data_size << std::endl;

        // Allocate space for the matrix to receive
        cusfloat* recv_data = generate_empty_vector<cusfloat>( data_size );

        // Receive data
        MPI_Recv(
                    recv_data,
                    data_size,
                    mpi_cusfloat,
                    0,
                    0,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                );

        // Print-out data to screen
        std::cout << "Data received: "; print_vector( data_size, recv_data, 0, 1 );

        // Deallocate recv_data object from heap memory
        mkl_free( recv_data );

    }

    // Finalize MPI environment
    MPI_Finalize( );

    return 0;
}