
// Include general usage libraries
#include "mpi.h"

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
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

    // Generate data to send
    int         np          = 4;
    cusfloat*   send_data   = nullptr;
    if ( proc_rank == 0 )
    {
        send_data = generate_empty_vector<cusfloat>( np );
        for ( int i=0; i<np; i++ )
        {
            send_data[i] = i;
        }
    }

    // Generate receiving data buffer
    int         recv_np     = 1;
    cusfloat*   recv_data   = generate_empty_vector<cusfloat>( recv_np );

    // Scatter data
    MPI_Scatter(
                    send_data,
                    1,
                    mpi_cusfloat,
                    recv_data,
                    recv_np,
                    mpi_cusfloat,
                    0,
                    MPI_COMM_WORLD
                );

    std::cout << "Processor rank: " << proc_rank << " - Received data: " << recv_data[0] << std::endl;

    // Process data at each processor
    recv_data[0] *= 3;

    // Gather information in the root process
    cusfloat* root_data = nullptr;
    if ( proc_rank == 0 )
    {
        root_data = generate_empty_vector<cusfloat>( np );
    }

    MPI_Gather(
                    recv_data,
                    recv_np,
                    mpi_cusfloat,
                    root_data,
                    recv_np,
                    mpi_cusfloat,
                    0,
                    MPI_COMM_WORLD
                );

    // Print out processed information at root process
    if( proc_rank == 0 )
    {
        std::cout << "Received data at root process: ";
        print_vector( np, root_data, 0, 1 );
    }

    // Delete local heap memory
    if ( proc_rank == 0 )
    {
        mkl_free( send_data );
        mkl_free( root_data );
    }
    mkl_free( recv_data );

    // Close MPI environment
    MPI_Finalize( );

    return 0;
}