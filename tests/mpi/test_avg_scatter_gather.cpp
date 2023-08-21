
// Include general usage libraries
#include <iostream>
#include "mpi.h"

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
#include "../../src/config.hpp"
#include "../../src/math/math_tools.hpp"


void ceil_div_int( int a, int b, int& c)
{
    c = ( a + b - 1 ) / b;
}


int main( void )
{
    // int np = 100;
    // int proc_total = 7;
    // int ipp = np / 7;
    // int ippf = np - (proc_total-1)*ipp;

    // std::cout << "np: " << np << " - proc_total: " << proc_total << " - ipp: " << ipp << " - ippf: " << ippf << std::endl;

    // Initialize MPI environment
    MPI_Init( NULL, NULL );

    // Get total number of processors
    int proc_total = 0;
    MPI_Comm_size(
                    MPI_COMM_WORLD,
                    &proc_total
                );

    // Get current process rank
    int proc_rank = 0;
    MPI_Comm_rank(
                    MPI_COMM_WORLD,
                    &proc_rank
                );

    // Generate data to average across all processes
    int         np      = 100;
    cusfloat*   data    = nullptr;
    if ( proc_rank == 0 )
    {
        data = generate_empty_vector<cusfloat>( np );
        for ( int i=0; i<np; i++ )
        {
            data[i] = i;
        }
    }

    // Calculate data to receive by each processor
    int send_np = np / proc_total;
    int recv_np = 0;
    if ( proc_rank < ( proc_total - 1 ) )
    {
        recv_np = send_np;
    }
    else
    {
        recv_np = np - ( proc_total - 1 )*send_np;
    }

    // Allocate space for the receiving buffer
    cusfloat* local_data = generate_empty_vector<cusfloat>( recv_np );

    // Scatter data
    MPI_Scatter(
                    data,
                    send_np,
                    mpi_cusfloat,
                    local_data,
                    recv_np,
                    mpi_cusfloat,
                    0,
                    MPI_COMM_WORLD
                );

    // Sum the local vector data
    cusfloat sum = 0.0;
    for ( int i=0; i<recv_np; i++ )
    {
        sum += local_data[i];
    }

    // Receive data from all processes
    cusfloat* gather_data = nullptr;
    if ( proc_rank == 0 )
    {
        gather_data = generate_empty_vector<cusfloat>( proc_total );
    }

    MPI_Gather(
                    &sum,
                    1,
                    mpi_cusfloat,
                    gather_data,
                    1,
                    mpi_cusfloat,
                    0,
                    MPI_COMM_WORLD
                );
    
    if ( proc_rank == 0 )
    {
        std::cout << "Gathered data: "; print_vector( proc_total, gather_data, 0, 6 );

        cusfloat gather_sum = 0.0;
        for ( int i=0; i<proc_total; i++ )
        {
            gather_sum += gather_data[i];
        }

        std::cout << "Data mean value: " << gather_sum / (cusfloat)(np) << std::endl;
    }

    // Deallocate heap memory
    if ( proc_rank == 0 )
    {
        mkl_free( data );
        mkl_free( gather_data );
    }
    mkl_free( local_data );

    // Close MPI environment
    MPI_Finalize( );

    return 0;
}