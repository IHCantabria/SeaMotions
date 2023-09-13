
// Include general usage libraries
#include <iostream>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "mpi_config.hpp"


void    MpiConfig::get_1d_bounds(
                                    int     np,
                                    int&    start_pos,
                                    int&    end_pos
                                )
{
    // Check if there is more processors than data
    if ( this->procs_total > np )
    {
        std::cerr << "There is more processes available than ";
        std::cerr << "data to distribute over them." << std::endl;
        throw std::runtime_error( "" );
    }

    // Divide data in chunks
    int chunk_size = 1 + ( ( np - 1 ) / this->procs_total );

    // Set interval bounds
    start_pos   = this->proc_rank * chunk_size;
    end_pos     = ( this->proc_rank + 1 ) * chunk_size;

    // Set limits to the upper bound
    end_pos     = ( end_pos > np ) ? np : end_pos;
}


bool    MpiConfig::is_root( void )
{
    return this->proc_root == this->proc_rank;
}


MpiConfig::MpiConfig(
                        int proc_rank_in,
                        int procs_total_in,
                        int proc_root_in,
                        MPI_Comm mpi_comm_in
                    )
{
    this->mpi_comm      = mpi_comm_in;
    this->proc_rank     = proc_rank_in;
    this->procs_total   = procs_total_in;
    this->proc_root     = proc_root_in;
}