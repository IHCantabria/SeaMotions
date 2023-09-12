
// Include local modules
#include "mpi_config.hpp"


bool MpiConfig::is_root( void )
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