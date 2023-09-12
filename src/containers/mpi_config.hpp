
#ifndef __mpi_config_hpp
#define __mpi_config_hpp

// Include general usage libraries
#include "mpi.h"


struct MpiConfig
{
public:
    // Define class attributes
    MPI_Comm    mpi_comm;
    int         proc_rank   = 0;
    int         procs_total = 0;
    int         proc_root   = 0;

    // Define class constructors and destructor
    MpiConfig(
                    int proc_rank_in,
                    int procs_total_in,
                    int proc_root_in,
                    MPI_Comm mpi_comm_in
                );

    // Define class methods
    bool is_root( void );

};

#endif