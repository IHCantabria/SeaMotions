
#ifndef __mpi_config_hpp
#define __mpi_config_hpp

// Include general usage libraries
#include "mpi.h"


struct MpiConfig
{
private:
    // Define class private attributes
    bool        _is_parallel    = true;

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
    void    get_1d_bounds(
                            int     np,
                            int&    start_pos,
                            int&    end_pos
                        );

    bool    is_root(    
                            void 
                    );

    void    set_parallel(
                            void
                        );

    void    set_serial( 
                            void 
                        );

};

#endif