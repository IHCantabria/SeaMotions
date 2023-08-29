
#ifndef __mpi_config_hpp
#define __mpi_config_hpp

struct MpiConfig
{
public:
    // Define class attributes
    int proc_rank   = 0;
    int procs_total = 0;
    int proc_root   = 0;

    // Define class methods
    bool is_root( void )
    {
        return this->proc_root == this->proc_rank;
    }

};

#endif