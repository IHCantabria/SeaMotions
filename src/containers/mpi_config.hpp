
#ifndef __mpi_config_hpp
#define __mpi_config_hpp

struct MpiConfig
{
public:
    int proc_rank   = 0;
    int procs_total = 0;
    int proc_root   = 0;

};

#endif