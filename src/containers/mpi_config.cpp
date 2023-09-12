
// Include local modules
#include "mpi_config.hpp"


bool MpiConfig::is_root( void )
{
    return this->proc_root == this->proc_rank;
}