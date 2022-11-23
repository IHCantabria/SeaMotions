
// Include external libraries
#include "mkl.h"

// Include local libraries
#include "config.hpp"
#include "tools.hpp"


cusfloat* generate_empty_vector(int size)
{
    cusfloat* new_vector = (cusfloat*)mkl_calloc( size, sizeof(cusfloat), FLOATING_PRECISION );

    return new_vector;
}