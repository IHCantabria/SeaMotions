
// Include external libraries
#include "mkl.h"

// Include local libraries
#include "config.hpp"
#include "tools.hpp"


cusfloat* generate_empty_vector(int size)
{
    cusfloat* new_vector = (cusfloat*)mkl_malloc( sizeof(cusfloat)*size, 64 );
    for (int i=0; i<size; i++)
    {
        new_vector[i] = 0.0;
    }

    return new_vector;
}