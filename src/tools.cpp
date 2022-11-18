
// Include external libraries
#include "mkl.h"

// Include local libraries
#include "tools.hpp"


double* generate_empty_vector(int size)
{
    double* new_vector = (double*)mkl_malloc( sizeof(double)*size, 64 );
    for (int i=0; i<size; i++)
    {
        new_vector[i] = 0.0;
    }

    return new_vector;
}