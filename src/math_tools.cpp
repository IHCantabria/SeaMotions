
// Include external libraries
#include "mkl.h"

// Include local libraries
#include "config.hpp"
#include "math_tools.hpp"


//////////////////////////////////////////////
/////// FUNCTION DEFINITION BLOCK ////////////
//////////////////////////////////////////////
signed long long factorial(int n)
{
    signed long long pk = 1;
    for (int i=1; i<=n; i++)
    {
        pk *= i;
    }

    return pk;
}

cusfloat* generate_empty_vector(int size)
{
    cusfloat* new_vector = (cusfloat*)mkl_calloc( size, sizeof(cusfloat), FLOATING_PRECISION );

    return new_vector;
}