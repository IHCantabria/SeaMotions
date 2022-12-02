
// Include external libraries
#include "mkl.h"

// Include local libraries
#include "config.hpp"
#include "math_tools.hpp"


//////////////////////////////////////////////
////// MATHEMATICAL CONSTANTS BLOCK //////////
//////////////////////////////////////////////
cusfloat PI = 3.141592653589793;


//////////////////////////////////////////////
/////// FUNCTION DEFINITION BLOCK ////////////
//////////////////////////////////////////////
cusfloat* generate_empty_vector(int size)
{
    cusfloat* new_vector = (cusfloat*)mkl_calloc( size, sizeof(cusfloat), FLOATING_PRECISION );

    return new_vector;
}