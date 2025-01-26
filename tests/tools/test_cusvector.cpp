
// Include general usage modules
#include <iostream>
#include <vector>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/math/cusvector.hpp"
#include "../../src/math/math_tools.hpp"


int main( void )
{
    // Define test properties
    int                     N       = 3;
    cusfloat                eps     = 1e-3;

    // Create Vector containers
    CusVector<cusfloat>     my_vector( N );
    std::vector<cusfloat>   vec( N );

    // Fill in vector containers
    for ( int i=0; i<N; i++ )
    {
        my_vector[i]    = i;
        vec[i]          = i;
    }

    for ( int i=0; i<N; i++ )
    {
        my_vector[i]    += i;
        vec[i]          += i;
    }

    // Check if the vectors are equal
    bool pass = true;
    for ( int i=0; i<N; i++ )
    {
        if ( std::abs( my_vector[i] - vec[i] ) > eps )
        {
            pass = false;
            break;
        }

    }
    
    if ( !pass )
    {
        std::cout << "CusVector and std::vector are not equal!" << std::endl;
        return 1;
    }

    return 0;
}