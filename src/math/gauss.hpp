
#ifndef __gauss_hpp
#define __gauss_hpp

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
#include "../config.hpp"
#include "math_tools.hpp"


/**************************************/
/********* Declare functions***********/
/**************************************/
void get_gauss_legendre(
                            int         num_points, 
                            cusfloat*   roots, 
                            cusfloat*   weights
                        );


/**************************************/
/********* Declare classes ************/
/**************************************/
struct GaussPoints
{
    // Define local variables
    int         np          = 0;
    cusfloat*   roots       = nullptr;
    cusfloat*   weights     = nullptr;

    // Define constructors and destructor
    GaussPoints( int np_in )
    {
        // Storage input data
        this->np        = np_in;

        // Allocate space for roots and weights
        this->roots     = generate_empty_vector<cusfloat>( this->np );
        this->weights   = generate_empty_vector<cusfloat>( this->np );

        // Get the required number of gauss points
        get_gauss_legendre(
                                this->np,
                                this->roots,
                                this->weights
                            );
        
    }

    ~GaussPoints( void )
    {
        mkl_free( this->roots );
        mkl_free( this->weights );
    }
    
};

#endif