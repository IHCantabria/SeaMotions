
// Include general usage libraries
#include <iostream>

// Include local modules
#include "../config.hpp"
#include "shape_functions.hpp"


int     dofs_rectangular_region(
                                    int poly_order
                                )
{
    int dofs_total = 0;
    if ( poly_order == 0 )
    {
        dofs_total = 1;
    }
    else if ( poly_order > 0 )
    {
        dofs_total = ( poly_order + 1 ) * ( poly_order + 1 );
    }
    else
    {
        std::cerr << "ERROR - dofs_rectangular_region" << std::endl;
        std::cerr << "Polynomial order: " << poly_order << " not valid." << std::endl;
        throw std::runtime_error( "" );
    }

    return dofs_total;
}


int     dofs_triangular_region(
                                    int poly_order
                                )
{
    int dofs_total = 0;
    if ( poly_order == 0 )
    {
        dofs_total = 1;
    }
    else if ( poly_order > 0 )
    {
        for ( int i=1; i<poly_order+2; i++ )
        {
            dofs_total += i;
        }
    }
    else
    {
        std::cerr << "ERROR - dofs_rectangular_region" << std::endl;
        std::cerr << "Polynomial order: " << poly_order << " not valid." << std::endl;
        throw std::runtime_error( "" );
    }

    return dofs_total;
}


cusfloat shape_functions(
                            int p_order,
                            int q_order,
                            cusfloat xi,
                            cusfloat eta
                        )
{
    cusfloat val = 0.0;
    if ( p_order == 0 && q_order == 0 )
    {
        val = 1.0;
    }
    else
    {
        std::cerr << "ERROR - shape_functions" << std::endl;
        std::cerr << "Not implemented yet!" << std::endl;
        throw std::runtime_error( "" );
    }

    return val;
}