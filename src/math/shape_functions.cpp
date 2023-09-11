
// Include local modules
#include "shape_functions.hpp"


int     dofs_rectangular_region(
                                    int poly_order
                                )
{
    return poly_order * poly_order;
}


int     dofs_triangular_region(
                                    int poly_order
                                )
{
    int dofs_total = 0;
    for ( int i=1; i<poly_order+1; i++ )
    {
        dofs_total += i;
    }

    return dofs_total;
}