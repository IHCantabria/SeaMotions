
// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "config.hpp"
#include "math_tools.hpp"


void cheby_roots_fun(int num_points, cusfloat* roots)
{
    for (int k=1; k<=num_points; k++)
    {
        roots[k-1] = -cos((2*k-1)*PI/(2*num_points));
    }
}