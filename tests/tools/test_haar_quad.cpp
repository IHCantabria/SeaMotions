
// Include general usage libraries
#include <functional>
#include <iostream>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/math/math_tools.hpp"


cusfloat quad_haar_1d(
                        std::function <cusfloat(cusfloat)> f_def,
                        cusfloat a,
                        cusfloat b,
                        int M
                        )
{
    // Calculate dependent parameters
    cusfloat dab = b-a;

    // Loop over interval to sum up contributions
    cusfloat sol = 0.0;
    for (int i=1; i<2*M+1; i++)
    {
        sol += f_def(a + dab*(i-0.5)/2/M);
    }
    sol *= dab/2/M;

    return sol;
}


int main(void)
{
    cusfloat int_value = quad_haar_1d(
                                        [](cusfloat x){return std::sin(x*x);},
                                        0.0,
                                        1.0,
                                        400
                                        );
    
    std::cout.precision(16);
    std::cout << "sin(x): " << int_value << std::endl;

    return 0;
}