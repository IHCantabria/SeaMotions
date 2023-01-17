
// Include external scientific libraries
#include <cmath>
#include "mkl.h"

// Include local libraries
#include "../config.hpp"
#include "math_tools.hpp"


//////////////////////////////////////////////
/////// FUNCTION DEFINITION BLOCK ////////////
//////////////////////////////////////////////
cusfloat bisection(std::function<cusfloat(cusfloat)> f_def, cusfloat a, cusfloat b, 
                    cusfloat abs_prec, int max_iter, bool verbose)
{
    /**
     * @brief Bisection method to solve non-linear monoparametric equations.
     * 
     * The bisection method solves non-linear equations dependent on one parameter. It is an 
     * iterative solver based on the interval sub-division. More information can be found 
     * at: Numerical Recipes The Art of Scientific Computing. W.Press S.Teukolsky
     */

    // Start the iterative method
    cusfloat fa = f_def(a);
    cusfloat c;
    cusfloat fc;

    // Start iterative loop
    int count_iter = 0;
    while (true)
    {
        // Calculate new value
        c = (a+b)/2.0;
        fc = f_def(c);

        // Print iterative process
        if (verbose)
        {
            std::cout << "Iter: " << count_iter << " - x: " << c << " - f(x): " << fc << std::endl;
        }

        // Check for convergence
        if (abs(fc)<=abs_prec)
        {
            break;
        }

        // Check for maximum iterations limit
        if (count_iter > max_iter)
        {
            std::cerr << "WARNING: Bisection method could not find the solution ";
            std::cerr << "with the accurancy requested. Residual Value: " << fc << std::endl;
            break;
        }
        count_iter++;

        // Update interval values for the next iteration
        if (fa*fc>0)
        {
            a = c;
            fa = fc;
        }
        else
        {
            b = c;
        }
    }

    return c;
}


signed long long factorial(int n)
{
    signed long long pk = 1;
    for (int i=1; i<=n; i++)
    {
        pk *= i;
    }

    return pk;
}