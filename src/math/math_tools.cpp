
// Include external scientific libraries
#include <cmath>
#include "mkl.h"

// Include local libraries
#include "../config.hpp"
#include "math_tools.hpp"


//////////////////////////////////////////////
/////// FUNCTION DEFINITION BLOCK ////////////
//////////////////////////////////////////////
void bisection(std::function<cusfloat(cusfloat)> f_def, cusfloat a, cusfloat b, 
                cusfloat abs_prec, cusfloat rel_prec, int max_iter, bool verbose,
                cusfloat &sol, int &info)
{
    /**
     * @brief Bisection method to solve non-linear monoparametric equations.
     * 
     * The bisection method solves non-linear equations dependent on one parameter. It is an 
     * iterative solver based on the interval sub-division. More information can be found 
     * at: Numerical Recipes The Art of Scientific Computing. W.Press S.Teukolsky
     * 
     * \param f_def Target function
     * \param a Lower bound of the interval in which is assumed to have the zero.
     * \param b Upper bound of the interval in which is assumed to have the zero.
     * \param abs_prec Absolute precision to stop the search of the zero.
     * \param rel_prec Relative precision to stop the search of the zero.
     * \param max_iter Maximum iterations allowed to find the zero of the target function.
     * \param verbose Bool flag to print out to screen the iterative path to find the solution.
     * \param sol Scalar parameter in which store the solution of the iterative process.
     * \param info This parameter describes the finish status of the iterative process.
     *              - info == 0 -> Solution found.
     *              - info == 1 -> Solution not found or convergence problems.
     */

    // Start the iterative method
    cusfloat fa = f_def(a);
    cusfloat c, cb=(a+b)/2.0+2*rel_prec+1.0;
    cusfloat fc;

    // Start iterative loop
    int count_iter = 0;
    info = 0;
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
        if ((abs(fc)<=abs_prec) || (abs(c-cb)<rel_prec))
        {
            break;
        }

        // Check for maximum iterations limit
        if (count_iter > max_iter)
        {
            std::cerr << "WARNING: Bisection method could not find the solution ";
            std::cerr << "with the accurancy requested. Residual Value: " << fc << std::endl;
            info = 1;
            break;
        }
        count_iter++;

        // Update interval values for the next iteration
        cb = c;
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

    // Store solution in the output variable
    sol = c;
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