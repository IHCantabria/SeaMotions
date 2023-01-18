
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
                cusfloat fabs_tol, cusfloat xrel_tol, int max_iter, bool verbose,
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
     * \param fabs_tol Absolute value tolerance of target function value to stop iterations: f_def(x) < fabs_tol
     * \param xrel_tol Relative tolerance of the independent variable to stop the iterations: dx < x*xrel_tol
     * \param max_iter Maximum iterations allowed to find the zero of the target function.
     * \param verbose Bool flag to print out to screen the iterative path to find the solution.
     * \param sol Scalar parameter in which store the solution of the iterative process.
     * \param info This parameter describes the finish status of the iterative process.
     *              - info == 0 -> Solution found.
     *              - info == 1 -> Solution not found or convergence problems.
     */

    // Start the iterative method
    cusfloat fa = f_def(a);
    cusfloat c, cb=(a+b)/2.0+2*xrel_tol+1.0;
    cusfloat fc;

    // Start iterative loop
    cusfloat abs_err = 0.0;
    int count_iter = 0;
    info = 0;
    cusfloat rel_err = 0.0;
    while (true)
    {
        // Calculate new value
        c = (a+b)/2.0;
        fc = f_def(c);
        
        // Calculate errors
        abs_err = abs(fc);
        rel_err = abs(c-cb);

        // Print iterative process
        if (verbose)
        {
            std::cout << "Iter: " << count_iter << " - x: " << c << " - f(x): " << fc << std::endl;
            std::cout << "  -> a: " << a << " - f(a): " << f_def(a) << std::endl;
            std::cout << "  -> b: " << b << " - f(b): " << f_def(b) << std::endl;
            std::cout << "      - Abs.Error: " << abs_err << " - Rel.Error: " << rel_err << std::endl;
        }

        // Check for convergence
        if ((abs(fc)<=fabs_tol) || (abs(c-cb)<abs(c*xrel_tol+1e-14)))
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


void newton_raphson(std::function<cusfloat(cusfloat)> f_def, std::function<cusfloat(cusfloat)> f_der_def,
    cusfloat x0, cusfloat fabs_tol, cusfloat xrel_tol, int max_iter, bool verbose, cusfloat &sol, 
    int &info)
{
    /**
     * @brief Newthon-Rapshon method to solve non-linear equations of the type x=T(x).
     * 
     * This function implements the Newthon-Rapshon method to solve non-linear equations. The
     * solver is prepared to solve monoparametric equations (1D). The solver is based on the 
     * clasical algorithm that can be found at: 
     *  - Numerical Recipes The Art of Scientific Computing. W.Press S.Teukolsky
     * 
     * \param f_def Target function definition
     * \param f_der_der Target function derivative definition
     * \param x0 First pont to start iterations
     * \param fabs_tol Absolute value tolerance of target function value to stop iterations: f_def(x) < fabs_tol
     * \param xrel_tol Relative tolerance of the independent variable to stop the iterations: dx < x*xrel_tol
     * \param max_iter Maximum iterations allowed
     * \param verbose Flag to activate the iterative process print out to screen.
     * \param sol Solution to the equation provided.
     * \param info This parameter describes the finish status of the iterative process.
     *              - info == 0 -> Solution found.
     *              - info == 1 -> Solution not found or convergence problems.
     * 
     */
    // Set info variable to sucess. It will be set to failed
    // if the convergence is poor or it fails.
    info = 0;

    // Perform iterative loop to look for the solution
    cusfloat x = x0, dx = 0.0;
    int count_iter = 0;
    while (true)
    {
        // Calculate new increment
        dx = f_def(x)/f_der_def(x);

        // Calculate new solution point
        x -= dx;

        // Print-out to screen the iterative process status
        if (verbose)
        {
            std::cout << "Iter: " << count_iter << " - fabs(x): " << f_def(x);
            std::cout << " - x: " << x << " - dx: " << dx << std::endl; 
        }

        // Check for convergence
        if ((std::abs(f_def(x))<fabs_tol || (std::abs(dx)<=std::abs(xrel_tol*x+1e-14))))
        {
            break;
        }

        // Update iterations and check for limits
        if (count_iter > max_iter)
        {
            std::cerr << "WARNING: Newton-Raphson method could not find the solution ";
            std::cerr << "with the accurancy requested. Residual Value: " << f_def(x) << std::endl;
            info = 1;
            break;
        }
        count_iter++;
    }

    // Storage result in output channel
    sol = x;
}