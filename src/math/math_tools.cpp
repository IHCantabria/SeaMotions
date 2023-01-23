
// Include general usage libraries
#include <functional>

// Include external scientific libraries
#include <cmath>
#include "mkl.h"

// Include local libraries
#include "../config.hpp"
#include "math_tools.hpp"


// Include namespaces
using namespace std::literals::complex_literals;


//////////////////////////////////////////////
/////// FUNCTION DEFINITION BLOCK ////////////
//////////////////////////////////////////////
int assert_complex_equality(cuscomplex u, cuscomplex v, cusfloat epsilon)
{
    int pass = 0;
    if (
        (std::abs(u.real()-v.real())<epsilon)
        &&
        (std::abs(u.imag()-v.imag())<epsilon)
        )
    {
        pass = 1;
    }

    return pass;
}


cuscomplex complex_integration(
                                std::function <cuscomplex(cuscomplex)> f_def,
                                cuscomplex a,
                                cuscomplex b,
                                cusfloat tol
                                )
/**
 * @brief Integrate complex function along the segment [a,b].
 * 
 * Perform integration in the complex plane along the segment [a,b]. The integration is performed 
 * using Romberg's method.
 * 
 * \param f_def Target complex function to integate.
 * \param a First point of the segment.
 * \param b Second point of the segment.
 * \param tol Tolerance of integration. The tolerance is measured using Rirchardson Extrapolation.
 */

{
    // Create auxiliary function to calculate the real and 
    // imaginary parts
    cuscomplex jac = b-a;
    auto coord_map = [a, b](cusfloat t)->cuscomplex{return (a+t*(b-a));};
    auto f_real = [f_def, coord_map, jac](cusfloat t)->cusfloat{return (jac*f_def(coord_map(t))).real();};
    auto f_imag = [f_def, coord_map, jac](cusfloat t)->cusfloat{return (jac*f_def(coord_map(t))).imag();};

    // Perform integration of the real and imaginary parts
    cusfloat real_part = romberg_quadrature(f_real, 0.0, 1.0, tol);
    cusfloat imag_part = romberg_quadrature(f_imag, 0.0, 1.0, tol);
    cuscomplex z0 = real_part + imag_part*1i;

    return z0;
}


cuscomplex complex_integration_path(
                                    std::function <cuscomplex(cuscomplex)> f_def,
                                    int num_way_points,
                                    cuscomplex* way_points,
                                    cusfloat tol,
                                    bool close_path,
                                    bool verbose
                                    )
/**
 * @brief Integrate function along a path defined by a list of way points.
 * 
 * The integration along the way points list is performed using the function 
 * complex_integration(). Therfore, the integration in between successive 
 * way points is treated as a straight segment in between them.
 * 
 * \param f_def Target complex function to integrate
 * \param num_way_points Number of way points specified
 * \param way_points Way points that specify the integration path
 * \param tol Tolerance for the integration.
 * \param close_path Flag that specifies if the way points should be treated as closed path.
 *                   If true the last way point is joint with the first one in the list.
 * \param verbose Flag to activate the print-out to screen. This print-out specifies the value
 *                  of the integral at each interval of way-points list.
 * 
 */
{
    // Declare solution storage variable
    cuscomplex ci = 0.0 + 0.0i;
    cuscomplex sol = 0.0 + 0.0i;

    // Calculate limits for loop integration
    int N = num_way_points-1;
    if (close_path)
    {
        N++;
    }

    // Loop over segments
    int j = 0;
    for (int i=0; i<N; i++)
    {
        // Calculate forward index
        j = (i+1)%num_way_points;

        // Calculate integral over the segment
        ci = complex_integration(f_def, way_points[i], way_points[j], tol);
        sol += ci;
        if (verbose)
        {
            std::cout << "Path Segment [" << i << "]: " << std::endl;
            std::cout << "  -> A: " << way_points[i] << std::endl;
            std::cout << "  -> B: " << way_points[j] << std::endl;
            std::cout << "  -> Int. Value: " << ci << std::endl;
            std::cout << "  -> Int. Value. Cum: " << sol << std::endl;
        }
    }

    return sol;
}


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