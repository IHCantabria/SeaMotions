
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

// Include general usage libraries
#include <functional>
#include <optional>

// Include external scientific libraries
#include <cmath>
#include "mkl.h"

// Include local libraries
#include "../config.hpp"
#include "math_tools.hpp"


//////////////////////////////////////////////
/////// FUNCTION DEFINITION BLOCK ////////////
//////////////////////////////////////////////
cusfloat    angfreq_to_freq( 
                                cusfloat angfreq 
                            )
{
    return angfreq / 2.0 / PI;
}


cusfloat    angfreq_to_period( 
                                cusfloat angfreq 
                            )
{
    return 2.0 * PI / angfreq;
}


int assert_complex_equality(
                                cuscomplex u, 
                                cuscomplex v, 
                                cusfloat epsilon
                            )
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


void bisection(
                    std::function<cusfloat(cusfloat)> f_def, 
                    cusfloat a, 
                    cusfloat b, 
                    cusfloat fabs_tol, 
                    cusfloat xrel_tol, 
                    int max_iter, 
                    bool verbose,
                    cusfloat &sol,
                    int &info
                )
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


cusfloat check_zero_eps( 
                            cusfloat value,
                            cusfloat eps
                        )
{
    if ( std::abs( value ) < eps )
    {
        value = eps;
    }
    
    return value;
}


void    conj_vector(
                            int         n,
                            cuscomplex* u,
                            cuscomplex* v
                        )
{
    for ( int i=0; i<n; i++ )
    {
        v[i] = std::conj( u[i] );
    }
}


cusfloat    cos_alpha(
                        int         alpha,
                        cusfloat    theta
                    )
{
	cusfloat val = 0.0;
    if ( alpha != 0 )
    {
        // Cast alpha value to cusfloat
        cusfloat af = static_cast<cusfloat>( alpha );

        // Evaluate primitive integral
        val = std::cos( theta * af ) / af;
    }
    else
    {
        val = 0.0;
    }

    return val;
}


cusfloat    cos3_int_0_2PI( 
                            int m,
                            int n,
                            int p
                        )
{
    // Define integral bounds
    cusfloat    th_0    = 0.0;
    cusfloat    th_1    = 2 * PI;

    // Calculate the integral of the first term
    int         alpha_0 = m - n - p;
    cusfloat    t0      = sin_alpha( alpha_0, th_1 ) - sin_alpha( alpha_0, th_0 );

    // Calculate the integral of the second term
    int         alpha_1 = m + n - p;
    cusfloat    t1      = sin_alpha( alpha_1, th_1 ) - sin_alpha( alpha_1, th_0 );

    // Calculate the integral of the third term
    int         alpha_2 = m - n + p;
    cusfloat    t2      = sin_alpha( alpha_2, th_1 ) - sin_alpha( alpha_2, th_0 );

    // Calculate the integral of the fourth term
    int         alpha_3 = m + n + p;
    cusfloat    t3      = sin_alpha( alpha_3, th_1 ) - sin_alpha( alpha_3, th_0 );

    return 0.25 * ( t0 + t1 + t2 + t3 );
}


cusfloat    cos2sin_int_0_2PI( 
                                int m,
                                int n,
                                int p
                            )
{
    // Define integral bounds
    cusfloat    th_0    = 0.0;
    cusfloat    th_1    = 2 * PI;

    // Calculate the integral of the first term
    int         alpha_0 = m - n - p;
    cusfloat    t0      = cos_alpha( alpha_0, th_1 ) - cos_alpha( alpha_0, th_0 );

    // Calculate the integral of the second term
    int         alpha_1 = m + n - p;
    cusfloat    t1      = cos_alpha( alpha_1, th_1 ) - cos_alpha( alpha_1, th_0 );

    // Calculate the integral of the third term
    int         alpha_2 = m - n + p;
    cusfloat    t2      = cos_alpha( alpha_2, th_1 ) - cos_alpha( alpha_2, th_0 );

    // Calculate the integral of the fourth term
    int         alpha_3 = m + n + p;
    cusfloat    t3      = cos_alpha( alpha_3, th_1 ) - cos_alpha( alpha_3, th_0 );

    return 0.25 * ( t0 + t1 - t2 - t3 );
}


cusfloat    cossin2_int_0_2PI( 
                                int m,
                                int n,
                                int p
                            )
{
    // Define integral bounds
    cusfloat    th_0    = 0.0;
    cusfloat    th_1    = 2 * PI;

    // Calculate the integral of the first term
    int         alpha_0 = m - n - p;
    cusfloat    t0      = sin_alpha( alpha_0, th_1 ) - sin_alpha( alpha_0, th_0 );

    // Calculate the integral of the second term
    int         alpha_1 = m + n - p;
    cusfloat    t1      = sin_alpha( alpha_1, th_1 ) - sin_alpha( alpha_1, th_0 );

    // Calculate the integral of the third term
    int         alpha_2 = m - n + p;
    cusfloat    t2      = sin_alpha( alpha_2, th_1 ) - sin_alpha( alpha_2, th_0 );

    // Calculate the integral of the fourth term
    int         alpha_3 = m + n + p;
    cusfloat    t3      = sin_alpha( alpha_3, th_1 ) - sin_alpha( alpha_3, th_0 );

    return 0.25 * ( -t0 + t1 + t2 - t3 );
}


cusfloat    deg_to_rad(  cusfloat deg )
{
    return deg*PI/180.0;
}


signed long long factorial( int n )
{
    signed long long pk = 1;
    for (int i=1; i<=n; i++)
    {
        pk *= i;
    }

    return pk;
}


cusfloat    freq_to_angfreq( cusfloat freq )
{
    return 2*PI*freq;
}


void        newton_raphson(
                               std::function<cusfloat(cusfloat)> f_def, 
                               std::function<cusfloat(cusfloat)> f_der_def,
                               cusfloat x0, 
                               cusfloat fabs_tol, 
                               cusfloat xrel_tol, 
                               int max_iter, 
                               bool verbose, 
                               cusfloat &sol, 
                               int &info
                           )
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


cusfloat    period_to_angfreq( cusfloat period )
{
    return 2*PI/period;
}


cusfloat    sin_alpha( 
                        int         alpha,
                        cusfloat    theta
                    )
{
    cusfloat val = 0.0;

    if ( alpha != 0 )
    {
        // Cast alpha integer value to cusfloat
        cusfloat af  = static_cast<cusfloat>( alpha );
        
        // Evaluate primitive integral 
        val = std::sin( theta * af ) / af;
    }
    else
    {
        // Evaluate primitive integral
        val = theta;
    }

    return val;
}


cusfloat    sin3_int_0_2PI( 
                                int m,
                                int n,
                                int p
                            )
{
    // Define integral bounds
    cusfloat    th_0    = 0.0;
    cusfloat    th_1    = 2 * PI;

    // Calculate the integral of the first term
    int         alpha_0 = m - n - p;
    cusfloat    t0      = cos_alpha( alpha_0, th_1 ) - cos_alpha( alpha_0, th_0 );

    // Calculate the integral of the second term
    int         alpha_1 = m + n - p;
    cusfloat    t1      = cos_alpha( alpha_1, th_1 ) - cos_alpha( alpha_1, th_0 );

    // Calculate the integral of the third term
    int         alpha_2 = m - n + p;
    cusfloat    t2      = cos_alpha( alpha_2, th_1 ) - cos_alpha( alpha_2, th_0 );

    // Calculate the integral of the fourth term
    int         alpha_3 = m + n + p;
    cusfloat    t3      = cos_alpha( alpha_3, th_1 ) - cos_alpha( alpha_3, th_0 );

    return 0.25 * ( t0 - t1 - t2 + t3 );
}


void triangle_geom_properties( 
                                    const   cusfloat*   x,
                                    const   cusfloat*   y,
                                    const   cusfloat*   z,
                                            cusfloat&   area,
                                            cusfloat*   centroid,
                                            cusfloat*   moments_fo,
                                            cusfloat*   moments_so,
                                    const   std::optional<cusfloat*>& ref_sys = std::nullopt
                                )
{
    // Calculate area using the cross product of two sides
    cusfloat v0[3] = { x[1]-x[0], y[1]-y[0], z[1]-z[0] };
    cusfloat v1[3] = { x[2]-x[0], y[2]-y[0], z[2]-z[0] };
    cusfloat v2[3] = { 0.0, 0.0, 0.0 };

    // Calculate cross product
    cross( v0, v1, v2 );

    // Calculate area of the triangle
    cusfloat v2_mod = std::sqrt( pow2s( v2[0] ) + pow2s( v2[1] ) + pow2s( v2[2] ) );
    area            = v2_mod / 2.0;

    // Calculate centroid of the triangle
    centroid[0] = ( x[0] + x[1] + x[2] ) / 3.0;
    centroid[1] = ( y[0] + y[1] + y[2] ) / 3.0;
    centroid[2] = ( z[0] + z[1] + z[2] ) / 3.0;

    // Choose reference system
    cusfloat ref_sys_local[3] = { 0.0, 0.0, 0.0 };
    if ( ref_sys.has_value() )
    {
        ref_sys_local[0] = ref_sys.value()[0];
        ref_sys_local[1] = ref_sys.value()[1];
        ref_sys_local[2] = ref_sys.value()[2];
    }
    else
    {
        ref_sys_local[0] = centroid[0];
        ref_sys_local[1] = centroid[1];
        ref_sys_local[2] = centroid[2];
    }

    // Shift vertices to reference system
    cusfloat vertex_cent[3][3];
    for ( int i=0; i<3; i++ )
    {
        vertex_cent[i][0] = x[i] - ref_sys_local[0];
        vertex_cent[i][1] = y[i] - ref_sys_local[1];
        vertex_cent[i][2] = z[i] - ref_sys_local[2];
    }

    // Calculate vertex summation
    cusfloat vertex_sum[3]  = { 0.0, 0.0, 0.0 };
    cusfloat vertex_sum2[3] = { 0.0, 0.0, 0.0 };
    for ( int i=0; i<3; i++ ) 
    {
        vertex_sum[0] += vertex_cent[i][0];
        vertex_sum[1] += vertex_cent[i][1];
        vertex_sum[2] += vertex_cent[i][2];

        vertex_sum2[0] += pow2s( vertex_cent[i][0] );
        vertex_sum2[1] += pow2s( vertex_cent[i][1] );
        vertex_sum2[2] += pow2s( vertex_cent[i][2] );
    }

    // Calculate first order moments
    for ( int i=0; i<3; i++ )
    {
        moments_fo[i] = area / 3.0 * vertex_sum[i];
    }

    // Calculate second order moments
    for ( int i=0; i<2; i++ )
    {
        int j = ( i + 1 ) % 2;
        moments_so[j]   = ( 
                                vertex_sum2[i]
                                +
                                vertex_cent[0][i] * vertex_cent[1][i]
                                +
                                vertex_cent[0][i] * vertex_cent[2][i]
                                +
                                vertex_cent[1][i] * vertex_cent[2][i]
                            ) * area / 6.0;
    }

    moments_so[2]   =   ( 
                            2.0 * vertex_cent[0][0] * vertex_cent[0][1]
                            +
                            vertex_cent[1][0] * vertex_cent[0][1]
                            +
                            vertex_cent[2][0] * vertex_cent[0][1]
                            +
                            vertex_cent[0][0] * vertex_cent[1][1]
                            +
                            2.0 * vertex_cent[1][0] * vertex_cent[1][1]
                            +
                            vertex_cent[2][0] * vertex_cent[1][1]
                            +
                            vertex_cent[0][0] * vertex_cent[2][1]
                            +
                            vertex_cent[1][0] * vertex_cent[2][1]
                            +
                            2.0 * vertex_cent[2][0] * vertex_cent[2][1]
                        ) * area / 12.0;
}