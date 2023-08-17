
// Include general usage libraries
#include <iostream>

// Include local modules
#include "integration.hpp"

// Include namespaces
using namespace std::literals::complex_literals;


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
