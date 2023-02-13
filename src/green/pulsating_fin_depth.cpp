
// Include general usage libraries
#include <cassert>
#include <iostream>
#include <tuple>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "../config.hpp"
#include "integrals_db.hpp"
#include "../math/math_tools.hpp"
#include "../math/special_math.hpp"
#include "pulsating_inf_depth.hpp"
#include "pulsating_fin_depth.hpp"
#include "../waves.hpp"

// Include namespaces
using namespace std;
using namespace std::literals::complex_literals;


cusfloat calculate_dr_dx(cusfloat R, cusfloat dx)
{
    /**
     * @brief Calcualte derivative of R with respect to X
     * 
     * dr_dx = -dx/R
     * where R:sqrt((x-xi)**2.0+(y-yi)**2.0)
     * 
     * \param R Euclidean distance between the field and the source points in the horizontal plane.
     * \param dx X distance in between the field and the source points in the horizontal plane
     */

    return -dx/R;
}


cusfloat calculate_dr_dy(cusfloat R, cusfloat dy)
{
    /**
     * @brief Calcualte derivative of R with respect to X
     * 
     * dr_dy = -dy/R
     * where R:sqrt((x-xi)**2.0+(y-yi)**2.0)
     * 
     * \param R Euclidean distance between the field and the source points in the horizontal plane.
     * \param dy Y distance in between the field and the source points in the horizontal plane
     */

    return -dy/R;
}


cusfloat calculate_r(
                    cusfloat x,
                    cusfloat y,
                    cusfloat xi,
                    cusfloat eta
                    )
{
    /**
     * @brief Calculate euclidean distance in between field and source point in the
     *        horizontal plane
     * 
     * \param x X Coordinate of the field point
     * \param y Y Coordinate of the field point
     * \param xi X Coordinate of the source point
     * \param eta Y Coordinate of the source point
    * 
    */

    return sqrt( pow2s(x-xi) + pow2s(y-eta) );
}


cusfloat G1(
            cusfloat A,
            cusfloat B,
            cusfloat H,
            IntegralsDb &idb
            )
{
    cusfloat sol = 0.0;
    if (H > 1)
    {
        cusfloat A2 = pow2s(A);
        sol = (
                idb.l3->get_value_abh(A, B, H)
                -
                2/sqrt(A2 + pow2s(2+B))
                -
                2/sqrt(A2 + pow2s(2-B))
                );
    }
    else
    {
        sol = idb.l1->get_value_abh(A, B, H) + idb.l2->int_1d;
    }

    return sol;
}


cusfloat G1_dA(
                cusfloat A,
                cusfloat B,
                cusfloat H,
                IntegralsDb &idb
                )
    {
    cusfloat sol = 0.0;
    if (H > 1)
    {
        cusfloat A2 = pow2s(A);
        sol = (
                idb.l3_da->get_value_abh(A, B, H)
                +
                2*A/pow(A2 + pow2s(2+B), 3.0/2.0)
                +
                2*A/pow(A2 + pow2s(2-B), 3.0/2.0)
                );
    }
    else
    {
        sol = idb.l1_da->get_value_abh(A, B, H);
    }

    return sol;
}


cusfloat G1_dB(
                cusfloat A,
                cusfloat B,
                cusfloat H,
                IntegralsDb &idb
                )
{
    cusfloat sol = 0.0;
    if (H > 1)
    {
        cusfloat A2 = pow2s(A);
        sol = (
                idb.l3_db->get_value_abh(A, B, H)
                +
                2*(2+B)/pow(A2 + pow2s(2+B), 3.0/2.0)
                -
                2*(2-B)/pow(A2 + pow2s(2-B), 3.0/2.0)
                );
    }
    else
    {
        sol = idb.l1_db->get_value_abh(A, B, H);
    }

    return sol;
}


cusfloat G2(cusfloat A,
            cusfloat B,
            cusfloat H,
            IntegralsDb &idb
            )
{
    cusfloat sol = 0.0;
    if (H > 1.0)
    {
        sol = (
                idb.m3->get_value_abh(A, B, H)
                -
                2/sqrt(pow2s(A) + pow2s(2+B))
                );
    }
    else
    {
        sol = idb.m1->get_value_abh(A, B, H) + idb.m2->int_1d;
    }

    return sol;
}


cusfloat G2_dA(
                cusfloat A,
                cusfloat B,
                cusfloat H,
                IntegralsDb &idb
                )
{
    cusfloat sol = 0.0;
    if (H > 1.0)
    {
        sol = (
                idb.m3_da->get_value_abh(A, B, H)
                +
                2*A/pow(pow2s(A) + pow2s(2+B), 3.0/2.0)
                );
    }
    else
    {
        sol = idb.m1_da->get_value_abh(A, B, H);
    }

    return sol;
}


cusfloat G2_dB(cusfloat A,
                cusfloat B,
                cusfloat H,
                IntegralsDb &idb
                )
{
    cusfloat sol = 0.0;
    if (H > 1.0)
    {
        sol = (
                idb.m3_db->get_value_abh(A, B, H)
                +
                2*(2+B)/pow(pow2s(A) + pow2s(2+B), 3.0/2.0)
                );
    }
    else
    {
        sol = idb.m1_db->get_value_abh(A, B, H);
    }

    return sol;
}


cuscomplex G_integral(
                        cusfloat R,
                        cusfloat z,
                        cusfloat zeta,
                        cusfloat h,
                        WaveDispersionData &wave_data,
                        IntegralsDb &idb
                    )
{
    /**
     * @brief Calculate finite water depth Green function: steady sources + wave term
     *      at the integral region (R/h>1.0)
     * 
     * \param R Eucledian distance in the horizontal plane
     * \param z Z Coordinate of the field point
     * \param zeta Z Coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \param idb Integrals Database 
     * 
     */

    // Calculate derivative properties
    cusfloat v0 = abs(z-zeta);
    cusfloat v1 = z+zeta+2*h;
    cusfloat v2 = abs(z+zeta);
    cusfloat v3 = z-zeta+2*h;
    cusfloat v4 = zeta-z+2*h;
    cusfloat v5 = z+zeta+4*h;
    cusfloat r0 = sqrt(pow2s(R) + pow2s(v0));
    cusfloat r1 = sqrt(pow2s(R) + pow2s(v1));
    cusfloat r2 = sqrt(pow2s(R) + pow2s(v2));
    cusfloat r3 = sqrt(pow2s(R) + pow2s(v3));
    cusfloat r4 = sqrt(pow2s(R) + pow2s(v4));
    cusfloat r5 = sqrt(pow2s(R) + pow2s(v5));

    // Calculate steady part of the Green function
    cuscomplex green_steady = (1/r0 + 1/r1 + 1/r2 + 1/r3 + 1/r4 + 1/r5) + 0.0i;

    // Add wave term
    cuscomplex green_wave = wave_term_fin_depth_integral(R, z, zeta, h, wave_data, idb);

    // Compound total solution
    cuscomplex green_total = green_steady + green_wave;

    return green_total;
}


cuscomplex G_integral_dr(
                        cusfloat R,
                        cusfloat z,
                        cusfloat zeta,
                        cusfloat h,
                        WaveDispersionData &wave_data,
                        IntegralsDb &idb
                        )
{
    /**
     * @brief Calculate the derivative with respect to R of finite water depth 
     *      Green function: steady sources + wave term at the integral region
     *      (R/h>1.0)
     * 
     * \param R Eucledian distance in the horizontal plane
     * \param z Z Coordinate of the field point
     * \param zeta Z Coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \param idb Integrals Database 
     * 
     */

    // Calculate derivative properties
    cusfloat v0 = abs(z-zeta);
    cusfloat v1 = z+zeta+2*h;
    cusfloat v2 = abs(z+zeta);
    cusfloat v3 = z-zeta+2*h;
    cusfloat v4 = zeta-z+2*h;
    cusfloat v5 = z+zeta+4*h;
    cusfloat r0 = sqrt(pow2s(R) + pow2s(v0));
    cusfloat r1 = sqrt(pow2s(R) + pow2s(v1));
    cusfloat r2 = sqrt(pow2s(R) + pow2s(v2));
    cusfloat r3 = sqrt(pow2s(R) + pow2s(v3));
    cusfloat r4 = sqrt(pow2s(R) + pow2s(v4));
    cusfloat r5 = sqrt(pow2s(R) + pow2s(v5));

    // Calculate steady part of the Green function
    cuscomplex green_steady = -R*(
                                    1/pow3s(r0)
                                    +
                                    1/pow3s(r1)
                                    +
                                    1/pow3s(r2)
                                    +
                                    1/pow3s(r3)
                                    +
                                    1/pow3s(r4)
                                    +
                                    1/pow3s(r5)
                                    ) + 0.0i;

    // Add wave term
    cuscomplex green_wave = wave_term_fin_depth_integral_dr(R, z, zeta, h, wave_data, idb);

    // Compound total solution
    cuscomplex green_total = green_steady + green_wave;

    return green_total;
}


cuscomplex G_integral_dz(
                        cusfloat R,
                        cusfloat z,
                        cusfloat zeta,
                        cusfloat h,
                        WaveDispersionData &wave_data,
                        IntegralsDb &idb
                        )
{
    /**
     * @brief Calculate the derivative with respect to Z of finite water depth 
     *      Green function: steady sources + wave term at the integral region
     *      (R/h>1.0)
     * 
     * \param R Eucledian distance in the horizontal plane
     * \param z Z Coordinate of the field point
     * \param zeta Z Coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \param idb Integrals Database 
     * 
     */

    // Calculate derivative properties
    cusfloat v0 = abs(z-zeta);
    cusfloat v1 = z+zeta+2*h;
    cusfloat v2 = abs(z+zeta);
    cusfloat v3 = z-zeta+2*h;
    cusfloat v4 = zeta-z+2*h;
    cusfloat v5 = z+zeta+4*h;
    cusfloat r0 = sqrt(pow2s(R) + pow2s(v0));
    cusfloat r1 = sqrt(pow2s(R) + pow2s(v1));
    cusfloat r2 = sqrt(pow2s(R) + pow2s(v2));
    cusfloat r3 = sqrt(pow2s(R) + pow2s(v3));
    cusfloat r4 = sqrt(pow2s(R) + pow2s(v4));
    cusfloat r5 = sqrt(pow2s(R) + pow2s(v5));

    // Calculate steady part of the Green function
    cuscomplex green_steady = -(
                                sign(z-zeta)*v0/pow3s(r0)
                                +
                                v1/pow3s(r1)
                                +
                                sign(z+zeta)*v2/pow3s(r2)
                                +
                                v3/pow3s(r3)
                                -
                                v4/pow3s(r4)
                                +
                                v5/pow3s(r5)
                                ) + 0.0i;

    // Add wave term
    cuscomplex green_wave = wave_term_fin_depth_integral_dz(R, z, zeta, h, wave_data, idb);

    // Compound total solution
    cuscomplex green_total = green_steady + green_wave;

    return green_total;
}


cuscomplex john_series(
                        cusfloat R,
                        cusfloat z, 
                        cusfloat zeta, 
                        cusfloat h,
                        WaveDispersionData &wave_data
                        )
{
    /**
     * @brief John series representation of the finite water depth Green function
     * 
     * This is an eigenfunction expasion generated by Fritz John in the article
     * "On the Motion of Floating Bodies II". In the current implementation the 
     * series has been modified to work withou hyperbolic cosines in order to 
     * reduce numerical problems with big number. It is explained at:
     * "Consistent expression for the free-surfae Green function in finite water
     * depth" - Ed Mackay.
     * 
     * \param R Horizontal distance between source and the field point
     * \param z Vertical coordinate of the field point
     * \param zeta Vertical coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \return value of the Jonh series
     */

    // Check if John's coeffcients were precalculated
    assert(wave_data.is_john);

    // Get some data from wave data to have shorter names
    cusfloat k0 = wave_data.k0;

    // Calculate vertical distances between the source and
    // the field point
    cusfloat v3 = z+zeta;
    cusfloat v4 = z-zeta+2*h;
    cusfloat v5 = zeta-z+2*h;
    cusfloat v6 = z+zeta+4*h;

    // Calcuate real root series part
    cusfloat k0r = k0*R;
    cusfloat expsum = (
                        + exp(k0*v3)
                        + exp(-k0*v4)
                        + exp(-k0*v5)
                        + exp(-k0*v6)
                    );
    cusfloat c_real = -wave_data.k0nu*expsum;
    cuscomplex sol = c_real*(bessely0(k0r)+besselj0(k0r)*1i);

    // Calculate imag root series part
    cusfloat c_imag = 0.0;
    cusfloat ci = 0.0;
    int count_k = 0;
    cusfloat kni = 0.0;
    cusfloat zetah = zeta+h;
    cusfloat zh = z+h;
    while (true)
    {
        // Calculate i term of the series
        kni = wave_data.kn[count_k];
        ci = 4*wave_data.knnu[count_k]*cos(kni*zh)*cos(kni*zetah)*besselk0(kni*R);
        c_imag += ci;

        // Check for convergence
        if (abs(ci)<EPS_PRECISION)
        {
            break;
        }

        // Check for the limit in imaginary roots
        if (count_k > (wave_data.num_kn-2))
        {
            std::cerr << "Jonh series could not converge up to the precision required with" << std::endl;
            std::cerr << "Input parameters:"  << std::endl;
            std::cerr << "  - R/h: " << R/h << std::endl;
            std::cerr << "  - R: " << R << std::endl;
            std::cerr << "  - z: " << z << std::endl;
            std::cerr << "  - zeta: " << zeta << std::endl;
            std::cerr << "  - h: " << h << std::endl;
            std::cerr << "  - nu: " << wave_data.nu << std::endl;
            std::cerr << "  - k0: " << k0 << std::endl;
            std::cerr << "  - num_kn: " << wave_data.num_kn << std::endl;
            std::cerr << "* Current term value is: " << ci << std::endl;
            throw std::runtime_error("Jonh series value could not converge. See log file for details.");
        }

        // Update counter
        count_k++;
    }
    sol += c_imag;

    return sol;
}


cuscomplex john_series(
                        cusfloat x,
                        cusfloat y,
                        cusfloat z,
                        cusfloat xi,
                        cusfloat eta,
                        cusfloat zeta,
                        cusfloat h,
                        WaveDispersionData &wave_data
                        )
{
    /**
     * @brief John series representation of the finite water depth Green function
     * 
     * This is an overload for the function @ref john_series that acts as an interface.
     * 
     * \param x X Coordinate of the field point
     * \param y Y Coordinate of the field point
     * \param z Z Coordinate of the field point
     * \param xi X Coordinate of the source point
     * \param eta Y Coordinate of the source point
     * \param zeta Z Coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \return value of the Jonh series
     */

    // Calculate dependent parameters
    cusfloat R = std::sqrt(pow2s(x-xi)+pow2s(y-eta));

    // Call John kernel function
    return john_series(R, z, zeta, h, wave_data);
}


tuple<cuscomplex, cuscomplex> john_series_dhoriz(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionData &wave_data
                                                )
{
    /**
     * @brief Derivative with respect to X and Y of the John series representation of the 
     *        finite water depth Green function
     * 
     * Check out @see john_series() for more information.
     * 
     * \param x X Coordinate of the field point
     * \param y Y Coordinate of the field point
     * \param z Z Coordinate of the field point
     * \param xi X Coordinate of the source point
     * \param eta Y Coordinate of the source point
     * \param zeta Z Coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \return value of the Jonh series derivative with respect to R
     */

    // Calculate dependent parameters
    cusfloat dx = x-xi;
    cusfloat dy = y-eta;
    cusfloat R = sqrt(pow2s(dx)+pow2s(dy));

    // Calculate derivative with respect to R
    cuscomplex dr = john_series_dr(R, z, zeta, h, wave_data);

    // Calculate John series derivative with respect to X
    cuscomplex dj_dx = dr * calculate_dr_dx(R, dx);

    // Calculate John series derivative with respect to Y
    cuscomplex dj_dy = dr * calculate_dr_dy(R, dy);

    return make_tuple(dj_dx, dj_dy);
}


cuscomplex john_series_dr(
                            cusfloat R,
                            cusfloat z,
                            cusfloat zeta,
                            cusfloat h,
                            WaveDispersionData &wave_data
                            )
{
    /**
     * @brief Derivative with respect to R of the John series representation of the 
     *        finite water depth Green function
     * 
     * Check out @see john_series() for more information.
     * 
     * \param R Horizontal distance between source and the field point
     * \param z Vertical coordinate of the field point
     * \param zeta Vertical coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \return value of the Jonh series derivative with respect to R
     */

    // Check if John's coeffcients were precalculated
    assert(wave_data.is_john);

    // Get some data from wave data to have shorter names
    cusfloat k0 = wave_data.k0;

    // Calculate vertical distances between the source and
    // the field point
    cusfloat v3 = z+zeta;
    cusfloat v4 = z-zeta+2*h;
    cusfloat v5 = zeta-z+2*h;
    cusfloat v6 = z+zeta+4*h;

    // Calcuate real root series part
    cusfloat k0r = k0*R;
    cusfloat expsum = (
                        + exp(k0*v3)
                        + exp(-k0*v4)
                        + exp(-k0*v5)
                        + exp(-k0*v6)
                    );
    cusfloat c_real = -wave_data.k0nu*expsum;
    cuscomplex sol = -c_real*k0*(bessely1(k0r)+besselj1(k0r)*1i);

    // Calculate imag root series part
    cusfloat c_imag = 0.0;
    cusfloat ci = 0.0;
    int count_k = 0;
    cusfloat kni = 0.0;
    cusfloat zetah = zeta+h;
    cusfloat zh = z+h;
    while (true)
    {
        // Calculate i term of the series
        kni = wave_data.kn[count_k];
        ci = -4*wave_data.knnu[count_k]*cos(kni*zh)*cos(kni*zetah)*kni*besselk1(kni*R);
        c_imag += ci;

        // Check for convergence
        if (abs(ci)<EPS_PRECISION)
        {
            break;
        }

        // Check for the limit in imaginary roots
        if (count_k > (wave_data.num_kn-2))
        {
            std::cerr << "Jonh series dG/dR could not converge up to the precision required with" << std::endl;
            std::cerr << "Input parameters:"  << std::endl;
            std::cerr << "  - R/h: " << R/h << std::endl;
            std::cerr << "  - R: " << R << std::endl;
            std::cerr << "  - z: " << z << std::endl;
            std::cerr << "  - zeta: " << zeta << std::endl;
            std::cerr << "  - h: " << h << std::endl;
            std::cerr << "  - nu: " << wave_data.nu << std::endl;
            std::cerr << "  - k0: " << k0 << std::endl;
            std::cerr << "  - num_kn: " << wave_data.num_kn << std::endl;
            std::cerr << "* Current term value is: " << ci << std::endl;
            throw std::runtime_error("Jonh series value could not converge. See log file for details.");
        }

        // Update counter
        count_k++;
    }
    sol += c_imag;

    return sol;
}


cuscomplex john_series_dx(
                            cusfloat x,
                            cusfloat y,
                            cusfloat z,
                            cusfloat xi,
                            cusfloat eta,
                            cusfloat zeta,
                            cusfloat h,
                            WaveDispersionData &wave_data
                            )
{
    /**
     * @brief Derivative with respect to X of the John series representation of the 
     *        finite water depth Green function
     * 
     * Check out @see john_series() for more information.
     * 
     * \param x X Coordinate of the field point
     * \param y Y Coordinate of the field point
     * \param z Z Coordinate of the field point
     * \param xi X Coordinate of the source point
     * \param eta Y Coordinate of the source point
     * \param zeta Z Coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \return value of the Jonh series derivative with respect to R
     */

    // Calculate dependent parameters
    cusfloat dx = x-xi;
    cusfloat dy = y-eta;
    cusfloat R = sqrt(pow2s(dx)+pow2s(dy));

    // Calculate John series derivative with respect to X
    cuscomplex dj_dx = john_series_dr(R, z, zeta, h, wave_data) * calculate_dr_dx(R, dx);

    return dj_dx;
}


cuscomplex john_series_dy(
                            cusfloat x,
                            cusfloat y,
                            cusfloat z,
                            cusfloat xi,
                            cusfloat eta,
                            cusfloat zeta,
                            cusfloat h,
                            WaveDispersionData &wave_data
                            )
{
    /**
     * @brief Derivative with respect to Y of the John series representation of the 
     *        finite water depth Green function
     * 
     * Check out @see john_series() for more information.
     * 
     * \param x X Coordinate of the field point
     * \param y Y Coordinate of the field point
     * \param z Z Coordinate of the field point
     * \param xi X Coordinate of the source point
     * \param eta Y Coordinate of the source point
     * \param zeta Z Coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \return value of the Jonh series derivative with respect to R
     */

    // Calculate dependent parameters
    cusfloat dx = x-xi;
    cusfloat dy = y-eta;
    cusfloat R = sqrt(pow2s(dx)+pow2s(dy));

    // Calculate John series derivative with respect to Y
    cuscomplex dj_dy = john_series_dr(R, z, zeta, h, wave_data) * calculate_dr_dy(R, dy);

    return dj_dy;
}


cuscomplex john_series_dz(
                            cusfloat R,
                            cusfloat z,
                            cusfloat zeta,
                            cusfloat h,
                            WaveDispersionData &wave_data
                            )
{
    /**
     * @brief Derivative with respect to Z of the John series representation of the 
     *        finite water depth Green function
     * 
     * This is an eigenfunction expasion generated by Fritz John in the article
     * "On the Motion of Floating Bodies II". In the current implementation the 
     * series has been modified to work withou hyperbolic cosines in order to 
     * reduce numerical problems with big number. It is explained at:
     * "Consistent expression for the free-surfae Green function in finite water
     * depth" - Ed Mackay.
     * 
     * \param R Horizontal distance between source and the field point
     * \param z Vertical coordinate of the field point
     * \param zeta Vertical coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \return value of the Jonh series derivative with repect to Z
     */

    // Check if John's coeffcients were precalculated
    assert(wave_data.is_john);

    // Get some data from wave data to have shorter names
    cusfloat k0 = wave_data.k0;

    // Calculate vertical distances between the source and
    // the field point
    cusfloat v3 = z+zeta;
    cusfloat v4 = z-zeta+2*h;
    cusfloat v5 = zeta-z+2*h;
    cusfloat v6 = z+zeta+4*h;

    // Calcuate real root series part
    cusfloat k0r = k0*R;
    cusfloat expsum = -k0*(
                        - exp(k0*v3)
                        + exp(-k0*v4)
                        - exp(-k0*v5)
                        + exp(-k0*v6)
                    );
    cusfloat c_real = -wave_data.k0nu*expsum;
    cuscomplex sol = c_real*(bessely0(k0r)+besselj0(k0r)*1i);

    // Calculate imag root series part
    cusfloat c_imag = 0.0;
    cusfloat ci = 0.0;
    int count_k = 0;
    cusfloat kni = 0.0;
    cusfloat zetah = zeta+h;
    cusfloat zh = z+h;
    while (true)
    {
        // Calculate i term of the series
        kni = wave_data.kn[count_k];
        ci = -4*wave_data.knnu[count_k]*kni*sin(kni*zh)*cos(kni*zetah)*besselk0(kni*R);
        c_imag += ci;

        // Check for convergence
        if (abs(ci)<EPS_PRECISION)
        {
            break;
        }

        // Check for the limit in imaginary roots
        if (count_k > (wave_data.num_kn-2))
        {
            std::cerr << "Jonh series dG/dz could not converge up to the precision required with" << std::endl;
            std::cerr << "Input parameters:"  << std::endl;
            std::cerr << "  - R/h: " << R/h << std::endl;
            std::cerr << "  - R: " << R << std::endl;
            std::cerr << "  - z: " << z << std::endl;
            std::cerr << "  - zeta: " << zeta << std::endl;
            std::cerr << "  - h: " << h << std::endl;
            std::cerr << "  - nu: " << wave_data.nu << std::endl;
            std::cerr << "  - k0: " << k0 << std::endl;
            std::cerr << "  - num_kn: " << wave_data.num_kn << std::endl;
            std::cerr << "* Current term value is: " << ci << std::endl;
            throw std::runtime_error("Jonh series value could not converge. See log file for details.");
        }

        // Update counter
        count_k++;
    }
    sol += c_imag;

    return sol;
}


cuscomplex john_series_dz(
                            cusfloat x,
                            cusfloat y,
                            cusfloat z,
                            cusfloat xi,
                            cusfloat eta,
                            cusfloat zeta,
                            cusfloat h,
                            WaveDispersionData &wave_data
                            )
{
    /**
     * @brief Derivative with respect to Z of the John series representation of the 
     *        finite water depth Green function
     * 
     * This is an eigenfunction expasion generated by Fritz John in the article
     * "On the Motion of Floating Bodies II". In the current implementation the 
     * series has been modified to work withou hyperbolic cosines in order to 
     * reduce numerical problems with big number. It is explained at:
     * "Consistent expression for the free-surfae Green function in finite water
     * depth" - Ed Mackay.
     * 
     * \param x X Coordinate of the field point
     * \param y Y Coordinate of the field point
     * \param z Z Coordinate of the field point
     * \param xi X Coordinate of the source point
     * \param eta Y Coordinate of the source point
     * \param zeta Z Coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \return value of the Jonh series derivative with repect to Z
     */

    // Calculate derivative paramerters
    cusfloat R = calculate_r(x, y, xi, eta);

    // Calculate Z derivative using kernel function
    return john_series_dz(R, z, zeta, h, wave_data);
}


cuscomplex wave_term_fin_depth(
                                cusfloat x,
                                cusfloat y,
                                cusfloat z,
                                cusfloat xi,
                                cusfloat eta,
                                cusfloat zeta,
                                cusfloat h,
                                WaveDispersionData &wave_data,
                                IntegralsDb &idb
                                )
{
    // Calculate R value
    cusfloat R = calculate_r(x, y, xi, eta);

    // Define horizontal distance to depth ratio
    cusfloat rh = R/h;

    // Use series or integral formulations depending on rh value
    cuscomplex sol = 0.0;
    if (rh > 1.0)
    {
        sol = john_series(R, z, zeta, h , wave_data);
    }
    else
    {
        sol = wave_term_fin_depth_integral(R, z, zeta, h, wave_data, idb);
    }

    return sol;
}


cuscomplex wave_term_fin_depth_integral(
                                        cusfloat R,
                                        cusfloat z,
                                        cusfloat zeta,
                                        cusfloat h,
                                        WaveDispersionData &wave_data,
                                        IntegralsDb &idb
                                        )
{
    /**
     * @brief Wave term for the finite water depth function (Integral representation).
     * 
     * The formulation used here (valid for R/h<1.0) is mainly taken from:
     * "Consistent expressions for the free surface function in
     *  finite water depth - Ed Mackay".
     * 
     * \param R Eucleadian distance in between field and source points in the horizontal plane
     * \param z Z Coordinate of the field point
     * \param zeta Z Coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \param idb Integrals Database
     * \return value of the integral representation
     */

    // Copy variables to the stack
    cusfloat k0 = wave_data.k0;
    cusfloat nu = wave_data.nu;

    // Calculate dependent parameters
    cusfloat v0 = abs(z-zeta);
    cusfloat v1 = z+zeta+2*h;
    cusfloat v2 = abs(z+zeta);
    cusfloat v3 = z-zeta+2*h;
    cusfloat v4 = zeta-z+2*h;
    cusfloat v5 = z+zeta+4*h;

    cusfloat A = R/h;
    cusfloat B0 = v0/h;
    cusfloat B1 = v1/h;
    cusfloat H = nu*h;

    // Check that B0 and  B1 is in between limits
    assert( ((B0>=0.0) && (B0<=1.0)) && "B0 is out of interval [0, 1]" );
    assert( ((B1>=0.0) && (B1<=2.0)) && "B1 is out of interval [0, 2]" );

    // std::cout << "A: " << A << " - B1: " << B1 << " - H: " << H << std::endl;

    // Calculate real part
    cusfloat G_real = 0.0;
    if (B1 <= 1)
    {
        G_real = (G1(A, B0, H, idb) + G1(A, B1, H, idb))/h;
    }
    else
    {
        // Calculate dependent param
        cusfloat X = nu*R;
        cusfloat Y = nu*abs(z+zeta);
        cusfloat g1 = G1(A, B0, H, idb);
        cusfloat g2 = G2(A, B1, H, idb);
        G_real = (g1 + g2)/h + nu*wave_term_inf_depth(X, Y);
    }

    // Calculate imaginary part
    cusfloat expsum = (
                        +exp(-k0*v2)
                        +exp(-k0*v3)
                        +exp(-k0*v4)
                        +exp(-k0*v5)
                        );
    cusfloat G_imag = -wave_data.k0nu*expsum*besselj0(k0*R);

    // Define G integral as a complex number
    cuscomplex G = G_real + G_imag*1i;

    return G;
}


cuscomplex wave_term_fin_depth_integral_dr(
                                            cusfloat R,
                                            cusfloat z,
                                            cusfloat zeta,
                                            cusfloat h,
                                            WaveDispersionData &wave_data,
                                            IntegralsDb &idb
                                            )
{
    /**
     * @brief R derivative of the wave term for the finite water depth 
     *        function (Integral representation).
     * 
     * The formulation used here (valid for R/h<1.0) is mainly taken from:
     * "Consistent expressions for the free surface function in
     *  finite water depth - Ed Mackay".
     * 
     * \param R Eucleadian distance in between field and source points in the horizontal plane
     * \param z Z Coordinate of the field point
     * \param zeta Z Coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \param idb Integrals Database
     * \return value of the integral representation
     */

    // Copy variables to the stack
    cusfloat k0 = wave_data.k0;
    cusfloat nu = wave_data.nu;

    // Calculate dependent parameters
    cusfloat v0 = abs(z-zeta);
    cusfloat v1 = z+zeta+2*h;
    cusfloat v2 = abs(z+zeta);
    cusfloat v3 = z-zeta+2*h;
    cusfloat v4 = zeta-z+2*h;
    cusfloat v5 = z+zeta+4*h;

    cusfloat A = R/h;
    cusfloat B0 = v0/h;
    cusfloat B1 = v1/h;
    cusfloat H = nu*h;

    // Check that B0 and  B1 is in between limits
    assert( ((B0>=0.0) && (B0<=1.0)) && "B0 is out of interval [0, 1]" );
    assert( ((B1>=0.0) && (B1<=2.0)) && "B1 is out of interval [0, 2]" );

    // std::cout << "A: " << A << " - B1: " << B1 << " - H: " << H << std::endl;

    // Calculate real part
    cusfloat G_real = 0.0;
    if (B1 <= 1)
    {
        G_real = (G1_dA(A, B0, H, idb) + G1_dA(A, B1, H, idb))/pow2s(h);
    }
    else
    {
        // Calculate dependent param
        cusfloat X = nu*R;
        cusfloat Y = nu*abs(z+zeta);
        cusfloat g1 = G1_dA(A, B0, H, idb);
        cusfloat g2 = G2_dA(A, B1, H, idb);
        G_real = (g1 + g2)/pow2s(h) + pow2s(nu)*wave_term_inf_depth_dx(X, Y);
    }

    // Calculate imaginary part
    cusfloat expsum = (
                        +exp(-k0*v2)
                        +exp(-k0*v3)
                        +exp(-k0*v4)
                        +exp(-k0*v5)
                        );
    cusfloat G_imag = wave_data.k0nu*expsum*besselj1(k0*R)*k0;

    // Define G integral as a complex number
    cuscomplex G = G_real + G_imag*1i;

    return G;
}


cuscomplex wave_term_fin_depth_integral_dz(
                                            cusfloat R,
                                            cusfloat z,
                                            cusfloat zeta,
                                            cusfloat h,
                                            WaveDispersionData &wave_data,
                                            IntegralsDb &idb
                                            )
{
    /**
     * @brief Z derivative of the wave term for the finite water depth 
     *        function (Integral representation).
     * 
     * The formulation used here (valid for R/h<1.0) is mainly taken from:
     * "Consistent expressions for the free surface function in
     *  finite water depth - Ed Mackay".
     * 
     * \param R Eucleadian distance in between field and source points in the horizontal plane
     * \param z Z Coordinate of the field point
     * \param zeta Z Coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \param idb Integrals Database
     * \return value of the integral representation
     */

    // Copy variables to the stack
    cusfloat k0 = wave_data.k0;
    cusfloat nu = wave_data.nu;

    // Calculate dependent parameters
    cusfloat v0 = abs(z-zeta);
    cusfloat v1 = z+zeta+2*h;
    cusfloat v2 = abs(z+zeta);
    cusfloat v3 = z-zeta+2*h;
    cusfloat v4 = zeta-z+2*h;
    cusfloat v5 = z+zeta+4*h;

    cusfloat A = R/h;
    cusfloat B0 = v0/h;
    cusfloat B1 = v1/h;
    cusfloat H = nu*h;

    // Check that B0 and  B1 is in between limits
    assert( ((B0>=0.0) && (B0<=1.0)) && "B0 is out of interval [0, 1]" );
    assert( ((B1>=0.0) && (B1<=2.0)) && "B1 is out of interval [0, 2]" );

    // std::cout << "A: " << A << " - B1: " << B1 << " - H: " << H << std::endl;

    // Calculate real part
    cusfloat G_real = 0.0;
    if (B1 <= 1)
    {
        G_real = (G1_dB(A, B0, H, idb)*sign(z-zeta) + G1_dB(A, B1, H, idb))/pow2s(h);
    }
    else
    {
        // Calculate dependent param
        cusfloat X = nu*R;
        cusfloat Y = nu*abs(z+zeta);
        cusfloat g1 = G1_dB(A, B0, H, idb)*sign(z-zeta);
        cusfloat g2 = G2_dB(A, B1, H, idb);
        G_real = (g1 + g2)/pow2s(h) + pow2s(nu)*wave_term_inf_depth_dy(X, Y)*sign(z+zeta);
    }

    // Calculate imaginary part
    cusfloat expsum = (
                        +exp(-k0*v2)*sign(z+zeta)
                        +exp(-k0*v3)
                        -exp(-k0*v4)
                        +exp(-k0*v5)
                        );
    cusfloat G_imag = k0*wave_data.k0nu*expsum*besselj0(k0*R);

    // Define G integral as a complex number
    cuscomplex G = G_real + G_imag*1i;

    return G;
}