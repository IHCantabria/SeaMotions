
// Include general usage libraries
#include <cassert>
#include <iostream>
#include <tuple>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "../config.hpp"
#include "integrals_db.hpp"
#include "../math/bessel_factory.hpp"
#include "../math/math_interface.hpp"
#include "../math/math_tools.hpp"
#include "../math/special_math.hpp"
#include "pulsating_inf_depth.hpp"
#include "pulsating_fin_depth_v2.hpp"
#include "pulsating_inf_depth_v2.hpp"
#include "chebyshev_evaluator_interface.hpp"

// Include namespaces
using namespace std;
using namespace std::literals::complex_literals;


template<std::size_t N>
void         G1_Hlt1(
                                                const std::size_t   n,
                                                cusfloat*           A,
                                                cusfloat*           B,
                                                cusfloat            H,
                                                cusfloat*           results,
                                                cusfloat*           results_da,
                                                cusfloat*           results_db
                )
{
    // Clear vector to avoid spurious remaining data
    clear_vector<cusfloat, N>( results );

    // Calculate L1 integral using chebyshev expansions
    L1CEV<N>::evaluate( n, A, B, results );
    L1_dACEV<N>::evaluate( n, A, B, results_da );
    L1_dBCEV<N>::evaluate( n, A, B, results_db );

    // Add L2 integral scalar contribution
    for ( std::size_t i=0; i<n; i++ )
    {
        results[i] += ChebyshevTraits<L2C>::coeffs;
    }

}


template<std::size_t N>
void         G1_Hgt1(
                                                const std::size_t   n,
                                                cusfloat*           A,
                                                cusfloat*           B,
                                                cusfloat            H,
                                                cusfloat*           results,
                                                cusfloat*           results_da,
                                                cusfloat*           results_db
                )
{
    // Clear vector to avoid spurious remaining data
    clear_vector<cusfloat, N>( results );

    // Calculate auxiliar values
    cusfloat A2[N];
    cusfloat BP[N];
    cusfloat BM[N];

    for ( std::size_t i=0; i<n; i++ )
    {
        A2[i] = POW2S( A[i] );
    }

    for ( std::size_t i=0; i<n; i++ )
    {
        BP[i] = A2[i] + POW2S( 2+B[i] )
    }

    for ( std::size_t i=0; i<n; i++ )
    {
        BM[i] = A2[i] + POW2S( 2-B[i] )
    }

    lv_sqrt<cusfloat>( n, BP, BP );
    lv_sqrt<cusfloat>( n, BM, BM );
    
    // Calculate L3 using chebyshev expansion
    L3CEV<N>::evaluate( n, A, B, results );

    // Substract auxiliar variables
    for ( std::size_t i=0; i<n; i++ )
    {
        results[i] -= 2.0 / BP[i];
    }

    for ( std::size_t i=0; i<n; i++ )
    {
        results[i] -= 2.0 / BM[i];
    }
}


template<std::size_t N>
cusfloat    G1_dA(
                                                cusfloat*   A,
                                                cusfloat*   B,
                                                cusfloat    H,
                                                cusfloat*   results
                )
    {
    // Clear vector to avoid spurious remaining data
    clear_vector<cusfloat, N>( results );

    // Calculate G1 integral values
    if (H > 1)
    {
        // Calculate auxiliar values
        cusfloat A2[N];
        cusfloat BP[N];
        cusfloat BM[N];

        for ( std::size_t i=0; i<N; i++ )
        {
            A2[i] = POW2S( A[i] );
        }

        for ( std::size_t i=0; i<N; i++ )
        {
            BP[i] = A2[i] + POW2S( 2+B[i] )
        }

        for ( std::size_t i=0; i<N; i++ )
        {
            BM[i] = A2[i] + POW2S( 2-B[i] )
        }

        lv_pow3o2<cusfloat>( N, BP, BP );
        lv_pow3o2<cusfloat>( N, BM, BM );
        
        // Calculate L3_dA chebyshev expansion
        L3_dACEV<N>::evaluate( A, B, results );

        // Substract auxiliar variables
        for ( std::size_t i=0; i<N; i++ )
        {
            results[i] += 2.0 * A[i] / BP[i];
        }

        for ( std::size_t i=0; i<N; i++ )
        {
            results[i] += 2.0 * A[i] / BM[i];
        }
    }
    else
    {
        // Calculate L1_dA integral using chebyshev expansions
        L1_dACEV<N>::evaluate( A, B, results );
    }

    return sol;
}


template<std::size_t N>
cusfloat    G1_dB(
                                                cusfloat*   A,
                                                cusfloat*   B,
                                                cusfloat    H,
                                                cusfloat*   results
                )
{
    // Clear vector to avoid spurious remaining data
    clear_vector<cusfloat, N>( results );

    // Calculate G1 integral values
    if (H > 1)
    {
        // Calculate auxiliar values
        cusfloat A2[N];
        cusfloat BP[N];
        cusfloat BM[N];

        for ( std::size_t i=0; i<N; i++ )
        {
            A2[i] = POW2S( A[i] );
        }

        for ( std::size_t i=0; i<N; i++ )
        {
            BP[i] = A2[i] + POW2S( 2+B[i] )
        }

        for ( std::size_t i=0; i<N; i++ )
        {
            BM[i] = A2[i] + POW2S( 2-B[i] )
        }

        lv_pow3o2<cusfloat>( N, BP, BP );
        lv_pow3o2<cusfloat>( N, BM, BM );
        
        // Calculate L3_dA chebyshev expansion
        L3_dACEV<N>::evaluate( A, B, results );

        // Substract auxiliar variables
        for ( std::size_t i=0; i<N; i++ )
        {
            results[i] += 2.0 * A[i] / BP[i];
        }

        for ( std::size_t i=0; i<N; i++ )
        {
            results[i] += 2.0 * A[i] / BM[i];
        }

        cusfloat A2 = pow2s(A);
        sol = (
                idb.l3_db->get_value_ab(A, B)
                +
                2*(2+B)/pow(A2 + pow2s(2+B), 3.0/2.0)
                -
                2*(2-B)/pow(A2 + pow2s(2-B), 3.0/2.0)
                );
    }
    else
    {
        sol = idb.l1_db->get_value_ab(A, B);
    }

    return sol;
}


cusfloat    G2(
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
                idb.m3->get_value_ab(A, B)
                -
                2/sqrt(pow2s(A) + pow2s(2+B))
                );
    }
    else
    {
        sol = idb.m1->get_value_ab(A, B) + idb.m2->int_1d;
    }

    return sol;
}


cusfloat    G2_dA(
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
                idb.m3_da->get_value_ab(A, B)
                +
                2*A/pow(pow2s(A) + pow2s(2+B), 3.0/2.0)
                );
    }
    else
    {
        sol = idb.m1_da->get_value_ab(A, B);
    }

    return sol;
}


cusfloat    G2_dB(
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
                idb.m3_db->get_value_ab(A, B)
                +
                2*(2+B)/pow(pow2s(A) + pow2s(2+B), 3.0/2.0)
                );
    }
    else
    {
        sol = idb.m1_db->get_value_ab(A, B);
    }

    return sol;
}


cuscomplex  G_integral(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
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


cuscomplex  G_integral_steady(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h
                            )
{
    /**
     * @brief Calculate finite water depth Green function: steady sources + wave term
     *      at the integral region (R/h>1.0)
     * 
     * \param x X Coordinate of the field point
     * \param y Y Coordinate of the field point
     * \param z Z Coordinate of the field point
     * \param xi X Coordinate of the source point
     * \param eta Y Coordinate of the source point
     * \param zeta Z Coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \param idb Integrals Database 
     * 
     */

    // Calculate derivative properties
    cusfloat R = calculate_r(x, y, xi, eta);
    
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

    return green_steady;
}


cuscomplex  G_integral_wave(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
                                                IntegralsDb &idb
                            )
{
    /**
     * @brief Calculate finite water depth Green function: wave term
     *      at the integral region (R/h>1.0)
     * 
     * \param x X Coordinate of the field point
     * \param y Y Coordinate of the field point
     * \param z Z Coordinate of the field point
     * \param xi X Coordinate of the source point
     * \param eta Y Coordinate of the source point
     * \param zeta Z Coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \param idb Integrals Database 
     * 
     */

    // Calculate derivative properties
    cusfloat R = calculate_r(x, y, xi, eta);
    

    return wave_term_fin_depth_integral(R, z, zeta, h, wave_data, idb);
}



cuscomplex  G_integral_dr(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
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


cuscomplex  G_integral_steady_dr(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h
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

    return green_steady;
}


cuscomplex G_integral_wave_dr(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
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

    // Add wave term
    cuscomplex green_wave = wave_term_fin_depth_integral_dr(R, z, zeta, h, wave_data, idb);

    return green_wave;
}


cuscomplex  G_integral_dz(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
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


cuscomplex  G_integral_steady_dz(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h
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

    return green_steady;
}


cuscomplex  G_integral_wave_dz(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
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

    // Add wave term
    cuscomplex green_wave = wave_term_fin_depth_integral_dz(R, z, zeta, h, wave_data, idb);

    return green_wave;
}




cuscomplex  john_series(
                                                cusfloat R,
                                                cusfloat z, 
                                                cusfloat zeta, 
                                                cusfloat h,
                                                WaveDispersionFO &wave_data
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
    cuscomplex sol = c_real*(bessely0(k0r)-besselj0(k0r)*1i);

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


cuscomplex  john_series(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data
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


tuplecc     john_series_dhoriz(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data
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
    cuscomplex dj_dy = dr * calculate_dr_dx(R, dy);

    return make_tuple(dj_dx, dj_dy);
}


cuscomplex  john_series_dr(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data
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
    cuscomplex sol = -c_real*k0*(bessely1(k0r)-besselj1(k0r)*1i);

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


cuscomplex  john_series_dx(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data
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


cuscomplex  john_series_dy(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data
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
    cuscomplex dj_dy = john_series_dr(R, z, zeta, h, wave_data) * calculate_dr_dx(R, dy);

    return dj_dy;
}


cuscomplex  john_series_dz(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data
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
    cuscomplex sol = c_real*(bessely0(k0r)-besselj0(k0r)*1i);

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


cuscomplex  john_series_dz(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data
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


cuscomplex  wave_term_fin_depth(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
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


template<std::size_t N>
void        wave_term_fin_depth_integral(
                                                cusfloat*           R,
                                                cusfloat*           z,
                                                cusfloat*           zeta,
                                                cusfloat            h,
                                                BesselFactoryVec<N> &bessel_factory,
                                                WaveDispersionFO    &wave_data,
                                                cuscomplex*         G,
                                                cuscomplex*         G_dr,
                                                cuscomplex*         G_dz
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
     * \param k0 Wave number of the target wave to calculate in finite depth
     * \param nu Wave number of the target wave to calculate in infinite depth 
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \param idb Integrals Database
     * \return value of the integral representation
     */

    // Copy variables to the stack
    cusfloat k0     = wave_data.k0;
    cusfloat nu     = wave_data.nu;
    cusfloat k0nu   = wave_data.k0nu;

    // Allocate space for local variables
    cusfloat A[N];
    cusfloat B0[N];
    cusfloat B1[N];
    cusfloat res_fcn0[N];
    cusfloat res_fcn0_da[N];
    cusfloat res_fcn0_db[N];
    cusfloat res_fcn1_blt1[N];
    cusfloat res_fcn1_blt1_da[N];
    cusfloat res_fcn1_blt1_db[N];
    cusfloat res_fcn1_bgt1[N];
    cusfloat res_fcn1_bgt1_da[N];
    cusfloat res_fcn1_bgt1_db[N];
    cusfloat res_fxy_bgt1[N];
    cusfloat res_fxy_bgt1_dx[N];
    cusfloat res_fxy_bgt1_dy[N];
    cusfloat v0[N];
    cusfloat v1[N];
    cusfloat v2[N];
    cusfloat v3[N];
    cusfloat v4[N];
    cusfloat v5[N];

    // Calculate exponential terms expansion parameters
    for ( std::size_t i=0; i<N; i++ )
    {
        // Calculate dependent parameters
        v0[i] = abs( z[i] - zeta[i] );
        v1[i] = z[i] + zeta[i] + 2*h;
        v2[i] = abs( z[i] + zeta[i] );
        v3[i] = z[i] - zeta[i] + 2*h;
        v4[i] = zeta[i] - z[i] + 2*h;
        v5[i] = z[i] + zeta[i] + 4*h;

    }

    // Calculate integrals expansion hyperparameters
    H = nu * h;
    for ( std::size_t i=0; i<N; i++ )
    {
        A[i]    = R[i] / h;
        B0[i]   = v0[i] / h;
        B1[i]   = v1[i] / h;
    }

    // Check that B0 and  B1 is in between limits
    // assert( ((B0>=0.0) && (B0<=1.0)) && "B0 is out of interval [0, 1]" );
    // assert( ((B1>=0.0) && (B1<=2.0)) && "B1 is out of interval [0, 2]" );

    for ( std::size_t i=0; i<N; i++ )
    {
        if ( ( B0[i] < 0.0 ) || ( B0[i] > 1.0 ) )
        {
            B0[i] = std::min( std::max( B0[i], 0.0 ), 1.0 );
        }
    
        if ( ( B1[i] < 0.0 ) || ( B1[i] > 2.0 ) )
        {
            B1[i] = std::min( std::max( B1[i], 0.0 ), 2.0 );
        }
    }

    // Divide data into groups G1: B1 < 1 and G2: B1 > 1.0
    std::size_t b1_lt1_count = 0;
    std::size_t b1_gt1_count = 0;

    cusfloat    A_lt1[N];
    cusfloat    A_gt1[N];
    cusfloat    b1_lt1[N];
    cusfloat    b1_gt1[N];
    std::size_t b1_lt1_pos[N];
    std::size_t b1_gt1_pos[N];
    cusfloat    X_gt1[N];
    cusfloat    Y_gt1[N];

    for ( std::size_t i=0; i<N; i++ )
    {
        if ( B1[i] > 1.0 )
        {
            A_gt1[b1_gt1_count]         = A[i];
            b1_gt1[b1_gt1_count]        = B1[i];
            b1_gt1_pos[b1_gt1_count]    = i;
            X_gt1[b1_gt1_count]         = nu * R[i];
            Y_gt1[b1_gt1_count]         = nu * v2[i];
            b1_gt1_count                += 1;
        }
        else
        {
            A_lt1[b1_lt1_count]         = A[i];
            b1_lt1[b1_lt1_count]        = B1[i];
            b1_lt1_pos[b1_lt1_count]    = i;
            b1_lt1_count                += 1;
        }
    }

    // Update Bessel Factory
    bessel_factory.calculate_cheby( b1_gt1_count, X_gt1 );

    // Calculate green function main parts depends on non-dimensional depth parameter
    if ( H > 1 )
    {
        G1_Hgt1<N>( N, A, B0, res_fcn0, res_fcn0_da, res_fcn0_db );
        
        if ( b1_gt1_count < 1 )
        {
            G1_Hgt1<N>( b1_lt1_count, A_lt1, b1_lt1, res_fcn1_blt1, res_fcn1_blt1_da, res_fcn1_blt1_db );
        }
        else if ( b1_lt1_count < 1 )
        {
            G2_Hgt1<N>( b1_gt1_count, A_gt1, b1_gt1, res_fcn1_bgt1, res_fcn1_bgt1_da, res_fcn1_bgt1_db );
        }
        else
        {
            G1_Hgt1<N>( b1_lt1_count, A_lt1, b1_lt1, res_fcn1_blt1, res_fcn1_blt1_da, res_fcn1_blt1_db );
            G2_Hgt1<N>( b1_gt1_count, A_gt1, b1_gt1, res_fcn1_bgt1, res_fcn1_bgt1_da, res_fcn1_bgt1_db );
        }
    }
    else
    {
        G1_Hlt1<N>( N, A, B0, res_fcn0, res_fcn0_da, res_fcn0_db );
        
        if ( b1_gt1_count < 1 )
        {
            G1_Hlt1<N>( b1_lt1_count, A_lt1, b1_lt1, res_fcn1_blt1, res_fcn1_blt1_da, res_fcn1_blt1_db );
        }
        else if ( b1_lt1_count < 1 )
        {
            G2_Hlt1<N>( b1_gt1_count, A_gt1, b1_gt1, res_fcn1_bgt1, res_fcn1_bgt1_da, res_fcn1_bgt1_db );
            Fxy<N>( b1_gt1_count, X_gt1, Y_gt1, res_fxy_bgt1, res_fxy_bgt1_dx, res_fxy_bgt1_dy );
        }
        else
        {
            G1_Hlt1<N>( b1_lt1_count, A_lt1, b1_lt1, res_fcn1_blt1, res_fcn1_blt1_da, res_fcn1_blt1_db );
            G2_Hlt1<N>( b1_gt1_count, A_gt1, b1_gt1, res_fcn1_bgt1, res_fcn1_bgt1_da, res_fcn1_bgt1_db );
            Fxy<N>( b1_gt1_count, X_gt1, Y_gt1, res_fxy_bgt1, res_fxy_bgt1_dx, res_fxy_bgt1_dy );
        }
    }

    // Calculate real part
    cusfloat G_real[N];
    cusfloat G_dr_real[N];
    cusfloat G_dz_real[N];

    for ( std::size_t i=0; i<N; i++ )
    {
        G_real[i]       = res_fcn0[i];
        G_real_dr[i]    = res_fcn0_da[i];
        G_real_dz[i]    = res_fcn0_db[i];
    }

    for ( std::size_t i=0; i<b1_lt1_count; i++ )
    {
        G_real[b1_lt1_pos[i]]       += res_fcn1_blt1[i];
        G_real_dr[b1_lt1_pos[i]]    += res_fcn1_blt1_da[i];
        G_real_dz[b1_lt1_pos[i]]    += res_fcn1_blt1_db[i];
    }

    for ( std::size_t i=0; i<b1_gt1_count; i++ )
    {
        G_real[b1_gt1_pos[i]]       += res_fcn1_bgt1[i];
        G_real_dr[b1_gt1_pos[i]]    += res_fcn1_bgt1_da[i];
        G_real_dz[b1_gt1_pos[i]]    += res_fcn1_bgt1_db[i];
    }

    cusfloat h2     = h * h;
    cusfloat nu2    = nu * nu;
    for ( std::size_t i=0; i<N; i++ )
    {
        G_real[i]       /= h;
        G_real_dr[i]    /= h2;
        G_real_dz[i]    /= h2;
    }

    for ( std::size_t i=0; i<b1_gt1_count; i++ )
    {
        G_real[b1_gt1_pos[i]]       += nu * res_fxy_bgt1[i];
        G_real_dr[b1_gt1_pos[i]]    += nu2 * res_fxy_bgt1_dx[i];
        G_real_dz[b1_gt1_pos[i]]    += nu2 * res_fxy_bgt1_dy[i];
    }

    // Calculate exponential terms for imaginary part
    cusfloat expsum[N];
    cusfloat expsum_dz[N];
    cusfloat c0 = 0.0, c1 = 0.0, c2 = 0.0, c3 = 0.0;
    cusfloat m0 = 0.0;
    for ( std::size_t i=0; i<N; i++ )
    {
        // Calculate coeffs to reuse them
        c0 = exp( -k0 * v2[i] );
        c1 = exp( -k0 * v3[i] );
        c2 = exp( -k0 * v4[i] );
        c3 = exp( -k0 * v5[i] );

        // Calculate exponential sumation
        m0              = c1 + c3;
        expsum[i]       = m0 + c0 + c2;
        expsum_dz[i]    = m0 + c0 * sign(z+zeta) - c2;
    }

    // Calculate beseel functions
    cusfloat j0_vec[N];
    cusfloat j1_vec[N];

    for ( std::size_t i=0; i<N; i++ )
    {
        j0_vec[i] = besselj0( k0 * R );
        j1_vec[i] = besselj1( k0 * R );
    }

    // Calculate green function values
    for ( std::size_t i=0; i<N; i++ )
    {
        G[i]    = cuscomplex( G_real[i], k0nu * expsum[i] * j0_vec[i] );
        G_dr[i] = cuscomplex( G_real_dr[i], -k0nu * expsum[i] * j1_vec[i] * k0 );
        G_dz[i] = cuscomplex( G_real_dz[i], -k0nu * expsum[i] * j0_vec[i] * k0 );
    }

    // // Calculate real part
    // cusfloat G_real = 0.0;
    // if (B1 <= 1)
    // {
    //     G_real = (G1(A, B0, H, idb) + G1(A, B1, H, idb))/h;
    // }
    // else
    // {
    //     // Calculate dependent param
    //     cusfloat X = nu*R;
    //     cusfloat Y = nu*abs(z+zeta);
    //     cusfloat g1 = G1(A, B0, H, idb);
    //     cusfloat g2 = G2(A, B1, H, idb);
    //     std::size_t count = 0;
    //     G_real = (g1 + g2)/h + nu*newcode::wave_term_inf_depth(X, Y, idb, count);
    // }

    // // Calculate imaginary part
    // cusfloat expsum = (
    //                     +exp(-k0*v2)
    //                     +exp(-k0*v3)
    //                     +exp(-k0*v4)
    //                     +exp(-k0*v5)
    //                     );
    // cusfloat G_imag = k0nu*expsum*besselj0(k0*R);

    // // Define G integral as a complex number
    // cuscomplex G = G_real + G_imag*1i;

}


cuscomplex  wave_term_fin_depth_integral_dr(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
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
    // assert( ((B0>=0.0) && (B0<=1.0)) && "B0 is out of interval [0, 1]" );
    // assert( ((B1>=0.0) && (B1<=2.0)) && "B1 is out of interval [0, 2]" );

    B0 = std::min( std::max( B0, 0.0 ), 1.0 );
    B1 = std::min( std::max( B1, 0.0 ), 2.0 );

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
        std::size_t count = 0;
        G_real = (g1 + g2)/pow2s(h) + pow2s(nu)*wave_term_inf_depth_dxndim(X, Y, idb, count);
        // G_real = (g1 + g2)/pow2s(h) + pow2s(nu)*wave_term_inf_depth_dxndim(X, Y);
    }

    // Calculate imaginary part
    cusfloat expsum = (
                        +exp(-k0*v2)
                        +exp(-k0*v3)
                        +exp(-k0*v4)
                        +exp(-k0*v5)
                        );
    cusfloat G_imag = -wave_data.k0nu*expsum*besselj1(k0*R)*k0;

    // Define G integral as a complex number
    cuscomplex G = G_real + G_imag*1i;

    return G;
}


cuscomplex  wave_term_fin_depth_integral_dz(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
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
    // assert( ((B0>=0.0) && (B0<=1.0)) && "B0 is out of interval [0, 1]" );
    // assert( ((B1>=0.0) && (B1<=2.0)) && "B1 is out of interval [0, 2]" );

    B0 = std::min( std::max( B0, 0.0 ), 1.0 );
    B1 = std::min( std::max( B1, 0.0 ), 2.0 );

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
        G_real = (g1 + g2)/pow2s(h) + pow2s(nu)*wave_term_inf_depth_dyndim(X, Y, idb)*sign(z+zeta);
    }

    // Calculate imaginary part
    cusfloat expsum = (
                        +exp(-k0*v2)*sign(z+zeta)
                        +exp(-k0*v3)
                        -exp(-k0*v4)
                        +exp(-k0*v5)
                        );
    cusfloat G_imag = -k0*wave_data.k0nu*expsum*besselj0(k0*R);

    // Define G integral as a complex number
    cuscomplex G = G_real + G_imag*1i;

    return G;
}