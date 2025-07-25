
// Include general usage libraries
#include <cassert>
#include <iostream>
#include <tuple>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "../config.hpp"
#include "../math/math_interface.hpp"
#include "../math/math_tools.hpp"
#include "../math/special_math.hpp"
#include "pulsating_fin_depth_v2.hpp"
#include "chebyshev_evaluator_interface.hpp"

// Include namespaces
using namespace std::literals::complex_literals;


template<std::size_t N, int mode_loop, int mode_f, int mode_dfdr, int mode_dfdz>
void         Fxy(
                                                const std::size_t           n,
                                                cusfloat*                   X,
                                                cusfloat*                   Y,
                                                BesselFactoryVecUpTo<N>&    bessel_factory,
                                                cusfloat*                   results,
                                                cusfloat*                   results_dx,
                                                cusfloat*                   results_dy
                )
{
    // Clear vector to avoid spurious remaining data
    STATIC_COND( ONLY_FCN,      STATIC_CLEAR( n, N, results     ) )
    STATIC_COND( ONLY_FCNDR,    STATIC_CLEAR( n, N, results_dx  ) )
    STATIC_COND( ONLY_FCNDZ,    STATIC_CLEAR( n, N, results_dy  ) )

    // Calculate Bessel functions
    bessel_factory.calculate_cheby( n, X );

    // Calculate auxiliar variables
    cusfloat SQRT[N];
    cusfloat EXPY[N];
    cusfloat XINV[N];

    STATIC_LOOP( n, N, SQRT[i] = pow2s( X[i] ) + pow2s( Y[i] ); )
    STATIC_LOOP( n, N, EXPY[i] = -Y[i]; )
    STATIC_LOOP( n, N, XINV[i] = 1.0 / X[i]; )

    lv_sqrt<cusfloat>( n, SQRT, SQRT );
    lv_exp<cusfloat>( n, EXPY, EXPY );

    // Calculate Chebyshev expansions
    STATIC_COND( ONLY_FCN or ONLY_FCNDZ,        R11CEV<N>::evaluate( n, X, Y, results );        )
    STATIC_COND( ONLY_FCNDR,                    R11_dXCEV<N>::evaluate( n, X, Y, results_dx );  )

    // Add Bessel functions contribution
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCN or ONLY_FCNDZ,    results[i]      *= -2.0;           ) )
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDR,                results_dx[i]   *= -2.0 * XINV[i]; ) )
    
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCN or ONLY_FCNDZ,    results[i]      -= PI * EXPY[i] * ( bessel_factory.struve0[i] + bessel_factory.y0[i] ); ) )
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDR,                results_dx[i]   += 2.0 * XINV[i] * Y[i] / SQRT[i];                                      ) )

    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDR,                results_dx[i]   += PI * EXPY[i] * ( bessel_factory.y1[i] + bessel_factory.struve1[i] - 2.0 / PI );   ) )
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDZ,                results_dy[i]    = results[i] - 2.0 / SQRT[i];                                                       ) )

}


template<std::size_t N, int mode_loop, int mode_f, int mode_dfdr, int mode_dfdz>
void         G1_Hlt1(
                                                const std::size_t   n,
                                                cusfloat*           A,
                                                cusfloat*           B,
                                                cusfloat*           results,
                                                cusfloat*           results_da,
                                                cusfloat*           results_db
                )
{
    // Clear vector to avoid spurious remaining data
    STATIC_COND( ONLY_FCN,      STATIC_CLEAR( n, N, results     ) )
    STATIC_COND( ONLY_FCNDR,    STATIC_CLEAR( n, N, results_da  ) )
    STATIC_COND( ONLY_FCNDZ,    STATIC_CLEAR( n, N, results_db  ) )

    // Calculate L1 integral using chebyshev expansions
    STATIC_COND( ONLY_FCN,      L1CEV<N>::evaluate( n, A, B, results );         )
    STATIC_COND( ONLY_FCNDR,    L1_dACEV<N>::evaluate( n, A, B, results_da );   )
    STATIC_COND( ONLY_FCNDZ,    L1_dBCEV<N>::evaluate( n, A, B, results_db );   )

    // Add L2 integral scalar contribution
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCN, results[i] += ChebyshevTraits<L2C>::coeffs; ) )

}


template<std::size_t N, int mode_loop, int mode_f, int mode_dfdr, int mode_dfdz>
void         G1_Hgt1(
                                                const std::size_t   n,
                                                cusfloat*           A,
                                                cusfloat*           B,
                                                cusfloat*           results,
                                                cusfloat*           results_da,
                                                cusfloat*           results_db
                )
{
    // Clear vector to avoid spurious remaining data
    STATIC_COND( ONLY_FCN,      STATIC_CLEAR( n, N, results     ) )
    STATIC_COND( ONLY_FCNDR,    STATIC_CLEAR( n, N, results_da  ) )
    STATIC_COND( ONLY_FCNDZ,    STATIC_CLEAR( n, N, results_db  ) )

    // Calculate auxiliar values
    cusfloat A2[N];
    cusfloat BP[N];
    cusfloat BM[N];
    cusfloat BPM[N];
    cusfloat BMM[N];

    STATIC_LOOP( n, N, A2[i] = POW2S( A[i] );           )
    STATIC_LOOP( n, N, BP[i] = A2[i] + POW2S( 2+B[i] ); )
    STATIC_LOOP( n, N, BM[i] = A2[i] + POW2S( 2-B[i] ); )

    lv_sqrt<cusfloat>( n, BP, BPM );
    lv_sqrt<cusfloat>( n, BM, BMM );
    
    // Calculate L3 using chebyshev expansion
    STATIC_COND( ONLY_FCN,      L3CEV<N>::evaluate( n, A, B, results );         )
    STATIC_COND( ONLY_FCNDR,    L3_dACEV<N>::evaluate( n, A, B, results_da );   )
    STATIC_COND( ONLY_FCNDZ,    L3_dBCEV<N>::evaluate( n, A, B, results_db );   )

    // Substract auxiliar variables
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCN,      results[i]      -= 2.0 / BPM[i];                            ) )
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDR,    results_da[i]   += 2.0 * A[i] / BPM[i] / BP[i];             ) )
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDZ,    results_db[i]   += 2.0 * ( 2.0 + B[i] ) / BPM[i] / BP[i];   ) )

    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCN,      results[i]      -= 2.0 / BMM[i];                            ) )
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDR,    results_da[i]   += 2.0 * A[i] / BMM[i] / BM[i];             ) )
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDZ,    results_db[i]   -= 2.0 * ( 2.0 - B[i] ) / BMM[i] / BM[i];   ) )
}


template<std::size_t N, int mode_loop, int mode_f, int mode_dfdr, int mode_dfdz>
void         G2_Hlt1(
                                                const std::size_t   n,
                                                cusfloat*           A,
                                                cusfloat*           B,
                                                cusfloat*           results,
                                                cusfloat*           results_da,
                                                cusfloat*           results_db
                )
{
    // Clear vector to avoid spurious remaining data
    STATIC_COND( ONLY_FCN,      STATIC_CLEAR( n, N, results     ) )
    STATIC_COND( ONLY_FCNDR,    STATIC_CLEAR( n, N, results_da  ) )
    STATIC_COND( ONLY_FCNDZ,    STATIC_CLEAR( n, N, results_db  ) )

    // Calculate M1 integral using chebyshev expansions
    STATIC_COND( ONLY_FCN,      M1CEV<N>::evaluate( n, A, B, results );         )
    STATIC_COND( ONLY_FCNDR,    M1_dACEV<N>::evaluate( n, A, B, results_da );   )
    STATIC_COND( ONLY_FCNDZ,    M1_dBCEV<N>::evaluate( n, A, B, results_db );   )

    // Add M2 integral scalar contribution
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCN,  results[i] += ChebyshevTraits<M2C>::coeffs; ) )
}


template<std::size_t N, int mode_loop, int mode_f, int mode_dfdr, int mode_dfdz>
void         G2_Hgt1(
                                                const std::size_t   n,
                                                cusfloat*           A,
                                                cusfloat*           B,
                                                cusfloat*           results,
                                                cusfloat*           results_da,
                                                cusfloat*           results_db
                )
{
    // Clear vector to avoid spurious remaining data
    STATIC_COND( ONLY_FCN,      STATIC_CLEAR( n, N, results     ) )
    STATIC_COND( ONLY_FCNDR,    STATIC_CLEAR( n, N, results_da  ) )
    STATIC_COND( ONLY_FCNDZ,    STATIC_CLEAR( n, N, results_db  ) )

    // Calculate auxiliar values
    cusfloat A2[N];
    cusfloat BP[N];
    cusfloat BPM[N];

    STATIC_LOOP( n, N, A2[i] = POW2S( A[i] );           )
    STATIC_LOOP( n, N, BP[i] = A2[i] + POW2S( 2+B[i] ); )

    lv_sqrt<cusfloat>( n, BP, BPM );
    
    // Calculate L3 using chebyshev expansion
    STATIC_COND( ONLY_FCN,      M3CEV<N>::evaluate( n, A, B, results );          )
    STATIC_COND( ONLY_FCNDR,    M3_dACEV<N>::evaluate( n, A, B, results_da );    )
    STATIC_COND( ONLY_FCNDZ,    M3_dBCEV<N>::evaluate( n, A, B, results_db );    )

    // Substract auxiliar variables
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCN,      results[i]      -= 2.0 / BPM[i];                            ) )
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDR,    results_da[i]   += 2.0 * A[i] / BPM[i] / BP[i];             ) )
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDZ,    results_db[i]   += 2.0 * ( 2.0 + B[i] ) / BPM[i] / BP[i];   ) )

}


template<std::size_t N, int mode_f, int mode_dfdr, int mode_dfdz>
void        john_series(
                                                cusfloat*               R,
                                                cusfloat*               z,
                                                cusfloat*               zeta,
                                                cusfloat                h,
                                                BesselFactoryVecUpTo<N> &bessel_factory,
                                                WaveDispersionFONK      &wave_data,
                                                cuscomplex*             G,
                                                cuscomplex*             G_dr,
                                                cuscomplex*             G_dz
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
    cusfloat ccz[N];
    cusfloat scz[N];
    cusfloat ci_g[N];
    cusfloat ci_g_dr[N];
    cusfloat ci_g_dz[N];
    cusfloat c_real[N];
    cusfloat c_real_dz[N];
    cusfloat c_g_imag[N];
    cusfloat c_g_dr_imag[N];
    cusfloat c_g_dz_imag[N];
    cusfloat expsum[N];
    cusfloat expsum_dz[N];
    cusfloat k0r[N];
    cusfloat knir[N];
    cusfloat knizetah[N];
    cusfloat knizh[N];
    cusfloat v3[N];
    cusfloat v4[N];
    cusfloat v5[N];
    cusfloat v6[N];
    cusfloat zh[N];
    cusfloat zetah[N];

    STATIC_COND( ONLY_FCN,      ( clear_vector<cusfloat,N>( c_g_imag )    );      )
    STATIC_COND( ONLY_FCNDR,    ( clear_vector<cusfloat,N>( c_g_dr_imag ) );      )
    STATIC_COND( ONLY_FCNDZ,    ( clear_vector<cusfloat,N>( c_g_dz_imag ) );      )

    for ( std::size_t i=0; i<N; i++ )
    {
        v3[i] = z[i] + zeta[i];
        v4[i] = z[i] - zeta[i] + 2.0 * h;
        v5[i] = zeta[i] - z[i] + 2.0 * h;
        v6[i] = z[i] + zeta[i] + 4.0 * h;
    }

    // Calcuate real root series part
    for ( std::size_t i=0; i<N; i++ )
    {
        k0r[i] = k0 * R[i];
    }

    bessel_factory.calculate_cheby( N, k0r );

    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( 
                        ONLY_FCN or ONLY_FCNDR,  
                        expsum[i]       = (
                                            + exp( k0 * v3[i] )
                                            + exp( -k0 * v4[i] )
                                            + exp( -k0 * v5[i] )
                                            + exp( -k0 * v6[i] )
                                        );
                    )
        
        STATIC_COND( 
                        ONLY_FCNDZ, 
                        expsum_dz[i]    = - (
                                            - exp( k0 * v3[i] )
                                            + exp( -k0 * v4[i] )
                                            - exp( -k0 * v5[i] )
                                            + exp( -k0 * v6[i] )
                                        ) * k0;
                    )
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( ONLY_FCN or ONLY_FCNDR,    c_real[i]       = - wave_data.k0nu * expsum[i]; )
        STATIC_COND( ONLY_FCNDZ,                c_real_dz[i]    = - wave_data.k0nu * expsum_dz[i]; )
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( ONLY_FCN,      G[i]    = + c_real[i] * ( bessel_factory.y0[i] - 1i * bessel_factory.j0[i] );       )
        STATIC_COND( ONLY_FCNDR,    G_dr[i] = - c_real[i] * k0 * ( bessel_factory.y1[i] - 1i * bessel_factory.j1[i] );  )
        STATIC_COND( ONLY_FCNDZ,    G_dz[i] = + c_real_dz[i] * ( bessel_factory.y0[i] - 1i * bessel_factory.j0[i] );    )
    }

    // Calculate imag root series part
    cusfloat    c_max       = 0.0;
    cusfloat    c_g_max     = 0.0;
    cusfloat    c_g_dr_max  = 0.0;
    cusfloat    c_g_dz_max  = 0.0;
    int         count_k     = 0;
    cusfloat    kni         = 0.0;

    while (true)
    {
        // Get local copy of th imaginary wave number
        kni = wave_data.kn[count_k];

        // Calculate local loop variables
        for ( std::size_t i=0; i<N; i++ )
        {
            zetah[i]    = zeta[i] + h;
            zh[i]       = z[i] + h;
        }

        for ( std::size_t i=0; i<N; i++ )
        {
            knir[i]     = kni * R[i];
            knizetah[i] = kni * zetah[i];
            knizh[i]    = kni * zh[i];
        }

        for ( std::size_t i=0; i<N; i++ )
        {
            ccz[i]     = 4.0 * wave_data.knnu[count_k] * cos( knizetah[i] );
            scz[i]     = ccz[i];
        }

        for ( std::size_t i=0; i<N; i++ )
        {
            STATIC_COND( ONLY_FCN or ONLY_FCNDR,    ccz[i]     = cos( knizh[i] ) * ccz[i]; )
            STATIC_COND( ONLY_FCNDZ,                scz[i]     = sin( knizh[i] ) * scz[i]; )
        }

        // Calculate Bessel function for kni*R
        bessel_factory.calculate_modified( N, knir );

        // Calculate i term of the series
        for ( std::size_t i=0; i<N; i++ )
        {
            STATIC_COND( ONLY_FCN,      ci_g[i]     = + ccz[i] * bessel_factory.k0[i];          )
            STATIC_COND( ONLY_FCNDR,    ci_g_dr[i]  = - ccz[i] * kni * bessel_factory.k1[i];    )
            STATIC_COND( ONLY_FCNDZ,    ci_g_dz[i]  = - scz[i] * kni * bessel_factory.k0[i];    )
        }

        for ( std::size_t i=0; i<N; i++ )
        {
            STATIC_COND( ONLY_FCN,      c_g_imag[i]     += ci_g[i];     )
            STATIC_COND( ONLY_FCNDR,    c_g_dr_imag[i]  += ci_g_dr[i];  )
            STATIC_COND( ONLY_FCNDZ,    c_g_dz_imag[i]  += ci_g_dz[i];  )
        }

        // std::cout << "knnu: " << wave_data.knnu[count_k] << " - knir: " << knir[0] << " - ccz: " << ccz[0] << " - scz: " << scz[0] << " - besselk0(kni*R): " << bessel_factory.k0[0] << " - G_imag: " << c_g_imag[0] << std::endl;

        // Check for convergence
        for ( std::size_t i=0; i<N; i++ )
        {
            STATIC_COND( ONLY_FCN,      ci_g[i]         = std::abs( ci_g[i] );      )
            STATIC_COND( ONLY_FCNDR,    ci_g_dr[i]      = std::abs( ci_g_dr[i] );   )
            STATIC_COND( ONLY_FCNDZ,    ci_g_dz[i]      = std::abs( ci_g_dz[i] );   )
        }

        c_g_max     = 0.0;
        c_g_dr_max  = 0.0;
        c_g_dz_max  = 0.0;
        for ( std::size_t i=0; i<N; i++ )
        {
            STATIC_COND( ONLY_FCN,      c_g_max         = std::max( c_g_max, ci_g[i] );         )
            STATIC_COND( ONLY_FCNDR,    c_g_dr_max      = std::max( c_g_dr_max, ci_g_dr[i] );   )
            STATIC_COND( ONLY_FCNDZ,    c_g_dz_max      = std::max( c_g_dz_max, ci_g_dz[i] );   )
        }

        c_max = std::max( c_g_max, std::max( c_g_dr_max, c_g_dz_max ) );
        if ( c_max < 1e-6 )
        {
            break;
        }

        // Check for the limit in imaginary roots
        #ifdef DEBUG_BUILD
        if (count_k > (wave_data.num_kn-2))
        {
            std::cerr << "Jonh series could not converge up to the precision required with" << std::endl;
            std::cerr << "Input parameters:"  << std::endl;
            std::cerr << "  - R/h: " << R[0]/h << std::endl;
            std::cerr << "  - R: " << R[0] << std::endl;
            std::cerr << "  - z: " << z[0] << std::endl;
            std::cerr << "  - zeta: " << zeta[0] << std::endl;
            std::cerr << "  - h: " << h << std::endl;
            std::cerr << "  - nu: " << wave_data.nu << std::endl;
            std::cerr << "  - k0: " << k0 << std::endl;
            std::cerr << "  - num_kn: " << wave_data.num_kn << std::endl;
            std::cerr << "* Current integral value is: " << c_g_imag[0] << " - " << c_g_dr_imag[0] << " - " << c_g_dz_imag[0] << std::endl;
            std::cerr << "* Current term value is: " << ci_g[0] << " - " << ci_g_dr[0] << " - " << ci_g_dz[0] << std::endl;
            throw std::runtime_error("Jonh series value could not converge. See log file for details.");
        }
        #endif

        // Update counter
        count_k++;
    }

    // Add series values
    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( ONLY_FCN,      G[i]    += c_g_imag[i];     )
        STATIC_COND( ONLY_FCNDR,    G_dr[i] += c_g_dr_imag[i];  )
        STATIC_COND( ONLY_FCNDZ,    G_dr[i] += c_g_dz_imag[i];  )
    }

}


template<std::size_t N>
void        wave_term_fin_depth(
                                                cusfloat*               R,
                                                cusfloat*               z,
                                                cusfloat*               zeta,
                                                cusfloat                h,
                                                BesselFactoryVecUpTo<N> &bessel_factory,
                                                WaveDispersionFONK      &wave_data,
                                                cuscomplex*             G,
                                                cuscomplex*             G_dr,
                                                cuscomplex*             G_dz
                                )
{
    // Divide data into groups rh < 1.0 (Chebyshev polynomials) and rh > 1.0 (John series)
    std::size_t rh_lt1_count = 0;
    std::size_t rh_gt1_count = 0;

    cuscomplex  G_lt1[N];
    cuscomplex  G_gt1[N];
    cuscomplex  G_dr_lt1[N];
    cuscomplex  G_dr_gt1[N];
    cuscomplex  G_dz_lt1[N];
    cuscomplex  G_dz_gt1[N];
    cusfloat    R_lt1[N];
    cusfloat    R_gt1[N];
    cusfloat    rh[N];
    std::size_t rh_lt1_pos[N];
    std::size_t rh_gt1_pos[N];
    cusfloat    z_lt1[N];
    cusfloat    z_gt1[N];
    cusfloat    zeta_lt1[N];
    cusfloat    zeta_gt1[N];

    for ( std::size_t i=0; i<N; i++ )
    {
        rh[i] = R[i] / h;
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        if ( rh[i] > 1.0 )
        {
            R_gt1[rh_gt1_count]         = R[i];
            z_gt1[rh_gt1_count]         = z[i];
            zeta_gt1[rh_gt1_count]      = zeta[i];
            rh_gt1_pos[rh_gt1_count]    = i;
            rh_gt1_count                += 1;
        }
        else
        {
            R_lt1[rh_lt1_count]         = R[i];
            z_lt1[rh_lt1_count]         = z[i];
            zeta_lt1[rh_lt1_count]      = zeta[i];
            rh_lt1_pos[rh_lt1_count]    = i;
            rh_lt1_count                += 1;
        }
    }

    // Calculate integrals
    if ( rh_gt1_count < 1 )
    {
        wave_term_fin_depth_integral<N>( N, R, z, zeta, h, bessel_factory, wave_data, G, G_dr, G_dz );
    }
    else if ( rh_lt1_count < 1 )
    {
        john_series<N>( N, R, z, zeta, h, bessel_factory, wave_data, G, G_dr, G_dz );
    }
    else
    {
        // Calculate cases where rh < 1.0
        wave_term_fin_depth_integral<N>( 
                                            rh_lt1_count, 
                                            R_lt1, 
                                            z_lt1, 
                                            zeta_lt1, 
                                            h, 
                                            bessel_factory, 
                                            wave_data, 
                                            G_lt1, 
                                            G_dr_lt1, 
                                            G_dz_lt1
                                        );
        
        // Calculate cases where rh > 1.0
        john_series<N>( 
                            rh_gt1_count, 
                            R_gt1, 
                            z_gt1, 
                            zeta_gt1, 
                            h, 
                            bessel_factory, 
                            wave_data, 
                            G_gt1, 
                            G_dr_gt1, 
                            G_dz_gt1
                        );

        // Convine results
        for ( std::size_t i=0; i<rh_lt1_count; i++ )
        {
            G[rh_lt1_pos[i]]       += G_lt1[i];
            G_dr[rh_lt1_pos[i]]    += G_dr_lt1[i];
            G_dz[rh_lt1_pos[i]]    += G_dz_lt1[i];
        }

        for ( std::size_t i=0; i<rh_gt1_count; i++ )
        {
            G[rh_gt1_pos[i]]       += G_gt1[i];
            G_dr[rh_gt1_pos[i]]    += G_dr_gt1[i];
            G_dz[rh_gt1_pos[i]]    += G_dz_gt1[i];
        }
    }
}


template<std::size_t N>
void        wave_term_fin_depth_integral(
                                                const std::size_t           n,
                                                cusfloat*                   R,
                                                cusfloat*                   z,
                                                cusfloat*                   zeta,
                                                cusfloat                    h,
                                                BesselFactoryVecUpTo<N>     &bessel_factory,
                                                WaveDispersionFONK          &wave_data,
                                                cuscomplex*                 G,
                                                cuscomplex*                 G_dr,
                                                cuscomplex*                 G_dz
                                        )
{
    /**
     * @brief Wave term for the finite water depth function (Integral representation).
     * 
     * The formulation used here (valid for R/h<1.0) is mainly taken from:
     * "Consistent expressions for the free surface function in
     *  finite water depth - Ed Mackay".
     * 
     * \param n Number of points to compute
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
    for ( std::size_t i=0; i<n; i++ )
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
    cusfloat H = nu * h;
    for ( std::size_t i=0; i<n; i++ )
    {
        A[i]    = R[i] / h;
        B0[i]   = v0[i] / h;
        B1[i]   = v1[i] / h;
    }

    // Check that B0 and  B1 is in between limits
    // assert( ((B0>=0.0) && (B0<=1.0)) && "B0 is out of interval [0, 1]" );
    // assert( ((B1>=0.0) && (B1<=2.0)) && "B1 is out of interval [0, 2]" );

    for ( std::size_t i=0; i<n; i++ )
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

    for ( std::size_t i=0; i<n; i++ )
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
        G1_Hgt1<N>( n, A, B0, res_fcn0, res_fcn0_da, res_fcn0_db );
        
        if ( b1_gt1_count < 1 )
        {
            G1_Hgt1<N>( n, A_lt1, b1_lt1, res_fcn1_blt1, res_fcn1_blt1_da, res_fcn1_blt1_db );
        }
        else if ( b1_lt1_count < 1 )
        {
            G2_Hgt1<N>( n, A_gt1, b1_gt1, res_fcn1_bgt1, res_fcn1_bgt1_da, res_fcn1_bgt1_db );
            Fxy<N>( n, X_gt1, Y_gt1, bessel_factory, res_fxy_bgt1, res_fxy_bgt1_dx, res_fxy_bgt1_dy );
        }
        else
        {
            G1_Hgt1<N>( b1_lt1_count, A_lt1, b1_lt1, res_fcn1_blt1, res_fcn1_blt1_da, res_fcn1_blt1_db );
            G2_Hgt1<N>( b1_gt1_count, A_gt1, b1_gt1, res_fcn1_bgt1, res_fcn1_bgt1_da, res_fcn1_bgt1_db );
            Fxy<N>( b1_gt1_count, X_gt1, Y_gt1, bessel_factory, res_fxy_bgt1, res_fxy_bgt1_dx, res_fxy_bgt1_dy );
        }
    }
    else
    {
        G1_Hlt1<N>( n, A, B0, res_fcn0, res_fcn0_da, res_fcn0_db );
        
        if ( b1_gt1_count < 1 )
        {
            G1_Hlt1<N>( n, A_lt1, b1_lt1, res_fcn1_blt1, res_fcn1_blt1_da, res_fcn1_blt1_db );
        }
        else if ( b1_lt1_count < 1 )
        {
            G2_Hlt1<N>( n, A_gt1, b1_gt1, res_fcn1_bgt1, res_fcn1_bgt1_da, res_fcn1_bgt1_db );
            Fxy<N>( n, X_gt1, Y_gt1, bessel_factory, res_fxy_bgt1, res_fxy_bgt1_dx, res_fxy_bgt1_dy );
        }
        else
        {
            G1_Hlt1<N>( b1_lt1_count, A_lt1, b1_lt1, res_fcn1_blt1, res_fcn1_blt1_da, res_fcn1_blt1_db );
            G2_Hlt1<N>( b1_gt1_count, A_gt1, b1_gt1, res_fcn1_bgt1, res_fcn1_bgt1_da, res_fcn1_bgt1_db );
            Fxy<N>( b1_gt1_count, X_gt1, Y_gt1, bessel_factory, res_fxy_bgt1, res_fxy_bgt1_dx, res_fxy_bgt1_dy );
        }
    }

    // Calculate real part
    cusfloat G_real[N];
    cusfloat G_dr_real[N];
    cusfloat G_dz_real[N];

    for ( std::size_t i=0; i<n; i++ )
    {
        G_real[i]       = res_fcn0[i];
        G_dr_real[i]    = res_fcn0_da[i];
        G_dz_real[i]    = res_fcn0_db[i];
    }

    for ( std::size_t i=0; i<b1_lt1_count; i++ )
    {
        G_real[b1_lt1_pos[i]]       += res_fcn1_blt1[i];
        G_dr_real[b1_lt1_pos[i]]    += res_fcn1_blt1_da[i];
        G_dz_real[b1_lt1_pos[i]]    += res_fcn1_blt1_db[i];
    }

    for ( std::size_t i=0; i<b1_gt1_count; i++ )
    {
        G_real[b1_gt1_pos[i]]       += res_fcn1_bgt1[i];
        G_dr_real[b1_gt1_pos[i]]    += res_fcn1_bgt1_da[i];
        G_dz_real[b1_gt1_pos[i]]    += res_fcn1_bgt1_db[i];
    }

    cusfloat h2     = h * h;
    cusfloat nu2    = nu * nu;
    for ( std::size_t i=0; i<n; i++ )
    {
        G_real[i]       /= h;
        G_dr_real[i]    /= h2;
        G_dz_real[i]    /= h2;
    }

    for ( std::size_t i=0; i<b1_gt1_count; i++ )
    {
        G_real[b1_gt1_pos[i]]       += nu * res_fxy_bgt1[i];
        G_dr_real[b1_gt1_pos[i]]    += nu2 * res_fxy_bgt1_dx[i];
        G_dz_real[b1_gt1_pos[i]]    += nu2 * res_fxy_bgt1_dy[i];
    }

    // Calculate exponential terms for imaginary part
    cusfloat expsum[N];
    cusfloat expsum_dz[N];
    cusfloat c0 = 0.0, c1 = 0.0, c2 = 0.0, c3 = 0.0;
    cusfloat m0 = 0.0;
    for ( std::size_t i=0; i<n; i++ )
    {
        // Calculate coeffs to reuse them
        c0 = exp( -k0 * v2[i] );
        c1 = exp( -k0 * v3[i] );
        c2 = exp( -k0 * v4[i] );
        c3 = exp( -k0 * v5[i] );

        // Calculate exponential sumation
        m0              = c1 + c3;
        expsum[i]       = m0 + c0 + c2;
        expsum_dz[i]    = m0 + c0 * sign( z[i] + zeta[i] ) - c2;
    }

    // Calculate beseel functions
    cusfloat j0_vec[N];
    cusfloat j1_vec[N];

    for ( std::size_t i=0; i<n; i++ )
    {
        j0_vec[i] = besselj0( k0 * R[i] );
        j1_vec[i] = besselj1( k0 * R[i] );
    }

    // Calculate green function values
    for ( std::size_t i=0; i<n; i++ )
    {
        G[i]    = cuscomplex( G_real[i], k0nu * expsum[i] * j0_vec[i] );
        G_dr[i] = cuscomplex( G_dr_real[i], -k0nu * expsum[i] * j1_vec[i] * k0 );
        G_dz[i] = cuscomplex( G_dz_real[i], -k0nu * expsum[i] * j0_vec[i] * k0 );
    }

}


template<std::size_t N, int mode_f, int mode_dfdr, int mode_dfdz>
void        wave_term_integral(
                                                cusfloat*                   R,
                                                cusfloat*                   z,
                                                cusfloat*                   zeta,
                                                cusfloat                    h,
                                                BesselFactoryVecUpTo<N>     &bessel_factory,
                                                WaveDispersionFONK          &wave_data,
                                                cuscomplex*                 G,
                                                cuscomplex*                 G_dr,
                                                cuscomplex*                 G_dz
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
    cusfloat H = nu * h;
    for ( std::size_t i=0; i<N; i++ )
    {
        A[i]    = R[i] / h;
        B0[i]   = v0[i] / h;
        B1[i]   = v1[i] / h;
    }

    // Check that B0 and  B1 is in between limits
    for ( std::size_t i=0; i<N; i++ )
    {
        // Bound B0 to the interval [0, 1]" );
        if ( ( B0[i] < 0.0 ) || ( B0[i] > 1.0 ) )
        {
            B0[i] = std::min( std::max( B0[i], 0.0 ), 1.0 );
        }
        
        // Bound B1 to the interval [0, 2]" );
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
        G1_Hgt1<N, STATIC_LOOP_ON, mode_f, mode_dfdr, mode_dfdz>( N, A, B0, res_fcn0, res_fcn0_da, res_fcn0_db );
        
        if ( b1_gt1_count < 1 )
        {
            G1_Hgt1<N, STATIC_LOOP_ON, mode_f, mode_dfdr, mode_dfdz>( N, A_lt1, b1_lt1, res_fcn1_blt1, res_fcn1_blt1_da, res_fcn1_blt1_db );
        }
        else if ( b1_lt1_count < 1 )
        {
            G2_Hgt1<N, STATIC_LOOP_ON, mode_f, mode_dfdr, mode_dfdz>( N, A_gt1, b1_gt1, res_fcn1_bgt1, res_fcn1_bgt1_da, res_fcn1_bgt1_db );
            Fxy<N, STATIC_LOOP_ON, mode_f, mode_dfdr, mode_dfdz>( N, X_gt1, Y_gt1, bessel_factory, res_fxy_bgt1, res_fxy_bgt1_dx, res_fxy_bgt1_dy );
        }
        else
        {
            G1_Hgt1<N, STATIC_LOOP_OFF, mode_f, mode_dfdr, mode_dfdz>( b1_lt1_count, A_lt1, b1_lt1, res_fcn1_blt1, res_fcn1_blt1_da, res_fcn1_blt1_db );
            G2_Hgt1<N, STATIC_LOOP_OFF, mode_f, mode_dfdr, mode_dfdz>( b1_gt1_count, A_gt1, b1_gt1, res_fcn1_bgt1, res_fcn1_bgt1_da, res_fcn1_bgt1_db );
            Fxy<N, STATIC_LOOP_OFF, mode_f, mode_dfdr, mode_dfdz>( b1_gt1_count, X_gt1, Y_gt1, bessel_factory, res_fxy_bgt1, res_fxy_bgt1_dx, res_fxy_bgt1_dy );
        }
    }
    else
    {
        G1_Hlt1<N, STATIC_LOOP_ON, mode_f, mode_dfdr, mode_dfdz>( N, A, B0, res_fcn0, res_fcn0_da, res_fcn0_db );
        
        if ( b1_gt1_count < 1 )
        {
            G1_Hlt1<N, STATIC_LOOP_ON, mode_f, mode_dfdr, mode_dfdz>( N, A_lt1, b1_lt1, res_fcn1_blt1, res_fcn1_blt1_da, res_fcn1_blt1_db );
        }
        else if ( b1_lt1_count < 1 )
        {
            G2_Hlt1<N, STATIC_LOOP_ON, mode_f, mode_dfdr, mode_dfdz>( N, A_gt1, b1_gt1, res_fcn1_bgt1, res_fcn1_bgt1_da, res_fcn1_bgt1_db );
            Fxy<N, STATIC_LOOP_ON, mode_f, mode_dfdr, mode_dfdz>( N, X_gt1, Y_gt1, bessel_factory, res_fxy_bgt1, res_fxy_bgt1_dx, res_fxy_bgt1_dy );
        }
        else
        {
            G1_Hlt1<N, STATIC_LOOP_OFF, mode_f, mode_dfdr, mode_dfdz>( b1_lt1_count, A_lt1, b1_lt1, res_fcn1_blt1, res_fcn1_blt1_da, res_fcn1_blt1_db );
            G2_Hlt1<N, STATIC_LOOP_OFF, mode_f, mode_dfdr, mode_dfdz>( b1_gt1_count, A_gt1, b1_gt1, res_fcn1_bgt1, res_fcn1_bgt1_da, res_fcn1_bgt1_db );
            Fxy<N, STATIC_LOOP_OFF, mode_f, mode_dfdr, mode_dfdz>( b1_gt1_count, X_gt1, Y_gt1, bessel_factory, res_fxy_bgt1, res_fxy_bgt1_dx, res_fxy_bgt1_dy );
        }
    }

    // Calculate real part
    cusfloat G_real[N];
    cusfloat G_dr_real[N];
    cusfloat G_dz_real[N];

    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( ONLY_FCN,      G_real[i]       = res_fcn0[i];      )
        STATIC_COND( ONLY_FCNDR,    G_dr_real[i]    = res_fcn0_da[i];   )
        STATIC_COND( ONLY_FCNDZ,    G_dz_real[i]    = res_fcn0_db[i];   )
    }

    for ( std::size_t i=0; i<b1_lt1_count; i++ )
    {
        STATIC_COND( ONLY_FCN,      G_real[b1_lt1_pos[i]]       += res_fcn1_blt1[i];    )
        STATIC_COND( ONLY_FCNDR,    G_dr_real[b1_lt1_pos[i]]    += res_fcn1_blt1_da[i]; )
        STATIC_COND( ONLY_FCNDZ,    G_dz_real[b1_lt1_pos[i]]    += res_fcn1_blt1_db[i]; )
    }

    for ( std::size_t i=0; i<b1_gt1_count; i++ )
    {
        STATIC_COND( ONLY_FCN,      G_real[b1_gt1_pos[i]]       += res_fcn1_bgt1[i];    )
        STATIC_COND( ONLY_FCNDR,    G_dr_real[b1_gt1_pos[i]]    += res_fcn1_bgt1_da[i]; )
        STATIC_COND( ONLY_FCNDZ,    G_dz_real[b1_gt1_pos[i]]    += res_fcn1_bgt1_db[i]; )
    }

    cusfloat h2  = h * h;
    cusfloat nu2 = nu * nu;
    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( ONLY_FCN,      G_real[i]       /= h;   )
        STATIC_COND( ONLY_FCNDR,    G_dr_real[i]    /= h2;  )
        STATIC_COND( ONLY_FCNDZ,    G_dz_real[i]    /= h2;  )
    }

    for ( std::size_t i=0; i<b1_gt1_count; i++ )
    {
        STATIC_COND( ONLY_FCN,      G_real[b1_gt1_pos[i]]       += nu * res_fxy_bgt1[i];        )
        STATIC_COND( ONLY_FCNDR,    G_dr_real[b1_gt1_pos[i]]    += nu2 * res_fxy_bgt1_dx[i];    )
        STATIC_COND( ONLY_FCNDZ,    G_dz_real[b1_gt1_pos[i]]    += nu2 * res_fxy_bgt1_dy[i];    )
    }

    // Calculate exponential terms for imaginary part
    cusfloat expsum[N], expsum_dz[N];
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
        m0  = c1 + c3;

        STATIC_COND( ONLY_FCN or ONLY_FCNDR,    expsum[i]       = m0 + c0 + c2;                             )
        STATIC_COND( ONLY_FCNDZ,                expsum_dz[i]    = m0 + c0 * sign( z[i] + zeta[i] ) - c2;    )
    }

    // Calculate beseel functions
    cusfloat j0_vec[N];
    cusfloat j1_vec[N];

    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( ONLY_FCN or ONLY_FCNDZ,    j0_vec[i] = besselj0( k0 * R[i] ); )
        STATIC_COND( ONLY_FCNDR,                j1_vec[i] = besselj1( k0 * R[i] ); )
    }

    // Calculate green function values
    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( ONLY_FCN, G[i]    = cuscomplex( G_real[i], k0nu * expsum[i] * j0_vec[i] ); )
        STATIC_COND( ONLY_FCNDR, G_dr[i] = cuscomplex( G_dr_real[i], -k0nu * expsum[i] * j1_vec[i] * k0 ); )
        STATIC_COND( ONLY_FCNDZ, G_dz[i] = cuscomplex( G_dz_real[i], -k0nu * expsum[i] * j0_vec[i] * k0 ); )
    }

}


template<std::size_t N>
void        custom_template(
                                cusfloat* R
                            )
{
    double a = 0.0;
}