
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
#include <cassert>
#include <iostream>
#include <tuple>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "../config.hpp"
#include "../static_tools.hpp"
#include "../math/math_interface.hpp"
#include "../math/math_tools.hpp"
#include "../math/special_math.hpp"
#include "pulsating_fin_depth.hpp"
#include "chebyshev_evaluator_interface.hpp"

// Include namespaces
using namespace std::literals::complex_literals;


template<std::size_t N, int mode_loop, int mode_f, int mode_dfdr, int mode_dfdz, int mode_fslid>
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

    // Fit in bounds X value
    cusfloat XB[N];
    cusfloat YB[N];
    STATIC_COPY( n, N, X, XB )
    STATIC_COPY( n, N, Y, YB )
    R11CEV<N, mode_loop>::check_boundaries_raw( n, XB, YB );

    // Calculate Bessel functions
    bessel_factory.calculate_cheby( n, XB );

    // Check for region location
    std::size_t r00_count = 0;
    std::size_t r00_pos[N];
    cusfloat    r00_results_dx[N];
    std::size_t r11_count = 0;
    std::size_t r11_pos[N];
    cusfloat    r11_results_dx[N];
    cusfloat    x_r00[N];
    cusfloat    x_r11[N];
    cusfloat    y_r00[N];
    cusfloat    y_r11[N];
    cusfloat    x_min = 1e-2;
    cusfloat    y_min = 10.0;
    for ( std::size_t i=0; i<N; i++ )
    {
        if ( XB[i] < x_min )
        {
            r00_pos[r00_count]  = i;
            x_r00[r00_count]    = XB[i];
            y_r00[r00_count]    = YB[i];
            r00_count++;
        }
        else
        {
            r11_pos[r11_count]  = i;
            x_r11[r11_count]    = XB[i];
            y_r11[r11_count]    = YB[i];
            r11_count++;
        }
    }

    // Calculate auxiliar variables
    cusfloat SQRT[N];
    cusfloat EXPY[N];
    cusfloat XINV[N];
    cusfloat LOG_SQRT[N];

    STATIC_LOOP( n, N, SQRT[i] = pow2s( XB[i] ) + pow2s( YB[i] ); )
    STATIC_LOOP( n, N, EXPY[i] = -YB[i]; )
    STATIC_LOOP( n, N, XINV[i] = 1.0 / XB[i]; )
    
    lv_sqrt<cusfloat>( n, SQRT, SQRT );
    lv_exp<cusfloat>( n, EXPY, EXPY );
    
    // Calculate auxiliar variables in case of lid panel
    STATIC_LOOP( n, N, STATIC_COND( IS_LID_PANEL, LOG_SQRT[i] = SQRT[i] + YB[i]; ) )

    STATIC_COND( IS_LID_PANEL, lv_log<cusfloat>( n, LOG_SQRT, LOG_SQRT ); );

    // Calculate FXY and dFdY functions
    STATIC_COND( ONLY_FCN or ONLY_FCNDZ, (R11CEV<N, mode_loop>::evaluate( n, XB, YB, results )); )

    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCN or ONLY_FCNDZ,    results[i]      *= - 2.0;                                                                    ) )
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCN or ONLY_FCNDZ,    results[i]      -= PI * EXPY[i] * ( bessel_factory.struve0[i] + bessel_factory.y0[i] );      ) )

    // Calculate dFdX
    if ( r11_count < 1 ) // All the points are in the R00 region
    {
        STATIC_COND( ONLY_FCNDR,  (R00_dXCEV<N, mode_loop>::evaluate( n, XB, YB, results_dx ));  )
    }
    else if ( r00_count < 1 ) // All the points are in the R11 region
    {
        STATIC_COND( ONLY_FCNDR,  (R11_dXCEV<N, mode_loop>::evaluate( n, XB, YB, results_dx ));  )

        STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDR,  results_dx[i]   *= -2.0 * XINV[i];                                                                     ) )
        STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDR,  results_dx[i]   += 2.0 * XINV[i] * YB[i] / SQRT[i];                                                    ) )
        STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDR,  results_dx[i]   += PI * EXPY[i] * ( bessel_factory.y1[i] + bessel_factory.struve1[i] - 2.0 / PI );     ) )
    }
    else // Partially some of the points are in the R00 region and the others int the R11 region
    {
        STATIC_COND( ONLY_FCNDR, (R00_dXCEV<N, STATIC_LOOP_OFF>::evaluate( r00_count, x_r00, y_r00, r00_results_dx )  ); )
        STATIC_COND( ONLY_FCNDR, (R11_dXCEV<N, STATIC_LOOP_OFF>::evaluate( r11_count, x_r11, y_r11, r11_results_dx )  ); )

        LOOP_DEF( r11_count, STATIC_COND( ONLY_FCNDR,  r11_results_dx[i]        *= -2.0 * XINV[r11_pos[i]];                                                                                     ) )
        LOOP_DEF( r11_count, STATIC_COND( ONLY_FCNDR,  r11_results_dx[i]        += 2.0 * XINV[r11_pos[i]] * YB[r11_pos[i]] / SQRT[r11_pos[i]];                                                  ) )
        LOOP_DEF( r11_count, STATIC_COND( ONLY_FCNDR,  r11_results_dx[i]        += PI * EXPY[r11_pos[i]] * ( bessel_factory.y1[r11_pos[i]] + bessel_factory.struve1[r11_pos[i]] - 2.0 / PI );   ) )

        LOOP_DEF( r00_count, STATIC_COND( ONLY_FCNDR,  results_dx[r00_pos[i]]   = r00_results_dx[i];                                                                                            ) )
        LOOP_DEF( r11_count, STATIC_COND( ONLY_FCNDR,  results_dx[r11_pos[i]]   = r11_results_dx[i];                                                                                            ) )

    }

    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDR, results_dx[i] = ( XB[i] < x_min && YB[i] > y_min ) ? 0.0 : results_dx[i]; ))
    
    // Remove singular value if panel LID
    STATIC_LOOP( n, N, STATIC_COND( IS_LID_PANEL,              results[i]      -= 2.0 * ( LOG2_GAMMA - LOG_SQRT[i] ); ) )
    
    // Check for extreme cases
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCN,   results[i]      = ( X[i] > 0.25 && Y[i] < 1e-5 ) ? -PI * EXPY[i] * ( bessel_factory.struve0[i] + bessel_factory.y0[i] ) : results[i]; ))
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDR, results_dx[i]   = ( X[i] > 0.25 && Y[i] < 1e-5 ) ? PI * EXPY[i] * ( bessel_factory.y1[i] + bessel_factory.struve1[i] - 2.0 / PI ) : results_dx[i]; ))
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDR, results_dx[i]   = ( X[i] < 1e-5 && Y[i] > 0.5 ) ? 0.0 : results_dx[i]; ))

    // Calculate dFdY
    STATIC_LOOP( n, N, STATIC_COND( ONLY_FCNDZ,                results_dy[i]    = - results[i] - 2.0 / SQRT[i];                                             ) )

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
    STATIC_COND( ONLY_FCN,      (L1CEV<N, mode_loop>::evaluate( n, A, B, results ));         )
    STATIC_COND( ONLY_FCNDR,    (L1_dACEV<N, mode_loop>::evaluate( n, A, B, results_da ));   )
    STATIC_COND( ONLY_FCNDZ,    (L1_dBCEV<N, mode_loop>::evaluate( n, A, B, results_db ));   )

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
    STATIC_COND( ONLY_FCN,      (L3CEV<N, mode_loop>::evaluate( n, A, B, results ));         )
    STATIC_COND( ONLY_FCNDR,    (L3_dACEV<N, mode_loop>::evaluate( n, A, B, results_da ));   )
    STATIC_COND( ONLY_FCNDZ,    (L3_dBCEV<N, mode_loop>::evaluate( n, A, B, results_db ));   )

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
    STATIC_COND( ONLY_FCN,      (M1CEV<N, mode_loop>::evaluate( n, A, B, results ));         )
    STATIC_COND( ONLY_FCNDR,    (M1_dACEV<N, mode_loop>::evaluate( n, A, B, results_da ));   )
    STATIC_COND( ONLY_FCNDZ,    (M1_dBCEV<N, mode_loop>::evaluate( n, A, B, results_db ));   )

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
    STATIC_COND( ONLY_FCN,      (M3CEV<N, mode_loop>::evaluate( n, A, B, results ));          )
    STATIC_COND( ONLY_FCNDR,    (M3_dACEV<N, mode_loop>::evaluate( n, A, B, results_da ));    )
    STATIC_COND( ONLY_FCNDZ,    (M3_dBCEV<N, mode_loop>::evaluate( n, A, B, results_db ));    )

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
                                                cuscomplex*             G_dz,
                                                cuscomplex*             G_dzeta
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
    cusfloat csz[N];
    cusfloat scz[N];
    cusfloat ci_g[N];
    cusfloat ci_g_dr[N];
    cusfloat ci_g_dz[N];
    cusfloat ci_g_dzeta[N];
    cusfloat c_real[N];
    cusfloat c_real_dz[N];
    cusfloat c_real_dzeta[N];
    cusfloat c_g_imag[N];
    cusfloat c_g_dr_imag[N];
    cusfloat c_g_dz_imag[N];
    cusfloat c_g_dzeta_imag[N];
    cusfloat expsum[N];
    cusfloat expsum_dz[N];
    cusfloat expsum_dzeta[N];
    cusfloat k0r[N];
    cusfloat knir[N];
    cusfloat knizetah[N];
    cusfloat knizh[N];
    cusfloat rank_f[N];
    cusfloat rank_dfdr[N];
    cusfloat rank_dfdz[N];
    cusfloat rank_dfdzeta[N];
    cusfloat R2[N];
    cusfloat r1_inv[N];
    cusfloat r2_inv[N];
    cusfloat r3_inv[N];
    cusfloat r4_inv[N];
    cusfloat r5_inv[N];
    cusfloat r6_inv[N];
    cusfloat r1_inv_2[N];
    cusfloat r2_inv_2[N];
    cusfloat r3_inv_2[N];
    cusfloat r4_inv_2[N];
    cusfloat r5_inv_2[N];
    cusfloat r6_inv_2[N];
    cusfloat v1[N];
    cusfloat v2[N];
    cusfloat v3[N];
    cusfloat v4[N];
    cusfloat v5[N];
    cusfloat v6[N];
    cusfloat zh[N];
    cusfloat zetah[N];

    STATIC_COND( ONLY_FCN,      ( clear_vector<cusfloat,N>( c_g_imag )    );      )
    STATIC_COND( ONLY_FCNDR,    ( clear_vector<cusfloat,N>( c_g_dr_imag ) );      )
    STATIC_COND( ONLY_FCNDZ,    ( clear_vector<cusfloat,N>( c_g_dz_imag ) );      )
    STATIC_COND( ONLY_FCNDZ,    ( clear_vector<cusfloat,N>( c_g_dzeta_imag ) );   )

    for ( std::size_t i=0; i<N; i++ )
    {
        v1[i] = std::abs( z[i] - zeta[i] );
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        v2[i] = z[i] + zeta[i] + 2.0 * h;
        v3[i] = z[i] + zeta[i];
        v4[i] = z[i] - zeta[i] + 2.0 * h;
        v5[i] = zeta[i] - z[i] + 2.0 * h;
        v6[i] = z[i] + zeta[i] + 4.0 * h;
    }

    // Calculate rankine terms to substract from Green function
    for ( std::size_t i=0; i<N; i++ )
    {
        R2[i] = pow2s( R[i] );
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        r1_inv_2[i] = 1.0 / ( R2[i] + pow2s( v1[i] ) );
        r2_inv_2[i] = 1.0 / ( R2[i] + pow2s( v2[i] ) );
        r3_inv_2[i] = 1.0 / ( R2[i] + pow2s( v3[i] ) );
        r4_inv_2[i] = 1.0 / ( R2[i] + pow2s( v4[i] ) );
        r5_inv_2[i] = 1.0 / ( R2[i] + pow2s( v5[i] ) );
        r6_inv_2[i] = 1.0 / ( R2[i] + pow2s( v6[i] ) );
    }

    lv_sqrt<cusfloat>( N, r1_inv_2, r1_inv );
    lv_sqrt<cusfloat>( N, r2_inv_2, r2_inv );
    lv_sqrt<cusfloat>( N, r3_inv_2, r3_inv );
    lv_sqrt<cusfloat>( N, r4_inv_2, r4_inv );
    lv_sqrt<cusfloat>( N, r5_inv_2, r5_inv );
    lv_sqrt<cusfloat>( N, r6_inv_2, r6_inv );

    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( ONLY_FCN, rank_f[i] = r1_inv[i] + r2_inv[i] + r3_inv[i] + r4_inv[i] + r5_inv[i] + r6_inv[i]; )
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( 
                        ONLY_FCNDR, 
                        rank_dfdr[i] = ( 
                                            r1_inv[i] * r1_inv_2[i]
                                            + 
                                            r2_inv[i] * r2_inv_2[i]
                                            + 
                                            r3_inv[i] * r3_inv_2[i]
                                            + 
                                            r4_inv[i] * r4_inv_2[i]
                                            + 
                                            r5_inv[i] * r5_inv_2[i]
                                            + 
                                            r6_inv[i] * r6_inv_2[i]
                                        ) * R[i];
                )
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( 
                        ONLY_FCNDZ,
                        rank_dfdz[i] = ( 
                                            r1_inv[i] * r1_inv_2[i] * v1[i] * sign( v1[i] )
                                            + 
                                            r2_inv[i] * r2_inv_2[i] * v2[i]
                                            + 
                                            r3_inv[i] * r3_inv_2[i] * v3[i]
                                            + 
                                            r4_inv[i] * r4_inv_2[i] * v4[i]
                                            -
                                            r5_inv[i] * r5_inv_2[i] * v5[i]
                                            + 
                                            r6_inv[i] * r6_inv_2[i] * v6[i]
                                        );
                    )
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( 
                        ONLY_FCNDZ,
                        rank_dfdzeta[i] = ( 
                                                r1_inv[i] * r1_inv_2[i] * v1[i] * sign( v1[i] )
                                                + 
                                                r2_inv[i] * r2_inv_2[i] * v2[i]
                                                + 
                                                r3_inv[i] * r3_inv_2[i] * v3[i]
                                                - 
                                                r4_inv[i] * r4_inv_2[i] * v4[i]
                                                +
                                                r5_inv[i] * r5_inv_2[i] * v5[i]
                                                + 
                                                r6_inv[i] * r6_inv_2[i] * v6[i]
                                            );
                    )
    }

    // Calcuate real root series part
    for ( std::size_t i=0; i<N; i++ )
    {
        k0r[i] = k0 * R[i];
    }

    bessel_factory.calculate_cheby( N, k0r );

    cusfloat c0, c1, c2, c3;
    for ( std::size_t i=0; i<N; i++ )
    {
        c0 = exp(  k0 * v3[i] );
        c1 = exp( -k0 * v4[i] );
        c2 = exp( -k0 * v5[i] );
        c3 = exp( -k0 * v6[i] );

        STATIC_COND( 
                        ONLY_FCN or ONLY_FCNDR,  
                        expsum[i]       = (
                                            + c0
                                            + c1
                                            + c2
                                            + c3
                                        );
                    )
        
        STATIC_COND( 
                        ONLY_FCNDZ, 
                        expsum_dz[i]    = - (
                                            - c0
                                            + c1
                                            - c2
                                            + c3
                                        ) * k0;
                    )

        STATIC_COND( 
                        ONLY_FCNDZ, 
                        expsum_dzeta[i] = - (
                                            - c0
                                            - c1
                                            + c2
                                            + c3
                                        ) * k0;
                    )
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( ONLY_FCN or ONLY_FCNDR,    c_real[i]       = - wave_data.k0nu * expsum[i];         )
        STATIC_COND( ONLY_FCNDZ,                c_real_dz[i]    = - wave_data.k0nu * expsum_dz[i];      )
        STATIC_COND( ONLY_FCNDZ,                c_real_dzeta[i] = - wave_data.k0nu * expsum_dzeta[i];   )
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( ONLY_FCN,                  G[i]            = + c_real[i] * ( bessel_factory.y0[i] - 1i * bessel_factory.j0[i] );       )
        STATIC_COND( ONLY_FCNDR,                G_dr[i]         = - c_real[i] * k0 * ( bessel_factory.y1[i] - 1i * bessel_factory.j1[i] );  )
        STATIC_COND( ONLY_FCNDZ,                G_dz[i]         = + c_real_dz[i] * ( bessel_factory.y0[i] - 1i * bessel_factory.j0[i] );    )
        STATIC_COND( ONLY_FCNDZ,                G_dzeta[i]      = + c_real_dzeta[i] * ( bessel_factory.y0[i] - 1i * bessel_factory.j0[i] ); )
    }

    // Pre-calculate static loop terms
    for ( std::size_t i=0; i<N; i++ )
    {
        zetah[i]    = zeta[i] + h;
        zh[i]       = z[i] + h;
    }

    // Calculate imag root series part
    cusfloat    c_max           = 0.0;
    cusfloat    c_g_max         = 0.0;
    cusfloat    c_g_dr_max      = 0.0;
    cusfloat    c_g_dz_max      = 0.0;
    cusfloat    c_g_dzeta_max   = 0.0;
    int         count_k         = 0;
    cusfloat    kni             = 0.0;

    while (true)
    {
        // Get local copy of th imaginary wave number
        kni = wave_data.kn[count_k];

        // Calculate local loop variables
        for ( std::size_t i=0; i<N; i++ )
        {
            knir[i]     = kni * R[i];
            knizetah[i] = kni * zetah[i];
            knizh[i]    = kni * zh[i];
        }

        for ( std::size_t i=0; i<N; i++ )
        {
            ccz[i]      = 4.0 * wave_data.knnu[count_k] * cos( knizetah[i] );
            scz[i]      = ccz[i];
            csz[i]      = 4.0 * wave_data.knnu[count_k] * cos( knizh[i] );
        }

        for ( std::size_t i=0; i<N; i++ )
        {
            STATIC_COND( ONLY_FCN or ONLY_FCNDR,    ccz[i]     = cos( knizh[i] )    * ccz[i]; )
            STATIC_COND( ONLY_FCNDZ,                scz[i]     = sin( knizh[i] )    * scz[i]; )
            STATIC_COND( ONLY_FCNDZ,                csz[i]     = sin( knizetah[i] ) * csz[i]; )
        }

        // Calculate Bessel function for kni*R
        bessel_factory.calculate_modified( N, knir );

        // Calculate i term of the series
        for ( std::size_t i=0; i<N; i++ )
        {
            STATIC_COND( ONLY_FCN,      ci_g[i]         = + ccz[i] * bessel_factory.k0[i];          )
            STATIC_COND( ONLY_FCNDR,    ci_g_dr[i]      = - ccz[i] * kni * bessel_factory.k1[i];    )
            STATIC_COND( ONLY_FCNDZ,    ci_g_dz[i]      = - scz[i] * kni * bessel_factory.k0[i];    )
            STATIC_COND( ONLY_FCNDZ,    ci_g_dzeta[i]   = - csz[i] * kni * bessel_factory.k0[i];    )
        }

        for ( std::size_t i=0; i<N; i++ )
        {
            STATIC_COND( ONLY_FCN,      c_g_imag[i]         += ci_g[i];         )
            STATIC_COND( ONLY_FCNDR,    c_g_dr_imag[i]      += ci_g_dr[i];      )
            STATIC_COND( ONLY_FCNDZ,    c_g_dz_imag[i]      += ci_g_dz[i];      )
            STATIC_COND( ONLY_FCNDZ,    c_g_dzeta_imag[i]   += ci_g_dzeta[i];   )
        }

        // Check for convergence
        for ( std::size_t i=0; i<N; i++ )
        {
            STATIC_COND( ONLY_FCN,      ci_g[i]         = std::abs( ci_g[i] );          )
            STATIC_COND( ONLY_FCNDR,    ci_g_dr[i]      = std::abs( ci_g_dr[i] );       )
            STATIC_COND( ONLY_FCNDZ,    ci_g_dz[i]      = std::abs( ci_g_dz[i] );       )
            STATIC_COND( ONLY_FCNDZ,    ci_g_dzeta[i]   = std::abs( ci_g_dzeta[i] );    )
        }

        c_g_max         = 0.0;
        c_g_dr_max      = 0.0;
        c_g_dz_max      = 0.0;
        c_g_dzeta_max   = 0.0;
        for ( std::size_t i=0; i<N; i++ )
        {
            STATIC_COND( ONLY_FCN,      c_g_max         = std::max( c_g_max, ci_g[i] );             )
            STATIC_COND( ONLY_FCNDR,    c_g_dr_max      = std::max( c_g_dr_max, ci_g_dr[i] );       )
            STATIC_COND( ONLY_FCNDZ,    c_g_dz_max      = std::max( c_g_dz_max, ci_g_dz[i] );       )
            STATIC_COND( ONLY_FCNDZ,    c_g_dzeta_max   = std::max( c_g_dzeta_max, ci_g_dzeta[i] ); )
        }

        c_max = std::max( c_g_max, std::max( c_g_dr_max, std::max( c_g_dz_max, c_g_dzeta_max ) ) );
        if ( c_max < 1e-6 )
        {
            break;
        }

        // Check for the limit in imaginary roots
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

        // Update counter
        count_k++;
    }

    // Add series values
    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( ONLY_FCN,      G[i]        += c_g_imag[i];         )
        STATIC_COND( ONLY_FCNDR,    G_dr[i]     += c_g_dr_imag[i];      )
        STATIC_COND( ONLY_FCNDZ,    G_dz[i]     += c_g_dz_imag[i];      )
        STATIC_COND( ONLY_FCNDZ,    G_dzeta[i]  += c_g_dzeta_imag[i];   )
    }

}


inline cuscomplex  john_series(
                                                cusfloat R,
                                                cusfloat z, 
                                                cusfloat zeta, 
                                                cusfloat h,
                                                WaveDispersionFONK &wave_data
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
    // assert(wave_data.is_john);

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
    cusfloat    c_real = -wave_data.k0nu*expsum;
    cuscomplex  sol( bessely0(k0r), -besselj0(k0r) );
                sol = c_real*sol;

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


template<std::size_t N, int mode_f, int mode_dfdr, int mode_dfdz, int mode_fslid>
void        wave_term_integral(
                                                cusfloat*                   Ri,
                                                cusfloat*                   zi,
                                                cusfloat*                   zeta,
                                                cusfloat                    h,
                                                BesselFactoryVecUpTo<N>     &bessel_factory,
                                                WaveDispersionFONK          &wave_data,
                                                cuscomplex*                 G,
                                                cuscomplex*                 G_dr,
                                                cuscomplex*                 G_dz,
                                                cuscomplex*                 G_dzeta
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

    // Check minimum distances from Source to field point
    cusfloat R[N];
    cusfloat z[N];
    for ( std::size_t i=0; i<N; i++ )
    {
        R[i] = ( Ri[i] < 1e-5 ) ? 1e-5: Ri[i];
        z[i] = ( std::abs( zi[i]-zeta[i] ) < 1e-5 ) ? zi[i]-1e-5: zi[i];
    }

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

    cusfloat sg_z_p_zeta[N];
    cusfloat sg_z_m_zeta[N];

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

    // Calculate sign( z - zeta ) for dB case
    for ( std::size_t i=0; i<N; i++ )
    {
        sg_z_p_zeta[i]  = sign( z[i] + zeta[i] );
        sg_z_m_zeta[i]  = sign( z[i] - zeta[i] );
    }

    // for ( std::size_t i=0; i<N; i++ )
    // {
    //     sg_z_p_zeta[i] = ( v2[i] < FS_SEL_THR ) ? 0.0 : sg_z_p_zeta[i];
    //     // sg_z_p_zeta[i] = ( v2[i] < FS_SEL_THR ) ? -1.0 : -1.0;
    //     sg_z_m_zeta[i] = ( v0[i] < FS_SEL_THR ) ? 0.0 : sg_z_m_zeta[i];
    // }

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

    // std::cout << "WI - R: " << R[0] << " - z: " << z[0] << " - zeta: " << zeta[0] << std::endl;

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
            Fxy<N, STATIC_LOOP_ON, mode_f, mode_dfdr, mode_dfdz, mode_fslid>( N, X_gt1, Y_gt1, bessel_factory, res_fxy_bgt1, res_fxy_bgt1_dx, res_fxy_bgt1_dy );
        }
        else
        {
            G1_Hgt1<N, STATIC_LOOP_OFF, mode_f, mode_dfdr, mode_dfdz>( b1_lt1_count, A_lt1, b1_lt1, res_fcn1_blt1, res_fcn1_blt1_da, res_fcn1_blt1_db );
            G2_Hgt1<N, STATIC_LOOP_OFF, mode_f, mode_dfdr, mode_dfdz>( b1_gt1_count, A_gt1, b1_gt1, res_fcn1_bgt1, res_fcn1_bgt1_da, res_fcn1_bgt1_db );
            Fxy<N, STATIC_LOOP_OFF, mode_f, mode_dfdr, mode_dfdz, mode_fslid>( b1_gt1_count, X_gt1, Y_gt1, bessel_factory, res_fxy_bgt1, res_fxy_bgt1_dx, res_fxy_bgt1_dy );
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
            Fxy<N, STATIC_LOOP_ON, mode_f, mode_dfdr, mode_dfdz, mode_fslid>( N, X_gt1, Y_gt1, bessel_factory, res_fxy_bgt1, res_fxy_bgt1_dx, res_fxy_bgt1_dy );
        }
        else
        {
            G1_Hlt1<N, STATIC_LOOP_OFF, mode_f, mode_dfdr, mode_dfdz>( b1_lt1_count, A_lt1, b1_lt1, res_fcn1_blt1, res_fcn1_blt1_da, res_fcn1_blt1_db );
            G2_Hlt1<N, STATIC_LOOP_OFF, mode_f, mode_dfdr, mode_dfdz>( b1_gt1_count, A_gt1, b1_gt1, res_fcn1_bgt1, res_fcn1_bgt1_da, res_fcn1_bgt1_db );
            Fxy<N, STATIC_LOOP_OFF, mode_f, mode_dfdr, mode_dfdz, mode_fslid>( b1_gt1_count, X_gt1, Y_gt1, bessel_factory, res_fxy_bgt1, res_fxy_bgt1_dx, res_fxy_bgt1_dy );
        }
    }

    // Calculate real part
    cusfloat G_real[N];
    cusfloat G_dr_real[N];
    cusfloat G_dz_real[N];
    cusfloat G_dzeta_real[N];

    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( ONLY_FCN,      G_real[i]       = res_fcn0[i];                          )
        STATIC_COND( ONLY_FCNDR,    G_dr_real[i]    = res_fcn0_da[i];                       )
        STATIC_COND( ONLY_FCNDZ,    G_dz_real[i]    = res_fcn0_db[i] * sg_z_m_zeta[i];      )
        STATIC_COND( ONLY_FCNDZ,    G_dzeta_real[i] = - res_fcn0_db[i] * sg_z_m_zeta[i];    )
    }

    

    for ( std::size_t i=0; i<b1_lt1_count; i++ )
    {
        STATIC_COND( ONLY_FCN,      G_real[b1_lt1_pos[i]]       += res_fcn1_blt1[i];    )
        STATIC_COND( ONLY_FCNDR,    G_dr_real[b1_lt1_pos[i]]    += res_fcn1_blt1_da[i]; )
        STATIC_COND( ONLY_FCNDZ,    G_dz_real[b1_lt1_pos[i]]    += res_fcn1_blt1_db[i]; )
        STATIC_COND( ONLY_FCNDZ,    G_dzeta_real[b1_lt1_pos[i]] += res_fcn1_blt1_db[i]; )
    }

    for ( std::size_t i=0; i<b1_gt1_count; i++ )
    {
        STATIC_COND( ONLY_FCN,      G_real[b1_gt1_pos[i]]       += res_fcn1_bgt1[i];    )
        STATIC_COND( ONLY_FCNDR,    G_dr_real[b1_gt1_pos[i]]    += res_fcn1_bgt1_da[i]; )
        STATIC_COND( ONLY_FCNDZ,    G_dz_real[b1_gt1_pos[i]]    += res_fcn1_bgt1_db[i]; )
        STATIC_COND( ONLY_FCNDZ,    G_dzeta_real[b1_gt1_pos[i]] += res_fcn1_bgt1_db[i]; )
    }

    cusfloat h2  = h * h;
    cusfloat nu2 = nu * nu;
    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( ONLY_FCN,      G_real[i]       /= h;   )
        STATIC_COND( ONLY_FCNDR,    G_dr_real[i]    /= h2;  )
        STATIC_COND( ONLY_FCNDZ,    G_dz_real[i]    /= h2;  )
        STATIC_COND( ONLY_FCNDZ,    G_dzeta_real[i] /= h2;  )
    }

    for ( std::size_t i=0; i<b1_gt1_count; i++ )
    {
        STATIC_COND( ONLY_FCN,      G_real[b1_gt1_pos[i]]       += nu * res_fxy_bgt1[i];                        )
        STATIC_COND( ONLY_FCNDR,    G_dr_real[b1_gt1_pos[i]]    += nu2 * res_fxy_bgt1_dx[i];                    )
        STATIC_COND( ONLY_FCNDZ,    G_dz_real[b1_gt1_pos[i]]    += nu2 * res_fxy_bgt1_dy[i] * sg_z_p_zeta[i];   )
        STATIC_COND( ONLY_FCNDZ,    G_dzeta_real[b1_gt1_pos[i]] += nu2 * res_fxy_bgt1_dy[i] * sg_z_p_zeta[i];   )
    }

    // Calculate exponential terms for imaginary part
    cusfloat expsum[N], expsum_dz[N], expsum_dzeta[N];
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

        STATIC_COND( ONLY_FCN or ONLY_FCNDR,    expsum[i]       = m0 + c0 + c2;                         )
        STATIC_COND( ONLY_FCNDZ,                expsum_dz[i]    = m0 + c0 * sg_z_p_zeta[i] - c2;        )
        STATIC_COND( ONLY_FCNDZ,                expsum_dzeta[i] = c0 * sg_z_p_zeta[i] - c1 + c2 + c3;   )
    }

    // Calculate beseel functions
    cusfloat j0_vec[N];
    cusfloat j1_vec[N];

    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( ONLY_FCN or ONLY_FCNDZ,    j0_vec[i] = besselj0( k0 * R[i] ); )
        STATIC_COND( ONLY_FCNDR,                j1_vec[i] = besselj1( k0 * R[i] ); )
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        STATIC_COND( ONLY_FCN,      G[i]       = cuscomplex( G_real[i],             k0nu * expsum[i] * j0_vec[i]               ); )
        STATIC_COND( ONLY_FCNDR,    G_dr[i]    = cuscomplex( G_dr_real[i],         -k0nu * expsum[i] * j1_vec[i] * k0          ); )
        STATIC_COND( ONLY_FCNDZ,    G_dz[i]    = cuscomplex( G_dz_real[i],         -k0nu * expsum_dz[i] * j0_vec[i] * k0       ); )
        STATIC_COND( ONLY_FCNDZ,    G_dzeta[i] = cuscomplex( G_dzeta_real[i],      -k0nu * expsum_dzeta[i] * j0_vec[i] * k0    ); )
    }

}