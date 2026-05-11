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

#pragma once

#include <cmath>
#include <algorithm>

#include "../config.hpp"
#include "../math/chebyshev_evaluation.hpp"


// ---------------------------------------------------------------------------
// 2D residual evaluator for time-domain types.
//
// TD must be a ChebyshevTraits specialisation for a time-domain 2D type
// (e.g., ChebyshevTraits<dGdtC>).  These types have:
//   - x: beta (linear, [x_min_global, x_max_global])
//   - y: log10(mu) (log-scale, [y_min_global, y_max_global])
//
// The block hash formula uses x/y-patch counts derived from the global domain
// extents and the minimum region size, NOT the intervals_np field (which is 1
// for all time types and encodes the max-refinement level within a patch, not
// the number of patches).
// ---------------------------------------------------------------------------
template<typename TD>
cusfloat eval_time_residual_2d( cusfloat beta, cusfloat log_mu )
{
    // Clamp inputs to domain
    const cusfloat x = std::max( std::min( beta,   TD::x_max_global ), TD::x_min_global );
    const cusfloat y = std::max( std::min( log_mu, TD::y_max_global ), TD::y_min_global );

    // Compute patch counts along each axis
    const int nx_blocks = static_cast<int>( std::round(
        ( TD::x_max_global - TD::x_min_global ) / TD::dx_min_region ) );
    const int ny_blocks = static_cast<int>( std::round(
        ( TD::y_max_global - TD::y_min_global ) / TD::dy_min_region ) );

    // Hash to block index
    int nx = static_cast<int>( std::floor( ( x - TD::x_min_global ) / TD::dx_min_region ) );
    int ny = static_cast<int>( std::floor( ( y - TD::y_min_global ) / TD::dy_min_region ) );
    nx     = std::min( std::max( nx, 0 ), nx_blocks - 1 );
    ny     = std::min( std::max( ny, 0 ), ny_blocks - 1 );
    const int nt = nx * ny_blocks + ny;

    // Block properties
    const std::size_t sp  = TD::blocks_start[nt];
    const std::size_t np  = TD::blocks_coeffs_np[nt];

    // Map to [-1, 1] using per-block domain info
    const cusfloat xm = static_cast<cusfloat>( 2.0 ) * ( x - TD::x_min_region[nt] ) / TD::dx_region[nt] - static_cast<cusfloat>( 1.0 );
    const cusfloat ym = static_cast<cusfloat>( 2.0 ) * ( y - TD::y_min_region[nt] ) / TD::dy_region[nt] - static_cast<cusfloat>( 1.0 );

    // Evaluate 2D Chebyshev: sum c_j * T_{ncx_j}(xm) * T_{ncy_j}(ym)
    cusfloat result = static_cast<cusfloat>( 0.0 );
    evaluate_chebyshev_polynomials_2d_t<TD>( sp, np, static_cast<std::size_t>( nt ), xm, ym, result );

    return result;
}


// ---------------------------------------------------------------------------
// 1D evaluator for time-domain A-coefficient types (dGdtA0C, dGdttA0C, etc.)
//
// AT must expose the 1D ChebyshevTraits fields:
//   x_min_global, x_max_global, dx_min_region, intervals_np,
//   blocks_start, blocks_coeffs_np, x_min_region, dx_region,
//   coeffs, ncx, max_cheby_order.
//
// Input: log_mu = log10(mu)
// ---------------------------------------------------------------------------
template<typename AT>
cusfloat eval_time_1d( cusfloat log_mu )
{
    // Clamp input to domain
    const cusfloat x = std::max( std::min( log_mu, AT::x_max_global ), AT::x_min_global );

    // Hash to block (linear in log-mu space)
    int nx = static_cast<int>( std::floor( ( x - AT::x_min_global ) / AT::dx_min_region ) );
    nx     = std::min( std::max( nx, 0 ), static_cast<int>( AT::intervals_np ) - 1 );

    const std::size_t sp = AT::blocks_start[nx];
    const std::size_t np = AT::blocks_coeffs_np[nx];

    // Map to [-1, 1]
    const cusfloat xm = static_cast<cusfloat>( 2.0 ) * ( x - AT::x_min_region[nx] ) / AT::dx_region[nx] - static_cast<cusfloat>( 1.0 );

    // Evaluate 1D Chebyshev: sum c_j * T_{ncx_j}(xm)
    cusfloat poly_x[ AT::max_cheby_order + 1 ];
    chebyshev_poly_upto_order( static_cast<std::size_t>( AT::blocks_max_cheby_order[nx] ), xm, poly_x );

    cusfloat result = static_cast<cusfloat>( 0.0 );
    for ( std::size_t i = 0; i < np; ++i )
    {
        result += AT::coeffs[sp + i] * poly_x[ AT::ncx[sp + i] ];
    }

    return result;
}


// ---------------------------------------------------------------------------
// G0 basis functions
//
// These are the analytic asymptotic basis functions used to represent the
// singular behaviour of the time-domain Green's function and its derivatives.
// Formulas derived from Wehausen & Laitone (1960) and their mu/beta derivatives.
//
// Convention:  x = beta^2/8 + alpha,  mu is the geometrical parameter z/r.
// ---------------------------------------------------------------------------

// dG/dt  G0 basis
inline cusfloat dGdt_G0_basis( cusfloat beta, cusfloat mu, cusfloat alpha )
{
    if ( beta < static_cast<cusfloat>( 1e-1 ) )
        return static_cast<cusfloat>( 0.0 );

    const cusfloat x2  = beta * beta / static_cast<cusfloat>( 8.0 ) + alpha;
    const cusfloat lt  = static_cast<cusfloat>( PI ) * beta * beta * beta
                         / ( static_cast<cusfloat>( 16.0 ) * static_cast<cusfloat>( std::sqrt( 2.0 ) ) )
                         * std::exp( -beta * beta * mu / static_cast<cusfloat>( 4.0 ) );

    cusfloat bt = std::cyl_bessel_j( static_cast<cusfloat>( 0.25 ),  x2 )
                * std::cyl_bessel_j( static_cast<cusfloat>( -0.25 ), x2 )
                + std::cyl_bessel_j( static_cast<cusfloat>( 0.75 ),  x2 )
                * std::cyl_bessel_j( static_cast<cusfloat>( -0.75 ), x2 );

    if ( std::isnan( bt ) )
        bt = static_cast<cusfloat>( 0.0 );

    return lt * bt;
}

// d^2G/(dt dx)  G0 basis  (x-derivative of dGdt_G0)
inline cusfloat dGdtx_G0_basis( cusfloat beta, cusfloat mu, cusfloat alpha )
{
    return ( -beta * beta / static_cast<cusfloat>( 4.0 ) ) * dGdt_G0_basis( beta, mu, alpha );
}

// d^3G/(dt dx^2)  G0 basis
inline cusfloat dGdtxx_G0_basis( cusfloat beta, cusfloat mu, cusfloat alpha )
{
    return ( -beta * beta / static_cast<cusfloat>( 4.0 ) ) * dGdtx_G0_basis( beta, mu, alpha );
}

// d^2G/dt^2  G0 basis
inline cusfloat dGdtt_G0_basis( cusfloat beta, cusfloat mu, cusfloat alpha )
{
    const cusfloat lt   = static_cast<cusfloat>( PI )
                          / ( static_cast<cusfloat>( 16.0 ) * static_cast<cusfloat>( std::sqrt( 2.0 ) ) );
    const cusfloat expt = std::exp( -beta * beta * mu / static_cast<cusfloat>( 4.0 ) );
    cusfloat x          = beta * beta / static_cast<cusfloat>( 8.0 ) + alpha;
    if ( x < static_cast<cusfloat>( 1e-6 ) )
        x = static_cast<cusfloat>( 0.0 );

    using std::cyl_bessel_j;

    // Term 1: contribution from the mu-exponential derivative
    cusfloat y0 = cyl_bessel_j( static_cast<cusfloat>(  0.25 ), x )
                * cyl_bessel_j( static_cast<cusfloat>( -0.25 ), x )
                + cyl_bessel_j( static_cast<cusfloat>(  0.75 ), x )
                * cyl_bessel_j( static_cast<cusfloat>( -0.75 ), x );
    if ( std::isnan( y0 ) ) y0 = static_cast<cusfloat>( 0.0 );
    const cusfloat T1 = -static_cast<cusfloat>( 0.5 ) * beta * beta * beta * y0 * beta * mu * expt;

    // Term 2: contribution from Bessel product differentiation
    const cusfloat a = static_cast<cusfloat>( 3.0 ) * beta * beta;
    const cusfloat b = beta * beta * beta * beta / static_cast<cusfloat>( 4.0 );
    const cusfloat c = b / static_cast<cusfloat>( 2.0 );

    cusfloat yt = a * ( cyl_bessel_j( static_cast<cusfloat>(  0.25 ), x )
                      * cyl_bessel_j( static_cast<cusfloat>( -0.25 ), x )
                      + cyl_bessel_j( static_cast<cusfloat>(  0.75 ), x )
                      * cyl_bessel_j( static_cast<cusfloat>( -0.75 ), x ) )
                + b * ( cyl_bessel_j( static_cast<cusfloat>( -0.75 ), x )
                      * cyl_bessel_j( static_cast<cusfloat>( -0.25 ), x )
                      - cyl_bessel_j( static_cast<cusfloat>(  0.75 ), x )
                      * cyl_bessel_j( static_cast<cusfloat>(  0.25 ), x ) )
                + c * ( -cyl_bessel_j( static_cast<cusfloat>(  1.25 ), x )
                       * cyl_bessel_j( static_cast<cusfloat>( -0.25 ), x )
                       + cyl_bessel_j( static_cast<cusfloat>( -1.25 ), x )
                       * cyl_bessel_j( static_cast<cusfloat>(  0.25 ), x ) )
                + c * ( -cyl_bessel_j( static_cast<cusfloat>(  1.75 ), x )
                       * cyl_bessel_j( static_cast<cusfloat>( -0.75 ), x )
                       + cyl_bessel_j( static_cast<cusfloat>( -1.75 ), x )
                       * cyl_bessel_j( static_cast<cusfloat>(  0.75 ), x ) );

    if ( std::isnan( yt ) ) yt = static_cast<cusfloat>( 0.0 );

    // Linear ramp for very small beta to avoid singular behaviour
    const cusfloat dt = static_cast<cusfloat>( 5e-2 );
    if ( beta < dt && dt > static_cast<cusfloat>( 0.0 ) )
    {
        // The Python code ramps yt linearly: yt[pos] = yt[dt_index] / dt * beta
        // Here we approximate the ramp as 0 since the function approaches 0 at beta=0
        yt = static_cast<cusfloat>( 0.0 );
    }

    const cusfloat T2 = yt * expt;

    return lt * ( T1 + T2 );
}

// d^3G/(dt^2 dx)  G0 basis
inline cusfloat dGdttx_G0_basis( cusfloat beta, cusfloat mu, cusfloat alpha )
{
    return ( -beta * beta / static_cast<cusfloat>( 4.0 ) ) * dGdtt_G0_basis( beta, mu, alpha );
}

// d^4G/(dt^2 dx^2)  G0 basis
inline cusfloat dGdttxx_G0_basis( cusfloat beta, cusfloat mu, cusfloat alpha )
{
    return ( -beta * beta / static_cast<cusfloat>( 4.0 ) ) * dGdttx_G0_basis( beta, mu, alpha );
}


// ---------------------------------------------------------------------------
// G0 contribution evaluator
//
// Evaluates:
//   G0(beta, mu) = sum_{k=0}^{K-1} A_k(log10(mu)) * basis(beta, mu, alpha_k)
//
// TD  : 2D residual traits type (carries G0_alpha_shift, G0_alpha_shift_np)
// A0T : 1D traits type for k=0 amplitude coefficient A_0(log10(mu))
// A1T : 1D traits type for k=1  (pass void_tag if not applicable)
// A2T : 1D traits type for k=2  (pass void_tag if not applicable)
// BasisFcn : function pointer type compatible with
//            cusfloat(cusfloat beta, cusfloat mu, cusfloat alpha)
// ---------------------------------------------------------------------------
template<typename TD, typename A0T, cusfloat BasisFcn( cusfloat, cusfloat, cusfloat )>
cusfloat eval_G0_1alpha( cusfloat beta, cusfloat mu )
{
    const cusfloat log_mu = std::log10( mu );
    const cusfloat A0     = eval_time_1d<A0T>( log_mu );
    return A0 * BasisFcn( beta, mu, TD::G0_alpha_shift[0] );
}

template<typename TD, typename A0T, typename A1T, typename A2T,
         cusfloat BasisFcn( cusfloat, cusfloat, cusfloat )>
cusfloat eval_G0_3alpha( cusfloat beta, cusfloat mu )
{
    const cusfloat log_mu = std::log10( mu );
    cusfloat result = static_cast<cusfloat>( 0.0 );
    result += eval_time_1d<A0T>( log_mu ) * BasisFcn( beta, mu, TD::G0_alpha_shift[0] );
    result += eval_time_1d<A1T>( log_mu ) * BasisFcn( beta, mu, TD::G0_alpha_shift[1] );
    result += eval_time_1d<A2T>( log_mu ) * BasisFcn( beta, mu, TD::G0_alpha_shift[2] );
    return result;
}


// ---------------------------------------------------------------------------
// Full integral evaluators  (G0 + residual)
//
// Notation:
//   dGdt   = d/dt of G        (time derivative)
//   dGdtx  = d^2/(dt dmu)     (time + mu derivative)
//   dGdtxx = d^3/(dt dmu^2)
//   dGdtt  = d^2/dt^2         (second time derivative)
//   dGdttx = d^3/(dt^2 dmu)
//   dGdttxx= d^4/(dt^2 dmu^2)
//
// Each function reconstructs:
//   G(beta, mu) = G0(beta, mu) + Residual(beta, log10(mu))
// ---------------------------------------------------------------------------

// Forward declarations of the coefficient traits specialisations.
// These are defined after including the traits headers — they need to be
// available at the call sites that use the full evaluators below.
#include "chebyshev_traits_time.hpp"

// -- dGdt --
inline cusfloat eval_dGdt( cusfloat beta, cusfloat mu )
{
    using TD  = ChebyshevTraits<dGdtC>;
    using A0T = ChebyshevTraits<dGdtA0C>;

    const cusfloat G0       = eval_G0_1alpha<TD, A0T, dGdt_G0_basis>( beta, mu );
    const cusfloat residual = eval_time_residual_2d<TD>( beta, std::log10( mu ) );
    return G0 + residual;
}

// -- dGdtx --
inline cusfloat eval_dGdtx( cusfloat beta, cusfloat mu )
{
    using TD  = ChebyshevTraits<dGdtxC>;
    using A0T = ChebyshevTraits<dGdtxA0C>;
    using A1T = ChebyshevTraits<dGdtxA1C>;
    using A2T = ChebyshevTraits<dGdtxA2C>;

    const cusfloat G0       = eval_G0_3alpha<TD, A0T, A1T, A2T, dGdtx_G0_basis>( beta, mu );
    const cusfloat residual = eval_time_residual_2d<TD>( beta, std::log10( mu ) );
    return G0 + residual;
}

// -- dGdtxx --
inline cusfloat eval_dGdtxx( cusfloat beta, cusfloat mu )
{
    using TD  = ChebyshevTraits<dGdtxxC>;
    using A0T = ChebyshevTraits<dGdtxxA0C>;
    using A1T = ChebyshevTraits<dGdtxxA1C>;
    using A2T = ChebyshevTraits<dGdtxxA2C>;

    const cusfloat G0       = eval_G0_3alpha<TD, A0T, A1T, A2T, dGdtxx_G0_basis>( beta, mu );
    const cusfloat residual = eval_time_residual_2d<TD>( beta, std::log10( mu ) );
    return G0 + residual;
}

// -- dGdtt --
inline cusfloat eval_dGdtt( cusfloat beta, cusfloat mu )
{
    using TD  = ChebyshevTraits<dGdttC>;
    using A0T = ChebyshevTraits<dGdttA0C>;
    using A1T = ChebyshevTraits<dGdttA1C>;
    using A2T = ChebyshevTraits<dGdttA2C>;

    const cusfloat G0       = eval_G0_3alpha<TD, A0T, A1T, A2T, dGdtt_G0_basis>( beta, mu );
    const cusfloat residual = eval_time_residual_2d<TD>( beta, std::log10( mu ) );
    return G0 + residual;
}

// -- dGdttx --
inline cusfloat eval_dGdttx( cusfloat beta, cusfloat mu )
{
    using TD  = ChebyshevTraits<dGdttxC>;
    using A0T = ChebyshevTraits<dGdttxA0C>;
    using A1T = ChebyshevTraits<dGdttxA1C>;
    using A2T = ChebyshevTraits<dGdttxA2C>;

    const cusfloat G0       = eval_G0_3alpha<TD, A0T, A1T, A2T, dGdttx_G0_basis>( beta, mu );
    const cusfloat residual = eval_time_residual_2d<TD>( beta, std::log10( mu ) );
    return G0 + residual;
}

// -- dGdttxx --
inline cusfloat eval_dGdttxx( cusfloat beta, cusfloat mu )
{
    using TD  = ChebyshevTraits<dGdttxxC>;
    using A0T = ChebyshevTraits<dGdttxxA0C>;
    using A1T = ChebyshevTraits<dGdttxxA1C>;
    using A2T = ChebyshevTraits<dGdttxxA2C>;

    const cusfloat G0       = eval_G0_3alpha<TD, A0T, A1T, A2T, dGdttxx_G0_basis>( beta, mu );
    const cusfloat residual = eval_time_residual_2d<TD>( beta, std::log10( mu ) );
    return G0 + residual;
}


// ---------------------------------------------------------------------------
// Component-only evaluators (for debugging and decomposed unit tests)
//
// *_G0      : returns only the analytic G0 term
// *_residual: returns only the Chebyshev-fitted residual term
// The full evaluator above equals the sum of these two.
// ---------------------------------------------------------------------------

// dGdt
inline cusfloat eval_dGdt_G0( cusfloat beta, cusfloat mu )
{
    using TD  = ChebyshevTraits<dGdtC>;
    using A0T = ChebyshevTraits<dGdtA0C>;
    return eval_G0_1alpha<TD, A0T, dGdt_G0_basis>( beta, mu );
}
inline cusfloat eval_dGdt_residual( cusfloat beta, cusfloat mu )
{
    using TD = ChebyshevTraits<dGdtC>;
    return eval_time_residual_2d<TD>( beta, std::log10( mu ) );
}

// dGdtx
inline cusfloat eval_dGdtx_G0( cusfloat beta, cusfloat mu )
{
    using TD  = ChebyshevTraits<dGdtxC>;
    using A0T = ChebyshevTraits<dGdtxA0C>;
    using A1T = ChebyshevTraits<dGdtxA1C>;
    using A2T = ChebyshevTraits<dGdtxA2C>;
    return eval_G0_3alpha<TD, A0T, A1T, A2T, dGdtx_G0_basis>( beta, mu );
}
inline cusfloat eval_dGdtx_residual( cusfloat beta, cusfloat mu )
{
    using TD = ChebyshevTraits<dGdtxC>;
    return eval_time_residual_2d<TD>( beta, std::log10( mu ) );
}

// dGdtxx
inline cusfloat eval_dGdtxx_G0( cusfloat beta, cusfloat mu )
{
    using TD  = ChebyshevTraits<dGdtxxC>;
    using A0T = ChebyshevTraits<dGdtxxA0C>;
    using A1T = ChebyshevTraits<dGdtxxA1C>;
    using A2T = ChebyshevTraits<dGdtxxA2C>;
    return eval_G0_3alpha<TD, A0T, A1T, A2T, dGdtxx_G0_basis>( beta, mu );
}
inline cusfloat eval_dGdtxx_residual( cusfloat beta, cusfloat mu )
{
    using TD = ChebyshevTraits<dGdtxxC>;
    return eval_time_residual_2d<TD>( beta, std::log10( mu ) );
}

// dGdtt
inline cusfloat eval_dGdtt_G0( cusfloat beta, cusfloat mu )
{
    using TD  = ChebyshevTraits<dGdttC>;
    using A0T = ChebyshevTraits<dGdttA0C>;
    using A1T = ChebyshevTraits<dGdttA1C>;
    using A2T = ChebyshevTraits<dGdttA2C>;
    return eval_G0_3alpha<TD, A0T, A1T, A2T, dGdtt_G0_basis>( beta, mu );
}
inline cusfloat eval_dGdtt_residual( cusfloat beta, cusfloat mu )
{
    using TD = ChebyshevTraits<dGdttC>;
    return eval_time_residual_2d<TD>( beta, std::log10( mu ) );
}

// dGdttx
inline cusfloat eval_dGdttx_G0( cusfloat beta, cusfloat mu )
{
    using TD  = ChebyshevTraits<dGdttxC>;
    using A0T = ChebyshevTraits<dGdttxA0C>;
    using A1T = ChebyshevTraits<dGdttxA1C>;
    using A2T = ChebyshevTraits<dGdttxA2C>;
    return eval_G0_3alpha<TD, A0T, A1T, A2T, dGdttx_G0_basis>( beta, mu );
}
inline cusfloat eval_dGdttx_residual( cusfloat beta, cusfloat mu )
{
    using TD = ChebyshevTraits<dGdttxC>;
    return eval_time_residual_2d<TD>( beta, std::log10( mu ) );
}

// dGdttxx
inline cusfloat eval_dGdttxx_G0( cusfloat beta, cusfloat mu )
{
    using TD  = ChebyshevTraits<dGdttxxC>;
    using A0T = ChebyshevTraits<dGdttxxA0C>;
    using A1T = ChebyshevTraits<dGdttxxA1C>;
    using A2T = ChebyshevTraits<dGdttxxA2C>;
    return eval_G0_3alpha<TD, A0T, A1T, A2T, dGdttxx_G0_basis>( beta, mu );
}
inline cusfloat eval_dGdttxx_residual( cusfloat beta, cusfloat mu )
{
    using TD = ChebyshevTraits<dGdttxxC>;
    return eval_time_residual_2d<TD>( beta, std::log10( mu ) );
}
