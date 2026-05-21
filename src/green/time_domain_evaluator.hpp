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

// d^3G/dt^3  G0 basis
//
// Derived by differentiating dGdtt_G02 = lt' * fj + lt * fj' once more
// via the product rule (matching the Python dGdttt_G0 implementation):
//
//   dGdttt_G0 = lt'' * fj + 2 * lt' * fj' + lt * fj''
//
// where:
//   lt    = π β³/(16√2) · exp(−μβ²/4)
//   lt'   = π β²/(16√2) · exp(−μβ²/4) · (3 − μβ²/2)
//   lt''  = π β /(16√2) · exp(−μβ²/4) · (6 − 7μβ²/2 + μ²β⁴/4)
//   fj    = J_{1/4}(x) J_{−1/4}(x) + J_{3/4}(x) J_{−3/4}(x),   x = β²/8
//   fj'   = −(β/4) · h,   h = J_{−1/4} J_{5/4} + 2 J_{3/4} J_{1/4} + J_{7/4} J_{−3/4}
//   fj''  = −(h/4 + β²/16 · dh/dx)
//           dh/dx = 1/2 J_{−5/4} J_{5/4} + 3/2 fj + 1/2 J_{7/4} J_{−7/4}
//                 − 3/2 J_{3/4}  J_{5/4} − 3/2 J_{7/4} J_{1/4}
//                 − 1/2 J_{−1/4} J_{9/4} − 1/2 J_{11/4} J_{−3/4}
inline cusfloat dGdttt_G0_basis( cusfloat beta, cusfloat mu, cusfloat alpha )
{
    if ( beta < static_cast<cusfloat>( 1e-1 ) )
        return static_cast<cusfloat>( 0.0 );

    const cusfloat scale = static_cast<cusfloat>( PI )
                           / ( static_cast<cusfloat>( 16.0 )
                               * static_cast<cusfloat>( std::sqrt( 2.0 ) ) );
    const cusfloat expt  = std::exp( -beta * beta * mu / static_cast<cusfloat>( 4.0 ) );
    cusfloat x = beta * beta / static_cast<cusfloat>( 8.0 ) + alpha;
    if ( x < static_cast<cusfloat>( 1e-6 ) )
        x = static_cast<cusfloat>( 0.0 );

    using std::cyl_bessel_j;

    // fj = J_{1/4} J_{-1/4} + J_{3/4} J_{-3/4}
    cusfloat fj = cyl_bessel_j( static_cast<cusfloat>(  0.25 ), x )
                * cyl_bessel_j( static_cast<cusfloat>( -0.25 ), x )
                + cyl_bessel_j( static_cast<cusfloat>(  0.75 ), x )
                * cyl_bessel_j( static_cast<cusfloat>( -0.75 ), x );
    if ( std::isnan( fj ) ) fj = static_cast<cusfloat>( 0.0 );

    // h = J_{-1/4} J_{5/4} + 2 J_{3/4} J_{1/4} + J_{7/4} J_{-3/4}
    cusfloat h =       cyl_bessel_j( static_cast<cusfloat>( -0.25 ), x )
                     * cyl_bessel_j( static_cast<cusfloat>(  1.25 ), x )
               + static_cast<cusfloat>( 2.0 )
                     * cyl_bessel_j( static_cast<cusfloat>(  0.75 ), x )
                     * cyl_bessel_j( static_cast<cusfloat>(  0.25 ), x )
               +       cyl_bessel_j( static_cast<cusfloat>(  1.75 ), x )
                     * cyl_bessel_j( static_cast<cusfloat>( -0.75 ), x );
    if ( std::isnan( h ) ) h = static_cast<cusfloat>( 0.0 );

    // fj' = -(β/4) * h
    const cusfloat fj_p = ( -beta / static_cast<cusfloat>( 4.0 ) ) * h;

    // dh/dx
    cusfloat dh =
          static_cast<cusfloat>( 0.5 ) * cyl_bessel_j( static_cast<cusfloat>( -1.25 ), x )
                                        * cyl_bessel_j( static_cast<cusfloat>(  1.25 ), x )
        + static_cast<cusfloat>( 1.5 ) * fj
        + static_cast<cusfloat>( 0.5 ) * cyl_bessel_j( static_cast<cusfloat>(  1.75 ), x )
                                        * cyl_bessel_j( static_cast<cusfloat>( -1.75 ), x )
        - static_cast<cusfloat>( 1.5 ) * cyl_bessel_j( static_cast<cusfloat>(  0.75 ), x )
                                        * cyl_bessel_j( static_cast<cusfloat>(  1.25 ), x )
        - static_cast<cusfloat>( 1.5 ) * cyl_bessel_j( static_cast<cusfloat>(  1.75 ), x )
                                        * cyl_bessel_j( static_cast<cusfloat>(  0.25 ), x )
        - static_cast<cusfloat>( 0.5 ) * cyl_bessel_j( static_cast<cusfloat>( -0.25 ), x )
                                        * cyl_bessel_j( static_cast<cusfloat>(  2.25 ), x )
        - static_cast<cusfloat>( 0.5 ) * cyl_bessel_j( static_cast<cusfloat>(  2.75 ), x )
                                        * cyl_bessel_j( static_cast<cusfloat>( -0.75 ), x );
    if ( std::isnan( dh ) ) dh = static_cast<cusfloat>( 0.0 );

    // fj'' = -(h/4 + β²/16 * dh)
    const cusfloat b2    = beta * beta;
    const cusfloat fj_pp = -( h / static_cast<cusfloat>( 4.0 )
                             + b2 / static_cast<cusfloat>( 16.0 ) * dh );

    // lt, lt', lt'' (exponential envelope and its β-derivatives)
    const cusfloat b3   = b2 * beta;
    const cusfloat b4   = b2 * b2;
    const cusfloat lt   = scale * b3 * expt;
    const cusfloat lt_p = scale * b2 * expt
                          * ( static_cast<cusfloat>( 3.0 )
                              - mu * b2 / static_cast<cusfloat>( 2.0 ) );
    const cusfloat lt_pp = scale * beta * expt
                           * ( static_cast<cusfloat>( 6.0 )
                               - static_cast<cusfloat>( 7.0 ) * mu * b2
                                 / static_cast<cusfloat>( 2.0 )
                               + mu * mu * b4 / static_cast<cusfloat>( 4.0 ) );

    return lt_pp * fj
         + static_cast<cusfloat>( 2.0 ) * lt_p * fj_p
         + lt * fj_pp;
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

// -- dGdttt --
inline cusfloat eval_dGdttt( cusfloat beta, cusfloat mu )
{
    using TD  = ChebyshevTraits<dGdtttC>;
    using A0T = ChebyshevTraits<dGdtttA0C>;

    const cusfloat G0       = eval_G0_1alpha<TD, A0T, dGdttt_G0_basis>( beta, mu );
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

// dGdttt
inline cusfloat eval_dGdttt_G0( cusfloat beta, cusfloat mu )
{
    using TD  = ChebyshevTraits<dGdtttC>;
    using A0T = ChebyshevTraits<dGdtttA0C>;
    return eval_G0_1alpha<TD, A0T, dGdttt_G0_basis>( beta, mu );
}
inline cusfloat eval_dGdttt_residual( cusfloat beta, cusfloat mu )
{
    using TD = ChebyshevTraits<dGdtttC>;
    return eval_time_residual_2d<TD>( beta, std::log10( mu ) );
}


// ---------------------------------------------------------------------------
// Vectorized 2D residual evaluator for time-domain types.
//
// Processes n ≤ N points in one call.  The block hash uses dx_min_region /
// dy_min_region (NOT intervals_np, which is 1 for all time-domain types).
//
// When all n points fall in the same 2D block the call is forwarded to
// evaluate_chebyshev_polynomials_2d_vector_t<TD,N,mode_loop>, which can
// exploit SIMD / loop-unrolling via the mode_loop template parameter.
// Otherwise each point is evaluated independently with the scalar path.
//
// Template parameters:
//   TD        – ChebyshevTraits specialisation for the desired function
//               (e.g. ChebyshevTraits<dGdtC>)
//   N         – compile-time maximum vector width (must be ≥ n at runtime)
//   mode_loop – 0 = runtime loop, 1 = compile-time unrolled loop (passed
//               through to evaluate_chebyshev_polynomials_2d_vector_t)
//
// Inputs (all length-n arrays):
//   beta   – beta values (linear, not log-scaled)
//   log_mu – log10(mu) values (already log-scaled, matching the y-axis)
//   result – output array of length n; written in place
// ---------------------------------------------------------------------------
template<typename TD, std::size_t N, int mode_loop = 0>
void eval_time_residual_2d_vec( std::size_t     n,
                                const cusfloat* beta,
                                const cusfloat* log_mu,
                                cusfloat*       result )
{
    // Precompute block counts (constexpr-friendly via static local)
    const int nx_blocks = static_cast<int>( std::round(
        ( TD::x_max_global - TD::x_min_global ) / TD::dx_min_region ) );
    const int ny_blocks = static_cast<int>( std::round(
        ( TD::y_max_global - TD::y_min_global ) / TD::dy_min_region ) );

    // Per-point: clamp → hash → map to [-1, 1]
    cusfloat    xs[N],  ys[N];
    cusfloat    xm[N],  ym[N];
    std::size_t sp[N],  np_b[N], nt[N];

    for ( std::size_t i = 0; i < n; ++i )
    {
        xs[i] = std::max( std::min( beta[i],   TD::x_max_global ), TD::x_min_global );
        ys[i] = std::max( std::min( log_mu[i], TD::y_max_global ), TD::y_min_global );

        int inx = static_cast<int>( std::floor( ( xs[i] - TD::x_min_global ) / TD::dx_min_region ) );
        int iny = static_cast<int>( std::floor( ( ys[i] - TD::y_min_global ) / TD::dy_min_region ) );
        inx = std::min( std::max( inx, 0 ), nx_blocks - 1 );
        iny = std::min( std::max( iny, 0 ), ny_blocks - 1 );

        nt[i]   = static_cast<std::size_t>( inx * ny_blocks + iny );
        sp[i]   = TD::blocks_start[nt[i]];
        np_b[i] = TD::blocks_coeffs_np[nt[i]];

        xm[i] = static_cast<cusfloat>( 2.0 )
                * ( xs[i] - TD::x_min_region[nt[i]] ) / TD::dx_region[nt[i]]
                - static_cast<cusfloat>( 1.0 );
        ym[i] = static_cast<cusfloat>( 2.0 )
                * ( ys[i] - TD::y_min_region[nt[i]] ) / TD::dy_region[nt[i]]
                - static_cast<cusfloat>( 1.0 );
    }

    // Fast path: all points in the same block → use vectorised inner kernel
    bool single_block = true;
    for ( std::size_t i = 1; i < n; ++i )
    {
        if ( nt[i] != nt[0] ) { single_block = false; break; }
    }

    if ( single_block )
    {
        evaluate_chebyshev_polynomials_2d_vector_t<TD, N, mode_loop>(
            sp[0], np_b[0], nt[0], n, xm, ym, result );
    }
    else
    {
        for ( std::size_t i = 0; i < n; ++i )
        {
            evaluate_chebyshev_polynomials_2d_t<TD>(
                sp[i], np_b[i], nt[i], xm[i], ym[i], result[i] );
        }
    }
}


// ---------------------------------------------------------------------------
// Vectorized full/component evaluators  (free-function wrappers)
//
// These mirror the scalar eval_dGdt / eval_dGdt_G0 / eval_dGdt_residual
// functions but process n ≤ N points in one call.
//
// Template parameters:
//   N         – compile-time maximum vector width (must be ≥ n at runtime)
//   mode_loop – forwarded to eval_time_residual_2d_vec / the vectorised
//               Chebyshev kernel (0 = runtime, 1 = compile-time unrolled)
//
// All functions take raw mu[] (not log-scaled); log10 is applied internally.
// The _residual_vec variants accept log_mu[] directly to avoid recomputing
// the logarithm when the caller already has it.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Helper: compute log10(mu[i]) into log_mu[i] for i < n
// ---------------------------------------------------------------------------
template<std::size_t N>
inline void compute_log_mu( std::size_t n, const cusfloat* mu, cusfloat* log_mu )
{
    for ( std::size_t i = 0; i < n; ++i )
        log_mu[i] = std::log10( mu[i] );
}


// -- dGdt --
template<std::size_t N, int mode_loop = 0>
void eval_dGdt_vec( std::size_t n, const cusfloat* beta, const cusfloat* mu, cusfloat* result )
{
    using TD  = ChebyshevTraits<dGdtC>;
    cusfloat log_mu[N];
    compute_log_mu<N>( n, mu, log_mu );
    eval_time_residual_2d_vec<TD, N, mode_loop>( n, beta, log_mu, result );
    for ( std::size_t i = 0; i < n; ++i )
        result[i] += eval_dGdt_G0( beta[i], mu[i] );
}

template<std::size_t N, int mode_loop = 0>
void eval_dGdt_G0_vec( std::size_t n, const cusfloat* beta, const cusfloat* mu, cusfloat* result )
{
    for ( std::size_t i = 0; i < n; ++i )
        result[i] = eval_dGdt_G0( beta[i], mu[i] );
}

template<std::size_t N, int mode_loop = 0>
void eval_dGdt_residual_vec( std::size_t n, const cusfloat* beta, const cusfloat* log_mu, cusfloat* result )
{
    eval_time_residual_2d_vec<ChebyshevTraits<dGdtC>, N, mode_loop>( n, beta, log_mu, result );
}


// -- dGdtx --
template<std::size_t N, int mode_loop = 0>
void eval_dGdtx_vec( std::size_t n, const cusfloat* beta, const cusfloat* mu, cusfloat* result )
{
    using TD  = ChebyshevTraits<dGdtxC>;
    cusfloat log_mu[N];
    compute_log_mu<N>( n, mu, log_mu );
    eval_time_residual_2d_vec<TD, N, mode_loop>( n, beta, log_mu, result );
    for ( std::size_t i = 0; i < n; ++i )
        result[i] += eval_dGdtx_G0( beta[i], mu[i] );
}

template<std::size_t N, int mode_loop = 0>
void eval_dGdtx_G0_vec( std::size_t n, const cusfloat* beta, const cusfloat* mu, cusfloat* result )
{
    for ( std::size_t i = 0; i < n; ++i )
        result[i] = eval_dGdtx_G0( beta[i], mu[i] );
}

template<std::size_t N, int mode_loop = 0>
void eval_dGdtx_residual_vec( std::size_t n, const cusfloat* beta, const cusfloat* log_mu, cusfloat* result )
{
    eval_time_residual_2d_vec<ChebyshevTraits<dGdtxC>, N, mode_loop>( n, beta, log_mu, result );
}


// -- dGdtxx --
template<std::size_t N, int mode_loop = 0>
void eval_dGdtxx_vec( std::size_t n, const cusfloat* beta, const cusfloat* mu, cusfloat* result )
{
    using TD  = ChebyshevTraits<dGdtxxC>;
    cusfloat log_mu[N];
    compute_log_mu<N>( n, mu, log_mu );
    eval_time_residual_2d_vec<TD, N, mode_loop>( n, beta, log_mu, result );
    for ( std::size_t i = 0; i < n; ++i )
        result[i] += eval_dGdtxx_G0( beta[i], mu[i] );
}

template<std::size_t N, int mode_loop = 0>
void eval_dGdtxx_G0_vec( std::size_t n, const cusfloat* beta, const cusfloat* mu, cusfloat* result )
{
    for ( std::size_t i = 0; i < n; ++i )
        result[i] = eval_dGdtxx_G0( beta[i], mu[i] );
}

template<std::size_t N, int mode_loop = 0>
void eval_dGdtxx_residual_vec( std::size_t n, const cusfloat* beta, const cusfloat* log_mu, cusfloat* result )
{
    eval_time_residual_2d_vec<ChebyshevTraits<dGdtxxC>, N, mode_loop>( n, beta, log_mu, result );
}


// -- dGdtt --
template<std::size_t N, int mode_loop = 0>
void eval_dGdtt_vec( std::size_t n, const cusfloat* beta, const cusfloat* mu, cusfloat* result )
{
    using TD  = ChebyshevTraits<dGdttC>;
    cusfloat log_mu[N];
    compute_log_mu<N>( n, mu, log_mu );
    eval_time_residual_2d_vec<TD, N, mode_loop>( n, beta, log_mu, result );
    for ( std::size_t i = 0; i < n; ++i )
        result[i] += eval_dGdtt_G0( beta[i], mu[i] );
}

template<std::size_t N, int mode_loop = 0>
void eval_dGdtt_G0_vec( std::size_t n, const cusfloat* beta, const cusfloat* mu, cusfloat* result )
{
    for ( std::size_t i = 0; i < n; ++i )
        result[i] = eval_dGdtt_G0( beta[i], mu[i] );
}

template<std::size_t N, int mode_loop = 0>
void eval_dGdtt_residual_vec( std::size_t n, const cusfloat* beta, const cusfloat* log_mu, cusfloat* result )
{
    eval_time_residual_2d_vec<ChebyshevTraits<dGdttC>, N, mode_loop>( n, beta, log_mu, result );
}


// -- dGdttx --
template<std::size_t N, int mode_loop = 0>
void eval_dGdttx_vec( std::size_t n, const cusfloat* beta, const cusfloat* mu, cusfloat* result )
{
    using TD  = ChebyshevTraits<dGdttxC>;
    cusfloat log_mu[N];
    compute_log_mu<N>( n, mu, log_mu );
    eval_time_residual_2d_vec<TD, N, mode_loop>( n, beta, log_mu, result );
    for ( std::size_t i = 0; i < n; ++i )
        result[i] += eval_dGdttx_G0( beta[i], mu[i] );
}

template<std::size_t N, int mode_loop = 0>
void eval_dGdttx_G0_vec( std::size_t n, const cusfloat* beta, const cusfloat* mu, cusfloat* result )
{
    for ( std::size_t i = 0; i < n; ++i )
        result[i] = eval_dGdttx_G0( beta[i], mu[i] );
}

template<std::size_t N, int mode_loop = 0>
void eval_dGdttx_residual_vec( std::size_t n, const cusfloat* beta, const cusfloat* log_mu, cusfloat* result )
{
    eval_time_residual_2d_vec<ChebyshevTraits<dGdttxC>, N, mode_loop>( n, beta, log_mu, result );
}


// -- dGdttxx --
template<std::size_t N, int mode_loop = 0>
void eval_dGdttxx_vec( std::size_t n, const cusfloat* beta, const cusfloat* mu, cusfloat* result )
{
    using TD  = ChebyshevTraits<dGdttxxC>;
    cusfloat log_mu[N];
    compute_log_mu<N>( n, mu, log_mu );
    eval_time_residual_2d_vec<TD, N, mode_loop>( n, beta, log_mu, result );
    for ( std::size_t i = 0; i < n; ++i )
        result[i] += eval_dGdttxx_G0( beta[i], mu[i] );
}

template<std::size_t N, int mode_loop = 0>
void eval_dGdttxx_G0_vec( std::size_t n, const cusfloat* beta, const cusfloat* mu, cusfloat* result )
{
    for ( std::size_t i = 0; i < n; ++i )
        result[i] = eval_dGdttxx_G0( beta[i], mu[i] );
}

template<std::size_t N, int mode_loop = 0>
void eval_dGdttxx_residual_vec( std::size_t n, const cusfloat* beta, const cusfloat* log_mu, cusfloat* result )
{
    eval_time_residual_2d_vec<ChebyshevTraits<dGdttxxC>, N, mode_loop>( n, beta, log_mu, result );
}


// -- dGdttt --
template<std::size_t N, int mode_loop = 0>
void eval_dGdttt_vec( std::size_t n, const cusfloat* beta, const cusfloat* mu, cusfloat* result )
{
    using TD  = ChebyshevTraits<dGdtttC>;
    cusfloat log_mu[N];
    compute_log_mu<N>( n, mu, log_mu );
    eval_time_residual_2d_vec<TD, N, mode_loop>( n, beta, log_mu, result );
    for ( std::size_t i = 0; i < n; ++i )
        result[i] += eval_dGdttt_G0( beta[i], mu[i] );
}

template<std::size_t N, int mode_loop = 0>
void eval_dGdttt_G0_vec( std::size_t n, const cusfloat* beta, const cusfloat* mu, cusfloat* result )
{
    for ( std::size_t i = 0; i < n; ++i )
        result[i] = eval_dGdttt_G0( beta[i], mu[i] );
}

template<std::size_t N, int mode_loop = 0>
void eval_dGdttt_residual_vec( std::size_t n, const cusfloat* beta, const cusfloat* log_mu, cusfloat* result )
{
    eval_time_residual_2d_vec<ChebyshevTraits<dGdtttC>, N, mode_loop>( n, beta, log_mu, result );
}
