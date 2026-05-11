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

// ---------------------------------------------------------------------------
// time_domain_evaluator_gauss.hpp
//
// Optimised time-domain Green's function evaluators for Duhamel integration
// at a predetermined set of Gauss quadrature points.
//
// Motivation
// ----------
// In the Duhamel integral the time argument tau is evaluated at a fixed set
// of Gauss quadrature nodes.  These nodes map to a fixed set of beta values
// (beta = sqrt(g/r) * tau) so it is wasteful to re-evaluate the beta
// direction of the approximation at every call.  Two precomputation steps
// are performed once at construction time:
//
//   1. Beta-folded residual  (BetaFoldedResidual<TD>)
//      The 2D Chebyshev expansion  R(beta, log_mu)  is collapsed along the
//      beta dimension for each Gauss point, producing a compact 1D expansion
//      in log_mu.  Run-time evaluation requires only a cheap 1D Chebyshev
//      sum.
//
//   2. G0 time-cache  (G0Gauss1AlphaEvaluator / G0Gauss3AlphaEvaluator)
//      The G0 basis functions involve expensive Bessel-function products that
//      depend only on beta (and the fixed alpha-shift constants).  These are
//      split into a pure "time" factor F0(beta, alpha) and an optional
//      "time × mu" factor F1(beta, alpha) such that:
//
//          BasisFcn(beta, mu, alpha) = [F0 + F1*mu] * exp(-beta²/4 * mu)
//
//      F0 and F1 are precomputed once per Gauss point.  Run-time evaluation
//      costs one exp(), a handful of 1D Chebyshev A_k(log_mu) evaluations,
//      and a few multiplications.
//
// Usage
// -----
//   std::vector<cusfloat> beta_gauss = { ... };   // Gauss time points
//
//   GaussEval_dGdt   eval_dGdt(  beta_gauss );
//   GaussEval_dGdtx  eval_dGdtx( beta_gauss );
//   // ...
//
//   // Inside Duhamel loop:
//   cusfloat val = eval_dGdt.evaluate( gauss_index, mu );
//
// The original (non-optimised) evaluators are kept intact in
// time_domain_evaluator.hpp, which this header depends on.
// ---------------------------------------------------------------------------

#include <cmath>
#include <vector>
#include <algorithm>

// Pull in eval_time_1d<AT>, the G0 basis functions, and ChebyshevTraits
// specialisations for all time-domain types.
#include "time_domain_evaluator.hpp"


// ===========================================================================
// BetaFoldedResidual<TD>
// ===========================================================================
// Precomputes the beta-direction of the 2D Chebyshev residual for each Gauss
// point, leaving a 1D expansion in log10(mu) per mu-patch that can be
// evaluated very cheaply at runtime.
// ===========================================================================
template<typename TD>
class BetaFoldedResidual
{
public:
    // Construct from a set of Gauss beta values.
    explicit BetaFoldedResidual( const std::vector<cusfloat>& beta_gauss )
    {
        const int ng = static_cast<int>( beta_gauss.size() );

        // Number of uniform patches in each direction
        const int nx_blocks = static_cast<int>( std::round(
            ( TD::x_max_global - TD::x_min_global ) / TD::dx_min_region ) );
        const int ny_blocks = static_cast<int>( std::round(
            ( TD::y_max_global - TD::y_min_global ) / TD::dy_min_region ) );

        num_gauss_ = ng;
        ny_blocks_ = ny_blocks;
        patches_.resize( static_cast<std::size_t>( ng * ny_blocks ) );

        for ( int g = 0; g < ng; ++g )
        {
            // Clamp beta to the valid domain and locate its x-patch
            const cusfloat xc = std::max( std::min( beta_gauss[g],
                                                     TD::x_max_global ),
                                           TD::x_min_global );
            int nx = static_cast<int>(
                std::floor( ( xc - TD::x_min_global ) / TD::dx_min_region ) );
            nx = std::min( std::max( nx, 0 ), nx_blocks - 1 );

            for ( int ny = 0; ny < ny_blocks; ++ny )
            {
                // 2D block index (row-major: x varies slowest)
                const int nt = nx * ny_blocks + ny;

                const std::size_t sp = TD::blocks_start[nt];
                const std::size_t np = TD::blocks_coeffs_np[nt];
                const int         mo = static_cast<int>( TD::blocks_max_cheby_order[nt] );

                // Map beta into the per-block normalised coordinate [-1, 1]
                const cusfloat xm =
                    static_cast<cusfloat>( 2.0 )
                    * ( xc - TD::x_min_region[nt] ) / TD::dx_region[nt]
                    - static_cast<cusfloat>( 1.0 );

                // Evaluate T_0(xm) .. T_mo(xm)
                cusfloat poly_x[ TD::max_cheby_order + 1 ];
                chebyshev_poly_upto_order(
                    static_cast<std::size_t>( mo ), xm, poly_x );

                // Fold the beta dimension:
                //   folded[k] = sum_{j : ncy[sp+j]==k} coeffs[sp+j] * T_{ncx[sp+j]}(xm)
                cusfloat folded[ TD::max_cheby_order + 1 ] = {};
                for ( std::size_t j = 0; j < np; ++j )
                {
                    folded[ TD::ncy[sp + j] ] +=
                        TD::coeffs[sp + j] * poly_x[ TD::ncx[sp + j] ];
                }

                // Store patch metadata (y = log_mu direction)
                PatchInfo& pi       = patches_[
                    static_cast<std::size_t>( g * ny_blocks + ny ) ];
                pi.x_min            = TD::y_min_region[nt];
                pi.dx               = TD::dy_region[nt];
                pi.start            = static_cast<int>( flat_coeffs_.size() );
                pi.num_terms        = mo + 1;   // dense, includes possible zeros

                // Append to flat coefficient storage
                for ( int k = 0; k <= mo; ++k )
                    flat_coeffs_.push_back( folded[k] );
            }
        }
    }

    // Evaluate the 1D folded residual for Gauss point g at log10(mu).
    cusfloat evaluate( int g, cusfloat log_mu ) const
    {
        const cusfloat yc = std::max( std::min( log_mu,
                                                 TD::y_max_global ),
                                       TD::y_min_global );

        // Locate the mu-patch
        int ny = static_cast<int>(
            std::floor( ( yc - TD::y_min_global ) / TD::dy_min_region ) );
        ny = std::min( std::max( ny, 0 ), ny_blocks_ - 1 );

        const PatchInfo& pi = patches_[
            static_cast<std::size_t>( g * ny_blocks_ + ny ) ];

        // Normalised mu coordinate
        const cusfloat ym =
            static_cast<cusfloat>( 2.0 )
            * ( yc - pi.x_min ) / pi.dx
            - static_cast<cusfloat>( 1.0 );

        // 1D Chebyshev evaluation
        cusfloat poly_y[ TD::max_cheby_order + 1 ];
        chebyshev_poly_upto_order(
            static_cast<std::size_t>( pi.num_terms - 1 ), ym, poly_y );

        const cusfloat* c  = flat_coeffs_.data() + pi.start;
        cusfloat        res = static_cast<cusfloat>( 0.0 );
        for ( int k = 0; k < pi.num_terms; ++k )
            res += c[k] * poly_y[k];

        return res;
    }

    int num_gauss()  const { return num_gauss_; }
    int ny_blocks()  const { return ny_blocks_; }

private:
    struct PatchInfo
    {
        cusfloat x_min;     ///< log_mu minimum of the patch
        cusfloat dx;        ///< log_mu width of the patch
        int      start;     ///< index into flat_coeffs_
        int      num_terms; ///< number of stored 1D Chebyshev coefficients
    };

    int                    num_gauss_;
    int                    ny_blocks_;
    std::vector<PatchInfo> patches_;      ///< [g * ny_blocks_ + ny]
    std::vector<cusfloat>  flat_coeffs_;  ///< densely packed folded coefficients
};


// ===========================================================================
// G0 time-only factored basis functions
// ===========================================================================
//
// Every G0 basis function has the structure
//
//     BasisFcn(beta, mu, alpha) = [ F0(beta, alpha) + F1(beta, alpha) * mu ]
//                                  * exp( -beta²/4 * mu )
//
// The functions below return the beta/alpha-only factors F0 and F1.  These
// are precomputed once per Gauss point; only exp() and the A_k(log_mu)
// Chebyshev evaluations remain at runtime.
//
// The x-derivative relations dGdtx = (-β²/4)·dGdt  and
// dGdttx = (-β²/4)·dGdtt  (and their second-order counterparts) allow the
// F0/F1 functions of derivative types to be expressed as simple scalings.
// ===========================================================================

// ---------------------------------------------------------------------------
// dGdt  (one alpha shift, F1 = 0)
// ---------------------------------------------------------------------------
inline cusfloat dGdt_G0_F0( cusfloat beta, cusfloat alpha )
{
    if ( beta < static_cast<cusfloat>( 1e-1 ) )
        return static_cast<cusfloat>( 0.0 );

    const cusfloat x2    = beta * beta / static_cast<cusfloat>( 8.0 ) + alpha;
    const cusfloat scale = static_cast<cusfloat>( PI )
                           * beta * beta * beta
                           / ( static_cast<cusfloat>( 16.0 )
                               * static_cast<cusfloat>( std::sqrt( 2.0 ) ) );

    using std::cyl_bessel_j;
    cusfloat bt = cyl_bessel_j( static_cast<cusfloat>(  0.25 ), x2 )
                * cyl_bessel_j( static_cast<cusfloat>( -0.25 ), x2 )
                + cyl_bessel_j( static_cast<cusfloat>(  0.75 ), x2 )
                * cyl_bessel_j( static_cast<cusfloat>( -0.75 ), x2 );
    if ( std::isnan( bt ) ) bt = static_cast<cusfloat>( 0.0 );
    return scale * bt;
}

inline cusfloat dGdt_G0_F1( cusfloat /*beta*/, cusfloat /*alpha*/ )
{
    return static_cast<cusfloat>( 0.0 );
}

// ---------------------------------------------------------------------------
// dGdtx = (-β²/4) * dGdt  (three alpha shifts, F1 = 0)
// ---------------------------------------------------------------------------
inline cusfloat dGdtx_G0_F0( cusfloat beta, cusfloat alpha )
{
    return ( -beta * beta / static_cast<cusfloat>( 4.0 ) )
           * dGdt_G0_F0( beta, alpha );
}

inline cusfloat dGdtx_G0_F1( cusfloat /*beta*/, cusfloat /*alpha*/ )
{
    return static_cast<cusfloat>( 0.0 );
}

// ---------------------------------------------------------------------------
// dGdtxx = (-β²/4) * dGdtx  (three alpha shifts, F1 = 0)
// ---------------------------------------------------------------------------
inline cusfloat dGdtxx_G0_F0( cusfloat beta, cusfloat alpha )
{
    return ( -beta * beta / static_cast<cusfloat>( 4.0 ) )
           * dGdtx_G0_F0( beta, alpha );
}

inline cusfloat dGdtxx_G0_F1( cusfloat /*beta*/, cusfloat /*alpha*/ )
{
    return static_cast<cusfloat>( 0.0 );
}

// ---------------------------------------------------------------------------
// dGdtt  (three alpha shifts; F1 ≠ 0 due to the mu * exp(…) term in T1)
//
// The original formula is:
//   lt     = π / (16 √2)
//   T1     = lt * (-0.5 * β⁴ * y0) * μ * exp(-β²μ/4)
//   T2     = lt * yt              * exp(-β²μ/4)
//   result = T1 + T2
//
// Therefore:
//   F0 = lt * yt          (coefficient of exp(-cμ))
//   F1 = lt * (-0.5β⁴y0)  (coefficient of μ·exp(-cμ))
// ---------------------------------------------------------------------------
inline cusfloat dGdtt_G0_F0( cusfloat beta, cusfloat alpha )
{
    const cusfloat lt = static_cast<cusfloat>( PI )
                        / ( static_cast<cusfloat>( 16.0 )
                            * static_cast<cusfloat>( std::sqrt( 2.0 ) ) );

    cusfloat x = beta * beta / static_cast<cusfloat>( 8.0 ) + alpha;
    if ( x < static_cast<cusfloat>( 1e-6 ) )
        x = static_cast<cusfloat>( 0.0 );

    using std::cyl_bessel_j;
    const cusfloat a      = static_cast<cusfloat>( 3.0 ) * beta * beta;
    const cusfloat b      = beta * beta * beta * beta / static_cast<cusfloat>( 4.0 );
    const cusfloat c_coef = b / static_cast<cusfloat>( 2.0 );

    cusfloat yt =
          a      * ( cyl_bessel_j( static_cast<cusfloat>(  0.25 ), x )
                   * cyl_bessel_j( static_cast<cusfloat>( -0.25 ), x )
                   + cyl_bessel_j( static_cast<cusfloat>(  0.75 ), x )
                   * cyl_bessel_j( static_cast<cusfloat>( -0.75 ), x ) )
        + b      * ( cyl_bessel_j( static_cast<cusfloat>( -0.75 ), x )
                   * cyl_bessel_j( static_cast<cusfloat>( -0.25 ), x )
                   - cyl_bessel_j( static_cast<cusfloat>(  0.75 ), x )
                   * cyl_bessel_j( static_cast<cusfloat>(  0.25 ), x ) )
        + c_coef * ( -cyl_bessel_j( static_cast<cusfloat>(  1.25 ), x )
                    * cyl_bessel_j( static_cast<cusfloat>( -0.25 ), x )
                    + cyl_bessel_j( static_cast<cusfloat>( -1.25 ), x )
                    * cyl_bessel_j( static_cast<cusfloat>(  0.25 ), x ) )
        + c_coef * ( -cyl_bessel_j( static_cast<cusfloat>(  1.75 ), x )
                    * cyl_bessel_j( static_cast<cusfloat>( -0.75 ), x )
                    + cyl_bessel_j( static_cast<cusfloat>( -1.75 ), x )
                    * cyl_bessel_j( static_cast<cusfloat>(  0.75 ), x ) );

    if ( std::isnan( yt ) ) yt = static_cast<cusfloat>( 0.0 );

    // Mirror the linear ramp from the original basis function:
    // the T2 (F0) contribution is zeroed for very small beta.
    if ( beta < static_cast<cusfloat>( 5e-2 ) )
        yt = static_cast<cusfloat>( 0.0 );

    return lt * yt;
}

inline cusfloat dGdtt_G0_F1( cusfloat beta, cusfloat alpha )
{
    const cusfloat lt = static_cast<cusfloat>( PI )
                        / ( static_cast<cusfloat>( 16.0 )
                            * static_cast<cusfloat>( std::sqrt( 2.0 ) ) );

    cusfloat x = beta * beta / static_cast<cusfloat>( 8.0 ) + alpha;
    if ( x < static_cast<cusfloat>( 1e-6 ) )
        x = static_cast<cusfloat>( 0.0 );

    using std::cyl_bessel_j;
    cusfloat y0 = cyl_bessel_j( static_cast<cusfloat>(  0.25 ), x )
                * cyl_bessel_j( static_cast<cusfloat>( -0.25 ), x )
                + cyl_bessel_j( static_cast<cusfloat>(  0.75 ), x )
                * cyl_bessel_j( static_cast<cusfloat>( -0.75 ), x );
    if ( std::isnan( y0 ) ) y0 = static_cast<cusfloat>( 0.0 );

    return lt
           * ( -static_cast<cusfloat>( 0.5 )
               * beta * beta * beta * beta * y0 );
}

// ---------------------------------------------------------------------------
// dGdttx = (-β²/4) * dGdtt  (three alpha shifts, F1 ≠ 0)
// ---------------------------------------------------------------------------
inline cusfloat dGdttx_G0_F0( cusfloat beta, cusfloat alpha )
{
    return ( -beta * beta / static_cast<cusfloat>( 4.0 ) )
           * dGdtt_G0_F0( beta, alpha );
}

inline cusfloat dGdttx_G0_F1( cusfloat beta, cusfloat alpha )
{
    return ( -beta * beta / static_cast<cusfloat>( 4.0 ) )
           * dGdtt_G0_F1( beta, alpha );
}

// ---------------------------------------------------------------------------
// dGdttxx = (-β²/4) * dGdttx  (three alpha shifts, F1 ≠ 0)
// ---------------------------------------------------------------------------
inline cusfloat dGdttxx_G0_F0( cusfloat beta, cusfloat alpha )
{
    return ( -beta * beta / static_cast<cusfloat>( 4.0 ) )
           * dGdttx_G0_F0( beta, alpha );
}

inline cusfloat dGdttxx_G0_F1( cusfloat beta, cusfloat alpha )
{
    return ( -beta * beta / static_cast<cusfloat>( 4.0 ) )
           * dGdttx_G0_F1( beta, alpha );
}


// ===========================================================================
// G0Gauss1AlphaEvaluator<TD, A0T>
// ===========================================================================
// Precomputed G0 evaluator for functions that use a single alpha shift
// (currently only dGdt).
//
// Precomputed per Gauss point:  c = β²/4,  F0,  F1
// Runtime evaluation:  A0(log_mu) * (F0 + F1*μ) * exp(-c*μ)
// ===========================================================================
template<typename TD, typename A0T>
class G0Gauss1AlphaEvaluator
{
public:
    using F0FcnPtr = cusfloat (*)( cusfloat, cusfloat );
    using F1FcnPtr = cusfloat (*)( cusfloat, cusfloat );

    G0Gauss1AlphaEvaluator( const std::vector<cusfloat>& beta_gauss,
                             F0FcnPtr                     f0_fcn,
                             F1FcnPtr                     f1_fcn )
    {
        cache_.resize( beta_gauss.size() );
        for ( std::size_t g = 0; g < beta_gauss.size(); ++g )
        {
            const cusfloat beta = beta_gauss[g];
            GaussData&     d    = cache_[g];
            d.c  = beta * beta / static_cast<cusfloat>( 4.0 );
            d.F0 = f0_fcn( beta, TD::G0_alpha_shift[0] );
            d.F1 = f1_fcn( beta, TD::G0_alpha_shift[0] );
        }
    }

    // Evaluate the G0 contribution at Gauss point g and geometrical param μ.
    cusfloat evaluate( int g, cusfloat mu ) const
    {
        const GaussData& d     = cache_[static_cast<std::size_t>( g )];
        const cusfloat   log_mu = std::log10( mu );
        const cusfloat   A0    = eval_time_1d<A0T>( log_mu );
        const cusfloat   e     = std::exp( -d.c * mu );
        return A0 * ( d.F0 * e + d.F1 * mu * e );
    }

private:
    struct GaussData { cusfloat c, F0, F1; };
    std::vector<GaussData> cache_;
};


// ===========================================================================
// G0Gauss3AlphaEvaluator<TD, A0T, A1T, A2T>
// ===========================================================================
// Precomputed G0 evaluator for functions that use three alpha shifts
// (dGdtx, dGdtxx, dGdtt, dGdttx, dGdttxx).
//
// Precomputed per Gauss point:  c = β²/4,  F0[3],  F1[3]
// Runtime evaluation:
//   sum_{k=0}^{2} A_k(log_mu) * (F0[k] + F1[k]*μ) * exp(-c*μ)
// ===========================================================================
template<typename TD, typename A0T, typename A1T, typename A2T>
class G0Gauss3AlphaEvaluator
{
public:
    using F0FcnPtr = cusfloat (*)( cusfloat, cusfloat );
    using F1FcnPtr = cusfloat (*)( cusfloat, cusfloat );

    G0Gauss3AlphaEvaluator( const std::vector<cusfloat>& beta_gauss,
                             F0FcnPtr                     f0_fcn,
                             F1FcnPtr                     f1_fcn )
    {
        cache_.resize( beta_gauss.size() );
        for ( std::size_t g = 0; g < beta_gauss.size(); ++g )
        {
            const cusfloat beta = beta_gauss[g];
            GaussData&     d    = cache_[g];
            d.c = beta * beta / static_cast<cusfloat>( 4.0 );
            for ( int k = 0; k < 3; ++k )
            {
                d.F0[k] = f0_fcn( beta, TD::G0_alpha_shift[k] );
                d.F1[k] = f1_fcn( beta, TD::G0_alpha_shift[k] );
            }
        }
    }

    // Evaluate the G0 contribution at Gauss point g and geometrical param μ.
    cusfloat evaluate( int g, cusfloat mu ) const
    {
        const GaussData& d      = cache_[static_cast<std::size_t>( g )];
        const cusfloat   log_mu = std::log10( mu );
        const cusfloat   A0     = eval_time_1d<A0T>( log_mu );
        const cusfloat   A1     = eval_time_1d<A1T>( log_mu );
        const cusfloat   A2     = eval_time_1d<A2T>( log_mu );
        const cusfloat   e      = std::exp( -d.c * mu );
        const cusfloat   me     = mu * e;
        return ( A0 * ( d.F0[0] * e + d.F1[0] * me )
               + A1 * ( d.F0[1] * e + d.F1[1] * me )
               + A2 * ( d.F0[2] * e + d.F1[2] * me ) );
    }

private:
    struct GaussData
    {
        cusfloat c;
        cusfloat F0[3];
        cusfloat F1[3];
    };
    std::vector<GaussData> cache_;
};


// ===========================================================================
// Concrete per-function-type Gauss evaluators
// ===========================================================================
//
// Each class wraps a BetaFoldedResidual and the appropriate G0 evaluator.
// The constructor takes the vector of Gauss beta values; after construction
// the object is fully precomputed and evaluate(g, mu) is cheap.
//
// Component accessors evaluate_G0() and evaluate_residual() mirror those in
// time_domain_evaluator.hpp and are useful for testing and decomposed output.
// ===========================================================================

// ---------------------------------------------------------------------------
// dGdt
// ---------------------------------------------------------------------------
class GaussEval_dGdt
{
    using TDT  = ChebyshevTraits<dGdtC>;
    using A0T  = ChebyshevTraits<dGdtA0C>;

public:
    explicit GaussEval_dGdt( const std::vector<cusfloat>& beta_gauss )
        : residual_( beta_gauss )
        , g0_( beta_gauss, dGdt_G0_F0, dGdt_G0_F1 )
    {}

    cusfloat evaluate( int g, cusfloat mu ) const
    {
        return g0_.evaluate( g, mu )
             + residual_.evaluate( g, std::log10( mu ) );
    }

    cusfloat evaluate_G0      ( int g, cusfloat mu ) const { return g0_.evaluate( g, mu ); }
    cusfloat evaluate_residual( int g, cusfloat mu ) const { return residual_.evaluate( g, std::log10( mu ) ); }

    int num_gauss() const { return residual_.num_gauss(); }

private:
    BetaFoldedResidual<TDT>           residual_;
    G0Gauss1AlphaEvaluator<TDT, A0T>  g0_;
};


// ---------------------------------------------------------------------------
// dGdtx
// ---------------------------------------------------------------------------
class GaussEval_dGdtx
{
    using TDT  = ChebyshevTraits<dGdtxC>;
    using A0T  = ChebyshevTraits<dGdtxA0C>;
    using A1T  = ChebyshevTraits<dGdtxA1C>;
    using A2T  = ChebyshevTraits<dGdtxA2C>;

public:
    explicit GaussEval_dGdtx( const std::vector<cusfloat>& beta_gauss )
        : residual_( beta_gauss )
        , g0_( beta_gauss, dGdtx_G0_F0, dGdtx_G0_F1 )
    {}

    cusfloat evaluate( int g, cusfloat mu ) const
    {
        return g0_.evaluate( g, mu )
             + residual_.evaluate( g, std::log10( mu ) );
    }

    cusfloat evaluate_G0      ( int g, cusfloat mu ) const { return g0_.evaluate( g, mu ); }
    cusfloat evaluate_residual( int g, cusfloat mu ) const { return residual_.evaluate( g, std::log10( mu ) ); }

    int num_gauss() const { return residual_.num_gauss(); }

private:
    BetaFoldedResidual<TDT>                    residual_;
    G0Gauss3AlphaEvaluator<TDT, A0T, A1T, A2T> g0_;
};


// ---------------------------------------------------------------------------
// dGdtxx
// ---------------------------------------------------------------------------
class GaussEval_dGdtxx
{
    using TDT  = ChebyshevTraits<dGdtxxC>;
    using A0T  = ChebyshevTraits<dGdtxxA0C>;
    using A1T  = ChebyshevTraits<dGdtxxA1C>;
    using A2T  = ChebyshevTraits<dGdtxxA2C>;

public:
    explicit GaussEval_dGdtxx( const std::vector<cusfloat>& beta_gauss )
        : residual_( beta_gauss )
        , g0_( beta_gauss, dGdtxx_G0_F0, dGdtxx_G0_F1 )
    {}

    cusfloat evaluate( int g, cusfloat mu ) const
    {
        return g0_.evaluate( g, mu )
             + residual_.evaluate( g, std::log10( mu ) );
    }

    cusfloat evaluate_G0      ( int g, cusfloat mu ) const { return g0_.evaluate( g, mu ); }
    cusfloat evaluate_residual( int g, cusfloat mu ) const { return residual_.evaluate( g, std::log10( mu ) ); }

    int num_gauss() const { return residual_.num_gauss(); }

private:
    BetaFoldedResidual<TDT>                    residual_;
    G0Gauss3AlphaEvaluator<TDT, A0T, A1T, A2T> g0_;
};


// ---------------------------------------------------------------------------
// dGdtt
// ---------------------------------------------------------------------------
class GaussEval_dGdtt
{
    using TDT  = ChebyshevTraits<dGdttC>;
    using A0T  = ChebyshevTraits<dGdttA0C>;
    using A1T  = ChebyshevTraits<dGdttA1C>;
    using A2T  = ChebyshevTraits<dGdttA2C>;

public:
    explicit GaussEval_dGdtt( const std::vector<cusfloat>& beta_gauss )
        : residual_( beta_gauss )
        , g0_( beta_gauss, dGdtt_G0_F0, dGdtt_G0_F1 )
    {}

    cusfloat evaluate( int g, cusfloat mu ) const
    {
        return g0_.evaluate( g, mu )
             + residual_.evaluate( g, std::log10( mu ) );
    }

    cusfloat evaluate_G0      ( int g, cusfloat mu ) const { return g0_.evaluate( g, mu ); }
    cusfloat evaluate_residual( int g, cusfloat mu ) const { return residual_.evaluate( g, std::log10( mu ) ); }

    int num_gauss() const { return residual_.num_gauss(); }

private:
    BetaFoldedResidual<TDT>                    residual_;
    G0Gauss3AlphaEvaluator<TDT, A0T, A1T, A2T> g0_;
};


// ---------------------------------------------------------------------------
// dGdttx
// ---------------------------------------------------------------------------
class GaussEval_dGdttx
{
    using TDT  = ChebyshevTraits<dGdttxC>;
    using A0T  = ChebyshevTraits<dGdttxA0C>;
    using A1T  = ChebyshevTraits<dGdttxA1C>;
    using A2T  = ChebyshevTraits<dGdttxA2C>;

public:
    explicit GaussEval_dGdttx( const std::vector<cusfloat>& beta_gauss )
        : residual_( beta_gauss )
        , g0_( beta_gauss, dGdttx_G0_F0, dGdttx_G0_F1 )
    {}

    cusfloat evaluate( int g, cusfloat mu ) const
    {
        return g0_.evaluate( g, mu )
             + residual_.evaluate( g, std::log10( mu ) );
    }

    cusfloat evaluate_G0      ( int g, cusfloat mu ) const { return g0_.evaluate( g, mu ); }
    cusfloat evaluate_residual( int g, cusfloat mu ) const { return residual_.evaluate( g, std::log10( mu ) ); }

    int num_gauss() const { return residual_.num_gauss(); }

private:
    BetaFoldedResidual<TDT>                    residual_;
    G0Gauss3AlphaEvaluator<TDT, A0T, A1T, A2T> g0_;
};


// ---------------------------------------------------------------------------
// dGdttxx
// ---------------------------------------------------------------------------
class GaussEval_dGdttxx
{
    using TDT  = ChebyshevTraits<dGdttxxC>;
    using A0T  = ChebyshevTraits<dGdttxxA0C>;
    using A1T  = ChebyshevTraits<dGdttxxA1C>;
    using A2T  = ChebyshevTraits<dGdttxxA2C>;

public:
    explicit GaussEval_dGdttxx( const std::vector<cusfloat>& beta_gauss )
        : residual_( beta_gauss )
        , g0_( beta_gauss, dGdttxx_G0_F0, dGdttxx_G0_F1 )
    {}

    cusfloat evaluate( int g, cusfloat mu ) const
    {
        return g0_.evaluate( g, mu )
             + residual_.evaluate( g, std::log10( mu ) );
    }

    cusfloat evaluate_G0      ( int g, cusfloat mu ) const { return g0_.evaluate( g, mu ); }
    cusfloat evaluate_residual( int g, cusfloat mu ) const { return residual_.evaluate( g, std::log10( mu ) ); }

    int num_gauss() const { return residual_.num_gauss(); }

private:
    BetaFoldedResidual<TDT>                    residual_;
    G0Gauss3AlphaEvaluator<TDT, A0T, A1T, A2T> g0_;
};
