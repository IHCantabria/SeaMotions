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
// Pull in fold_time_residual_x and ChebyshevTraits infrastructure.
#include "integrals_database.hpp"


// ===========================================================================
// FoldedResidualSnapshot<NS>
// ===========================================================================
// Value-type snapshot of the mutable fold-X storage held in
// ChebyshevTraits<NS> immediately after a fold_time_residual_x<NS>(beta)
// call.  Storing one snapshot per Gauss point allows BetaFoldedResidual to
// cache all Gauss-point fold results independently of the global static
// storage, which can only hold one beta value at a time.
//
// NS must be a raw-data type instantiated with CHEBYSHEV_TIME_2DF_TRAITS.
// ===========================================================================
template<typename NS>
struct FoldedResidualSnapshot
{
    using CT = ChebyshevTraits<NS>;

    static constexpr std::size_t NY  = CT::ny_blocks_f;
    static constexpr std::size_t NFC = CT::max_f_coeffs;

    cusfloat    coeffs_f[NFC]          = {};
    std::size_t ncy_f[NFC]             = {};
    std::size_t blocks_start_f[NY]     = {};
    std::size_t blocks_np_f[NY]        = {};
    std::size_t blocks_max_order_f[NY] = {};
    cusfloat    y_min_region_f[NY]     = {};
    cusfloat    dy_region_f[NY]        = {};

    // Capture the current global-static fold result from ChebyshevTraits<NS>.
    // Call immediately after fold_time_residual_x<NS>(beta).
    void capture()
    {
        for ( std::size_t i = 0; i < NFC; ++i )
        {
            coeffs_f[i] = CT::coeffs_f[i];
            ncy_f[i]    = CT::ncy_f[i];
        }
        for ( std::size_t i = 0; i < NY; ++i )
        {
            blocks_start_f[i]     = CT::blocks_start_f[i];
            blocks_np_f[i]        = CT::blocks_np_f[i];
            blocks_max_order_f[i] = CT::blocks_max_order_f[i];
            y_min_region_f[i]     = CT::y_min_region_f[i];
            dy_region_f[i]        = CT::dy_region_f[i];
        }
    }

    // Evaluate the 1-D (Y = log_mu) Chebyshev expansion at mu.
    // Applies log-scale if NS::y_log_scale is true (which it is for all
    // standard time-domain types where y = log10(mu)).
    cusfloat evaluate( cusfloat mu ) const
    {
        cusfloat ys = ( NS::y_log_scale ) ? std::log10( mu ) : mu;
        ys = std::max( std::min( ys, NS::y_max_global ), NS::y_min_global );

        int ny = static_cast<int>(
            std::floor( ( ys - NS::y_min_global ) / NS::dy_min_region ) );
        ny = std::max( 0, std::min( ny, static_cast<int>( NY ) - 1 ) );

        const std::size_t sp        = blocks_start_f[ny];
        const std::size_t np        = blocks_np_f[ny];
        const std::size_t max_order = blocks_max_order_f[ny];
        const cusfloat    y_min     = y_min_region_f[ny];
        const cusfloat    dy        = dy_region_f[ny];
        const cusfloat    ym        = 2.0f * ( ys - y_min ) / dy - 1.0f;

        cusfloat poly_y[ NS::max_cheby_order + 1 ];
        chebyshev_poly_upto_order( max_order, ym, poly_y );

        cusfloat result = 0.0f;
        for ( std::size_t j = 0; j < np; ++j )
            result += coeffs_f[sp + j] * poly_y[ ncy_f[sp + j] ];

        return result;
    }
};


// ===========================================================================
// BetaFoldedResidual<NS>
// ===========================================================================
// Precomputes the beta direction of the 2D Chebyshev residual for each
// Gauss point using fold_time_residual_x<NS>, then captures the resulting
// 1-D (Y = log_mu) state into a per-Gauss FoldedResidualSnapshot<NS>.
// Run-time evaluation is a cheap 1-D Chebyshev sum (no 2-D lookup).
//
// NS must be a raw-data type instantiated with CHEBYSHEV_TIME_2DF_TRAITS.
//
// Note: fold_time_residual_x writes to GLOBAL-STATIC storage, so calling
// this constructor from multiple threads for the same NS type concurrently
// is NOT safe.  Constructors for different NS types may run concurrently.
// ===========================================================================
template<typename NS>
class BetaFoldedResidual
{
public:
    explicit BetaFoldedResidual( const std::vector<cusfloat>& beta_gauss )
    {
        snapshots_.resize( beta_gauss.size() );
        for ( std::size_t g = 0; g < beta_gauss.size(); ++g )
        {
            // Fold the beta dimension into the global-static storage of
            // ChebyshevTraits<NS>, then capture a private copy.
            fold_time_residual_x<NS>( beta_gauss[g] );
            snapshots_[g].capture();
        }
    }

    // Evaluate the 1-D folded residual for Gauss point g at viscosity mu.
    cusfloat evaluate( int g, cusfloat mu ) const
    {
        return snapshots_[ static_cast<std::size_t>( g ) ].evaluate( mu );
    }

    // Evaluate using log10(mu) directly (convenience overload used by tests).
    cusfloat evaluate_log( int g, cusfloat log_mu ) const
    {
        // FoldedResidualSnapshot::evaluate applies log-scale internally, so
        // we must pass the raw mu value; compute it here from log_mu.
        return snapshots_[ static_cast<std::size_t>( g ) ]
               .evaluate( std::pow( 10.0f, log_mu ) );
    }

    int num_gauss() const { return static_cast<int>( snapshots_.size() ); }

private:
    std::vector<FoldedResidualSnapshot<NS>> snapshots_;
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
             + residual_.evaluate( g, mu );
    }

    cusfloat evaluate_G0      ( int g, cusfloat mu ) const { return g0_.evaluate( g, mu ); }
    cusfloat evaluate_residual( int g, cusfloat mu ) const { return residual_.evaluate( g, mu ); }

    int num_gauss() const { return residual_.num_gauss(); }

private:
    BetaFoldedResidual<dGdtC>         residual_;
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
             + residual_.evaluate( g, mu );
    }

    cusfloat evaluate_G0      ( int g, cusfloat mu ) const { return g0_.evaluate( g, mu ); }
    cusfloat evaluate_residual( int g, cusfloat mu ) const { return residual_.evaluate( g, mu ); }

    int num_gauss() const { return residual_.num_gauss(); }

private:
    BetaFoldedResidual<dGdtxC>                 residual_;
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
             + residual_.evaluate( g, mu );
    }

    cusfloat evaluate_G0      ( int g, cusfloat mu ) const { return g0_.evaluate( g, mu ); }
    cusfloat evaluate_residual( int g, cusfloat mu ) const { return residual_.evaluate( g, mu ); }

    int num_gauss() const { return residual_.num_gauss(); }

private:
    BetaFoldedResidual<dGdtxxC>                residual_;
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
             + residual_.evaluate( g, mu );
    }

    cusfloat evaluate_G0      ( int g, cusfloat mu ) const { return g0_.evaluate( g, mu ); }
    cusfloat evaluate_residual( int g, cusfloat mu ) const { return residual_.evaluate( g, mu ); }

    int num_gauss() const { return residual_.num_gauss(); }

private:
    BetaFoldedResidual<dGdttC>                 residual_;
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
             + residual_.evaluate( g, mu );
    }

    cusfloat evaluate_G0      ( int g, cusfloat mu ) const { return g0_.evaluate( g, mu ); }
    cusfloat evaluate_residual( int g, cusfloat mu ) const { return residual_.evaluate( g, mu ); }

    int num_gauss() const { return residual_.num_gauss(); }

private:
    BetaFoldedResidual<dGdttxC>                residual_;
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
             + residual_.evaluate( g, mu );
    }

    cusfloat evaluate_G0      ( int g, cusfloat mu ) const { return g0_.evaluate( g, mu ); }
    cusfloat evaluate_residual( int g, cusfloat mu ) const { return residual_.evaluate( g, mu ); }

    int num_gauss() const { return residual_.num_gauss(); }

private:
    BetaFoldedResidual<dGdttxxC>               residual_;
    G0Gauss3AlphaEvaluator<TDT, A0T, A1T, A2T> g0_;
};
