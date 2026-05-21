// Test suite for the Gauss-point-optimised time-domain evaluators implemented
// in src/green/time_domain_evaluator_gauss.hpp.
//
// Structural tests (no data files required):
//   1. test_gauss_construction    – GaussEval_* objects construct without
//      throwing and report the correct Gauss count.
//   2. test_residual_fold_finite  – BetaFoldedResidual produces finite values
//      at every (g, mu) in a representative grid.
//   3. test_G0_F0_F1_finite       – dGdtt_G0_F0/F1 are finite and match the
//      factored form of dGdtt_G0_basis (cross-check of the factorisation).
//   4. test_G0_factor_chain       – F0/F1 of the x-derivative types satisfy
//      F0_x = (-β²/4)·F0_parent (same chain rule as the original basis).
//   5. test_gauss_vs_original     – GaussEval_*::evaluate(g, mu) agrees with
//      the original free function (eval_dGdt, …) to within floating-point
//      round-off; this is the primary correctness check.
//   6. test_component_sum         – evaluate_G0 + evaluate_residual ==
//      evaluate for every evaluator type.
//   7. test_all_evaluators_finite – all six GaussEval_* types produce finite
//      results across a representative (beta, mu) grid.
//
// Reference-data tests (require the same 18 .dat files used by
// test_time_domain_evaluator, passed via argv[1..18]):
//   8. test_gauss_against_reference – compare each GaussEval_*::evaluate
//      against the precomputed reference values with the same tolerances as
//      the original test.  This verifies that the precomputed optimisation
//      does not degrade accuracy.

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/config.hpp"
#include "../../src/green/time_domain_evaluator_gauss.hpp"


// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
static void check( bool cond, const char* msg )
{
    if ( !cond )
    {
        std::cerr << "FAIL " << msg << std::endl;
        throw std::runtime_error( msg );
    }
}

template<typename T>
static bool near( T a, T b,
                  T rel_tol = static_cast<T>( 1e-4 ),
                  T abs_tol = static_cast<T>( 1e-10 ) )
{
    return std::abs( a - b ) <= abs_tol + rel_tol * ( std::abs( a ) + std::abs( b ) );
}

// Representative set of Gauss beta values used throughout the tests.
// These span the full beta domain [0, 30] of the dGdtC database.
static std::vector<cusfloat> make_test_betas()
{
    return { 0.5f, 2.0f, 5.0f, 8.0f, 12.0f, 17.0f, 22.0f, 27.0f };
}


// ---------------------------------------------------------------------------
// Test 1 – construction and Gauss count
// ---------------------------------------------------------------------------
static void test_gauss_construction()
{
    const std::vector<cusfloat> betas = make_test_betas();
    const int ng = static_cast<int>( betas.size() );

    GaussEval_dGdt    ev_dGdt   ( betas );
    GaussEval_dGdtx   ev_dGdtx  ( betas );
    GaussEval_dGdtxx  ev_dGdtxx ( betas );
    GaussEval_dGdtt   ev_dGdtt  ( betas );
    GaussEval_dGdttx  ev_dGdttx ( betas );
    GaussEval_dGdttxx ev_dGdttxx( betas );
    GaussEval_dGdttt  ev_dGdttt ( betas );

    check( ev_dGdt.num_gauss()    == ng, "test_gauss_construction: dGdt    wrong Gauss count" );
    check( ev_dGdtx.num_gauss()   == ng, "test_gauss_construction: dGdtx   wrong Gauss count" );
    check( ev_dGdtxx.num_gauss()  == ng, "test_gauss_construction: dGdtxx  wrong Gauss count" );
    check( ev_dGdtt.num_gauss()   == ng, "test_gauss_construction: dGdtt   wrong Gauss count" );
    check( ev_dGdttx.num_gauss()  == ng, "test_gauss_construction: dGdttx  wrong Gauss count" );
    check( ev_dGdttxx.num_gauss() == ng, "test_gauss_construction: dGdttxx wrong Gauss count" );
    check( ev_dGdttt.num_gauss()  == ng, "test_gauss_construction: dGdttt  wrong Gauss count" );

    std::cout << "PASS test_gauss_construction" << std::endl;
}


// ---------------------------------------------------------------------------
// Test 2 – BetaFoldedResidual produces finite values
// ---------------------------------------------------------------------------
static void test_residual_fold_finite()
{
    const std::vector<cusfloat> betas = make_test_betas();
    BetaFoldedResidual<dGdtC> res( betas );

    const cusfloat mus[] = { 1e-4f, 1e-3f, 1e-2f, 5e-2f, 2e-1f };
    const int nm = 5;
    const int ng = static_cast<int>( betas.size() );

    for ( int g = 0; g < ng; ++g )
    {
        for ( int im = 0; im < nm; ++im )
        {
            const cusfloat val    = res.evaluate( g, mus[im] );
            check( std::isfinite( val ),
                   "test_residual_fold_finite: non-finite folded residual" );
        }
    }

    std::cout << "PASS test_residual_fold_finite" << std::endl;
}


// ---------------------------------------------------------------------------
// Test 3 – G0 factorisation cross-check for dGdtt
//
// Verify that for each (beta, alpha, mu):
//   dGdtt_G0_F0(beta, alpha) * exp(-beta^2/4 * mu)
//   + dGdtt_G0_F1(beta, alpha) * mu * exp(-beta^2/4 * mu)
//   == dGdtt_G0_basis(beta, mu, alpha)
//
// This confirms that the F0/F1 factorisation is correct for the only basis
// type where F1 ≠ 0 in the 1-alpha family (dGdt F1=0, so that case is
// trivially correct).
// ---------------------------------------------------------------------------
static void test_G0_F0_F1_finite()
{
    const cusfloat betas[]  = { 0.5f, 2.0f, 5.0f, 10.0f, 18.0f };
    const cusfloat mus[]    = { 1e-4f, 1e-2f, 0.3f };
    const cusfloat alphas[] = { 0.0f, 0.5f, 2.0f };
    constexpr int nb = 5, nm = 3, na = 3;
    constexpr cusfloat tol = static_cast<cusfloat>( 1e-5 );

    for ( int ib = 0; ib < nb; ++ib )
    for ( int im = 0; im < nm; ++im )
    for ( int ia = 0; ia < na; ++ia )
    {
        const cusfloat beta  = betas[ib];
        const cusfloat mu    = mus[im];
        const cusfloat alpha = alphas[ia];

        const cusfloat c    = beta * beta / static_cast<cusfloat>( 4.0 );
        const cusfloat e    = std::exp( -c * mu );
        const cusfloat F0   = dGdtt_G0_F0( beta, alpha );
        const cusfloat F1   = dGdtt_G0_F1( beta, alpha );

        check( std::isfinite( F0 ), "test_G0_F0_F1_finite: F0 is not finite" );
        check( std::isfinite( F1 ), "test_G0_F0_F1_finite: F1 is not finite" );

        // Reconstructed value from factored form
        const cusfloat reconstructed = ( F0 + F1 * mu ) * e;

        // Reference from original basis function
        const cusfloat reference = dGdtt_G0_basis( beta, mu, alpha );

        check( near( reconstructed, reference, tol ),
               "test_G0_F0_F1_finite: F0/F1 factorisation does not match dGdtt_G0_basis" );
    }

    std::cout << "PASS test_G0_F0_F1_finite" << std::endl;
}


// ---------------------------------------------------------------------------
// Test 4 – x-derivative chain for F0 and F1
//   F0_dGdtx(beta, alpha)   == (-beta^2/4) * F0_dGdt(beta, alpha)
//   F0_dGdttx(beta, alpha)  == (-beta^2/4) * F0_dGdtt(beta, alpha)
//   F0_dGdttxx(beta, alpha) == (-beta^2/4) * F0_dGdttx(beta, alpha)
//   (same for F1)
// ---------------------------------------------------------------------------
static void test_G0_factor_chain()
{
    const cusfloat betas[]  = { 0.5f, 2.0f, 5.0f, 15.0f };
    const cusfloat alphas[] = { 0.0f, 0.5f, 1.5f };
    constexpr int nb = 4, na = 3;
    constexpr cusfloat tol = static_cast<cusfloat>( 1e-6 );

    for ( int ib = 0; ib < nb; ++ib )
    for ( int ia = 0; ia < na; ++ia )
    {
        const cusfloat beta   = betas[ib];
        const cusfloat alpha  = alphas[ia];
        const cusfloat factor = -beta * beta / static_cast<cusfloat>( 4.0 );

        // dGdt → dGdtx chain
        check( near( dGdtx_G0_F0( beta, alpha ), factor * dGdt_G0_F0( beta, alpha ), tol ),
               "test_G0_factor_chain: F0_dGdtx != (-β²/4) * F0_dGdt" );
        check( near( dGdtx_G0_F1( beta, alpha ), factor * dGdt_G0_F1( beta, alpha ), tol ),
               "test_G0_factor_chain: F1_dGdtx != (-β²/4) * F1_dGdt" );

        // dGdtx → dGdtxx chain
        check( near( dGdtxx_G0_F0( beta, alpha ), factor * dGdtx_G0_F0( beta, alpha ), tol ),
               "test_G0_factor_chain: F0_dGdtxx != (-β²/4) * F0_dGdtx" );
        check( near( dGdtxx_G0_F1( beta, alpha ), factor * dGdtx_G0_F1( beta, alpha ), tol ),
               "test_G0_factor_chain: F1_dGdtxx != (-β²/4) * F1_dGdtx" );

        // dGdtt → dGdttx chain
        check( near( dGdttx_G0_F0( beta, alpha ), factor * dGdtt_G0_F0( beta, alpha ), tol ),
               "test_G0_factor_chain: F0_dGdttx != (-β²/4) * F0_dGdtt" );
        check( near( dGdttx_G0_F1( beta, alpha ), factor * dGdtt_G0_F1( beta, alpha ), tol ),
               "test_G0_factor_chain: F1_dGdttx != (-β²/4) * F1_dGdtt" );

        // dGdttx → dGdttxx chain
        check( near( dGdttxx_G0_F0( beta, alpha ), factor * dGdttx_G0_F0( beta, alpha ), tol ),
               "test_G0_factor_chain: F0_dGdttxx != (-β²/4) * F0_dGdttx" );
        check( near( dGdttxx_G0_F1( beta, alpha ), factor * dGdttx_G0_F1( beta, alpha ), tol ),
               "test_G0_factor_chain: F1_dGdttxx != (-β²/4) * F1_dGdttx" );
    }

    std::cout << "PASS test_G0_factor_chain" << std::endl;
}


// ---------------------------------------------------------------------------
// Test 4c – G0 three-term (quadratic μ) factorisation for dGdttt
//
// Verify that for each (beta, mu):
//   (F0 + F1*μ + F2*μ²) * exp(-β²μ/4)
//   == dGdttt_G0_basis(beta, mu, alpha=0)
//
// This confirms that dGdttt_G0_F0/F1/F2 correctly factor out the μ
// dependence, which is the assumption that G0Gauss1Alpha3TermEvaluator
// relies on.
// ---------------------------------------------------------------------------
static void test_G0_3term_factorisation()
{
    const cusfloat betas[]  = { 0.5f, 2.0f, 5.0f, 10.0f, 18.0f };
    const cusfloat mus[]    = { 1e-4f, 1e-2f, 0.3f, 0.9f };
    constexpr int nb = 5, nm = 4;
    constexpr cusfloat tol = static_cast<cusfloat>( 1e-5 );

    // dGdttt has a single alpha shift == 0
    const cusfloat alpha = static_cast<cusfloat>( 0.0 );

    for ( int ib = 0; ib < nb; ++ib )
    for ( int im = 0; im < nm; ++im )
    {
        const cusfloat beta = betas[ib];
        const cusfloat mu   = mus[im];

        const cusfloat c    = beta * beta / static_cast<cusfloat>( 4.0 );
        const cusfloat e    = std::exp( -c * mu );
        const cusfloat F0   = dGdttt_G0_F0( beta, alpha );
        const cusfloat F1   = dGdttt_G0_F1( beta, alpha );
        const cusfloat F2   = dGdttt_G0_F2( beta, alpha );

        check( std::isfinite( F0 ), "test_G0_3term_factorisation: F0 is not finite" );
        check( std::isfinite( F1 ), "test_G0_3term_factorisation: F1 is not finite" );
        check( std::isfinite( F2 ), "test_G0_3term_factorisation: F2 is not finite" );

        // Reconstruct from factored form using Horner
        const cusfloat reconstructed = ( F0 + mu * ( F1 + mu * F2 ) ) * e;

        // Reference from original basis function
        const cusfloat reference = dGdttt_G0_basis( beta, mu, alpha );

        check( near( reconstructed, reference, tol ),
               "test_G0_3term_factorisation: F0/F1/F2 factorisation does not match dGdttt_G0_basis" );
    }

    std::cout << "PASS test_G0_3term_factorisation" << std::endl;
}


// ---------------------------------------------------------------------------
// Test 5 – GaussEval_*::evaluate agrees with the original free functions
//
// GaussEval_*::evaluate(g, mu) must produce exactly (up to floating-point
// round-off) the same result as eval_dGdt(beta_g, mu), etc., because both
// compute the same G0 + Residual decomposition.  We accept a tight relative
// tolerance of 1e-4 to allow for accumulated single-precision rounding in the
// folded coefficient storage.
// ---------------------------------------------------------------------------
static void test_gauss_vs_original()
{
    const std::vector<cusfloat> betas = make_test_betas();
    const int ng = static_cast<int>( betas.size() );

    GaussEval_dGdt    ev_dGdt   ( betas );
    GaussEval_dGdtx   ev_dGdtx  ( betas );
    GaussEval_dGdtxx  ev_dGdtxx ( betas );
    GaussEval_dGdtt   ev_dGdtt  ( betas );
    GaussEval_dGdttx  ev_dGdttx ( betas );
    GaussEval_dGdttxx ev_dGdttxx( betas );
    GaussEval_dGdttt  ev_dGdttt ( betas );

    const cusfloat mus[] = { 1e-4f, 1e-3f, 1e-2f, 0.1f };
    constexpr int nm = 4;

    // Tight tolerance: the two code paths differ only in how they accumulate
    // the Chebyshev sums, so agreement should be at the level of floating-
    // point round-off.  1e-4 relative is generous for single precision.
    constexpr cusfloat rel_tol = static_cast<cusfloat>( 1e-4 );
    constexpr cusfloat abs_tol = static_cast<cusfloat>( 1e-10 );

    for ( int g = 0; g < ng; ++g )
    {
        const cusfloat beta = betas[static_cast<std::size_t>( g )];

        for ( int im = 0; im < nm; ++im )
        {
            const cusfloat mu = mus[im];

            const cusfloat ref_dGdt    = eval_dGdt   ( beta, mu );
            const cusfloat ref_dGdtx   = eval_dGdtx  ( beta, mu );
            const cusfloat ref_dGdtxx  = eval_dGdtxx ( beta, mu );
            const cusfloat ref_dGdtt   = eval_dGdtt  ( beta, mu );
            const cusfloat ref_dGdttx  = eval_dGdttx ( beta, mu );
            const cusfloat ref_dGdttxx = eval_dGdttxx( beta, mu );

            check( near( ev_dGdt.evaluate   ( g, mu ), ref_dGdt,    rel_tol, abs_tol ),
                   "test_gauss_vs_original: dGdt mismatch" );
            check( near( ev_dGdtx.evaluate  ( g, mu ), ref_dGdtx,   rel_tol, abs_tol ),
                   "test_gauss_vs_original: dGdtx mismatch" );
            check( near( ev_dGdtxx.evaluate ( g, mu ), ref_dGdtxx,  rel_tol, abs_tol ),
                   "test_gauss_vs_original: dGdtxx mismatch" );
            check( near( ev_dGdtt.evaluate  ( g, mu ), ref_dGdtt,   rel_tol, abs_tol ),
                   "test_gauss_vs_original: dGdtt mismatch" );
            check( near( ev_dGdttx.evaluate ( g, mu ), ref_dGdttx,  rel_tol, abs_tol ),
                   "test_gauss_vs_original: dGdttx mismatch" );
            check( near( ev_dGdttxx.evaluate( g, mu ), ref_dGdttxx, rel_tol, abs_tol ),
                   "test_gauss_vs_original: dGdttxx mismatch" );

            const cusfloat ref_dGdttt = eval_dGdttt( beta, mu );
            check( near( ev_dGdttt.evaluate( g, mu ), ref_dGdttt, rel_tol, abs_tol ),
                   "test_gauss_vs_original: dGdttt mismatch" );
        }
    }

    std::cout << "PASS test_gauss_vs_original" << std::endl;
}


// ---------------------------------------------------------------------------
// Test 6 – evaluate_G0 + evaluate_residual == evaluate
// ---------------------------------------------------------------------------
static void test_component_sum()
{
    const std::vector<cusfloat> betas = make_test_betas();
    const int ng = static_cast<int>( betas.size() );

    GaussEval_dGdt    ev_dGdt   ( betas );
    GaussEval_dGdtx   ev_dGdtx  ( betas );
    GaussEval_dGdtxx  ev_dGdtxx ( betas );
    GaussEval_dGdtt   ev_dGdtt  ( betas );
    GaussEval_dGdttx  ev_dGdttx ( betas );
    GaussEval_dGdttxx ev_dGdttxx( betas );
    GaussEval_dGdttt  ev_dGdttt ( betas );

    const cusfloat mus[] = { 1e-3f, 1e-2f, 0.1f };
    constexpr int nm = 3;
    constexpr cusfloat abs_tol = static_cast<cusfloat>( 1e-10 );

    for ( int g = 0; g < ng; ++g )
    {
        for ( int im = 0; im < nm; ++im )
        {
            const cusfloat mu = mus[im];

            // G0 + residual must equal the total to within machine epsilon
            auto check_sum = [&]( auto& ev, const char* name )
            {
                const cusfloat total    = ev.evaluate          ( g, mu );
                const cusfloat g0       = ev.evaluate_G0       ( g, mu );
                const cusfloat residual = ev.evaluate_residual ( g, mu );
                check( std::abs( total - ( g0 + residual ) ) <= abs_tol,
                       name );
            };

            check_sum( ev_dGdt,    "test_component_sum: dGdt    G0+residual != total" );
            check_sum( ev_dGdtx,   "test_component_sum: dGdtx   G0+residual != total" );
            check_sum( ev_dGdtxx,  "test_component_sum: dGdtxx  G0+residual != total" );
            check_sum( ev_dGdtt,   "test_component_sum: dGdtt   G0+residual != total" );
            check_sum( ev_dGdttx,  "test_component_sum: dGdttx  G0+residual != total" );
            check_sum( ev_dGdttxx, "test_component_sum: dGdttxx G0+residual != total" );
            check_sum( ev_dGdttt,  "test_component_sum: dGdttt  G0+residual != total" );
        }
    }

    std::cout << "PASS test_component_sum" << std::endl;
}


// ---------------------------------------------------------------------------
// Test 7 – all evaluators produce finite values across a wide grid
// ---------------------------------------------------------------------------
static void test_all_evaluators_finite()
{
    // All evaluators except dGdttt can use the full beta range [0, 30]
    const std::vector<cusfloat> betas = { 0.5f, 2.0f, 5.0f, 10.0f, 15.0f, 22.0f, 27.0f };

    GaussEval_dGdt    ev_dGdt   ( betas );
    GaussEval_dGdtx   ev_dGdtx  ( betas );
    GaussEval_dGdtxx  ev_dGdtxx ( betas );
    GaussEval_dGdtt   ev_dGdtt  ( betas );
    GaussEval_dGdttx  ev_dGdttx ( betas );
    GaussEval_dGdttxx ev_dGdttxx( betas );

    // dGdttt has beta domain [0, 19] — use a separate, smaller set
    const std::vector<cusfloat> betas_ttt = { 0.5f, 2.0f, 5.0f, 8.0f, 12.0f, 17.0f };
    GaussEval_dGdttt  ev_dGdttt ( betas_ttt );

    const cusfloat mus[] = { 1e-4f, 1e-3f, 1e-2f, 5e-2f, 2e-1f };
    constexpr int nm = 5;
    const int ng     = static_cast<int>( betas.size() );
    const int ng_ttt = static_cast<int>( betas_ttt.size() );

    for ( int g = 0; g < ng; ++g )
    {
        for ( int im = 0; im < nm; ++im )
        {
            const cusfloat mu = mus[im];

            check( std::isfinite( ev_dGdt.evaluate   ( g, mu ) ), "test_all_evaluators_finite: dGdt"    );
            check( std::isfinite( ev_dGdtx.evaluate  ( g, mu ) ), "test_all_evaluators_finite: dGdtx"   );
            check( std::isfinite( ev_dGdtxx.evaluate ( g, mu ) ), "test_all_evaluators_finite: dGdtxx"  );
            check( std::isfinite( ev_dGdtt.evaluate  ( g, mu ) ), "test_all_evaluators_finite: dGdtt"   );
            check( std::isfinite( ev_dGdttx.evaluate ( g, mu ) ), "test_all_evaluators_finite: dGdttx"  );
            check( std::isfinite( ev_dGdttxx.evaluate( g, mu ) ), "test_all_evaluators_finite: dGdttxx" );
        }
    }
    for ( int g = 0; g < ng_ttt; ++g )
    {
        for ( int im = 0; im < nm; ++im )
        {
            const cusfloat mu = mus[im];
            check( std::isfinite( ev_dGdttt.evaluate( g, mu ) ), "test_all_evaluators_finite: dGdttt" );
        }
    }

    std::cout << "PASS test_all_evaluators_finite" << std::endl;
}


// ---------------------------------------------------------------------------
// Test 8 – compare GaussEval_* against the same reference .dat files used by
//           test_time_domain_evaluator.
//
// Because the Gauss evaluators must reproduce the original evaluators exactly
// (up to round-off), the same reference tolerances apply.
//
// The function accepts a method pointer for GaussEval so that the same helper
// can test the full, G0-only, and residual-only variants.
// ---------------------------------------------------------------------------
struct RefTolerance { double abs_tol; double rel_tol; };

static constexpr RefTolerance TOL_dGdt    = { 5e-4,  5e-4 };
static constexpr RefTolerance TOL_dGdtx   = { 5e-2,  5e-4 };
static constexpr RefTolerance TOL_dGdtxx  = { 5e+1,  5e-4 };
static constexpr RefTolerance TOL_dGdtt   = { 2e-1,  5e-4 };
static constexpr RefTolerance TOL_dGdttx  = { 5e+0,  5e-4 };
static constexpr RefTolerance TOL_dGdttxx = { 5e+1,  5e-3 };

static constexpr RefTolerance TOL_dGdt_residual    = { 5e-4,  1e+1 };
static constexpr RefTolerance TOL_dGdtx_residual   = { 5e-2,  1e+3 };
static constexpr RefTolerance TOL_dGdtxx_residual  = { 5e+1,  1e+6 };
static constexpr RefTolerance TOL_dGdtt_residual   = { 2e-1,  1e+2 };
static constexpr RefTolerance TOL_dGdttx_residual  = { 5e+0,  1e+2 };
static constexpr RefTolerance TOL_dGdttxx_residual = { 5e+1,  1e+6 };
// dGdttt: same domain as dGdtt (beta_max=19), use matching tolerances
static constexpr RefTolerance TOL_dGdttt           = { 5e+0,  5e-4 };  // 3rd t-deriv: larger fit error than dGdtt
static constexpr RefTolerance TOL_dGdttt_residual  = { 5e+0,  1e+2 };

// Generic Gauss-evaluator test: the evaluator instance is already constructed
// with the beta values read from the reference file; we build a new evaluator
// per test call to avoid polluting state.
template<typename GaussEvalT, typename EvalMethod>
static void test_gauss_against_reference(
    const std::string& filepath,
    const std::string& name,
    EvalMethod         eval_method,
    const RefTolerance tol,
    const std::string& diff_dir = "" )
{
    std::cout << "test_gauss_against_reference: " << name << " ..." << std::endl;

    std::ifstream fin( filepath );
    if ( !fin.is_open() )
        throw std::runtime_error( "cannot open reference file: " + filepath );

    int n_points = 0;
    fin >> n_points;
    if ( n_points <= 0 )
        throw std::runtime_error( "invalid n_points in reference file: " + filepath );

    // Read all (beta, mu, expected) triples
    std::vector<cusfloat> file_betas( static_cast<std::size_t>( n_points ) );
    std::vector<cusfloat> file_mus(   static_cast<std::size_t>( n_points ) );
    std::vector<cusfloat> file_exp(   static_cast<std::size_t>( n_points ) );

    for ( int i = 0; i < n_points; ++i )
    {
        double b, m, e;
        fin >> b >> m >> e;
        file_betas[static_cast<std::size_t>( i )] = static_cast<cusfloat>( b );
        file_mus  [static_cast<std::size_t>( i )] = static_cast<cusfloat>( m );
        file_exp  [static_cast<std::size_t>( i )] = static_cast<cusfloat>( e );
    }
    fin.close();

    // Collect unique sorted beta values to use as the Gauss set
    std::vector<cusfloat> unique_betas = file_betas;
    std::sort( unique_betas.begin(), unique_betas.end() );
    unique_betas.erase(
        std::unique( unique_betas.begin(), unique_betas.end() ),
        unique_betas.end() );

    // Construct the evaluator with these beta values
    GaussEvalT ev( unique_betas );

    // Build a beta -> Gauss-index lookup
    // (unique_betas is sorted so lower_bound works)
    auto gauss_index = [&]( cusfloat beta ) -> int {
        auto it = std::lower_bound( unique_betas.begin(), unique_betas.end(), beta );
        return static_cast<int>( std::distance( unique_betas.begin(), it ) );
    };

    // Open optional CSV diff file
    std::ofstream fcsv;
    if ( !diff_dir.empty() )
    {
        std::filesystem::create_directories( diff_dir );
        const std::string csv_path = diff_dir + "/" + name + "_gauss_diff.csv";
        fcsv.open( csv_path );
        if ( fcsv.is_open() )
            fcsv << "beta,mu,expected,computed,abs_err,rel_err,pass\n";
    }

    int    n_fail   = 0;
    double max_err  = 0.0;
    double max_rerr = 0.0;

    for ( int i = 0; i < n_points; ++i )
    {
        const std::size_t idx      = static_cast<std::size_t>( i );
        const cusfloat    beta     = file_betas[idx];
        const cusfloat    mu       = file_mus  [idx];
        const cusfloat    expected = file_exp  [idx];

        const int      g        = gauss_index( beta );
        const cusfloat computed = ( ev.*eval_method )( g, mu );

        const double abs_err  = std::abs( static_cast<double>( computed )
                                        - static_cast<double>( expected ) );
        const double rel_err  = abs_err / ( std::abs( static_cast<double>( expected ) ) + 1e-30 );
        const double threshold = tol.abs_tol + tol.rel_tol * std::abs( static_cast<double>( expected ) );
        const bool   passed   = ( abs_err <= threshold );

        if ( !passed )
        {
            ++n_fail;
            if ( n_fail <= 5 )
            {
                std::cerr << "  FAIL point " << i
                          << ": beta=" << static_cast<double>( beta )
                          << " mu="   << static_cast<double>( mu )
                          << " expected=" << static_cast<double>( expected )
                          << " computed=" << static_cast<double>( computed )
                          << " abs_err=" << abs_err
                          << " threshold=" << threshold << std::endl;
            }
        }
        max_err  = std::max( max_err,  abs_err );
        max_rerr = std::max( max_rerr, rel_err );

        if ( fcsv.is_open() )
        {
            fcsv << std::scientific << std::setprecision( 12 )
                 << static_cast<double>( beta ) << ","
                 << static_cast<double>( mu )   << ","
                 << static_cast<double>( expected ) << ","
                 << static_cast<double>( computed ) << ","
                 << abs_err << "," << rel_err << ","
                 << ( passed ? 1 : 0 ) << "\n";
        }
    }
    if ( fcsv.is_open() ) fcsv.close();

    const int max_allowed = std::max( 1, n_points / 100 );
    if ( n_fail > max_allowed )
    {
        std::cerr << "FAIL " << name
                  << ": " << n_fail << "/" << n_points
                  << " points exceed tolerance (max_abs_err=" << max_err
                  << " max_rel_err=" << max_rerr << ")" << std::endl;
        throw std::runtime_error( "test_gauss_against_reference failed for " + name );
    }

    std::cout << "PASS test_gauss_against_reference " << name
              << " (" << n_fail << " outliers / " << n_points
              << " points, max_abs_err=" << max_err
              << " max_rel_err=" << max_rerr << ")" << std::endl;
}


// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    // --- Structural tests (no data files) ---
    test_gauss_construction();
    test_residual_fold_finite();
    test_G0_F0_F1_finite();
    test_G0_factor_chain();
    test_G0_3term_factorisation();
    test_gauss_vs_original();
    test_component_sum();
    test_all_evaluators_finite();

    // Optional --diff-output <dir>: write per-test CSV diff files
    std::string diff_dir;
    for ( int i = 1; i < argc - 1; ++i )
    {
        if ( std::string( argv[i] ) == "--diff-output" )
        {
            diff_dir = argv[i + 1];
            std::cout << "NOTE: diff CSV files will be written to: " << diff_dir << std::endl;
            break;
        }
    }

    // --- Reference-data tests ---
    // argv[1..6]   = total (G0 + residual) reference files
    // argv[7..12]  = G0-only reference files
    // argv[13..18] = residual-only reference files
    // (same file set as test_time_domain_evaluator)

    // --- Reference-data tests ---
    // argv[1..7]   = total (G0 + residual) reference files
    // argv[8..14]  = G0-only reference files
    // argv[15..21] = residual-only reference files

    if ( argc < 8 )
    {
        std::cout << "NOTE: No reference data files provided (need 7 paths). "
                     "Skipping reference comparison tests." << std::endl;
    }
    else
    {
        test_gauss_against_reference<GaussEval_dGdt   >( argv[1], "dGdt",    &GaussEval_dGdt::evaluate,    TOL_dGdt,    diff_dir );
        test_gauss_against_reference<GaussEval_dGdtx  >( argv[2], "dGdtx",   &GaussEval_dGdtx::evaluate,   TOL_dGdtx,   diff_dir );
        test_gauss_against_reference<GaussEval_dGdtxx >( argv[3], "dGdtxx",  &GaussEval_dGdtxx::evaluate,  TOL_dGdtxx,  diff_dir );
        test_gauss_against_reference<GaussEval_dGdtt  >( argv[4], "dGdtt",   &GaussEval_dGdtt::evaluate,   TOL_dGdtt,   diff_dir );
        test_gauss_against_reference<GaussEval_dGdttx >( argv[5], "dGdttx",  &GaussEval_dGdttx::evaluate,  TOL_dGdttx,  diff_dir );
        test_gauss_against_reference<GaussEval_dGdttxx>( argv[6], "dGdttxx", &GaussEval_dGdttxx::evaluate, TOL_dGdttxx, diff_dir );
        test_gauss_against_reference<GaussEval_dGdttt >( argv[7], "dGdttt",  &GaussEval_dGdttt::evaluate,  TOL_dGdttt,  diff_dir );
    }

    if ( argc >= 15 )
    {
        // G0-only
        test_gauss_against_reference<GaussEval_dGdt   >( argv[8],  "dGdt_G0",    &GaussEval_dGdt::evaluate_G0,    TOL_dGdt,    diff_dir );
        test_gauss_against_reference<GaussEval_dGdtx  >( argv[9],  "dGdtx_G0",   &GaussEval_dGdtx::evaluate_G0,   TOL_dGdtx,   diff_dir );
        test_gauss_against_reference<GaussEval_dGdtxx >( argv[10], "dGdtxx_G0",  &GaussEval_dGdtxx::evaluate_G0,  TOL_dGdtxx,  diff_dir );
        test_gauss_against_reference<GaussEval_dGdtt  >( argv[11], "dGdtt_G0",   &GaussEval_dGdtt::evaluate_G0,   TOL_dGdtt,   diff_dir );
        test_gauss_against_reference<GaussEval_dGdttx >( argv[12], "dGdttx_G0",  &GaussEval_dGdttx::evaluate_G0,  TOL_dGdttx,  diff_dir );
        test_gauss_against_reference<GaussEval_dGdttxx>( argv[13], "dGdttxx_G0", &GaussEval_dGdttxx::evaluate_G0, TOL_dGdttxx, diff_dir );
        test_gauss_against_reference<GaussEval_dGdttt >( argv[14], "dGdttt_G0",  &GaussEval_dGdttt::evaluate_G0,  TOL_dGdttt,  diff_dir );
    }

    if ( argc >= 22 )
    {
        // Residual-only
        test_gauss_against_reference<GaussEval_dGdt   >( argv[15], "dGdt_residual",    &GaussEval_dGdt::evaluate_residual,    TOL_dGdt_residual,    diff_dir );
        test_gauss_against_reference<GaussEval_dGdtx  >( argv[16], "dGdtx_residual",   &GaussEval_dGdtx::evaluate_residual,   TOL_dGdtx_residual,   diff_dir );
        test_gauss_against_reference<GaussEval_dGdtxx >( argv[17], "dGdtxx_residual",  &GaussEval_dGdtxx::evaluate_residual,  TOL_dGdtxx_residual,  diff_dir );
        test_gauss_against_reference<GaussEval_dGdtt  >( argv[18], "dGdtt_residual",   &GaussEval_dGdtt::evaluate_residual,   TOL_dGdtt_residual,   diff_dir );
        test_gauss_against_reference<GaussEval_dGdttx >( argv[19], "dGdttx_residual",  &GaussEval_dGdttx::evaluate_residual,  TOL_dGdttx_residual,  diff_dir );
        test_gauss_against_reference<GaussEval_dGdttxx>( argv[20], "dGdttxx_residual", &GaussEval_dGdttxx::evaluate_residual, TOL_dGdttxx_residual, diff_dir );
        test_gauss_against_reference<GaussEval_dGdttt >( argv[21], "dGdttt_residual",  &GaussEval_dGdttt::evaluate_residual,  TOL_dGdttt_residual,  diff_dir );
    }

    std::cout << "All tests PASSED." << std::endl;
    return 0;
}
