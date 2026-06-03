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

// Standard library
#include <cmath>
#include <iomanip>
#include <iostream>

// Local modules
#include "../src/config.hpp"
#include "../src/math/gauss_kronrod_t.hpp"
#include "../src/math/math_constants.hpp"

// =============================================================================
// Tolerance: conservative enough for both float (SIMPLE_PREC) and double.
//   TOL_SUM  — weight sum and smooth-function integration
//   TOL_POLY — polynomial exactness (accounts for float node truncation at
//               high degree: ~1e-5 relative; 1e-3 is very conservative)
//   TOL_SYM  — node anti-symmetry / weight symmetry
// =============================================================================
static constexpr cusfloat TOL_SUM  = static_cast<cusfloat>( 1e-4 );
static constexpr cusfloat TOL_POLY = static_cast<cusfloat>( 1e-10 );
static constexpr cusfloat TOL_SYM  = static_cast<cusfloat>( 1e-5 );


// =============================================================================
// Helper: integrate f over [a, b] using the NGK-point Kronrod rule.
//   t(x) = mid + half*(x),  dt = half*dx
// =============================================================================
template<int NGK, typename Fcn>
cusfloat gk_integrate( Fcn f, cusfloat a, cusfloat b )
{
    const cusfloat half   = static_cast<cusfloat>( 0.5 ) * ( b - a );
    const cusfloat mid    = static_cast<cusfloat>( 0.5 ) * ( a + b );
    cusfloat       result = static_cast<cusfloat>( 0.0 );
    for ( int i = 0; i < NGK; ++i )
    {
        result += GaussKronrodT<NGK>::weights_x[i]
                * f( mid + half * GaussKronrodT<NGK>::roots_x[i] );
    }
    return result * half;
}

// Helper: x^n evaluated via double then cast back to cusfloat.
// Avoids float-specific std::pow overload resolution issues.
static cusfloat pow_cf( cusfloat x, int n )
{
    return static_cast<cusfloat>( std::pow( static_cast<double>( x ), n ) );
}


// =============================================================================
// Test 1 — sum of Kronrod weights equals 2  (= ∫_{-1}^{1} 1 dx)
// =============================================================================
bool test_weight_sum( )
{
    bool pass = true;

    auto check = [&]( const char* name, int n, const cusfloat* w )
    {
        cusfloat s = static_cast<cusfloat>( 0.0 );
        for ( int i = 0; i < n; ++i ) s += w[i];
        const cusfloat err = std::abs( s - static_cast<cusfloat>( 2.0 ) );
        if ( err > TOL_SUM )
        {
            std::cerr << "  [" << name << "] weight sum = " << s
                      << "  (expected 2.0, err = " << err << ")\n";
            pass = false;
        }
    };

    check( "G7K15",  15, GaussKronrodT<15>::weights_x );
    check( "G10K21", 21, GaussKronrodT<21>::weights_x );

    return pass;
}


// =============================================================================
// Test 2 — node anti-symmetry and weight symmetry
//   roots_x[i]   == -roots_x[N-1-i]   (odd function on [-1,1])
//   weights_x[i] ==  weights_x[N-1-i] (even function on [-1,1])
// =============================================================================
template<int NGK>
bool test_symmetry_rule( const char* name )
{
    bool pass = true;
    const cusfloat* x = GaussKronrodT<NGK>::roots_x;
    const cusfloat* w = GaussKronrodT<NGK>::weights_x;

    for ( int i = 0; i < NGK / 2; ++i )
    {
        const int j = NGK - 1 - i;

        const cusfloat node_err   = std::abs( x[i] + x[j] );
        const cusfloat weight_err = std::abs( w[i] - w[j] );

        if ( node_err > TOL_SYM )
        {
            std::cerr << "  [" << name << "] node anti-symmetry failed at i=" << i
                      << ":  x[i]+x[j] = " << ( x[i] + x[j] ) << "\n";
            pass = false;
        }
        if ( weight_err > TOL_SYM )
        {
            std::cerr << "  [" << name << "] weight symmetry failed at i=" << i
                      << ":  w[i]-w[j] = " << ( w[i] - w[j] ) << "\n";
            pass = false;
        }
    }
    return pass;
}

bool test_symmetry( )
{
    bool pass = true;
    pass &= test_symmetry_rule<15>( "G7K15"  );
    pass &= test_symmetry_rule<21>( "G10K21" );
    return pass;
}


// =============================================================================
// Test 3 — polynomial exactness
//
// GaussKronrodT<N> is exact for polynomials of degree ≤ (2N-1)/2 + N - 1,
// specifically:
//   G3K7   → exact for degree ≤ 13
//   G7K15  → exact for degree ≤ 29
//   G10K21 → exact for degree ≤ 41
//
// Exact values:  ∫_{-1}^{1} x^(2k) dx = 2 / (2k+1)
//                ∫_{-1}^{1} x^(2k+1) dx = 0  (odd integrand)
// =============================================================================
bool test_polynomial_exactness( )
{
    bool pass = true;

    auto check = [&]( const char* name, cusfloat result, cusfloat exact, int deg )
    {
        const cusfloat err = std::abs( result - exact );
        const cusfloat ref = std::abs( exact ) > static_cast<cusfloat>( 1e-30 )
                           ? std::abs( exact )
                           : static_cast<cusfloat>( 1.0 );
        if ( err / ref > TOL_POLY )
        {
            std::cerr << "  [" << name << "] x^" << deg
                      << "  result=" << result << "  exact=" << exact
                      << "  relerr=" << err / ref << "\n";
            pass = false;
        }
    };

    // --- G7K15 (exact ≤ degree 29) ------------------------------------------
    {
        check( "G7K15",
               gk_integrate<15>( []( cusfloat x ){ return pow_cf( x, 14 ); },
                                  static_cast<cusfloat>( -1 ), static_cast<cusfloat>( 1 ) ),
               static_cast<cusfloat>( 2.0 / 15.0 ), 14 );

        check( "G7K15",
               gk_integrate<15>( []( cusfloat x ){ return pow_cf( x, 15 ); },
                                  static_cast<cusfloat>( -1 ), static_cast<cusfloat>( 1 ) ),
               static_cast<cusfloat>( 0.0 ), 15 );
    }

    // --- G10K21 (exact ≤ degree 41) -----------------------------------------
    {
        check( "G10K21",
               gk_integrate<21>( []( cusfloat x ){ return pow_cf( x, 20 ); },
                                  static_cast<cusfloat>( -1 ), static_cast<cusfloat>( 1 ) ),
               static_cast<cusfloat>( 2.0 / 21.0 ), 20 );

        check( "G10K21",
               gk_integrate<21>( []( cusfloat x ){ return pow_cf( x, 21 ); },
                                  static_cast<cusfloat>( -1 ), static_cast<cusfloat>( 1 ) ),
               static_cast<cusfloat>( 0.0 ), 21 );
    }

    return pass;
}


// =============================================================================
// Test 4 — smooth-function integration using G7K15 on general intervals
//
//   (a)  ∫_0^1    exp(x) dx   = e - 1     ≈ 1.7182818284590452
//   (b)  ∫_0^π    sin(x) dx   = 2
//   (c)  ∫_1^2    1/x    dx   = ln(2)     ≈ 0.6931471805599453
//   (d)  ∫_0^π/2  cos(x) dx   = 1
// =============================================================================
bool test_smooth_functions( )
{
    bool pass = true;

    auto check = [&]( const char* label, cusfloat result, cusfloat exact )
    {
        const cusfloat ref = std::abs( exact ) > static_cast<cusfloat>( 1e-30 )
                           ? std::abs( exact )
                           : static_cast<cusfloat>( 1.0 );
        const cusfloat relerr = std::abs( result - exact ) / ref;
        if ( relerr > TOL_SUM )
        {
            std::cerr << "  [G7K15] " << label
                      << "  result=" << result << "  exact=" << exact
                      << "  relerr=" << relerr << "\n";
            pass = false;
        }
    };

    // (a) ∫_0^1 exp(x) dx = e - 1
    check( "exp(x) on [0,1]",
           gk_integrate<15>( []( cusfloat x ){ return static_cast<cusfloat>( std::exp( static_cast<double>( x ) ) ); },
                              static_cast<cusfloat>( 0 ), static_cast<cusfloat>( 1 ) ),
           static_cast<cusfloat>( std::exp( 1.0 ) - 1.0 ) );

    // (b) ∫_0^π sin(x) dx = 2
    check( "sin(x) on [0,pi]",
           gk_integrate<15>( []( cusfloat x ){ return static_cast<cusfloat>( std::sin( static_cast<double>( x ) ) ); },
                              static_cast<cusfloat>( 0 ), PI ),
           static_cast<cusfloat>( 2.0 ) );

    // (c) ∫_1^2 1/x dx = ln(2)
    check( "1/x on [1,2]",
           gk_integrate<15>( []( cusfloat x ){ return static_cast<cusfloat>( 1.0 / static_cast<double>( x ) ); },
                              static_cast<cusfloat>( 1 ), static_cast<cusfloat>( 2 ) ),
           static_cast<cusfloat>( std::log( 2.0 ) ) );

    // (d) ∫_0^π/2 cos(x) dx = 1
    check( "cos(x) on [0,pi/2]",
           gk_integrate<15>( []( cusfloat x ){ return static_cast<cusfloat>( std::cos( static_cast<double>( x ) ) ); },
                              static_cast<cusfloat>( 0 ), PI / static_cast<cusfloat>( 2 ) ),
           static_cast<cusfloat>( 1.0 ) );

    return pass;
}


// =============================================================================
// main
// =============================================================================
int main( )
{
    std::cout << std::scientific << std::setprecision( 8 );
    std::cout << "--- Gauss-Kronrod quadrature tests ---\n";

    int rc = 0;

    bool pass;

    pass = test_weight_sum( );
    std::cout << "test_weight_sum           : " << ( pass ? "PASS" : "FAIL" ) << "\n";
    if ( !pass ) rc = 1;

    pass = test_symmetry( );
    std::cout << "test_symmetry             : " << ( pass ? "PASS" : "FAIL" ) << "\n";
    if ( !pass ) rc = 1;

    pass = test_polynomial_exactness( );
    std::cout << "test_polynomial_exactness : " << ( pass ? "PASS" : "FAIL" ) << "\n";
    if ( !pass ) rc = 1;

    pass = test_smooth_functions( );
    std::cout << "test_smooth_functions     : " << ( pass ? "PASS" : "FAIL" ) << "\n";
    if ( !pass ) rc = 1;

    std::cout << "--- " << ( rc == 0 ? "ALL PASS" : "SOME FAILURES" ) << " ---\n";
    return rc;
}
