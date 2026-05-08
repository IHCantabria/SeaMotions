
// Include general usage libraries
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/green/integrals_database.hpp"
#include "../../src/green/chebyshev_traits.hpp"
#include "../../src/math/chebyshev.hpp"


// =============================================================================
//  Synthetic 2D database for fold_database_2d unit test
// =============================================================================
//
//  Polynomial:
//    f(xm, ym) = A*T0(xm)*T0(ym)  +  B*T0(xm)*T1(ym)  +  C*T0(xm)*T2(ym)
//              + D*T1(xm)*T0(ym)  +  E*T1(xm)*T1(ym)
//              + F*T2(xm)*T0(ym)
//
//    A = 3.0, B = 2.0, C = 0.5, D = 1.5, E = 0.75, F = 1.0
//
//  Domain:  x ∈ [0, 1],  y ∈ [0, 1],  single block,  intervals_np = 1
//  Mapping: xm = 2*x - 1,  ym = 2*y - 1
//
//  ncx ordering (sorted, each group consecutive):  {0, 0, 0, 1, 1, 2}
//  ncy ordering (varies within each ncx group):    {0, 1, 2, 0, 1, 0}
//
//  After folding y at ym0 = 2*y0 - 1 the folded 1D polynomial is:
//    g(xm; ym0) = [A + B*ym0 + C*(2*ym0^2-1)] * T0(xm)
//               + [D + E*ym0]                  * T1(xm)
//               + F                             * T2(xm)
// =============================================================================

struct TestDB2D
{
    static constexpr int         max_ref_level               = 1;
    static constexpr int         intervals_np                = 1;
    static constexpr int         max_cheby_order             = 2;
    static constexpr int         blocks_np                   = 1;
    static constexpr int         blocks_np_f                 = 1;

    static constexpr std::size_t blocks_start[1]             = {0};
    static constexpr std::size_t blocks_coeffs_np[1]         = {6};
    static constexpr std::size_t blocks_max_cheby_order[1]   = {2};

    static constexpr bool        x_log_scale                 = false;
    static constexpr cusfloat    x_min_global                = 0.0;
    static constexpr cusfloat    x_max_global                = 1.0;
    static constexpr cusfloat    dx_min_region               = 1.0;
    static constexpr cusfloat    x_min_region[1]             = {0.0};
    static constexpr cusfloat    x_max_region[1]             = {1.0};
    static constexpr cusfloat    dx_region[1]                = {1.0};

    static constexpr bool        y_log_scale                 = false;
    static constexpr cusfloat    y_min_global                = 0.0;
    static constexpr cusfloat    y_max_global                = 1.0;
    static constexpr cusfloat    dy_min_region               = 1.0;
    static constexpr cusfloat    y_min_region[1]             = {0.0};
    static constexpr cusfloat    y_max_region[1]             = {1.0};
    static constexpr cusfloat    dy_region[1]                = {1.0};

    static constexpr bool        fcn_log_scale               = false;
    static constexpr std::size_t num_c                       = 6;
    static constexpr std::size_t num_cf                      = 6;

    // Polynomial coefficients [A, B, C, D, E, F]
    static constexpr cusfloat    c[6]    = {3.0, 2.0, 0.5, 1.5, 0.75, 1.0};
    // ncx sorted ascending; consecutive entries sharing the same ncx form a group
    static constexpr std::size_t ncx[6]  = {0, 0, 0, 1, 1, 2};
    // ncy varies within each ncx group
    static constexpr std::size_t ncy[6]  = {0, 1, 2, 0, 1, 0};

    // Mutable storage written by fold_database_2d
    inline static std::size_t   blocks_start_f[1];
    inline static std::size_t   blocks_coeffs_np_f[1];
    inline static std::size_t   blocks_max_cheby_order_f[1];
    inline static cusfloat      x_min_region_f[1];
    inline static cusfloat      x_max_region_f[1];
    inline static cusfloat      dx_region_f[1];
    inline static cusfloat      cf[6];
    inline static std::size_t   ncxf[6];
};

// ChebyshevTraits specialization for TestDB2D
// Must expose the same interface that fold_database_2d writes to.
template<>
struct ChebyshevTraits<TestDB2D>
{
    static constexpr int        intervals_np            = TestDB2D::intervals_np;
    static constexpr int        max_cheby_order         = TestDB2D::max_cheby_order;
    static constexpr bool       x_log_scale             = TestDB2D::x_log_scale;
    static constexpr cusfloat   x_min_global            = TestDB2D::x_min_global;
    static constexpr cusfloat   x_max_global            = TestDB2D::x_max_global;
    static constexpr cusfloat   dx_min_region           = TestDB2D::dx_min_region;

    inline static std::size_t*  blocks_start            = TestDB2D::blocks_start_f;
    inline static std::size_t*  blocks_coeffs_np        = TestDB2D::blocks_coeffs_np_f;
    inline static std::size_t*  blocks_max_cheby_order  = TestDB2D::blocks_max_cheby_order_f;
    inline static cusfloat*     x_min_region            = TestDB2D::x_min_region_f;
    inline static cusfloat*     x_max_region            = TestDB2D::x_max_region_f;
    inline static cusfloat*     dx_region               = TestDB2D::dx_region_f;
    inline static cusfloat*     coeffs                  = TestDB2D::cf;
    inline static std::size_t*  ncx                     = TestDB2D::ncxf;
};


// =============================================================================
//  Helper: direct evaluation of the L2C (1D) polynomial at water depth H
// =============================================================================
static cusfloat eval_L2C_direct( cusfloat H )
{
    cusfloat Hreg = std::log10( H );
    Hreg = std::min( std::max( Hreg, L2C::x_min_global + (cusfloat)1e-6 ),
                                     L2C::x_max_global - (cusfloat)1e-6 );

    std::size_t blk = 0;
    for ( std::size_t i = 0; i < static_cast<std::size_t>( L2C::intervals_np ); i++ )
    {
        if ( L2C::x_min_region[i] < Hreg && L2C::x_max_region[i] > Hreg )
        {
            blk = i;
            break;
        }
    }

    cusfloat h_map = 2.0 * ( Hreg - L2C::x_min_region[blk] ) / L2C::dx_region[blk] - 1.0;

    cusfloat poly[ L2C::max_cheby_order + 1 ];
    chebyshev_poly_upto_order( L2C::blocks_max_cheby_order[blk], h_map, poly );

    cusfloat result = 0.0;
    std::size_t bs = L2C::blocks_start[blk];
    std::size_t nt = L2C::blocks_coeffs_np[blk];
    for ( std::size_t j = bs; j < bs + nt; j++ )
    {
        result += L2C::c[j] * poly[ L2C::ncx[j] ];
    }
    return result;
}


// =============================================================================
//  Helper: direct evaluation of the L1C (3D) polynomial at (x, y, H)
// =============================================================================
static cusfloat eval_L1C_direct( cusfloat x, cusfloat y, cusfloat H )
{
    cusfloat Hreg = std::log10( H );
    Hreg = std::min( std::max( Hreg, L1C::z_min_global + (cusfloat)1e-6 ),
                                     L1C::z_max_global - (cusfloat)1e-6 );

    int blk = -1;
    for ( int b = 0; b < L1C::blocks_np; b++ )
    {
        if ( x    > L1C::x_min_region[b] && x    < L1C::x_max_region[b] &&
             y    > L1C::y_min_region[b] && y    < L1C::y_max_region[b] &&
             Hreg > L1C::z_min_region[b] && Hreg < L1C::z_max_region[b] )
        {
            blk = b;
            break;
        }
    }
    if ( blk < 0 )
        throw std::runtime_error( "eval_L1C_direct: no 3-D block contains the test point" );

    cusfloat xm = 2.0 * ( x    - L1C::x_min_region[blk] ) / L1C::dx_region[blk] - 1.0;
    cusfloat ym = 2.0 * ( y    - L1C::y_min_region[blk] ) / L1C::dy_region[blk] - 1.0;
    cusfloat zm = 2.0 * ( Hreg - L1C::z_min_region[blk] ) / L1C::dz_region[blk] - 1.0;

    cusfloat poly_x[ L1C::max_cheby_order + 1 ];
    cusfloat poly_y[ L1C::max_cheby_order + 1 ];
    cusfloat poly_z[ L1C::max_cheby_order + 1 ];
    std::size_t ord = L1C::blocks_max_cheby_order[blk];
    chebyshev_poly_upto_order( ord, xm, poly_x );
    chebyshev_poly_upto_order( ord, ym, poly_y );
    chebyshev_poly_upto_order( ord, zm, poly_z );

    std::size_t bs = L1C::blocks_start[blk];
    std::size_t nt = L1C::blocks_coeffs_np[blk];

    cusfloat result = 0.0;
    for ( std::size_t j = bs; j < bs + nt; j++ )
    {
        result += L1C::c[j]
                * poly_x[ L1C::ncx[j] ]
                * poly_y[ L1C::ncy[j] ]
                * poly_z[ L1C::ncz[j] ];
    }
    return result;
}


// =============================================================================
//  Helper: evaluate the 2-D folded polynomial stored in ChebyshevTraits<IDB>
//          at point (x, y).  Requires IDB to be a 3D database type whose
//          ChebyshevTraits has been populated by fold_database_3d.
// =============================================================================
template<typename IDB>
static cusfloat eval_folded_2d( cusfloat x, cusfloat y )
{
    constexpr cusfloat dx = ( IDB::x_max_global - IDB::x_min_global ) / IDB::intervals_np;
    constexpr cusfloat dy = ( IDB::y_max_global - IDB::y_min_global ) / IDB::intervals_np;

    auto clamp = []( std::size_t v, std::size_t lo, std::size_t hi ) -> std::size_t
    {
        return v < lo ? lo : ( v > hi ? hi : v );
    };

    std::size_t nx  = clamp(
        static_cast<std::size_t>( std::floor( ( x - IDB::x_min_global ) / dx ) ),
        0, static_cast<std::size_t>( IDB::intervals_np ) - 1 );
    std::size_t ny  = clamp(
        static_cast<std::size_t>( std::floor( ( y - IDB::y_min_global ) / dy ) ),
        0, static_cast<std::size_t>( IDB::intervals_np ) - 1 );
    std::size_t nt  = nx * static_cast<std::size_t>( IDB::intervals_np ) + ny;

    std::size_t sp  = ChebyshevTraits<IDB>::blocks_start[nt];
    std::size_t np  = ChebyshevTraits<IDB>::blocks_coeffs_np[nt];
    std::size_t ord = ChebyshevTraits<IDB>::blocks_max_cheby_order[nt];

    cusfloat xm = 2.0 * ( x - ChebyshevTraits<IDB>::x_min_region[nt] )
                            / ChebyshevTraits<IDB>::dx_region[nt] - 1.0;
    cusfloat ym = 2.0 * ( y - ChebyshevTraits<IDB>::y_min_region[nt] )
                            / ChebyshevTraits<IDB>::dy_region[nt] - 1.0;

    cusfloat poly_x[ ChebyshevTraits<IDB>::max_cheby_order + 1 ];
    cusfloat poly_y[ ChebyshevTraits<IDB>::max_cheby_order + 1 ];
    chebyshev_poly_upto_order( ord, xm, poly_x );
    chebyshev_poly_upto_order( ord, ym, poly_y );

    cusfloat result = 0.0;
    for ( std::size_t j = 0; j < np; j++ )
    {
        result += ChebyshevTraits<IDB>::coeffs[sp + j]
                * poly_x[ ChebyshevTraits<IDB>::ncx[sp + j] ]
                * poly_y[ ChebyshevTraits<IDB>::ncy[sp + j] ];
    }
    return result;
}


// =============================================================================
//  Helper: evaluate the 1-D folded polynomial stored in ChebyshevTraits<IDB>
//          (populated by fold_database_2d) at point x.
//          Requires IDB::intervals_np == 1.
// =============================================================================
template<typename IDB>
static cusfloat eval_folded_1d( cusfloat x )
{
    constexpr std::size_t nt = 0; // single interval

    std::size_t sp  = ChebyshevTraits<IDB>::blocks_start[nt];
    std::size_t np  = ChebyshevTraits<IDB>::blocks_coeffs_np[nt];
    std::size_t ord = ChebyshevTraits<IDB>::blocks_max_cheby_order[nt];

    cusfloat xm = 2.0 * ( x - ChebyshevTraits<IDB>::x_min_region[nt] )
                            / ChebyshevTraits<IDB>::dx_region[nt] - 1.0;

    cusfloat poly[ ChebyshevTraits<IDB>::max_cheby_order + 1 ];
    chebyshev_poly_upto_order( ord, xm, poly );

    cusfloat result = 0.0;
    for ( std::size_t j = 0; j < np; j++ )
    {
        result += ChebyshevTraits<IDB>::coeffs[sp + j]
                * poly[ ChebyshevTraits<IDB>::ncx[sp + j] ];
    }
    return result;
}


// =============================================================================
//  Helper: direct evaluation of the TestDB2D polynomial at (x, y)
// =============================================================================
static cusfloat eval_TestDB2D_direct( cusfloat x, cusfloat y )
{
    cusfloat xm = 2.0 * ( x - TestDB2D::x_min_region[0] ) / TestDB2D::dx_region[0] - 1.0;
    cusfloat ym = 2.0 * ( y - TestDB2D::y_min_region[0] ) / TestDB2D::dy_region[0] - 1.0;

    cusfloat poly_x[ TestDB2D::max_cheby_order + 1 ];
    cusfloat poly_y[ TestDB2D::max_cheby_order + 1 ];
    chebyshev_poly_upto_order( TestDB2D::max_cheby_order, xm, poly_x );
    chebyshev_poly_upto_order( TestDB2D::max_cheby_order, ym, poly_y );

    cusfloat result = 0.0;
    for ( std::size_t j = 0; j < TestDB2D::num_c; j++ )
    {
        result += TestDB2D::c[j] * poly_x[ TestDB2D::ncx[j] ] * poly_y[ TestDB2D::ncy[j] ];
    }
    return result;
}


// =============================================================================
//  Test 1 – fold_database_1d<L2C>
//
//  Strategy: for several water-depth values H, call fold_database_1d<L2C>(H)
//  and verify that the resulting scalar ChebyshevTraits<L2C>::coeffs equals
//  a direct evaluation of the L2C Chebyshev polynomial at log10(H).
// =============================================================================
static void test_fold_1d()
{
    constexpr cusfloat tol = 1000 * std::numeric_limits<cusfloat>::epsilon();

    // H values that map to each of the two L2C z-intervals:
    //   log10(0.1)   = -1.0  -> region [-2.5,  0.0)
    //   log10(0.01)  = -2.0  -> region [-2.5,  0.0)
    //   log10(0.001) = -3.0  -> region [-5.0, -2.5)
    //   log10(1e-4)  = -4.0  -> region [-5.0, -2.5)
    const cusfloat H_values[] = { 0.1, 0.01, 0.001, 1e-4 };
    constexpr int  n_test     = 4;

    for ( int t = 0; t < n_test; t++ )
    {
        cusfloat H = H_values[t];

        fold_database_1d<L2C>( H );

        cusfloat ref  = eval_L2C_direct( H );
        cusfloat diff = std::abs( ChebyshevTraits<L2C>::coeffs - ref );

        if ( diff > tol )
        {
            std::cerr << "FAIL test_fold_1d at H = " << H << std::endl;
            std::cerr << "  fold result = " << ChebyshevTraits<L2C>::coeffs << std::endl;
            std::cerr << "  reference   = " << ref  << std::endl;
            std::cerr << "  abs diff    = " << diff << "  (tol = " << tol << ")" << std::endl;
            throw std::runtime_error( "test_fold_1d failed" );
        }
    }

    std::cout << "PASS test_fold_1d" << std::endl;
}


// =============================================================================
//  Test 2 – fold_database_3d<L1C>
//
//  Strategy: for several (x, y, H) triples, call fold_database_3d<L1C>(H)
//  and evaluate the resulting 2-D folded polynomial at (x, y).  Compare the
//  value against a direct evaluation of the original 3-D L1C polynomial at
//  (x, y, H).  Both evaluations must agree to within floating-point roundoff.
// =============================================================================
static void test_fold_3d()
{
    constexpr cusfloat tol = 1000 * std::numeric_limits<cusfloat>::epsilon();

    // (x, y, H) test points – one per 2D block (2×2 grid) and two H values
    // covering both z-intervals of L1C (z_min=-5, z_max=0, log-scaled).
    const cusfloat xs[] = { 0.25, 0.75, 0.25, 0.75 };
    const cusfloat ys[] = { 0.25, 0.75, 0.75, 0.25 };
    const cusfloat Hs[] = { 0.1,  0.1,  0.001, 0.001 };
    constexpr int  n_test = 4;

    for ( int t = 0; t < n_test; t++ )
    {
        cusfloat x = xs[t];
        cusfloat y = ys[t];
        cusfloat H = Hs[t];

        fold_database_3d<L1C>( H );

        cusfloat folded = eval_folded_2d<L1C>( x, y );
        cusfloat ref    = eval_L1C_direct( x, y, H );
        cusfloat diff   = std::abs( folded - ref );

        if ( diff > tol * ( std::abs( ref ) + 1.0 ) )
        {
            std::cerr << "FAIL test_fold_3d at (x=" << x << ", y=" << y << ", H=" << H << ")" << std::endl;
            std::cerr << "  folded result = " << folded << std::endl;
            std::cerr << "  reference     = " << ref    << std::endl;
            std::cerr << "  abs diff      = " << diff   << "  (tol = " << tol << ")" << std::endl;
            throw std::runtime_error( "test_fold_3d failed" );
        }
    }

    std::cout << "PASS test_fold_3d" << std::endl;
}


// =============================================================================
//  Test 3 – fold_database_2d<TestDB2D>
//
//  Strategy: for two y-fold values and three x-evaluation points, call
//  fold_database_2d<TestDB2D>(y0) and evaluate the resulting 1-D folded
//  polynomial at x.  Compare against a direct evaluation of the full 2-D
//  TestDB2D polynomial at (x, y0).
// =============================================================================
static void test_fold_2d()
{
    constexpr cusfloat tol = 1000 * std::numeric_limits<cusfloat>::epsilon();

    const cusfloat y_fold_values[] = { 0.3, 0.7 };
    const cusfloat x_eval_values[] = { 0.2, 0.5, 0.8 };
    constexpr int  n_y = 2;
    constexpr int  n_x = 3;

    for ( int iy = 0; iy < n_y; iy++ )
    {
        cusfloat y0 = y_fold_values[iy];

        fold_database_2d<TestDB2D>( y0 );

        for ( int ix = 0; ix < n_x; ix++ )
        {
            cusfloat x  = x_eval_values[ix];

            cusfloat folded = eval_folded_1d<TestDB2D>( x );
            cusfloat ref    = eval_TestDB2D_direct( x, y0 );
            cusfloat diff   = std::abs( folded - ref );

            if ( diff > tol * ( std::abs( ref ) + 1.0 ) )
            {
                std::cerr << "FAIL test_fold_2d at (x=" << x << ", y0=" << y0 << ")" << std::endl;
                std::cerr << "  folded result = " << folded << std::endl;
                std::cerr << "  reference     = " << ref    << std::endl;
                std::cerr << "  abs diff      = " << diff   << "  (tol = " << tol << ")" << std::endl;
                throw std::runtime_error( "test_fold_2d failed" );
            }
        }
    }

    std::cout << "PASS test_fold_2d" << std::endl;
}


// =============================================================================
//  Entry point
// =============================================================================
int main()
{
    test_fold_1d();
    test_fold_3d();
    test_fold_2d();

    return 0;
}
