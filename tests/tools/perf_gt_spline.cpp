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

/**
 * Performance benchmark for 2D spline interpolation of Gt(beta, mu).
 *
 * Input (argv[1]):
 *   Path to aux_data/0_integrals_database/1_time_domain/Gt.h5
 *   Datasets used: beta (6000,), mu (100,), fcn (6000x100).
 *
 * Coordinate transform:
 *   The mu axis spans ~4 decades; the spline is built in log10(mu) space
 *   so the y resolution is uniform.  beta is uniformly spaced (step=0.01).
 *
 * Methods compared
 *   A. Bilinear              – no precomputation, O(log N) lookup
 *   B. Bicubic Catmull-Rom   – no precomputation, 4x4 kernel, C1 continuity
 *   C. Natural cubic spline  – precomputed second-derivative table (Thomas
 *                              algorithm, one per beta row in log_mu direction),
 *                              cubic in log_mu + Catmull-Rom in beta, C2 in mu
 *
 * Benchmarks reported
 *   1. Construction time (std::chrono, averaged over N_BUILD_REPEAT builds)
 *   2. Batch evaluation time for N = 1 … 100 000
 *   3. Single-point evaluation time (tight loop, N_SINGLE_CALLS iterations)
 *   4. On-grid accuracy   (max / mean |error| vs tabulated values)
 *   5. Off-grid accuracy  (50 000 random interior points, ref = bilinear)
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <functional>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "hdf5.h"


// ---------------------------------------------------------------------------
// Types and aliases
// ---------------------------------------------------------------------------

using real_t = double;
using Clock  = std::chrono::high_resolution_clock;
using Dur    = std::chrono::duration<double>;


// ---------------------------------------------------------------------------
// HDF5 data loading
// ---------------------------------------------------------------------------

struct GridData
{
    std::vector<real_t> beta;       // size N_BETA
    std::vector<real_t> log_mu;     // size N_MU,  = log10(mu)
    std::vector<real_t> fcn;        // size N_BETA * N_MU, row-major
    int n_beta = 0;
    int n_mu   = 0;
};

static std::vector<real_t> read_h5_dataset_1d( hid_t file, const char* name )
{
    hid_t dset  = H5Dopen( file, name, H5P_DEFAULT );
    hid_t space = H5Dget_space( dset );

    hsize_t dim = 0;
    H5Sget_simple_extent_dims( space, &dim, nullptr );

    std::vector<real_t> buf( static_cast<std::size_t>( dim ) );
    H5Dread( dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data() );

    H5Sclose( space );
    H5Dclose( dset );
    return buf;
}

static std::vector<real_t> read_h5_dataset_2d( hid_t file, const char* name,
                                                int& rows, int& cols )
{
    hid_t dset  = H5Dopen( file, name, H5P_DEFAULT );
    hid_t space = H5Dget_space( dset );

    hsize_t dims[2] = {};
    H5Sget_simple_extent_dims( space, dims, nullptr );
    rows = static_cast<int>( dims[0] );
    cols = static_cast<int>( dims[1] );

    std::vector<real_t> buf( static_cast<std::size_t>( rows * cols ) );
    H5Dread( dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data() );

    H5Sclose( space );
    H5Dclose( dset );
    return buf;
}

GridData load_data( const std::string& h5_path )
{
    hid_t file = H5Fopen( h5_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
    if ( file < 0 )
        throw std::runtime_error( "Cannot open HDF5 file: " + h5_path );

    GridData gd;

    // Read beta and mu axes
    gd.beta = read_h5_dataset_1d( file, "beta" );
    std::vector<real_t> mu = read_h5_dataset_1d( file, "mu" );

    // Transform mu -> log10(mu)
    gd.log_mu.resize( mu.size() );
    for ( std::size_t i = 0; i < mu.size(); ++i )
        gd.log_mu[i] = std::log10( mu[i] );

    // Read fcn (row-major: beta × mu)
    gd.fcn = read_h5_dataset_2d( file, "fcn", gd.n_beta, gd.n_mu );

    H5Fclose( file );

    if ( gd.n_beta != static_cast<int>( gd.beta.size() ) ||
         gd.n_mu   != static_cast<int>( gd.log_mu.size() ) )
        throw std::runtime_error( "HDF5 dimension mismatch" );

    return gd;
}


// ---------------------------------------------------------------------------
// Index lookup helpers
// ---------------------------------------------------------------------------

// Uniform grid: direct index from value (faster than binary search).
inline int find_cell_uniform( real_t x, real_t x0, real_t dx, int n_cells )
{
    int ix = static_cast<int>( ( x - x0 ) / dx );
    return std::max( 0, std::min( ix, n_cells - 1 ) );
}

// Non-uniform grid: binary search.
inline int find_cell( const real_t* grid, int n, real_t x )
{
    // Returns i such that grid[i] <= x < grid[i+1]  (clamped to [0, n-2])
    const real_t* it = std::upper_bound( grid, grid + n, x );
    int i = static_cast<int>( it - grid ) - 1;
    return std::max( 0, std::min( i, n - 2 ) );
}


// ---------------------------------------------------------------------------
// A. Bilinear interpolator
// ---------------------------------------------------------------------------

struct BilinearInterp
{
    const GridData& g;
    real_t beta0;
    real_t dbeta;   // uniform beta spacing

    explicit BilinearInterp( const GridData& gd )
        : g( gd )
        , beta0( gd.beta.front() )
        , dbeta( gd.beta[1] - gd.beta[0] )
    {}

    void build() const {}   // no precomputation

    real_t eval( real_t beta, real_t log_mu ) const
    {
        const int nb = g.n_beta;
        const int nm = g.n_mu;
        const int ix = find_cell_uniform( beta,   beta0, dbeta, nb - 1 );
        const int iy = find_cell( g.log_mu.data(), nm, log_mu );

        const real_t tx = ( beta   - g.beta[ix]   ) / ( g.beta[ix+1]   - g.beta[ix]   );
        const real_t ty = ( log_mu - g.log_mu[iy] ) / ( g.log_mu[iy+1] - g.log_mu[iy] );

        const real_t f00 = g.fcn[ ix    * nm + iy   ];
        const real_t f10 = g.fcn[(ix+1) * nm + iy   ];
        const real_t f01 = g.fcn[ ix    * nm + iy+1 ];
        const real_t f11 = g.fcn[(ix+1) * nm + iy+1 ];

        return ( 1.0-tx ) * ( ( 1.0-ty ) * f00 + ty * f01 )
             +      tx    * ( ( 1.0-ty ) * f10 + ty * f11 );
    }
};


// ---------------------------------------------------------------------------
// B. Bicubic Catmull-Rom interpolator
//    No precomputation; uses the 4x4 neighbourhood at each query point.
//    Boundary: clamp index to valid range (repeat edge value).
// ---------------------------------------------------------------------------

namespace detail {

// Standard Catmull-Rom basis weights for parameter t in [0, 1].
// Corresponds to the Keys cubic convolution kernel with a = -0.5.
inline void catmull_rom_weights( real_t t, real_t& w0, real_t& w1,
                                  real_t& w2, real_t& w3 )
{
    const real_t t2 = t  * t;
    const real_t t3 = t2 * t;
    w0 = -0.5*t3 + 1.0*t2 - 0.5*t;
    w1 =  1.5*t3 - 2.5*t2 + 1.0;
    w2 = -1.5*t3 + 2.0*t2 + 0.5*t;
    w3 =  0.5*t3 - 0.5*t2;
}

inline int clamp( int v, int lo, int hi ) { return v < lo ? lo : ( v > hi ? hi : v ); }

} // namespace detail


struct BicubicCatmullRomInterp
{
    const GridData& g;
    real_t beta0;
    real_t dbeta;

    explicit BicubicCatmullRomInterp( const GridData& gd )
        : g( gd )
        , beta0( gd.beta.front() )
        , dbeta( gd.beta[1] - gd.beta[0] )
    {}

    void build() const {}   // no precomputation

    real_t eval( real_t beta, real_t log_mu ) const
    {
        const int nb = g.n_beta;
        const int nm = g.n_mu;
        const int ix = find_cell_uniform( beta,   beta0, dbeta, nb - 1 );
        const int iy = find_cell( g.log_mu.data(), nm, log_mu );

        const real_t tx = ( beta   - g.beta[ix]   ) / dbeta;
        const real_t ty = ( log_mu - g.log_mu[iy] ) / ( g.log_mu[iy+1] - g.log_mu[iy] );

        real_t wx[4], wy[4];
        detail::catmull_rom_weights( tx, wx[0], wx[1], wx[2], wx[3] );
        detail::catmull_rom_weights( ty, wy[0], wy[1], wy[2], wy[3] );

        // Gather 4x4 neighbourhood (with boundary clamping)
        real_t val = 0.0;
        for ( int di = 0; di < 4; ++di )
        {
            const int ri = detail::clamp( ix + di - 1, 0, nb - 1 );
            for ( int dj = 0; dj < 4; ++dj )
            {
                const int rj = detail::clamp( iy + dj - 1, 0, nm - 1 );
                val += wx[di] * wy[dj] * g.fcn[ ri * nm + rj ];
            }
        }
        return val;
    }
};


// ---------------------------------------------------------------------------
// C. Natural cubic spline (log_mu direction, Catmull-Rom in beta)
//
//    Precompute: for each beta row i, solve the natural cubic spline
//    tridiagonal system in log_mu direction using the Thomas algorithm.
//    The solution gives second derivatives m[i][j] at each knot.
//
//    Evaluate: cubic spline formula in log_mu for 4 bracketing beta rows,
//    then Catmull-Rom in beta across those 4 values.
// ---------------------------------------------------------------------------

// Thomas algorithm for a symmetric positive-definite tridiagonal system
//   diag_main[k] * x[k] + off_diag[k-1] * x[k-1] + off_diag[k] * x[k+1] = rhs[k]
// in-place solution stored in x[].
// n: size of system, off_diag has size n-1.
static void thomas_solve( const real_t* diag_main, const real_t* off_diag,
                           const real_t* rhs, real_t* x, int n )
{
    // Forward sweep (LU decomposition into vectors)
    std::vector<real_t> c_prime( n - 1 );
    std::vector<real_t> d_prime( n );

    c_prime[0] = off_diag[0] / diag_main[0];
    d_prime[0] = rhs[0]       / diag_main[0];

    for ( int k = 1; k < n; ++k )
    {
        const real_t w = diag_main[k] - off_diag[k-1] * c_prime[k-1];
        if ( k < n - 1 )
            c_prime[k] = off_diag[k] / w;
        d_prime[k] = ( rhs[k] - off_diag[k-1] * d_prime[k-1] ) / w;
    }

    // Back substitution
    x[n-1] = d_prime[n-1];
    for ( int k = n - 2; k >= 0; --k )
        x[k] = d_prime[k] - c_prime[k] * x[k+1];
}


struct NaturalCubicMuInterp
{
    const GridData& g;
    real_t beta0;
    real_t dbeta;

    // m2[i * n_mu + j] = second derivative of the log_mu spline for row i at knot j.
    std::vector<real_t> m2;

    explicit NaturalCubicMuInterp( const GridData& gd )
        : g( gd )
        , beta0( gd.beta.front() )
        , dbeta( gd.beta[1] - gd.beta[0] )
    {}

    void build()
    {
        const int nb = g.n_beta;
        const int nm = g.n_mu;
        const real_t* lmu = g.log_mu.data();

        // Pre-build the tridiagonal structure (same for every beta row)
        // Interior knots: j = 1 .. nm-2  (system size n = nm-2)
        const int n = nm - 2;

        std::vector<real_t> diag( n ), off( n - 1 ), rhs( n ), sol( n );

        // h[j] = lmu[j+1] - lmu[j]
        std::vector<real_t> h( nm - 1 );
        for ( int j = 0; j < nm - 1; ++j )
            h[j] = lmu[j+1] - lmu[j];

        // Build tridiagonal (same structure for all beta rows)
        for ( int k = 0; k < n; ++k )
        {
            diag[k] = 2.0 * ( h[k] + h[k+1] );
            if ( k < n - 1 ) off[k] = h[k+1];
        }

        // Allocate second-derivative table
        m2.assign( static_cast<std::size_t>( nb * nm ), 0.0 );

        // Solve for each beta row
        for ( int i = 0; i < nb; ++i )
        {
            const real_t* fi = g.fcn.data() + i * nm;

            // Build right-hand side: 6 * ((f[k+2]-f[k+1])/h[k+1] - (f[k+1]-f[k])/h[k])
            for ( int k = 0; k < n; ++k )
            {
                const int j = k + 1;   // interior knot index
                rhs[k] = 6.0 * ( ( fi[j+1] - fi[j] ) / h[j]
                                - ( fi[j]   - fi[j-1] ) / h[j-1] );
            }

            // Solve tridiagonal (natural BC: m[0] = m[nm-1] = 0)
            thomas_solve( diag.data(), off.data(), rhs.data(), sol.data(), n );

            // Store into m2 table; boundary knots remain 0
            real_t* mi = m2.data() + i * nm;
            mi[0]    = 0.0;
            mi[nm-1] = 0.0;
            for ( int k = 0; k < n; ++k )
                mi[k+1] = sol[k];
        }
    }

    // Evaluate the natural cubic spline in log_mu for beta row ri at point lmu_query.
    real_t eval_mu_spline( int ri, real_t lmu_query ) const
    {
        const int nm    = g.n_mu;
        const int iy    = find_cell( g.log_mu.data(), nm, lmu_query );
        const real_t h  = g.log_mu[iy+1] - g.log_mu[iy];
        const real_t s  = lmu_query - g.log_mu[iy];

        const real_t* fi = g.fcn.data() + ri * nm;
        const real_t* mi = m2.data()    + ri * nm;

        // Standard cubic spline formula:
        //   f(x) = a + b*s + c*s^2 + d*s^3
        // where:
        //   a = fi[iy]
        //   b = (fi[iy+1] - fi[iy]) / h  -  h * (2*mi[iy] + mi[iy+1]) / 6
        //   c = mi[iy] / 2
        //   d = (mi[iy+1] - mi[iy]) / (6*h)
        const real_t a = fi[iy];
        const real_t b = ( fi[iy+1] - fi[iy] ) / h  -  h * ( 2.0*mi[iy] + mi[iy+1] ) / 6.0;
        const real_t c = mi[iy] / 2.0;
        const real_t d = ( mi[iy+1] - mi[iy] ) / ( 6.0 * h );

        return a + s * ( b + s * ( c + s * d ) );
    }

    real_t eval( real_t beta, real_t log_mu ) const
    {
        const int nb = g.n_beta;
        const int ix = find_cell_uniform( beta, beta0, dbeta, nb - 1 );

        // t in [0, 1] within the beta cell
        const real_t t = ( beta - g.beta[ix] ) / dbeta;

        // Evaluate mu-spline for 4 bracketing beta rows (Catmull-Rom in beta)
        const int r0 = detail::clamp( ix - 1, 0, nb - 1 );
        const int r1 = ix;
        const int r2 = detail::clamp( ix + 1, 0, nb - 1 );
        const int r3 = detail::clamp( ix + 2, 0, nb - 1 );

        const real_t v0 = eval_mu_spline( r0, log_mu );
        const real_t v1 = eval_mu_spline( r1, log_mu );
        const real_t v2 = eval_mu_spline( r2, log_mu );
        const real_t v3 = eval_mu_spline( r3, log_mu );

        // Catmull-Rom in beta using the 4 values
        real_t wx[4];
        detail::catmull_rom_weights( t, wx[0], wx[1], wx[2], wx[3] );

        return wx[0]*v0 + wx[1]*v1 + wx[2]*v2 + wx[3]*v3;
    }
};


// ---------------------------------------------------------------------------
// Printing helpers
// ---------------------------------------------------------------------------

static std::string fmt_time( double s )
{
    char buf[64];
    if      ( s >= 1.0   ) std::snprintf( buf, sizeof(buf), "%.3f s",  s );
    else if ( s >= 1e-3  ) std::snprintf( buf, sizeof(buf), "%.3f ms", s * 1e3 );
    else if ( s >= 1e-6  ) std::snprintf( buf, sizeof(buf), "%.3f us", s * 1e6 );
    else                   std::snprintf( buf, sizeof(buf), "%.1f ns", s * 1e9 );
    return buf;
}

static std::string fmt_e( double v )
{
    char buf[32];
    std::snprintf( buf, sizeof(buf), "%.3e", v );
    return buf;
}

static void section( const std::string& title )
{
    const std::string bar( 72, '=' );
    std::cout << "\n" << bar << "\n  " << title << "\n" << bar << "\n";
}

static void row( const std::string& label, const std::vector<std::string>& vals,
                 int lw = 38 )
{
    std::cout << "  " << std::left << std::setw(lw) << label;
    for ( const auto& v : vals )
        std::cout << "  " << std::right << std::setw(15) << v;
    std::cout << "\n";
}


// ---------------------------------------------------------------------------
// Benchmark helpers
// ---------------------------------------------------------------------------

static const std::vector<int> BATCH_SIZES = { 1, 10, 100, 1'000, 10'000, 100'000 };
static constexpr int N_SINGLE_CALLS = 100'000;
static constexpr int N_BUILD_REPEAT = 5;
static constexpr int N_EVAL_REPEAT  = 5;
static constexpr int N_OFF_GRID_PTS = 50'000;

// Generate N random (beta, log_mu) points inside the domain
static std::vector<real_t> random_pts( const GridData& g, int n, uint64_t seed = 42 )
{
    std::mt19937_64 rng( seed );
    std::uniform_real_distribution<real_t> db( g.beta.front(),   g.beta.back()   );
    std::uniform_real_distribution<real_t> dm( g.log_mu.front(), g.log_mu.back() );

    std::vector<real_t> pts( 2 * static_cast<std::size_t>(n) );
    for ( int i = 0; i < n; ++i )
    {
        pts[2*i  ] = db( rng );
        pts[2*i+1] = dm( rng );
    }
    return pts;
}

// Batch evaluation using a uniform signature: eval_fn receives all points
// and writes results to out[].
template<typename Interp>
static void eval_batch( const Interp& interp,
                        const real_t* pts, int n, real_t* out )
{
    for ( int i = 0; i < n; ++i )
        out[i] = interp.eval( pts[2*i], pts[2*i+1] );
}

template<typename Interp>
static double time_batch( const Interp& interp,
                          const real_t* pts_pool, int pool_size,
                          int batch_size )
{
    // Pick batch_size points from the pool (wrap if needed)
    std::vector<real_t> batch( 2 * static_cast<std::size_t>(batch_size) );
    for ( int i = 0; i < batch_size; ++i )
    {
        const int idx = i % pool_size;
        batch[2*i  ] = pts_pool[2*idx  ];
        batch[2*i+1] = pts_pool[2*idx+1];
    }
    std::vector<real_t> out( static_cast<std::size_t>(batch_size) );

    double total = 0.0;
    for ( int r = 0; r < N_EVAL_REPEAT; ++r )
    {
        auto t0 = Clock::now();
        eval_batch( interp, batch.data(), batch_size, out.data() );
        total += Dur( Clock::now() - t0 ).count();
    }
    return total / N_EVAL_REPEAT;
}

template<typename Interp>
static double time_single( const Interp& interp, real_t beta, real_t log_mu )
{
    volatile real_t sink = 0.0;   // prevent dead-code elimination
    auto t0 = Clock::now();
    for ( int i = 0; i < N_SINGLE_CALLS; ++i )
        sink += interp.eval( beta, log_mu );
    double total = Dur( Clock::now() - t0 ).count();
    ( void ) sink;
    return total / N_SINGLE_CALLS;
}


// ---------------------------------------------------------------------------
// Accuracy helpers
// ---------------------------------------------------------------------------

struct ErrStats { double max_abs, mean_abs; };

template<typename Interp>
static ErrStats on_grid_accuracy( const Interp& interp, const GridData& g )
{
    double max_abs = 0.0, sum_abs = 0.0;
    const std::size_t total = static_cast<std::size_t>( g.n_beta * g.n_mu );
    for ( int i = 0; i < g.n_beta; ++i )
    {
        for ( int j = 0; j < g.n_mu; ++j )
        {
            const real_t pred = interp.eval( g.beta[i], g.log_mu[j] );
            const real_t ref  = g.fcn[ i * g.n_mu + j ];
            const real_t diff = std::abs( pred - ref );
            max_abs  = std::max( max_abs, diff );
            sum_abs += diff;
        }
    }
    return { max_abs, sum_abs / static_cast<double>( total ) };
}

template<typename Interp>
static ErrStats off_grid_accuracy( const Interp& interp,
                                   const BilinearInterp& ref,
                                   const GridData& g )
{
    const auto pts = random_pts( g, N_OFF_GRID_PTS, 7 );
    double max_abs = 0.0, sum_abs = 0.0;
    for ( int i = 0; i < N_OFF_GRID_PTS; ++i )
    {
        const real_t b  = pts[2*i  ];
        const real_t lm = pts[2*i+1];
        const real_t pred = interp.eval( b, lm );
        const real_t rval = ref.eval(    b, lm );
        const real_t diff = std::abs( pred - rval );
        max_abs  = std::max( max_abs, diff );
        sum_abs += diff;
    }
    return { max_abs, sum_abs / static_cast<double>( N_OFF_GRID_PTS ) };
}


// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
    if ( argc < 2 )
    {
        std::cerr << "Usage: " << argv[0] << " <path/to/Gt.h5>\n";
        return 1;
    }

    // ---- Load data ---------------------------------------------------------
    section( "Loading data" );
    GridData g;
    {
        auto t0 = Clock::now();
        g = load_data( argv[1] );
        double dt = Dur( Clock::now() - t0 ).count();

        std::cout << "  beta    : " << g.n_beta << " pts   ["
                  << g.beta.front() << ", " << g.beta.back() << "]   step="
                  << ( g.beta[1] - g.beta[0] ) << "\n";
        std::cout << "  log_mu  : " << g.n_mu << " pts   ["
                  << g.log_mu.front() << ", " << g.log_mu.back() << "]\n";
        std::cout << "  fcn     : " << g.n_beta << " x " << g.n_mu << "\n";
        std::cout << "  Load time : " << fmt_time( dt ) << "\n";
    }

    // ---- Build interpolators -----------------------------------------------
    BilinearInterp          bilin( g );
    BicubicCatmullRomInterp catrom( g );
    NaturalCubicMuInterp    natcub( g );

    struct Method {
        std::string name;
        std::function<void()>            build_fn;
        std::function<real_t(real_t,real_t)> eval_fn;
    };

    std::vector<Method> methods = {
        { "Bilinear",
          [&]{ bilin.build(); },
          [&]( real_t b, real_t lm ){ return bilin.eval( b, lm ); } },

        { "Bicubic Catmull-Rom (4x4)",
          [&]{ catrom.build(); },
          [&]( real_t b, real_t lm ){ return catrom.eval( b, lm ); } },

        { "Natural cubic spline (log_mu)",
          [&]{ natcub.build(); },
          [&]( real_t b, real_t lm ){ return natcub.eval( b, lm ); } },
    };

    // =========================================================================
    // 1.  Construction cost
    // =========================================================================
    section( "Construction cost  (average over " + std::to_string(N_BUILD_REPEAT) + " builds)" );
    row( "Method", { "Build time (avg)" } );
    std::cout << "  " << std::string( 56, '-' ) << "\n";

    for ( auto& m : methods )
    {
        double total = 0.0;
        for ( int r = 0; r < N_BUILD_REPEAT; ++r )
        {
            auto t0 = Clock::now();
            m.build_fn();
            total += Dur( Clock::now() - t0 ).count();
        }
        m.build_fn();   // ensure the object is in a valid built state
        row( m.name, { fmt_time( total / N_BUILD_REPEAT ) } );
    }

    // =========================================================================
    // 2.  Batch evaluation time
    // =========================================================================
    section( "Batch evaluation time  (mean of " + std::to_string(N_EVAL_REPEAT) + " runs)" );
    {
        const int pool_size = std::min( 100'000, g.n_beta * g.n_mu );
        const auto pts_pool = random_pts( g, pool_size, 42 );

        std::vector<std::string> header = { "Method" };
        for ( int bs : BATCH_SIZES )
            header.push_back( "N=" + std::to_string(bs) );
        row( header[0], std::vector<std::string>( header.begin()+1, header.end() ) );
        std::cout << "  " << std::string( 72, '-' ) << "\n";

        for ( auto& m : methods )
        {
            std::vector<std::string> vals;
            for ( int bs : BATCH_SIZES )
            {
                // Wrap pool into a batch of exactly bs points
                std::vector<real_t> batch( 2 * static_cast<std::size_t>(bs) );
                for ( int i = 0; i < bs; ++i )
                {
                    const int idx = i % pool_size;
                    batch[2*i  ] = pts_pool[2*idx  ];
                    batch[2*i+1] = pts_pool[2*idx+1];
                }
                std::vector<real_t> out( static_cast<std::size_t>(bs) );

                double ttotal = 0.0;
                for ( int r = 0; r < N_EVAL_REPEAT; ++r )
                {
                    auto t0 = Clock::now();
                    for ( int i = 0; i < bs; ++i )
                        out[i] = m.eval_fn( batch[2*i], batch[2*i+1] );
                    ttotal += Dur( Clock::now() - t0 ).count();
                }
                vals.push_back( fmt_time( ttotal / N_EVAL_REPEAT ) );
            }
            row( m.name, vals );
        }
    }

    // =========================================================================
    // 3.  Throughput at N = 100 000
    // =========================================================================
    section( "Throughput at N = 100 000" );
    row( "Method", { "Total time", "Evals/sec" } );
    std::cout << "  " << std::string( 72, '-' ) << "\n";
    {
        const int bs = 100'000;
        const auto pts = random_pts( g, bs, 43 );
        std::vector<real_t> out( static_cast<std::size_t>(bs) );

        for ( auto& m : methods )
        {
            double ttotal = 0.0;
            for ( int r = 0; r < N_EVAL_REPEAT; ++r )
            {
                auto t0 = Clock::now();
                for ( int i = 0; i < bs; ++i )
                    out[i] = m.eval_fn( pts[2*i], pts[2*i+1] );
                ttotal += Dur( Clock::now() - t0 ).count();
            }
            const double t_mean = ttotal / N_EVAL_REPEAT;
            char thr[32];
            std::snprintf( thr, sizeof(thr), "%.0f", static_cast<double>(bs) / t_mean );
            row( m.name, { fmt_time( t_mean ), thr } );
        }
    }

    // =========================================================================
    // 4.  Single-point evaluation time
    // =========================================================================
    section( "Single-point evaluation time  (" +
             std::to_string(N_SINGLE_CALLS) + " calls, tight loop)" );
    row( "Method", { "Mean per-call", "Calls/sec" } );
    std::cout << "  " << std::string( 72, '-' ) << "\n";
    {
        const real_t b_q  = 0.5 * ( g.beta.front()   + g.beta.back()   );
        const real_t lm_q = 0.5 * ( g.log_mu.front() + g.log_mu.back() );
        for ( auto& m : methods )
        {
            volatile real_t sink = 0.0;
            auto t0 = Clock::now();
            for ( int i = 0; i < N_SINGLE_CALLS; ++i )
                sink += m.eval_fn( b_q, lm_q );
            const double total = Dur( Clock::now() - t0 ).count();
            ( void ) sink;
            const double per_call = total / N_SINGLE_CALLS;
            char thr[32];
            std::snprintf( thr, sizeof(thr), "%.0f", 1.0 / per_call );
            row( m.name, { fmt_time( per_call ), thr } );
        }
    }

    // =========================================================================
    // 5.  On-grid accuracy
    // =========================================================================
    section( "On-grid accuracy  (" +
             std::to_string(g.n_beta) + " x " + std::to_string(g.n_mu) +
             " = " + std::to_string(g.n_beta * g.n_mu) + " nodes)" );
    row( "Method", { "max|abs|", "mean|abs|" } );
    std::cout << "  " << std::string( 72, '-' ) << "\n";
    {
        for ( auto& m : methods )
        {
            double max_abs = 0.0, sum_abs = 0.0;
            const std::size_t total = static_cast<std::size_t>( g.n_beta * g.n_mu );
            for ( int i = 0; i < g.n_beta; ++i )
                for ( int j = 0; j < g.n_mu; ++j )
                {
                    const real_t pred = m.eval_fn( g.beta[i], g.log_mu[j] );
                    const real_t ref  = g.fcn[ i * g.n_mu + j ];
                    const real_t diff = std::abs( pred - ref );
                    max_abs  = std::max( max_abs, diff );
                    sum_abs += diff;
                }
            row( m.name, { fmt_e( max_abs ), fmt_e( sum_abs / static_cast<double>(total) ) } );
        }
    }

    // =========================================================================
    // 6.  Off-grid accuracy  (vs bilinear reference)
    // =========================================================================
    section( "Off-grid accuracy  (" + std::to_string(N_OFF_GRID_PTS) +
             " random interior points, ref = Bilinear)" );
    std::cout << "  Note: differences reflect higher-order vs bilinear, not true error.\n";
    row( "Method", { "max|abs|", "mean|abs|" } );
    std::cout << "  " << std::string( 72, '-' ) << "\n";
    {
        const auto off_pts = random_pts( g, N_OFF_GRID_PTS, 7 );
        for ( auto& m : methods )
        {
            if ( m.name == "Bilinear" )
            {
                row( m.name, { "(reference)", "---" } );
                continue;
            }
            double max_abs = 0.0, sum_abs = 0.0;
            for ( int i = 0; i < N_OFF_GRID_PTS; ++i )
            {
                const real_t b  = off_pts[2*i  ];
                const real_t lm = off_pts[2*i+1];
                const real_t pred  = m.eval_fn(    b, lm );
                const real_t rval  = bilin.eval( b, lm );
                const real_t diff  = std::abs( pred - rval );
                max_abs  = std::max( max_abs, diff );
                sum_abs += diff;
            }
            row( m.name, { fmt_e( max_abs ),
                           fmt_e( sum_abs / static_cast<double>(N_OFF_GRID_PTS) ) } );
        }
    }

    std::cout << "\nDone.\n\n";
    return 0;
}
