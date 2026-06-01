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

/**
 * Header-only 2D grid interpolation for the time-domain integral databases.
 *
 * All three classes share the same interface:
 *
 *   Constructor
 *     InterpolatorType(n_beta, beta[], n_mu, mu[], fcn[])
 *
 *     beta[]  – n_beta uniformly-spaced tabulation points (axis 0)
 *     mu[]    – n_mu  tabulation points in LINEAR scale, e.g. [0.0001, 0.9999].
 *               Internally stored as log10(mu) to achieve uniform coverage
 *               across the ~4 decades of mu.
 *     fcn[]   – function values, row-major: fcn[i*n_mu + j] = f(beta[i], mu[j])
 *
 *   Evaluation
 *     T eval(T beta, T log_mu) const
 *
 *     Callers supply log_mu = log10(mu) directly.  This avoids a log10()
 *     call at every evaluation site when the caller already works in log space
 *     (which matches the existing eval_time_residual_2d interface).
 *
 * Methods
 *   Bilinear2D<T>             – bilinear interpolation, O(log N) cell lookup,
 *                               no precomputation; exact at grid nodes.
 *
 *   BicubicCatmullRom2D<T>    – bicubic Catmull-Rom convolution (4×4 kernel),
 *                               no precomputation, C1 globally; near-exact at nodes.
 *
 *   NaturalCubicSpline2D<T>   – precomputed natural cubic spline in the log_mu
 *                               direction (Thomas algorithm, one tridiagonal
 *                               solve per beta row) combined with Catmull-Rom
 *                               across beta rows; C2 in mu, C1 in beta.
 *
 * Template parameter
 *   T – floating-point type (float or double; use cusfloat for project code).
 *
 * Dependencies
 *   Standard library only: <algorithm>, <cmath>, <stdexcept>, <vector>.
 *   No MKL, no HDF5, no project headers.
 */

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>


// ---------------------------------------------------------------------------
// Internal helpers (not part of public API)
// ---------------------------------------------------------------------------

namespace interp2d_detail {

// Locate the interval i such that axis[i] <= x < axis[i+1], clamped to
// [0, n-2].  Used for the non-uniform log_mu axis.
template<typename T>
inline int find_cell_nonuniform( const T* axis, int n, T x ) noexcept
{
    const T* it = std::upper_bound( axis, axis + n, x );
    int i = static_cast<int>( it - axis ) - 1;
    if ( i < 0 )     i = 0;
    if ( i > n - 2 ) i = n - 2;
    return i;
}

// Locate the cell for a uniform axis given the first point and step.
// Returns index clamped to [0, n_cells-1].
template<typename T>
inline int find_cell_uniform( T x, T x0, T dx, int n_cells ) noexcept
{
    int ix = static_cast<int>( ( x - x0 ) / dx );
    if ( ix < 0 )        ix = 0;
    if ( ix > n_cells-1) ix = n_cells - 1;
    return ix;
}

// Catmull-Rom basis weights for parameter t in [0, 1].
// Returns (w0, w1, w2, w3) for control points at indices -1, 0, 1, 2.
template<typename T>
inline void catmull_rom_weights( T t, T& w0, T& w1, T& w2, T& w3 ) noexcept
{
    const T t2 = t  * t;
    const T t3 = t2 * t;
    w0 = static_cast<T>(-0.5)*t3 + static_cast<T>( 1.0)*t2 + static_cast<T>(-0.5)*t;
    w1 = static_cast<T>( 1.5)*t3 + static_cast<T>(-2.5)*t2 + static_cast<T>( 1.0);
    w2 = static_cast<T>(-1.5)*t3 + static_cast<T>( 2.0)*t2 + static_cast<T>( 0.5)*t;
    w3 = static_cast<T>( 0.5)*t3 + static_cast<T>(-0.5)*t2;
}

// Integer clamping helper.
inline int clamp_idx( int v, int lo, int hi ) noexcept
{
    return v < lo ? lo : ( v > hi ? hi : v );
}

// Thomas algorithm: symmetric tridiagonal solve Ax = rhs.
// diag_main[k], off[k] for k=0..n-1 and k=0..n-2 respectively.
// Result written to x[].
template<typename T>
inline void thomas_solve( const T* diag_main, const T* off,
                          const T* rhs, T* x, int n )
{
    std::vector<T> c( static_cast<std::size_t>(n - 1) );
    std::vector<T> d( static_cast<std::size_t>(n)     );

    c[0] = off[0]  / diag_main[0];
    d[0] = rhs[0]  / diag_main[0];

    for ( int k = 1; k < n; ++k )
    {
        const T w = diag_main[k] - off[k-1] * c[k-1];
        if ( k < n - 1 ) c[k] = off[k] / w;
        d[k] = ( rhs[k] - off[k-1] * d[k-1] ) / w;
    }

    x[n-1] = d[n-1];
    for ( int k = n - 2; k >= 0; --k )
        x[k] = d[k] - c[k] * x[k+1];
}

} // namespace interp2d_detail


// ---------------------------------------------------------------------------
// A. Bilinear2D
//
// Piecewise bilinear interpolation on a 2D tensor-product grid.
// Beta axis is uniform; log_mu axis is non-uniform (but monotone).
// Exact at every grid node.  C0 globally (derivative discontinuous at
// cell boundaries).
// ---------------------------------------------------------------------------

template<typename T>
struct Bilinear2D
{
    int  _n_beta;
    int  _n_mu;
    T    _beta0;    // beta[0]
    T    _dbeta;    // uniform beta step

    std::vector<T> _log_mu;  // stored log10(mu) axis
    std::vector<T> _fcn;     // row-major: _fcn[i*n_mu + j] = f(beta[i], mu[j])

    /**
     * @param n_beta  Number of beta tabulation points.
     * @param beta    Beta values (uniform spacing assumed).
     * @param n_mu    Number of mu tabulation points.
     * @param mu      mu values in LINEAR scale; internally converted to log10.
     * @param fcn     Function values, row-major (beta × mu).
     */
    Bilinear2D( int n_beta, const T* beta,
                int n_mu,   const T* mu,
                const T* fcn )
        : _n_beta( n_beta )
        , _n_mu( n_mu )
        , _beta0( beta[0] )
        , _dbeta( ( n_beta > 1 ) ? ( beta[1] - beta[0] ) : static_cast<T>(1) )
        , _log_mu( static_cast<std::size_t>(n_mu)   )
        , _fcn   ( static_cast<std::size_t>(n_beta * n_mu) )
    {
        if ( n_beta < 2 || n_mu < 2 )
            throw std::invalid_argument( "Bilinear2D: grid must have at least 2 points per axis" );

        for ( int j = 0; j < n_mu; ++j )
            _log_mu[j] = std::log10( mu[j] );

        std::copy( fcn, fcn + n_beta * n_mu, _fcn.begin() );
    }

    /**
     * Evaluate the interpolant at (beta, log_mu).
     * @param beta    Query beta value; clamped to [beta[0], beta[n_beta-1]].
     * @param log_mu  Query log10(mu); clamped to [log_mu[0], log_mu[n_mu-1]].
     */
    T eval( T beta, T log_mu ) const noexcept
    {
        const int nb = _n_beta;
        const int nm = _n_mu;
        const int ix = interp2d_detail::find_cell_uniform( beta,   _beta0,   _dbeta, nb - 1 );
        const int iy = interp2d_detail::find_cell_nonuniform( _log_mu.data(), nm, log_mu );

        const T beta_lo = _beta0 + static_cast<T>(ix) * _dbeta;
        const T tx2     = ( beta   - beta_lo ) / _dbeta;
        const T ty      = ( log_mu - _log_mu[iy] ) / ( _log_mu[iy+1] - _log_mu[iy] );

        const T f00 = _fcn[  ix    * nm + iy   ];
        const T f10 = _fcn[ (ix+1) * nm + iy   ];
        const T f01 = _fcn[  ix    * nm + iy+1 ];
        const T f11 = _fcn[ (ix+1) * nm + iy+1 ];

        const T one = static_cast<T>(1);
        return ( one - tx2 ) * ( ( one - ty ) * f00 + ty * f01 )
             +        tx2    * ( ( one - ty ) * f10 + ty * f11 );
    }
};


// ---------------------------------------------------------------------------
// B. BicubicCatmullRom2D
//
// Bicubic Catmull-Rom convolution over a 4×4 neighbourhood.
// No precomputation.  C1 globally.  Recovers the tabulated value at grid
// nodes to within floating-point round-off.
// Boundary treatment: clamped (edge value repeated).
// ---------------------------------------------------------------------------

template<typename T>
struct BicubicCatmullRom2D
{
    int  _n_beta;
    int  _n_mu;
    T    _beta0;
    T    _dbeta;

    std::vector<T> _log_mu;
    std::vector<T> _fcn;

    BicubicCatmullRom2D( int n_beta, const T* beta,
                         int n_mu,   const T* mu,
                         const T* fcn )
        : _n_beta( n_beta )
        , _n_mu( n_mu )
        , _beta0( beta[0] )
        , _dbeta( ( n_beta > 1 ) ? ( beta[1] - beta[0] ) : static_cast<T>(1) )
        , _log_mu( static_cast<std::size_t>(n_mu)   )
        , _fcn   ( static_cast<std::size_t>(n_beta * n_mu) )
    {
        if ( n_beta < 4 || n_mu < 4 )
            throw std::invalid_argument( "BicubicCatmullRom2D: grid must have at least 4 points per axis" );

        for ( int j = 0; j < n_mu; ++j )
            _log_mu[j] = std::log10( mu[j] );

        std::copy( fcn, fcn + n_beta * n_mu, _fcn.begin() );
    }

    T eval( T beta, T log_mu ) const noexcept
    {
        const int nb = _n_beta;
        const int nm = _n_mu;
        const int ix = interp2d_detail::find_cell_uniform( beta,   _beta0, _dbeta, nb - 1 );
        const int iy = interp2d_detail::find_cell_nonuniform( _log_mu.data(), nm, log_mu );

        const T beta_lo = _beta0 + static_cast<T>(ix) * _dbeta;
        const T tx = ( beta   - beta_lo ) / _dbeta;
        const T ty = ( log_mu - _log_mu[iy] ) / ( _log_mu[iy+1] - _log_mu[iy] );

        T wx[4], wy[4];
        interp2d_detail::catmull_rom_weights( tx, wx[0], wx[1], wx[2], wx[3] );
        interp2d_detail::catmull_rom_weights( ty, wy[0], wy[1], wy[2], wy[3] );

        T val = static_cast<T>(0);
        for ( int di = 0; di < 4; ++di )
        {
            const int ri = interp2d_detail::clamp_idx( ix + di - 1, 0, nb - 1 );
            T col = static_cast<T>(0);
            for ( int dj = 0; dj < 4; ++dj )
            {
                const int rj = interp2d_detail::clamp_idx( iy + dj - 1, 0, nm - 1 );
                col += wy[dj] * _fcn[ ri * nm + rj ];
            }
            val += wx[di] * col;
        }
        return val;
    }
};


// ---------------------------------------------------------------------------
// C. NaturalCubicSpline2D
//
// Natural cubic spline in the log_mu direction (C2) combined with
// Catmull-Rom interpolation across beta rows (C1 in beta).
//
// Precomputation (constructor):
//   For each beta row i, solve the natural-spline tridiagonal system for the
//   second derivatives m[i][j] at each log_mu knot using the Thomas algorithm.
//   Storage: m2[i * n_mu + j].
//
// Evaluation:
//   1. Identify the 4 nearest beta rows (ix-1, ix, ix+1, ix+2).
//   2. Evaluate the natural cubic spline in log_mu for each of those rows.
//   3. Combine the 4 values with Catmull-Rom weights in beta.
// ---------------------------------------------------------------------------

template<typename T>
struct NaturalCubicSpline2D
{
    int  _n_beta;
    int  _n_mu;
    T    _beta0;
    T    _dbeta;

    std::vector<T> _log_mu;
    std::vector<T> _fcn;
    std::vector<T> _m2;    // second derivatives: _m2[i * n_mu + j]

    NaturalCubicSpline2D( int n_beta, const T* beta,
                          int n_mu,   const T* mu,
                          const T* fcn )
        : _n_beta( n_beta )
        , _n_mu( n_mu )
        , _beta0( beta[0] )
        , _dbeta( ( n_beta > 1 ) ? ( beta[1] - beta[0] ) : static_cast<T>(1) )
        , _log_mu( static_cast<std::size_t>(n_mu)   )
        , _fcn   ( static_cast<std::size_t>(n_beta * n_mu) )
        , _m2    ( static_cast<std::size_t>(n_beta * n_mu), static_cast<T>(0) )
    {
        if ( n_beta < 4 || n_mu < 3 )
            throw std::invalid_argument( "NaturalCubicSpline2D: need ≥4 beta and ≥3 mu points" );

        // Transform and copy axes / function values
        for ( int j = 0; j < n_mu; ++j )
            _log_mu[j] = std::log10( mu[j] );

        std::copy( fcn, fcn + n_beta * n_mu, _fcn.begin() );

        // ---- Precompute second derivatives (one tridiagonal per beta row) ----
        //
        // Natural cubic spline conditions: m[0] = m[n_mu-1] = 0.
        // Interior system for j = 1 .. n_mu-2  (n = n_mu-2 equations):
        //
        //   h[j-1]*m[j-1] + 2*(h[j-1]+h[j])*m[j] + h[j]*m[j+1]
        //     = 6 * [ (f[j+1]-f[j])/h[j] - (f[j]-f[j-1])/h[j-1] ]
        //
        const int nm  = n_mu;
        const int n   = nm - 2;        // interior knots

        // Precompute interval widths h[j] = log_mu[j+1] - log_mu[j]
        std::vector<T> h( static_cast<std::size_t>(nm - 1) );
        for ( int j = 0; j < nm - 1; ++j )
            h[j] = _log_mu[j+1] - _log_mu[j];

        // Build the (constant) tridiagonal structure
        std::vector<T> diag( static_cast<std::size_t>(n) );
        std::vector<T> off ( static_cast<std::size_t>(n - 1) );
        for ( int k = 0; k < n; ++k )
        {
            diag[k] = static_cast<T>(2) * ( h[k] + h[k+1] );
            if ( k < n - 1 ) off[k] = h[k+1];
        }

        // Allocate temporary RHS and solution vectors
        std::vector<T> rhs( static_cast<std::size_t>(n) );
        std::vector<T> sol( static_cast<std::size_t>(n) );

        // One tridiagonal solve per beta row
        for ( int i = 0; i < n_beta; ++i )
        {
            const T* fi = _fcn.data() + i * nm;

            for ( int k = 0; k < n; ++k )
            {
                const int j = k + 1;   // interior knot index
                rhs[k] = static_cast<T>(6) * (
                            ( fi[j+1] - fi[j]   ) / h[j]
                          - ( fi[j]   - fi[j-1] ) / h[j-1]
                         );
            }

            interp2d_detail::thomas_solve(
                diag.data(), off.data(), rhs.data(), sol.data(), n );

            // Natural BC: boundary second derivatives are 0 (already initialised)
            T* mi = _m2.data() + i * nm;
            for ( int k = 0; k < n; ++k )
                mi[k+1] = sol[k];
        }
    }

    // Evaluate the natural cubic spline in log_mu for row ri.
    T _eval_row( int ri, T log_mu ) const noexcept
    {
        const int nm = _n_mu;
        const int iy = interp2d_detail::find_cell_nonuniform( _log_mu.data(), nm, log_mu );

        const T h  = _log_mu[iy+1] - _log_mu[iy];
        const T s  = log_mu - _log_mu[iy];

        const T* fi = _fcn.data() + ri * nm;
        const T* mi = _m2.data()  + ri * nm;

        // Horner form: a + s*(b + s*(c + s*d))
        const T a = fi[iy];
        const T b = ( fi[iy+1] - fi[iy] ) / h
                  - h * ( static_cast<T>(2) * mi[iy] + mi[iy+1] ) / static_cast<T>(6);
        const T c = mi[iy] / static_cast<T>(2);
        const T d = ( mi[iy+1] - mi[iy] ) / ( static_cast<T>(6) * h );

        return a + s * ( b + s * ( c + s * d ) );
    }

    T eval( T beta, T log_mu ) const noexcept
    {
        const int nb = _n_beta;
        const int ix = interp2d_detail::find_cell_uniform( beta, _beta0, _dbeta, nb - 1 );

        const T beta_lo = _beta0 + static_cast<T>(ix) * _dbeta;
        const T t = ( beta - beta_lo ) / _dbeta;

        // Catmull-Rom across 4 surrounding beta rows
        const int r0 = interp2d_detail::clamp_idx( ix - 1, 0, nb - 1 );
        const int r1 = ix;
        const int r2 = interp2d_detail::clamp_idx( ix + 1, 0, nb - 1 );
        const int r3 = interp2d_detail::clamp_idx( ix + 2, 0, nb - 1 );

        const T v0 = _eval_row( r0, log_mu );
        const T v1 = _eval_row( r1, log_mu );
        const T v2 = _eval_row( r2, log_mu );
        const T v3 = _eval_row( r3, log_mu );

        T wx[4];
        interp2d_detail::catmull_rom_weights( t, wx[0], wx[1], wx[2], wx[3] );

        return wx[0]*v0 + wx[1]*v1 + wx[2]*v2 + wx[3]*v3;
    }
};
