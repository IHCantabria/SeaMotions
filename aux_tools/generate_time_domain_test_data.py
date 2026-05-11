#
# Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
#
# This file is part of SeaMotions Software.
#
# SeaMotions is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SeaMotions is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#

"""
Generate reduced test data for the time-domain Green's function evaluators.

This script reads the full-resolution integral databases located in
    aux_data/0_integrals_database/1_time_domain/
and writes a small subset of (beta, mu, value) test points to plain-text
files under
    tests/tests_data/time_domain/

The test points are chosen so that:
  * 60 beta values are uniformly distributed inside the C++ evaluator domain
    [beta_min, beta_max-0.5], avoiding boundary singularities.
  * 10 mu values are quasi-logarithmically spaced to cover all four y-blocks
    used by the Chebyshev fit (log10(mu) in [-4, 0]).
  * Each point is obtained by bilinear interpolation on the original grid via
    scipy.interpolate.RegularGridInterpolator.

File format  (for each function):
    n_points
    beta_1  mu_1  value_1
    beta_2  mu_2  value_2
    ...

Mapping from HDF5 file to C++ evaluator function:
    Gt.h5    →  eval_dGdt    (alpha=[0.0])
    Gtx.h5   →  eval_dGdtx   (alpha=[1e-3, 1e-2, 1e-1])
    Gtxx.h5  →  eval_dGdtxx  (alpha=[1e-3, 1e-2, 1e-1])
    Gtt.h5   →  eval_dGdtt   (alpha=[1e-3, 1e-2, 1e-1]; beta_max=19)
    Gttx.h5  →  eval_dGdttx  (alpha=[1e-3, 1e-2, 1e-1])
    Gttxx.h5 →  eval_dGdttxx (alpha=[1e-3, 1e-2, 1e-1])
"""

import os
import re
import sys

import h5py
import numpy as np
import scipy.interpolate
import scipy.special


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

def _root_path() -> str:
    """Return the project root directory (parent of aux_tools/)."""
    return os.path.abspath( os.path.join( os.path.dirname( __file__ ), ".." ) )


def _db_path() -> str:
    return os.path.join( _root_path(), "aux_data", "0_integrals_database", "1_time_domain" )


def _out_path() -> str:
    return os.path.join( _root_path(), "tests", "tests_data", "time_domain" )


# ---------------------------------------------------------------------------
# Test-point selection
# ---------------------------------------------------------------------------

# 10 mu values quasi-logarithmically spaced to cover all four y-blocks
# log10(mu) blocks: [-4,-3], [-3,-2], [-2,-1], [-1,~0]
# 2-3 values per block, away from block boundaries.
MU_TEST = np.array( [
    1.5e-4,   # log10 ≈ -3.82  (block 0)
    6.0e-4,   # log10 ≈ -3.22  (block 0)
    2.0e-3,   # log10 ≈ -2.70  (block 1)
    7.0e-3,   # log10 ≈ -2.15  (block 1)
    2.0e-2,   # log10 ≈ -1.70  (block 2)
    7.0e-2,   # log10 ≈ -1.15  (block 2)
    0.20,     # log10 ≈ -0.70  (block 3)
    0.50,     # log10 ≈ -0.30  (block 3)
    0.80,     # log10 ≈ -0.097 (block 3)
    0.95,     # log10 ≈ -0.022 (block 3)
] )

N_BETA  = 60  # number of beta test points per function
N_MU    = len( MU_TEST )  # 10

# C++ evaluator beta_max for each function
CPP_BETA_MAX = {
    "dGdt"   : 30.0,
    "dGdtx"  : 30.0,
    "dGdtxx" : 30.0,
    "dGdtt"  : 19.0,
    "dGdttx" : 30.0,
    "dGdttxx": 30.0,
}

# HDF5 file → function name mapping
DB_MAP = [
    ( "Gt.h5",    "dGdt"    ),
    ( "Gtx.h5",   "dGdtx"   ),
    ( "Gtxx.h5",  "dGdtxx"  ),
    ( "Gtt.h5",   "dGdtt"   ),
    ( "Gttx.h5",  "dGdttx"  ),
    ( "Gttxx.h5", "dGdttxx" ),
]


# ---------------------------------------------------------------------------
# 1D Chebyshev evaluator  (mirrors C++ eval_time_1d)
# ---------------------------------------------------------------------------

def _parse_1d_header( hpp_path: str ) -> dict:
    """
    Parse a 1D Chebyshev coefficient header (A-type, e.g. dGdttxA0.hpp) and
    return a dict with all data needed by _eval_time_1d.
    """
    txt = open( hpp_path ).read()

    def get_int( name ):
        m = re.search( r'\b' + name + r'\b\s*=\s*(\d+)', txt )
        return int( m.group(1) ) if m else None

    def get_float( name ):
        m = re.search( r'\b' + name + r'\b\s*=\s*([+-]?[\d.]+(?:[Ee][+-]?\d+)?)', txt )
        return float( m.group(1) ) if m else None

    def get_array_int( name ):
        m = re.search( r'\b' + name + r'\[(\d+)\]\s*=\s*\{([^}]+)\}', txt, re.DOTALL )
        if not m:
            return None
        # Strip // comments then split on whitespace-or-comma
        body = re.sub( r'//[^\n]*', '', m.group(2) )
        return [ int(v) for v in re.split( r'[\s,]+', body.strip() ) if v ]

    def get_array_float( name ):
        m = re.search( r'\b' + name + r'\[(\d+)\]\s*=\s*\{([^}]+)\}', txt, re.DOTALL )
        if not m:
            return None
        body = re.sub( r'//[^\n]*', '', m.group(2) )
        vals = []
        for v in re.split( r'[\s,]+', body.strip() ):
            if v:
                try:
                    vals.append( float( v ) )
                except ValueError:
                    pass
        return vals

    # In the header, Chebyshev coefficients are stored as `c[N]` and orders as
    # `ncx[N]`.  The traits wrapper exposes them as `coeffs` / `ncx`.
    return {
        'intervals_np'   : get_int( 'intervals_np' ),
        'x_min_global'   : get_float( 'x_min_global' ),
        'dx_min_region'  : get_float( 'dx_min_region' ),
        'blocks_start'   : get_array_int( 'blocks_start' ),
        'blocks_coeffs_np': get_array_int( 'blocks_coeffs_np' ),
        'x_min_region'   : get_array_float( 'x_min_region' ),
        'dx_region'      : get_array_float( 'dx_region' ),
        'c'              : get_array_float( 'c' ),          # Chebyshev coeffs
        'ncx'            : get_array_int( 'ncx' ),          # Chebyshev orders
        'max_cheby_order': get_int( 'max_cheby_order' ),
        'x_max_global'   : get_float( 'x_max_global' ),
    }


def _eval_time_1d( hdr: dict, log_mu: float ) -> float:
    """
    Evaluate a 1D piecewise Chebyshev function at log_mu.
    Mirrors C++ eval_time_1d<AT>(log_mu).
    """
    x_min = hdr['x_min_global']
    x_max = hdr['x_max_global']
    dx    = hdr['dx_min_region']
    n_int = hdr['intervals_np']

    x  = max( min( log_mu, x_max ), x_min )
    nx = int( np.floor( ( x - x_min ) / dx ) )
    nx = min( max( nx, 0 ), n_int - 1 )

    sp = hdr['blocks_start'][nx]
    np_ = hdr['blocks_coeffs_np'][nx]
    xm = 2.0 * ( x - hdr['x_min_region'][nx] ) / hdr['dx_region'][nx] - 1.0

    # Evaluate Chebyshev polynomials up to required order
    max_ord = hdr['max_cheby_order']
    T = [0.0] * ( max_ord + 1 )
    T[0] = 1.0
    if max_ord >= 1:
        T[1] = xm
    for k in range( 2, max_ord + 1 ):
        T[k] = 2.0 * xm * T[k-1] - T[k-2]

    result = 0.0
    for i in range( np_ ):
        result += hdr['c'][sp + i] * T[ hdr['ncx'][sp + i] ]

    return result


def _make_1d_evaluator( hpp_path: str ):
    """Parse a header file and return a callable A_k(log_mu)."""
    hdr = _parse_1d_header( hpp_path )
    return lambda log_mu: _eval_time_1d( hdr, log_mu )


# ---------------------------------------------------------------------------
# Locate coefficient header files
# ---------------------------------------------------------------------------

def _coeffs_dir() -> str:
    return os.path.join( _root_path(), "src", "green", "inf_depth_coeffs", "time" )


def _get_A_evaluators( fcn_name: str ):
    """
    Return a list of callables [A0(log_mu), ...] for the given function.
    For functions with constant coefficients (dGdt, dGdtx, dGdtxx, dGdtt)
    the callables simply return the constant value.
    For dGdttx and dGdttxx the callables evaluate the Chebyshev polynomial.
    """
    d = _coeffs_dir()
    # Header file names follow the pattern <fcn_name>A0.hpp, A1, A2
    a_files = {
        "dGdt"    : ["dGdtA0.hpp"],
        "dGdtx"   : ["dGdtxA0.hpp",  "dGdtxA1.hpp",  "dGdtxA2.hpp"],
        "dGdtxx"  : ["dGdtxxA0.hpp", "dGdtxxA1.hpp", "dGdtxxA2.hpp"],
        "dGdtt"   : ["dGdttA0.hpp",  "dGdttA1.hpp",  "dGdttA2.hpp"],
        "dGdttx"  : ["dGdttxA0.hpp", "dGdttxA1.hpp", "dGdttxA2.hpp"],
        "dGdttxx" : ["dGdttxxA0.hpp","dGdttxxA1.hpp","dGdttxxA2.hpp"],
    }
    return [ _make_1d_evaluator( os.path.join( d, fn ) ) for fn in a_files[fcn_name] ]


# ---------------------------------------------------------------------------
# G0 basis functions  (scalar, matches C++ time_domain_evaluator.hpp)
#
# These reproduce exactly the inline functions in time_domain_evaluator.hpp
# so that Python can compute G0 = sum_k A_k * basis(beta, mu, alpha_k).
# ---------------------------------------------------------------------------

_JV = scipy.special.jv  # fractional Bessel J

def _dGdt_G0_basis( beta: float, mu: float, alpha: float = 0.0 ) -> float:
    if beta < 1e-1:
        return 0.0
    x2  = beta * beta / 8.0 + alpha
    lt  = np.pi * beta**3 / ( 16.0 * np.sqrt( 2.0 ) ) * np.exp( -beta * beta * mu / 4.0 )
    bt  = ( _JV(  0.25, x2 ) * _JV( -0.25, x2 )
          + _JV(  0.75, x2 ) * _JV( -0.75, x2 ) )
    if np.isnan( bt ):
        bt = 0.0
    return float( lt * bt )

def _dGdtx_G0_basis( beta: float, mu: float, alpha: float = 0.0 ) -> float:
    return ( -beta * beta / 4.0 ) * _dGdt_G0_basis( beta, mu, alpha )

def _dGdtxx_G0_basis( beta: float, mu: float, alpha: float = 0.0 ) -> float:
    return ( -beta * beta / 4.0 ) * _dGdtx_G0_basis( beta, mu, alpha )

def _dGdtt_G0_basis( beta: float, mu: float, alpha: float = 0.0 ) -> float:
    lt   = np.pi / ( 16.0 * np.sqrt( 2.0 ) )
    expt = np.exp( -beta * beta * mu / 4.0 )
    x    = beta * beta / 8.0 + alpha
    if x < 1e-6:
        x = 0.0

    # Term T1: from differentiating the exponential envelope
    y0 = ( _JV(  0.25, x ) * _JV( -0.25, x )
         + _JV(  0.75, x ) * _JV( -0.75, x ) )
    if np.isnan( y0 ):
        y0 = 0.0
    T1 = -0.5 * beta**3 * float( y0 ) * beta * mu * expt

    # Term T2: from differentiating the Bessel products
    if beta < 5e-2:
        yt = 0.0
    else:
        a  = 3.0 * beta * beta
        b  = beta**4 / 4.0
        c  = b / 2.0
        yt = ( a * ( _JV(  0.25, x ) * _JV( -0.25, x ) + _JV(  0.75, x ) * _JV( -0.75, x ) )
             + b * ( _JV( -0.75, x ) * _JV( -0.25, x ) - _JV(  0.75, x ) * _JV(  0.25, x ) )
             + c * ( -_JV(  1.25, x ) * _JV( -0.25, x ) + _JV( -1.25, x ) * _JV(  0.25, x ) )
             + c * ( -_JV(  1.75, x ) * _JV( -0.75, x ) + _JV( -1.75, x ) * _JV(  0.75, x ) ) )
        yt = float( yt )
        if np.isnan( yt ):
            yt = 0.0

    return lt * ( T1 + yt * expt )

def _dGdttx_G0_basis( beta: float, mu: float, alpha: float = 0.0 ) -> float:
    return ( -beta * beta / 4.0 ) * _dGdtt_G0_basis( beta, mu, alpha )

def _dGdttxx_G0_basis( beta: float, mu: float, alpha: float = 0.0 ) -> float:
    return ( -beta * beta / 4.0 ) * _dGdttx_G0_basis( beta, mu, alpha )


# ---------------------------------------------------------------------------
# G0 amplitude functions — loaded lazily from coefficient header files
# ---------------------------------------------------------------------------

_G0_ALPHA_SHIFTS = {
    "dGdt"    : [ 0.0 ],
    "dGdtx"   : [ 1e-3, 1e-2, 1e-1 ],
    "dGdtxx"  : [ 1e-3, 1e-2, 1e-1 ],
    "dGdtt"   : [ 1e-3, 1e-2, 1e-1 ],
    "dGdttx"  : [ 1e-3, 1e-2, 1e-1 ],
    "dGdttxx" : [ 1e-3, 1e-2, 1e-1 ],
}

_G0_BASIS_FCN = {
    "dGdt"    : _dGdt_G0_basis,
    "dGdtx"   : _dGdtx_G0_basis,
    "dGdtxx"  : _dGdtxx_G0_basis,
    "dGdtt"   : _dGdtt_G0_basis,
    "dGdttx"  : _dGdttx_G0_basis,
    "dGdttxx" : _dGdttxx_G0_basis,
}

_A_EVALUATORS_CACHE: dict = {}


def _get_cached_A_evaluators( fcn_name: str ):
    if fcn_name not in _A_EVALUATORS_CACHE:
        _A_EVALUATORS_CACHE[fcn_name] = _get_A_evaluators( fcn_name )
    return _A_EVALUATORS_CACHE[fcn_name]


def _compute_G0( function_name: str, beta: float, mu: float ) -> float:
    """Evaluate the G0 analytic term for one (beta, mu) point.
    
    A_k coefficients are evaluated from the actual C++ Chebyshev coefficient
    headers via eval_time_1d, so results match the C++ eval_dGdt_G0 etc.
    exactly.
    """
    alphas   = _G0_ALPHA_SHIFTS[function_name]
    basis    = _G0_BASIS_FCN[function_name]
    A_evals  = _get_cached_A_evaluators( function_name )
    log_mu   = float( np.log10( mu ) )
    g0 = 0.0
    for A_fcn, alpha_k in zip( A_evals, alphas ):
        g0 += A_fcn( log_mu ) * basis( beta, mu, alpha_k )
    return g0


# ---------------------------------------------------------------------------
# Core generator
# ---------------------------------------------------------------------------

def generate_test_data_for_function(
    hdf5_filepath: str,
    function_name: str,
    output_dir: str,
    n_beta: int = N_BETA,
    mu_test: np.ndarray = MU_TEST,
) -> str:
    """
    Read an HDF5 database and write a reduced test-data text file.

    Parameters
    ----------
    hdf5_filepath : str
        Path to the HDF5 integral database file.
    function_name : str
        Name label for the output file (e.g. "dGdt").
    output_dir : str
        Directory in which to create the output .dat file.
    n_beta : int
        Number of beta test points.
    mu_test : np.ndarray
        Array of mu values at which to sample.

    Returns
    -------
    str
        Path to the written output file.
    """
    beta_max_cpp = CPP_BETA_MAX[ function_name ]

    # -- Read raw data ---------------------------------------------------
    with h5py.File( hdf5_filepath, "r" ) as fid:
        beta_db = fid[ "beta" ][:]
        mu_db   = fid[ "mu"   ][:]
        fcn_db  = fid[ "fcn"  ][:]

    # -- Build interpolator -----------------------------------------------
    # fcn_db has shape (n_beta_db, n_mu_db)
    interp = scipy.interpolate.RegularGridInterpolator(
        ( beta_db, mu_db ),
        fcn_db,
        method="cubic",
        bounds_error=False,
        fill_value=None,  # extrapolate if needed (shouldn't happen)
    )

    # -- Select beta test points ------------------------------------------
    # Use values within [0.5, beta_max_cpp - 0.5] to stay well inside the
    # C++ evaluator domain.
    beta_min_test = 0.5
    beta_max_test = beta_max_cpp - 0.5
    beta_test = np.linspace( beta_min_test, beta_max_test, n_beta )

    # -- Filter mu_test to database range ---------------------------------
    mu_lo = float( mu_db[0]  ) * 1.001
    mu_hi = float( mu_db[-1] ) * 0.999
    mu_valid = mu_test[ ( mu_test >= mu_lo ) & ( mu_test <= mu_hi ) ]

    if mu_valid.shape[0] == 0:
        raise ValueError( f"No valid mu values for function {function_name}" )

    # -- Evaluate interpolator at test grid -------------------------------
    beta_grid, mu_grid = np.meshgrid( beta_test, mu_valid, indexing="ij" )
    pts = np.stack( [ beta_grid.ravel(), mu_grid.ravel() ], axis=1 )
    values = interp( pts )

    # -- Write output file ------------------------------------------------
    os.makedirs( output_dir, exist_ok=True )
    out_file = os.path.join( output_dir, f"{function_name}.dat" )

    n_points = pts.shape[0]
    with open( out_file, "w" ) as fid:
        fid.write( f"{n_points}\n" )
        for i in range( n_points ):
            fid.write(
                f"{pts[i,0]:.16e}  {pts[i,1]:.16e}  {values[i]:.16e}\n"
            )

    print( f"  Written {n_points} test points → {out_file}" )
    return out_file


# ---------------------------------------------------------------------------
# Decomposed generator (G0 + residual)
# ---------------------------------------------------------------------------

def generate_decomposed_data_for_function(
    hdf5_filepath: str,
    function_name: str,
    output_dir: str,
    n_beta: int = N_BETA,
    mu_test: np.ndarray = MU_TEST,
) -> tuple:
    """
    Write `<function>_G0.dat` and `<function>_residual.dat` at the same test
    points used by generate_test_data_for_function.

    The G0 reference is computed analytically using _compute_G0().
    The residual reference is fcn_db - G0_ref.

    Returns
    -------
    (g0_path, residual_path) : tuple of str
    """
    beta_max_cpp = CPP_BETA_MAX[function_name]

    with h5py.File( hdf5_filepath, "r" ) as fid:
        beta_db = fid["beta"][:]
        mu_db   = fid["mu"][:]
        fcn_db  = fid["fcn"][:]

    interp = scipy.interpolate.RegularGridInterpolator(
        ( beta_db, mu_db ),
        fcn_db,
        method="cubic",
        bounds_error=False,
        fill_value=None,
    )

    beta_min_test = 0.5
    beta_max_test = beta_max_cpp - 0.5
    beta_test     = np.linspace( beta_min_test, beta_max_test, n_beta )

    mu_lo    = float( mu_db[0]  ) * 1.001
    mu_hi    = float( mu_db[-1] ) * 0.999
    mu_valid = mu_test[ ( mu_test >= mu_lo ) & ( mu_test <= mu_hi ) ]

    if mu_valid.shape[0] == 0:
        raise ValueError( f"No valid mu values for function {function_name}" )

    beta_grid, mu_grid = np.meshgrid( beta_test, mu_valid, indexing="ij" )
    pts           = np.stack( [ beta_grid.ravel(), mu_grid.ravel() ], axis=1 )
    total_values  = interp( pts )

    print( f"  Computing G0 at {pts.shape[0]} points for {function_name} ..." )
    g0_values = np.array( [
        _compute_G0( function_name, float( pts[i, 0] ), float( pts[i, 1] ) )
        for i in range( pts.shape[0] )
    ] )
    residual_values = total_values - g0_values

    os.makedirs( output_dir, exist_ok=True )
    n_points = pts.shape[0]

    g0_file  = os.path.join( output_dir, f"{function_name}_G0.dat" )
    res_file = os.path.join( output_dir, f"{function_name}_residual.dat" )

    for filepath, values in [ ( g0_file, g0_values ), ( res_file, residual_values ) ]:
        with open( filepath, "w" ) as fid:
            fid.write( f"{n_points}\n" )
            for i in range( n_points ):
                fid.write(
                    f"{pts[i,0]:.16e}  {pts[i,1]:.16e}  {values[i]:.16e}\n"
                )
        print( f"  Written {n_points} points → {filepath}" )

    return g0_file, res_file


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    db_dir  = _db_path()
    out_dir = _out_path()

    print( f"Reading databases from : {db_dir}" )
    print( f"Writing test data to   : {out_dir}" )
    print()

    for hdf5_file, fcn_name in DB_MAP:
        hdf5_path = os.path.join( db_dir, hdf5_file )
        if not os.path.isfile( hdf5_path ):
            print( f"  WARNING: {hdf5_path} not found – skipping {fcn_name}" )
            continue

        print( f"Processing {hdf5_file} → {fcn_name} ..." )
        generate_test_data_for_function(
            hdf5_filepath = hdf5_path,
            function_name = fcn_name,
            output_dir    = out_dir,
        )
        generate_decomposed_data_for_function(
            hdf5_filepath = hdf5_path,
            function_name = fcn_name,
            output_dir    = out_dir,
        )
        print()

        print()

    print( "Done." )


if __name__ == "__main__":
    main()
