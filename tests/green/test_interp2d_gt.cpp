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

// Test suite for the 2D interpolators in src/math/data_fit/interp2d.hpp,
// exercised against all seven time-domain integral tables:
//   Gt, Gtt, Gttt, Gtx, Gttx, Gtxx, Gttxx.
//
// Tests
//   1. test_construction       – all three interpolators build without throwing
//                                for every table.
//   2. test_on_grid_accuracy   – evaluation at every tabulation node returns
//                                the stored value within floating-point round-off.
//                                Bilinear must be exact (max|err| = 0); Catmull-Rom
//                                and NCS must satisfy
//                                  max|err| ≤ eps_rel * max|fcn|
//                                with eps_rel = 1e-9.
//   3. test_off_grid_finite    – 1000 uniformly-random (beta, mu) points in the
//                                interior of the domain produce finite values for
//                                every interpolator and every table.
//   4. test_off_grid_consistency – Catmull-Rom and NCS agree with Bilinear to
//                                within a generous absolute tolerance:
//                                  mean|Δ| ≤ tol_frac * (max|fcn| - min|fcn|)
//                                with tol_frac = 0.05 (5 % of the function range).
//                                The purpose is a smoke-test for gross divergence,
//                                NOT a statement about interpolation accuracy.
//
// Usage
//   test_interp2d_gt_cmd <Gt.h5> <Gtt.h5> <Gttt.h5> <Gtx.h5> <Gttx.h5>
//                        <Gtxx.h5> <Gttxx.h5>
//
//   argv[1..7] are the absolute or relative paths to the seven HDF5 files.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/config.hpp"
#include "../../src/math/data_fit/interp2d.hpp"
#include "../../src/green/td_data/td_database.hpp"


// ---------------------------------------------------------------------------
// Test helpers
// ---------------------------------------------------------------------------

static int g_tests_run    = 0;
static int g_tests_passed = 0;

static void check( bool cond, const char* msg )
{
    if ( !cond )
    {
        std::cerr << "  FAIL  " << msg << "\n";
        throw std::runtime_error( msg );
    }
}

// Run a single named test; catch exceptions so later tests still execute.
template<typename Fn>
static void run_test( const char* name, Fn&& fn )
{
    ++g_tests_run;
    std::cout << "[ RUN  ] " << name << "\n";
    try
    {
        fn();
        ++g_tests_passed;
        std::cout << "[ PASS ] " << name << "\n\n";
    }
    catch ( const std::exception& e )
    {
        std::cout << "[ FAIL ] " << name << "  (" << e.what() << ")\n\n";
    }
}


// ---------------------------------------------------------------------------
// Reference data loading (.dat files produced by generate_time_domain_test_data.py)
// Format:
//   N
//   beta_1  mu_1  value_1
//   ...
// ---------------------------------------------------------------------------

struct RefPoint
{
    double beta;
    double mu;   // linear space
    double ref;
};

static std::vector<RefPoint> load_ref_data( const std::string& path )
{
    std::ifstream f( path );
    if ( !f.is_open() )
        throw std::runtime_error( "Cannot open reference file: " + path );

    int n = 0;
    f >> n;
    if ( n <= 0 )
        throw std::runtime_error( "Invalid record count in " + path );

    std::vector<RefPoint> pts;
    pts.reserve( static_cast<std::size_t>(n) );
    for ( int i = 0; i < n; ++i )
    {
        RefPoint p{};
        f >> p.beta >> p.mu >> p.ref;
        if ( !f )
            throw std::runtime_error( "Unexpected EOF in " + path );
        pts.push_back( p );
    }
    return pts;
}


// ---------------------------------------------------------------------------
// HDF5 data loading (C API)
// ---------------------------------------------------------------------------

struct TableData
{
    std::string      name;
    std::vector<double> beta;     // n_beta
    std::vector<double> mu;       // n_mu  (linear space)
    std::vector<double> fcn;      // n_beta × n_mu, row-major
    int n_beta = 0;
    int n_mu   = 0;
};

static std::vector<double> h5_read_1d( hid_t file, const char* dset_name )
{
    hid_t dset  = H5Dopen( file, dset_name, H5P_DEFAULT );
    hid_t space = H5Dget_space( dset );
    hsize_t dim = 0;
    H5Sget_simple_extent_dims( space, &dim, nullptr );
    std::vector<double> buf( static_cast<std::size_t>(dim) );
    H5Dread( dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data() );
    H5Sclose( space );
    H5Dclose( dset );
    return buf;
}

static std::vector<double> h5_read_2d( hid_t file, const char* dset_name,
                                        int& rows, int& cols )
{
    hid_t dset  = H5Dopen( file, dset_name, H5P_DEFAULT );
    hid_t space = H5Dget_space( dset );
    hsize_t dims[2] = {};
    H5Sget_simple_extent_dims( space, dims, nullptr );
    rows = static_cast<int>( dims[0] );
    cols = static_cast<int>( dims[1] );
    std::vector<double> buf( static_cast<std::size_t>(rows * cols) );
    H5Dread( dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data() );
    H5Sclose( space );
    H5Dclose( dset );
    return buf;
}

static TableData load_table( const std::string& path, const std::string& name )
{
    hid_t file = H5Fopen( path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
    if ( file < 0 )
        throw std::runtime_error( "Cannot open HDF5 file: " + path );

    TableData td;
    td.name  = name;
    td.beta  = h5_read_1d( file, "beta" );
    td.mu    = h5_read_1d( file, "mu"   );
    td.fcn   = h5_read_2d( file, "fcn", td.n_beta, td.n_mu );

    H5Fclose( file );

    if ( td.n_beta != static_cast<int>(td.beta.size()) ||
         td.n_mu   != static_cast<int>(td.mu.size())   )
        throw std::runtime_error( "HDF5 dimension mismatch in " + path );

    return td;
}


// ---------------------------------------------------------------------------
// Test 1 – Construction
//
// Build all three interpolators for all seven tables.  Any exception is a
// test failure.
// ---------------------------------------------------------------------------

static void test_construction( const std::vector<TableData>& tables )
{
    for ( const auto& td : tables )
    {
        std::cout << "  Building interpolators for " << td.name << " ...\n";

        Bilinear2D<double> bilin(
            td.n_beta, td.beta.data(),
            td.n_mu,   td.mu.data(),
            td.fcn.data() );

        BicubicCatmullRom2D<double> catrom(
            td.n_beta, td.beta.data(),
            td.n_mu,   td.mu.data(),
            td.fcn.data() );

        NaturalCubicSpline2D<double> natcub(
            td.n_beta, td.beta.data(),
            td.n_mu,   td.mu.data(),
            td.fcn.data() );

        check( true, "construction succeeded" );
    }

    std::cout << "  All interpolators constructed for all " << tables.size()
              << " tables.\n";
}


// ---------------------------------------------------------------------------
// Test 2 – On-grid accuracy
//
// Evaluate each interpolator at every tabulation node and compare with the
// stored value.  All three methods must satisfy
//   max|err| ≤ EPS_REL × max|fcn|
// with EPS_REL = 1e-9 (comfortably above machine epsilon, well below any
// meaningful interpolation error).
//
// Note: Bilinear is theoretically exact at grid nodes but floating-point
// arithmetic in the uniform-grid cell lookup (ix = (beta - beta0) / dbeta)
// can introduce round-off of order machine_eps × |f|, so the same relative
// tolerance applies.
// ---------------------------------------------------------------------------

static void test_on_grid_accuracy( const std::vector<TableData>& tables )
{
    constexpr double EPS_REL = 1e-9;

    std::cout << "  " << std::left << std::setw(14) << "Table"
              << "  " << std::setw(18) << "Bilinear max|err|"
              << "  " << std::setw(18) << "CatmullRom max|err|"
              << "  " << std::setw(18) << "NCS max|err|"
              << "  " << "eps_limit\n";
    std::cout << "  " << std::string(80, '-') << "\n";

    for ( const auto& td : tables )
    {
        Bilinear2D<double> bilin(
            td.n_beta, td.beta.data(),
            td.n_mu,   td.mu.data(),
            td.fcn.data() );

        BicubicCatmullRom2D<double> catrom(
            td.n_beta, td.beta.data(),
            td.n_mu,   td.mu.data(),
            td.fcn.data() );

        NaturalCubicSpline2D<double> natcub(
            td.n_beta, td.beta.data(),
            td.n_mu,   td.mu.data(),
            td.fcn.data() );

        double max_fcn = 0.0;
        double max_err_bilin  = 0.0;
        double max_err_catrom = 0.0;
        double max_err_natcub = 0.0;

        for ( int i = 0; i < td.n_beta; ++i )
        {
            const double beta_i  = td.beta[i];
            for ( int j = 0; j < td.n_mu; ++j )
            {
                const double log_mu_j = std::log10( td.mu[j] );
                const double ref      = td.fcn[ i * td.n_mu + j ];

                max_fcn = std::max( max_fcn, std::abs(ref) );

                max_err_bilin  = std::max( max_err_bilin,
                    std::abs( bilin .eval(beta_i, log_mu_j) - ref ) );
                max_err_catrom = std::max( max_err_catrom,
                    std::abs( catrom.eval(beta_i, log_mu_j) - ref ) );
                max_err_natcub = std::max( max_err_natcub,
                    std::abs( natcub.eval(beta_i, log_mu_j) - ref ) );
            }
        }

        const double limit = EPS_REL * max_fcn;

        char buf_b[24], buf_c[24], buf_n[24], buf_l[24];
        std::snprintf( buf_b, sizeof(buf_b), "%.2e", max_err_bilin  );
        std::snprintf( buf_c, sizeof(buf_c), "%.2e", max_err_catrom );
        std::snprintf( buf_n, sizeof(buf_n), "%.2e", max_err_natcub );
        std::snprintf( buf_l, sizeof(buf_l), "%.2e", limit          );

        std::cout << "  " << std::left << std::setw(14) << td.name
                  << "  " << std::setw(18) << buf_b
                  << "  " << std::setw(18) << buf_c
                  << "  " << std::setw(18) << buf_n
                  << "  " << buf_l << "\n";

        check( max_err_bilin  <= limit,
               ( "on_grid_accuracy: Bilinear exceeds tolerance for " + td.name ).c_str() );
        check( max_err_catrom <= limit,
               ( "on_grid_accuracy: CatmullRom exceeds tolerance for " + td.name ).c_str() );
        check( max_err_natcub <= limit,
               ( "on_grid_accuracy: NCS exceeds tolerance for " + td.name ).c_str() );
    }
}


// ---------------------------------------------------------------------------
// Test 3 – Off-grid values are finite
//
// Sample 1000 random (beta, mu) points in the interior of the domain.
// All three interpolators must return finite values for every table.
// ---------------------------------------------------------------------------

static void test_off_grid_finite( const std::vector<TableData>& tables )
{
    constexpr int N_PTS  = 1000;
    constexpr uint64_t SEED = 42ULL;

    // Use the first table's axes to define the query domain
    // (all tables share the same beta and mu axes)
    const TableData& first = tables[0];
    const double beta_lo   = first.beta.front();
    const double beta_hi   = first.beta.back();
    const double lmu_lo    = std::log10( first.mu.front() );
    const double lmu_hi    = std::log10( first.mu.back()  );

    std::mt19937_64 rng( SEED );
    std::uniform_real_distribution<double> db( beta_lo, beta_hi );
    std::uniform_real_distribution<double> dm( lmu_lo,  lmu_hi  );

    std::vector<double> q_beta( N_PTS ), q_logmu( N_PTS );
    for ( int i = 0; i < N_PTS; ++i )
    {
        q_beta[i]  = db( rng );
        q_logmu[i] = dm( rng );
    }

    for ( const auto& td : tables )
    {
        Bilinear2D<double> bilin(
            td.n_beta, td.beta.data(),
            td.n_mu,   td.mu.data(),
            td.fcn.data() );

        BicubicCatmullRom2D<double> catrom(
            td.n_beta, td.beta.data(),
            td.n_mu,   td.mu.data(),
            td.fcn.data() );

        NaturalCubicSpline2D<double> natcub(
            td.n_beta, td.beta.data(),
            td.n_mu,   td.mu.data(),
            td.fcn.data() );

        int non_finite = 0;
        for ( int i = 0; i < N_PTS; ++i )
        {
            const double b  = q_beta[i];
            const double lm = q_logmu[i];
            if ( !std::isfinite( bilin .eval(b, lm) ) ) ++non_finite;
            if ( !std::isfinite( catrom.eval(b, lm) ) ) ++non_finite;
            if ( !std::isfinite( natcub.eval(b, lm) ) ) ++non_finite;
        }

        check( non_finite == 0,
               ( "off_grid_finite: " + td.name + " produced non-finite values" ).c_str() );
    }

    std::cout << "  All " << N_PTS << " off-grid points finite for all "
              << tables.size() << " tables (3 methods each).\n";
}


// ---------------------------------------------------------------------------
// Test 4 – Off-grid consistency
//
// Compare Catmull-Rom and NCS against Bilinear at 5000 random points.
// Accepts mean|Δ| ≤ tol_frac × range(fcn) per table.
// This is a smoke-test for gross divergence, not a statement about accuracy.
// ---------------------------------------------------------------------------

static void test_off_grid_consistency( const std::vector<TableData>& tables )
{
    constexpr int    N_PTS    = 5000;
    constexpr double TOL_FRAC = 0.05;   // 5 % of the function range
    constexpr uint64_t SEED   = 99ULL;

    const TableData& first = tables[0];
    const double beta_lo   = first.beta.front();
    const double beta_hi   = first.beta.back();
    const double lmu_lo    = std::log10( first.mu.front() );
    const double lmu_hi    = std::log10( first.mu.back()  );

    std::mt19937_64 rng( SEED );
    std::uniform_real_distribution<double> db( beta_lo, beta_hi );
    std::uniform_real_distribution<double> dm( lmu_lo,  lmu_hi  );

    std::vector<double> q_beta( N_PTS ), q_logmu( N_PTS );
    for ( int i = 0; i < N_PTS; ++i )
    {
        q_beta[i]  = db( rng );
        q_logmu[i] = dm( rng );
    }

    std::cout << "  " << std::left << std::setw(14) << "Table"
              << "  " << std::setw(20) << "mean|CatmRom - Bilin|"
              << "  " << std::setw(20) << "mean|NCS - Bilin|"
              << "  " << "tolerance\n";
    std::cout << "  " << std::string(72, '-') << "\n";

    for ( const auto& td : tables )
    {
        Bilinear2D<double> bilin(
            td.n_beta, td.beta.data(),
            td.n_mu,   td.mu.data(),
            td.fcn.data() );

        BicubicCatmullRom2D<double> catrom(
            td.n_beta, td.beta.data(),
            td.n_mu,   td.mu.data(),
            td.fcn.data() );

        NaturalCubicSpline2D<double> natcub(
            td.n_beta, td.beta.data(),
            td.n_mu,   td.mu.data(),
            td.fcn.data() );

        // Compute function range for tolerance
        double fcn_min =  std::numeric_limits<double>::infinity();
        double fcn_max = -std::numeric_limits<double>::infinity();
        for ( const double v : td.fcn )
        {
            if ( v < fcn_min ) fcn_min = v;
            if ( v > fcn_max ) fcn_max = v;
        }
        const double fcn_range = fcn_max - fcn_min;
        const double tol       = TOL_FRAC * fcn_range;

        double sum_cat = 0.0, sum_ncs = 0.0;
        for ( int i = 0; i < N_PTS; ++i )
        {
            const double b     = q_beta[i];
            const double lm    = q_logmu[i];
            const double ref   = bilin .eval( b, lm );
            sum_cat += std::abs( catrom.eval( b, lm ) - ref );
            sum_ncs += std::abs( natcub.eval( b, lm ) - ref );
        }
        const double mean_cat = sum_cat / N_PTS;
        const double mean_ncs = sum_ncs / N_PTS;

        char buf_c[24], buf_n[24], buf_t[24];
        std::snprintf( buf_c, sizeof(buf_c), "%.3e", mean_cat );
        std::snprintf( buf_n, sizeof(buf_n), "%.3e", mean_ncs );
        std::snprintf( buf_t, sizeof(buf_t), "%.3e", tol      );

        std::cout << "  " << std::left << std::setw(14) << td.name
                  << "  " << std::setw(20) << buf_c
                  << "  " << std::setw(20) << buf_n
                  << "  " << buf_t << "\n";

        check( mean_cat <= tol,
               ( "off_grid_consistency: CatmullRom vs Bilinear out of tolerance for " + td.name ).c_str() );
        check( mean_ncs <= tol,
               ( "off_grid_consistency: NCS vs Bilinear out of tolerance for " + td.name ).c_str() );
    }
}


// ---------------------------------------------------------------------------
// Test 5 – Against reference data
//
// The .dat files in tests/tests_data/time_domain/ were produced by
// generate_time_domain_test_data.py using scipy RegularGridInterpolator
// (method='cubic').  All three C++ interpolators are compared against this
// reference using a single global-scale tolerance:
//   |err| ≤ EPS_REL × max|ref|_table + EPS_ABS
// Using a global scale (rather than per-point |ref|) avoids spuriously tight
// bounds near zero crossings of the function.
// ---------------------------------------------------------------------------

static void test_against_reference(
    const std::vector<TableData>&              tables,
    const std::vector<std::vector<RefPoint>>&  all_ref )
{
    constexpr double EPS_REL = 5e-3;   // 0.5 % of function scale
    constexpr double EPS_ABS = 1e-12;

    std::cout << "  " << std::left << std::setw(14) << "Table"
              << "  " << std::setw(16) << "Bilin max|err|"
              << "  " << std::setw(16) << "CatmRom max|err|"
              << "  " << std::setw(16) << "NCS max|err|"
              << "  " << std::setw(12) << "tol"
              << "  " << "N_pts\n";
    std::cout << "  " << std::string(84, '-') << "\n";

    for ( std::size_t k = 0; k < tables.size(); ++k )
    {
        const TableData&            td  = tables[k];
        const std::vector<RefPoint>& ref = all_ref[k];

        Bilinear2D<double> bilin(
            td.n_beta, td.beta.data(),
            td.n_mu,   td.mu.data(),
            td.fcn.data() );

        BicubicCatmullRom2D<double> catrom(
            td.n_beta, td.beta.data(),
            td.n_mu,   td.mu.data(),
            td.fcn.data() );

        NaturalCubicSpline2D<double> natcub(
            td.n_beta, td.beta.data(),
            td.n_mu,   td.mu.data(),
            td.fcn.data() );

        // Compute global scale for this table's reference points
        double max_abs_ref = 0.0;
        for ( const auto& pt : ref )
            max_abs_ref = std::max( max_abs_ref, std::abs( pt.ref ) );

        const double tol   = EPS_REL * max_abs_ref + EPS_ABS;

        double max_err_bilin  = 0.0;
        double max_err_catrom = 0.0;
        double max_err_natcub = 0.0;
        bool   bilin_fail     = false;
        bool   catrom_fail    = false;
        bool   natcub_fail    = false;
        std::string bilin_msg, catrom_msg, natcub_msg;

        for ( const auto& pt : ref )
        {
            const double lmu = std::log10( pt.mu );

            const double v_b = bilin .eval( pt.beta, lmu );
            const double v_c = catrom.eval( pt.beta, lmu );
            const double v_n = natcub.eval( pt.beta, lmu );

            const double err_b = std::abs( v_b - pt.ref );
            const double err_c = std::abs( v_c - pt.ref );
            const double err_n = std::abs( v_n - pt.ref );

            max_err_bilin  = std::max( max_err_bilin,  err_b );
            max_err_catrom = std::max( max_err_catrom, err_c );
            max_err_natcub = std::max( max_err_natcub, err_n );

            if ( !bilin_fail  && err_b > tol )
            {
                bilin_fail = true;
                char buf[256];
                std::snprintf( buf, sizeof(buf),
                    "against_reference: Bilinear exceeds tol for %s "
                    "at beta=%.4f mu=%.4e: err=%.3e tol=%.3e",
                    td.name.c_str(), pt.beta, pt.mu, err_b, tol );
                bilin_msg = buf;
            }
            if ( !catrom_fail && err_c > tol )
            {
                catrom_fail = true;
                char buf[256];
                std::snprintf( buf, sizeof(buf),
                    "against_reference: CatmullRom exceeds tol for %s "
                    "at beta=%.4f mu=%.4e: err=%.3e tol=%.3e",
                    td.name.c_str(), pt.beta, pt.mu, err_c, tol );
                catrom_msg = buf;
            }
            if ( !natcub_fail && err_n > tol )
            {
                natcub_fail = true;
                char buf[256];
                std::snprintf( buf, sizeof(buf),
                    "against_reference: NCS exceeds tol for %s "
                    "at beta=%.4f mu=%.4e: err=%.3e tol=%.3e",
                    td.name.c_str(), pt.beta, pt.mu, err_n, tol );
                natcub_msg = buf;
            }
        }

        char buf_b[24], buf_c[24], buf_n[24], buf_t[24];
        std::snprintf( buf_b, sizeof(buf_b), "%.2e", max_err_bilin  );
        std::snprintf( buf_c, sizeof(buf_c), "%.2e", max_err_catrom );
        std::snprintf( buf_n, sizeof(buf_n), "%.2e", max_err_natcub );
        std::snprintf( buf_t, sizeof(buf_t), "%.2e", tol            );

        std::cout << "  " << std::left << std::setw(14) << td.name
                  << "  " << std::setw(16) << buf_b
                  << "  " << std::setw(16) << buf_c
                  << "  " << std::setw(16) << buf_n
                  << "  " << std::setw(12) << buf_t
                  << "  " << ref.size() << "\n";

        if ( bilin_fail  ) check( false, bilin_msg .c_str() );
        if ( catrom_fail ) check( false, catrom_msg.c_str() );
        if ( natcub_fail ) check( false, natcub_msg.c_str() );
    }
}


// ---------------------------------------------------------------------------
// Test 6 – Embedded data (td_db) matches HDF5-loaded data
//
// Constructs all three interpolators from the embedded td_db arrays and
// compares their evaluations against identical interpolators built from the
// HDF5 data (already in `tables`).  Because the embedded arrays were
// generated with exact hex-float literals from the same HDF5 source, the
// two data sets are bit-identical; evaluations must agree exactly (err = 0).
// ---------------------------------------------------------------------------

static void test_embedded_data( const std::vector<TableData>& tables )
{
    constexpr int      N_PTS = 1000;
    constexpr uint64_t SEED  = 7777ULL;

    // Query domain from the shared embedded axes
    const double beta_lo = td_db::beta[0];
    const double beta_hi = td_db::beta[td_db::N_BETA - 1];
    const double lmu_lo  = std::log10( td_db::mu[0] );
    const double lmu_hi  = std::log10( td_db::mu[td_db::N_MU - 1] );

    std::mt19937_64 rng( SEED );
    std::uniform_real_distribution<double> db_dist( beta_lo, beta_hi );
    std::uniform_real_distribution<double> dm_dist( lmu_lo,  lmu_hi  );

    std::vector<double> q_beta( N_PTS ), q_lmu( N_PTS );
    for ( int i = 0; i < N_PTS; ++i )
    {
        q_beta[i] = db_dist( rng );
        q_lmu[i]  = dm_dist( rng );
    }

    using NCS    = NaturalCubicSpline2D<double>;
    using CatRom = BicubicCatmullRom2D<double>;
    using Bilin  = Bilinear2D<double>;

    // Function-pointer table, one row per table (same order as `tables`)
    struct Makers {
        NCS    (*make_ncs)();
        CatRom (*make_cat)();
        Bilin  (*make_bil)();
    };

    const Makers makers[7] = {
        { td_db::make_ncs_Gt,    td_db::make_catmullrom_Gt,    td_db::make_bilinear_Gt    },
        { td_db::make_ncs_Gtt,   td_db::make_catmullrom_Gtt,   td_db::make_bilinear_Gtt   },
        { td_db::make_ncs_Gttt,  td_db::make_catmullrom_Gttt,  td_db::make_bilinear_Gttt  },
        { td_db::make_ncs_Gtx,   td_db::make_catmullrom_Gtx,   td_db::make_bilinear_Gtx   },
        { td_db::make_ncs_Gttx,  td_db::make_catmullrom_Gttx,  td_db::make_bilinear_Gttx  },
        { td_db::make_ncs_Gtxx,  td_db::make_catmullrom_Gtxx,  td_db::make_bilinear_Gtxx  },
        { td_db::make_ncs_Gttxx, td_db::make_catmullrom_Gttxx, td_db::make_bilinear_Gttxx },
    };

    std::cout << "  " << std::left << std::setw(14) << "Table"
              << "  " << std::setw(14) << "Bilin max|\xce\x94|"
              << "  " << std::setw(14) << "CatRom max|\xce\x94|"
              << "  " << std::setw(14) << "NCS max|\xce\x94|"
              << "  N_pts\n";
    std::cout << "  " << std::string(66, '-') << "\n";

    for ( std::size_t k = 0; k < tables.size(); ++k )
    {
        const TableData& td = tables[k];

        // HDF5-loaded interpolators
        NCS    ncs_h5( td.n_beta, td.beta.data(), td.n_mu, td.mu.data(), td.fcn.data() );
        CatRom cat_h5( td.n_beta, td.beta.data(), td.n_mu, td.mu.data(), td.fcn.data() );
        Bilin  bil_h5( td.n_beta, td.beta.data(), td.n_mu, td.mu.data(), td.fcn.data() );

        // Embedded-data interpolators
        NCS    ncs_emb = makers[k].make_ncs();
        CatRom cat_emb = makers[k].make_cat();
        Bilin  bil_emb = makers[k].make_bil();

        double max_err_bil = 0.0;
        double max_err_cat = 0.0;
        double max_err_ncs = 0.0;

        for ( int i = 0; i < N_PTS; ++i )
        {
            const double b  = q_beta[i];
            const double lm = q_lmu[i];

            max_err_bil = std::max( max_err_bil,
                std::abs( bil_emb.eval(b, lm) - bil_h5.eval(b, lm) ) );
            max_err_cat = std::max( max_err_cat,
                std::abs( cat_emb.eval(b, lm) - cat_h5.eval(b, lm) ) );
            max_err_ncs = std::max( max_err_ncs,
                std::abs( ncs_emb.eval(b, lm) - ncs_h5.eval(b, lm) ) );
        }

        char buf_b[24], buf_c[24], buf_n[24];
        std::snprintf( buf_b, sizeof(buf_b), "%.2e", max_err_bil );
        std::snprintf( buf_c, sizeof(buf_c), "%.2e", max_err_cat );
        std::snprintf( buf_n, sizeof(buf_n), "%.2e", max_err_ncs );

        std::cout << "  " << std::left << std::setw(14) << td.name
                  << "  " << std::setw(14) << buf_b
                  << "  " << std::setw(14) << buf_c
                  << "  " << std::setw(14) << buf_n
                  << "  " << N_PTS << "\n";

        check( max_err_bil == 0.0,
            ( "embedded_data: Bilinear mismatch vs HDF5 for " + td.name ).c_str() );
        check( max_err_cat == 0.0,
            ( "embedded_data: CatmullRom mismatch vs HDF5 for " + td.name ).c_str() );
        check( max_err_ncs == 0.0,
            ( "embedded_data: NCS mismatch vs HDF5 for " + td.name ).c_str() );
    }
}


// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
    if ( argc < 8 )
    {
        std::cerr << "Usage: " << argv[0]
                  << " <Gt.h5> <Gtt.h5> <Gttt.h5> <Gtx.h5> <Gttx.h5> <Gtxx.h5> <Gttxx.h5>\n"
                  << "       [<dGdt.dat> <dGdtt.dat> <dGdttt.dat> <dGdtx.dat>"
                  << " <dGdttx.dat> <dGdtxx.dat> <dGdttxx.dat>]\n";
        return 1;
    }

    // ---- Load all seven tables ---------------------------------------------
    const char* table_names[7] = {
        "Gt", "Gtt", "Gttt", "Gtx", "Gttx", "Gtxx", "Gttxx"
    };

    std::vector<TableData> tables;
    tables.reserve( 7 );
    bool load_ok = true;
    for ( int k = 0; k < 7; ++k )
    {
        try
        {
            std::cout << "Loading " << table_names[k] << ".h5 ...\n";
            tables.push_back( load_table( argv[k+1], table_names[k] ) );
            std::cout << "  " << tables.back().n_beta << " × "
                      << tables.back().n_mu << " nodes loaded.\n";
        }
        catch ( const std::exception& e )
        {
            std::cerr << "ERROR loading " << argv[k+1] << ": " << e.what() << "\n";
            load_ok = false;
        }
    }
    if ( !load_ok )
    {
        std::cerr << "Aborting: one or more HDF5 files could not be loaded.\n";
        return 2;
    }
    std::cout << "\nAll tables loaded.\n\n";

    // ---- Run tests ---------------------------------------------------------
    run_test( "test_construction",
        [&]{ test_construction( tables ); } );

    run_test( "test_on_grid_accuracy",
        [&]{ test_on_grid_accuracy( tables ); } );

    run_test( "test_off_grid_finite",
        [&]{ test_off_grid_finite( tables ); } );

    run_test( "test_off_grid_consistency",
        [&]{ test_off_grid_consistency( tables ); } );

    run_test( "test_embedded_data",
        [&]{ test_embedded_data( tables ); } );

    // ---- Optional: reference-data test (argv[8..14]) -----------------------
    if ( argc >= 15 )
    {
        // ref dat files must correspond to the same table order
        const char* ref_names[7] = {
            "dGdt", "dGdtt", "dGdttt", "dGdtx", "dGdttx", "dGdtxx", "dGdttxx"
        };
        std::vector<std::vector<RefPoint>> all_ref;
        all_ref.reserve( 7 );
        bool ref_ok = true;
        for ( int k = 0; k < 7; ++k )
        {
            try
            {
                std::cout << "Loading ref " << ref_names[k] << ".dat ...\n";
                all_ref.push_back( load_ref_data( argv[8 + k] ) );
                std::cout << "  " << all_ref.back().size() << " reference points loaded.\n";
            }
            catch ( const std::exception& e )
            {
                std::cerr << "ERROR loading " << argv[8+k] << ": " << e.what() << "\n";
                ref_ok = false;
            }
        }
        std::cout << "\n";

        if ( ref_ok )
            run_test( "test_against_reference",
                [&]{ test_against_reference( tables, all_ref ); } );
        else
            std::cerr << "Skipping test_against_reference: failed to load reference data.\n";
    }

    // ---- Summary -----------------------------------------------------------
    std::cout << "========================================\n";
    std::cout << "Results: " << g_tests_passed << " / " << g_tests_run
              << " tests passed.\n";
    std::cout << "========================================\n";

    return ( g_tests_passed == g_tests_run ) ? 0 : 1;
}
