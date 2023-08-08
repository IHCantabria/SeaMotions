
// Include general usage libraries
#include <fstream>
#include <iostream>

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
#include "../../src/green/integrals_db.hpp"
#include "../../src/green/pulsating_inf_depth_series.hpp"
#include "../../src/math/chebyshev.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/math/special_math.hpp"


struct RefData
{
public:
    double* cheby_coeffs    = nullptr;
    int*    cheby_order_x   = nullptr;
    int*    cheby_order_y   = nullptr;
    int     coeffs_np       = 0;
    double* values          = nullptr;
    int     values_np       = 0;
    double* x               = nullptr;
    double* xe              = nullptr;
    double* y               = nullptr;
    double* ye              = nullptr;

    RefData( std::string file_path )
    {
        this->load_data( file_path );
    }

    ~RefData( void )
    {
        // Deallocate polynomial data
        mkl_free( this->cheby_coeffs );
        mkl_free( this->cheby_order_x );
        mkl_free( this->cheby_order_y );

        // Deallocate reference values
        mkl_free( this->xe );
        mkl_free( this->ye );
        mkl_free( this->values );

    }

    void load_data( std::string file_path )
    {
        // Open file unit
        std::fstream infile( file_path );

        if ( !infile.good( ) )
        {
            std::cerr << "IOError: Input file could not be opened or found." << std::endl;
            std::cerr << "  - File: " << file_path << std::endl;
            throw std::runtime_error( "" );
        }

        // Read number of coefficients
        infile >> this->coeffs_np;

        // Allocate space for the coefficient matrixes
        this->cheby_coeffs  = generate_empty_vector<cusfloat>( this->coeffs_np );
        this->cheby_order_x = generate_empty_vector<int>( this->coeffs_np );
        this->cheby_order_y = generate_empty_vector<int>( this->coeffs_np );
        
        for ( int i=0; i<this->coeffs_np; i++ )
        {
            infile >> this->cheby_coeffs[i];
            infile >> this->cheby_order_x[i];
            infile >> this->cheby_order_y[i];
        }

        // Read number of solution points
        infile >> this->values_np;

        // Allocate space for the solution variables
        this->values    = generate_empty_vector<cusfloat>( this->values_np );
        this->x         = generate_empty_vector<cusfloat>( this->values_np );
        this->xe        = generate_empty_vector<cusfloat>( this->values_np );
        this->y         = generate_empty_vector<cusfloat>( this->values_np );
        this->ye        = generate_empty_vector<cusfloat>( this->values_np );

        // Loop over value points to get the data from disk
        for ( int i=0; i<this->values_np; i++ )
        {
            infile >> this->x[i];
            infile >> this->y[i];
            infile >> this->xe[i];
            infile >> this->ye[i];
            infile >> this->values[i];
        }

        // Close file unit
        infile.close( );

    }

};


void launch_real_test( RefData& ref_data )
{
    // Load integrals database
    IntegralsDb idb = IntegralsDb();
    build_integrals_db(idb);

    // Loop over values and compare with the reference data
    cusfloat abs_tol    = 1e-3;
    cusfloat c0         = 0.0;
    cusfloat diff       = 0.0;
    cusfloat diff_num   = 0.0;
    cusfloat R          = 0.0;
    cusfloat sol        = 0.0;
    cusfloat sol_cmp    = 0.0;
    cusfloat sol_num    = 0.0;
    cusfloat xf         = 0.0;
    cusfloat yf         = 0.0;
    for ( int i=0; i<ref_data.values_np; i++ )
    {
        // Evaluate chebyshev polynomials to get the interpolated value
        // Fit input values to boundaries
        xf = idb.r11b_dx->fit_boundary_x( ref_data.x[i] );
        yf = idb.r11b_dx->fit_boundary_y( ref_data.y[i] );

        // Calculate chebyshev part
        // std::cout << "xf: " << xf << " - yf: " << yf << std::endl;
        sol = idb.r11b_dx->calculate_xy_cheby( xf, yf );
        // std::cout << "sol: " << sol << std::endl;

        R  = sqrt(pow2s(xf) + pow2s(yf));
        c0 = (
                +2*besselj1(xf)*log(R+yf)
                -2*besselj0(xf)*xf/R/(R+yf)
                +PI*bessely1(xf)
                -2*besselj1(xf)*log(xf)
                +2*besselj0(xf)/xf
            )*exp(-yf);
        sol_cmp  = exp(-yf)*sol + c0;
        // std::cout << std::setprecision(16) << "t1: " << besselj1(xf) << "- " << log(R+yf) << std::endl;
        // std::cout << std::setprecision(16) << "t2: " << besselj0(xf) << " - " << xf/R/(R+yf) << std::endl;
        // std::cout << std::setprecision(16) << "t3: " << bessely1(xf) << std::endl;
        // std::cout << std::setprecision(16) << "t4: " << besselj1(xf) << " - " << log(xf) << std::endl;
        // std::cout << std::setprecision(16) << "t5: " << besselj0(xf) << "- " << 1/xf << std::endl;
        // std::cout << "sol: " << sol << std::endl;
        // std::cout << "c0: " << c0 << std::endl;
        // std::cout << "exp: " << exp(-yf) << std::endl;

        // Calculate numerical solution
        cusfloat int_value = expint_inf_depth_num_dxtndim(xf, yf);
        // std::cout << std::setprecision(16) << "int_value: " << int_value + 1e-12 << std::endl;
        sol_num = (
                        -2.0/xf*expint_inf_depth_num_dxtndim(xf, yf)
                        + 2.0*yf/(xf*std::sqrt(pow2s(xf)+pow2s(yf)))
                        + PI*std::exp(-yf)*(bessely1(xf)+struve1(xf)-2.0/PI)
                    );
        // std::cout << "sol_num: " << sol_num << std::endl;

        // Compare the interpolated value with the reference one
        diff = std::abs( sol - ref_data.values[i] );
        diff_num = std::abs( sol_cmp - sol_num );

        if ( ( diff > abs_tol ) || ( diff_num > abs_tol ) )
        {
            std::cerr << "ERROR: tools/test_chebyshev_eval/real_space failed!" << std::endl;
            std::cerr << "Xi: " << ref_data.x[i] << " Yi: " << ref_data.y[i] << std::endl;
            std::cerr << "Value: " << ref_data.values[i] << " - Diff: " << diff;
            std::cerr << " - Diff.Num: " << diff_num << std::endl;
            throw std::runtime_error( "" );
        }
    }

}


void launch_unitary_test( RefData& ref_data )
{
    // Loop over values and compare with the reference data
    cusfloat abs_tol    = 1e-12;
    cusfloat diff       = 0.0;
    cusfloat sol        = 0.0;
    for ( int i=0; i<ref_data.values_np; i++ )
    {
        // Evaluate chebyshev polynomials to get the interpolated value
        sol = 0.0;
        for ( int j=0; j<ref_data.coeffs_np; j++ )
        {
            sol += ( 
                        ref_data.cheby_coeffs[j]
                        *
                        chebyshev_poly_raw(ref_data.cheby_order_x[j], ref_data.xe[i])
                        *
                        chebyshev_poly_raw(ref_data.cheby_order_y[j], ref_data.ye[i])
                    );
        }

        // Compare the interpolated value with the reference one
        diff = std::abs( sol - ref_data.values[i] );

        if ( diff > abs_tol )
        {
            std::cerr << "ERROR: tools/test_chebyshev_eval/unitary_space failed!" << std::endl;
            std::cerr << "Xi: " << ref_data.xe[i] << " Yi: " << ref_data.ye[i] << std::endl;
            std::cerr << "Value: " << ref_data.values[i] << " - Diff: " << diff << std::endl;
            throw std::runtime_error( "" );
        }
    }

}


int main( int argc, char* argv[] )
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 1))
    {
        return 1;
    }

    std::string file_path( argv[1] );

    // Load data
    RefData ref_data( file_path );

    // Perform test in the unitary data space
    launch_unitary_test( ref_data );

    // Perform test in the real data space
    launch_real_test( ref_data );

    return 0;
}