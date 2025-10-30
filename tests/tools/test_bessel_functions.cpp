
// Include general usage libraries
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "../containers.hpp"
#include "../../src/config.hpp"
#include "../../src/math/bessel_factory.hpp"
#include "../../src/math/special_math.hpp"
#include "../../src/tools.hpp"


// Set relative precision for Bessel functions.
cusfloat        EPS_BESSEL  = 1e-6;
constexpr int   N           = 16;


bool launch_test(cusfloat (*f)(cusfloat), std::string file_path, cusfloat precision, bool rel_flag)
{
    // Read reference data
    DataRef data_ref;
    data_ref.read_single_channel(file_path);

    // Loop over reference data to check the solution of the 
    // function
    bool pass = true;
    cusfloat diff = 0.0;
    for (int i=0; i<data_ref.num_points; i++)
    {
        if (rel_flag)
        {
            diff = std::abs((f(data_ref.x[i]) - data_ref.y[i])/data_ref.y[i]);
        }
        else
        {
            diff = std::abs(f(data_ref.x[i]) - data_ref.y[i]);
        }

        if (diff > precision)
        {
            std::cout << std::setprecision(15) << "X: " << data_ref.x[i] << " - Fi(x): " << f(data_ref.x[i]);
            std::cout << std::setprecision(15) << " - Y_ref(x): " << data_ref.y[i] << " - Diff: " << diff << std::endl; 
            pass = false;
            break;
        }
    }

    return pass;
}


bool launch_test_factory( void )
{
    // Declare value for test pass
    bool pass = true;

    // Generate a uniform distribution to get values from
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen( rd( ) ); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dist( 1e-5, 1000 );

    // Declare local variables
    cusfloat x  = 0.0;
    cusfloat j0 = 0.0;
    cusfloat j1 = 0.0;
    cusfloat y0 = 0.0;
    cusfloat y1 = 0.0;
    cusfloat h0 = 0.0;
    cusfloat h1 = 0.0;
    cusfloat i0 = 0.0;
    cusfloat i1 = 0.0;
    cusfloat k0 = 0.0;
    cusfloat k1 = 0.0;

    BesselFactory bessel_factory;

    for ( int i=0; i<100; i++ )
    {
        // Get new point to calculate
        x  = dist( gen );

        // Calculate standard Bessel functions
        j0 = besselj0( x );
        j1 = besselj1( x );
        y0 = bessely0( x );
        y1 = bessely1( x );
        h0 = struve0( x );
        h1 = struve1( x );

        // Calculate modified Bessel functions
        i0 = besseli0( x );
        i1 = besseli1( x );
        k0 = besselk0( x );
        k1 = besselk1( x );

        // Calculate values using Bessel factory
        bessel_factory.calculate_series( x );

        // Compare BesselFactory values with the reference ones
        if ( !assert_scalar_equality( j0, bessel_factory.j0, EPS_BESSEL ) ) break;
        if ( !assert_scalar_equality( j1, bessel_factory.j1, EPS_BESSEL ) ) break;
        if ( !assert_scalar_equality( y0, bessel_factory.y0, EPS_BESSEL ) ) break;
        if ( !assert_scalar_equality( y1, bessel_factory.y1, EPS_BESSEL ) ) break;
        if ( !assert_scalar_equality( h0, bessel_factory.struve0, EPS_BESSEL ) ) break;
        if ( !assert_scalar_equality( h1, bessel_factory.struve1, EPS_BESSEL ) ) break;
        if ( !assert_scalar_equality( i0, bessel_factory.i0, EPS_BESSEL ) ) break;
        if ( !assert_scalar_equality( i1, bessel_factory.i1, EPS_BESSEL ) ) break;
        if ( !assert_scalar_equality( k0, bessel_factory.k0, EPS_BESSEL ) ) break;
        if ( !assert_scalar_equality( k1, bessel_factory.k1, EPS_BESSEL ) ) break;
    }

    return pass;
}


template<typename BesselFamily>
bool launch_test_factory_vec( void )
{
    // Declare value for test pass
    bool pass = true;

    // Generate a uniform distribution to get values from
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen( rd( ) ); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dist( 1e-5, 6.0 );

    // Declare local variables
    BesselFamily bessel_factory;
    cusfloat x[N];
    cusfloat j0[N];
    cusfloat j1[N];
    cusfloat y0[N];
    cusfloat y1[N];
    cusfloat h0[N];
    cusfloat h1[N];
    cusfloat i0[N];
    cusfloat i1[N];
    cusfloat k0[N];
    cusfloat k1[N];

    for ( int i=0; i<1000; i++ )
    {
        for ( int j=0; j<N; j++ )
        {
            // Get new point to calculate
            x[j]  = dist( gen );

            // Calculate standard Bessel functions
            j0[j] = besselj0( x[j] );
            j1[j] = besselj1( x[j] );
            y0[j] = bessely0( x[j] );
            y1[j] = bessely1( x[j] );
            h0[j] = struve0( x[j] );
            h1[j] = struve1( x[j] );
    
            // Calculate modified Bessel functions
            i0[j] = besseli0( x[j] );
            i1[j] = besseli1( x[j] );
            k0[j] = besselk0( x[j] );
            k1[j] = besselk1( x[j] );
        }


        // Calculate values using Bessel factory
        bessel_factory.calculate_series( x );

        // Compare BesselFactory values with the reference ones
        for ( int j=0; j<N; j++ )
        {
            if ( !assert_scalar_equality( j0[j], bessel_factory.j0[j], EPS_BESSEL ) ) break;
            if ( !assert_scalar_equality( j1[j], bessel_factory.j1[j], EPS_BESSEL ) ) break;
            if ( !assert_scalar_equality( y0[j], bessel_factory.y0[j], EPS_BESSEL ) ) break;
            if ( !assert_scalar_equality( y1[j], bessel_factory.y1[j], EPS_BESSEL ) ) break;
            if ( !assert_scalar_equality( h0[j], bessel_factory.struve0[j], EPS_BESSEL ) ) break;
            if ( !assert_scalar_equality( h1[j], bessel_factory.struve1[j], EPS_BESSEL ) ) break;
            if ( !assert_scalar_equality( i0[j], bessel_factory.i0[j], EPS_BESSEL ) ) break;
            if ( !assert_scalar_equality( i1[j], bessel_factory.i1[j], EPS_BESSEL ) ) break;
            if ( !assert_scalar_equality( k0[j], bessel_factory.k0[j], EPS_BESSEL ) ) break;
            if ( !assert_scalar_equality( k1[j], bessel_factory.k1[j], EPS_BESSEL ) ) break;
        }
    }

    return pass;
}


int main(int argc, char* argv[])
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 10))
    {
        return 1;
    }

    std::string file_path_besselj0(argv[1]);
    std::string file_path_besselj1(argv[2]);
    std::string file_path_bessely0(argv[3]);
    std::string file_path_bessely1(argv[4]);
    std::string file_path_bessels0(argv[5]);
    std::string file_path_bessels1(argv[6]);
    std::string file_path_besseli0(argv[7]);
    std::string file_path_besseli1(argv[8]);
    std::string file_path_besselk0(argv[9]);
    std::string file_path_besselk1(argv[10]);

    // Declare local variables
    bool pass = false;

    // Test first kind first order Bessel function
    pass = launch_test(besselj0, file_path_besselj0, EPS_BESSEL, false);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/besselj0 failed!" << std::endl;
        return 1;
    }

    // Test first kind second order Bessel function
    pass = launch_test(besselj1, file_path_besselj1, EPS_BESSEL, false);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/besselj1 failed!" << std::endl;
        return 1;
    }

    // Test second kind first order Bessel function
    pass = launch_test(bessely0, file_path_bessely0, EPS_BESSEL, false);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/bessely0 failed!" << std::endl;
        return 1;
    }

    // Test second kind second order Bessel function
    pass = launch_test(bessely1, file_path_bessely1, EPS_BESSEL, false);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/bessely1 failed!" << std::endl;
        return 1;
    }

    // Test first order struve function
    pass = launch_test(struve0, file_path_bessels0, EPS_BESSEL, false);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/struve0 failed!" << std::endl;
        return 1;
    }

    // Test second order struve function
    pass = launch_test(struve1, file_path_bessels1, EPS_BESSEL, false);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/struve1 failed!" << std::endl;
        return 1;
    }

    // Test mofified first kind first order Bessel function
    pass = launch_test(besseli0, file_path_besseli0, EPS_BESSEL, true);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/besseli0 failed!" << std::endl;
        return 1;
    }

    // Test modified first kind second order Bessel function
    pass = launch_test(besseli1, file_path_besseli1, EPS_BESSEL, true);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/besseli1 failed!" << std::endl;
        return 1;
    }

    // Test modified second kind first order Bessel function
    pass = launch_test(besselk0, file_path_besselk0, EPS_BESSEL, true);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/besselk0 failed!" << std::endl;
        return 1;
    }

    // Test modified second kind second order Bessel function
    pass = launch_test(besselk1, file_path_besselk1, EPS_BESSEL, true);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/besselk1 failed!" << std::endl;
        return 1;
    }

    // Test BesselFactory object
    pass = launch_test_factory(  );
    if ( !pass )
    {
        std::cerr << "test_bessel_factory failed!" << std::endl;
        return 1;
    }

    // Test BesselFactoryVec object
    pass = launch_test_factory_vec<BesselFactoryVec<N>>( );
    if ( !pass )
    {
        std::cerr << "test_bessel_factory_vec failed!" << std::endl;
        return 1;
    }

    pass = launch_test_factory_vec<BesselFactoryBranch<N>>( );
    if ( !pass )
    {
        std::cerr << "test_bessel_factory_branch failed!" << std::endl;
        return 1;
    }

    std::cout << "Program Finished!" << std::endl;

    return 0;
}