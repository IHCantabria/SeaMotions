
// Include general usage libraries
#include <iostream>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/containers/panel_geom.hpp"
#include "../../src/green/source.hpp"


cusfloat ABS_TOL = 1e-3;


bool compare_formulations( PanelGeom* panel, cusfloat* dir_vec  )
{
    // Define velocity and potential variables
    cusfloat    velocity[3]     = { 0.0, 0.0, 0.0 };
    cusfloat    velocity_m[3]   = { 0.0, 0.0, 0.0 };
    cusfloat    phi             = 0.0;
    cusfloat    phi_m           = 0.0;

    // Loop over different radial distances to check
    // multipole expansions
    cusfloat field_point[3] = { 0.0, 0.0, 0.0 };
    bool     phi_pass       = false;
    cusfloat r              = 0.0;
    bool     vel_pass       = false;

    for ( int i=0; i<30; i++ )
    {
        // Define radial distance
        r = panel->length * 2.5 + cusfloat( i ) ;

        // Define field point
        field_point[0] = r * dir_vec[0] + panel->center[0];
        field_point[1] = r * dir_vec[1] + panel->center[1];
        field_point[2] = r * dir_vec[2] + panel->center[2];

        // Calculate potential and velocity using full integration
        calculate_source_newman_t<1,1,1>(
                                            panel,
                                            field_point,
                                            0,
                                            0,
                                            velocity,
                                            phi
                                        );

        calculate_source_newman_t<1,1,1>(
                                            panel,
                                            field_point,
                                            0,
                                            1,
                                            velocity_m,
                                            phi_m
                                        );
        
        // Check results
        phi_pass    = assert_scalar_equality( phi, phi_m, ABS_TOL, 0.01 );
        vel_pass    = assert_vector_equality( 3, velocity, velocity_m, ABS_TOL, 0.01 );

        //  Calculate errors
        cusfloat phi_abs_err = std::abs( phi - phi_m );
        cusfloat phi_rel_err = std::abs( phi - phi_m ) / phi;
        cusfloat vel_abs_err[3];
        cusfloat vel_rel_err[3];
        for ( int j=0; j<3; j++ )
        {
            vel_abs_err[j] = std::abs( velocity[j] - velocity_m[j] );
            vel_rel_err[j] = std::abs( velocity[j] - velocity_m[j] ) / std::abs( velocity[j] );
        }

        // Print out results to compare in between 
        // full integration and multipole expansion
        if ( false )
        {
            std::cout << "r/panel.length = " << r/panel->length << std::endl;
            std::cout << "    Full integration phi      = " << phi << std::endl;
            std::cout << "    Full integration phi      = " << phi_m << std::endl;
            std::cout << "    Error Abs phi             = " << phi_abs_err << std::endl;
            std::cout << "    Error Rel phi             = " << phi_rel_err << std::endl;
            std::cout << "    Full integration velocity = (" << velocity[0] << ", " << velocity[1] << ", " << velocity[2] << ")" << std::endl;
            std::cout << "    Multipole velocity        = (" << velocity_m[0] << ", " << velocity_m[1] << ", " << velocity_m[2] << ")" << std::endl;
            std::cout << "    Error Abs velocity        = (" << vel_abs_err[0] << ", " << vel_abs_err[1] << ", " << vel_abs_err[2] << ")" << std::endl;
            std::cout << "    Error Rel velocity        = (" << vel_rel_err[0] << ", " << vel_rel_err[1] << ", " << vel_rel_err[2] << ")" << std::endl;
            std::cout << "----------------------------------------" << std::endl;
        }
        
        if ( !phi_pass )
        {
            std::cerr << "Potential test failed at r/panel.length = " << r/panel->length << std::endl;
            std::cerr << "    Full integration phi = " << phi << ", Multipole phi = " << phi_m << std::endl;
            std::cerr << "    Absolute error = " << phi_abs_err << ", Relative error = " << phi_rel_err << std::endl;
            std::cerr << "    Direction vector = (" << dir_vec[0] << ", " << dir_vec[1] << ", " << dir_vec[2] << ")" << std::endl;
            return false;
        }

        if ( !vel_pass )
        {
            std::cerr << "Velocity test failed at r/panel.length = " << r/panel->length << std::endl;
            std::cerr << "    Full integration velocity = (" << velocity[0] << ", " << velocity[1] << ", " << velocity[2] << ")" << std::endl;
            std::cerr << "    Multipole velocity = (" << velocity_m[0] << ", " << velocity_m[1] << ", " << velocity_m[2] << ")" << std::endl;
            return false;
        }
        
    }

    return true;
}


bool test_triangle_0( void )
{
    // Define panel
    constexpr int N = 3;
    cusfloat cog[3] = { 0.0, 0.0, 0.0 };
    cusfloat  x[N]  = { -3.0, 3.0, 3.0 };
    cusfloat  y[N]  = { 0.0, 0.0, 5.0 };
    cusfloat  z[N]  = { -1.0, -3.0, -3.0 };
    PanelGeom panel( 
                        N,
                        x,
                        y,
                        z,
                        true,
                        0,
                        cog
                    );

    // Define directions to test
    cusfloat s2 = std::sqrt( 2.0 );
    cusfloat s3 = std::sqrt( 3.0 );
    cusfloat dir_vec[5][3];
    dir_vec[0][0] = 1.0;        dir_vec[0][1] = 0.0;        dir_vec[0][2] = 0.0;
    dir_vec[1][0] = 0.0;        dir_vec[1][1] = 1.0;        dir_vec[1][2] = 0.0;
    dir_vec[2][0] = 0.0;        dir_vec[2][1] = 0.0;        dir_vec[2][2] = -1.0;
    dir_vec[3][0] = -1/s2;      dir_vec[3][1] = -1/s2;      dir_vec[3][2] = 0.0;
    dir_vec[4][0] = -1/s3;      dir_vec[4][1] = -1/s3;      dir_vec[4][2] = -1/s3;

    // Test formulations in different directions
    bool pass = true;
    for ( int i=0; i<5; i++ )
    {
        pass *= compare_formulations( &panel, dir_vec[i] );
    }   

    return pass;
}


bool test_quadrilateral_0( void )
{
    // Define panel
    constexpr int N = 4;
    cusfloat cog[3] = { 0.0, 0.0, 0.0 };
    cusfloat  x[N]  = { -3.0, 3.0, 3.0, -3.0 };
    cusfloat  y[N]  = { 0.0, 0.0, 5.0, 5.0 };
    cusfloat  z[N]  = { -1.0, -3.0, -3.0, -1.0 };
    PanelGeom panel( 
                        N,
                        x,
                        y,
                        z,
                        true,
                        0,
                        cog
                    );

    // Define directions to test
    cusfloat s2 = std::sqrt( 2.0 );
    cusfloat s3 = std::sqrt( 3.0 );
    cusfloat dir_vec[5][3];
    dir_vec[0][0] = 1.0;        dir_vec[0][1] = 0.0;        dir_vec[0][2] = 0.0;
    dir_vec[1][0] = 0.0;        dir_vec[1][1] = 1.0;        dir_vec[1][2] = 0.0;
    dir_vec[2][0] = 0.0;        dir_vec[2][1] = 0.0;        dir_vec[2][2] = -1.0;
    dir_vec[3][0] = -1/s2;      dir_vec[3][1] = -1/s2;      dir_vec[3][2] = 0.0;
    dir_vec[4][0] = -1/s3;      dir_vec[4][1] = -1/s3;      dir_vec[4][2] = -1/s3;

    // Test formulations in different directions
    bool pass = true;
    for ( int i=0; i<5; i++ )
    {
        pass *= compare_formulations( &panel, dir_vec[i] );
    }   

    return pass;
}


int main( void )
{
    // Define local auxiliar variables
    bool pass = false;

    // Run tests
    pass = test_triangle_0( );
    if ( !pass )
    {
        std::cerr << "test_triangle_0 failed." << std::endl;
        return 1;
    }

    pass = test_quadrilateral_0( );
    if ( !pass )
    {
        std::cerr << "test_quadrilateral_0 failed." << std::endl;
        return 1;
    }

    return 0;

}