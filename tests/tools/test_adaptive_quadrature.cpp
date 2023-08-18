
// Include general usage libraries
#include <iostream>

// Include local modules
#include "../../src/containers/panel_geom.hpp"
#include "../../src/math/integration.hpp"
#include "../../src/math/math_tools.hpp"

// Include namespaces
using namespace std::literals::complex_literals;

// Define global module variables
cusfloat EPS = 1e-6;


void sub_test_square_1( void )
{
    // Define square panel
    PanelGeom panel;

    panel.num_nodes = 4;

    panel.x[0] = -8.0;
    panel.x[1] =  8.0;
    panel.x[2] =  8.0;
    panel.x[3] = -8.0;

    panel.y[0] = -8.0;
    panel.y[1] = -8.0;
    panel.y[2] =  8.0;
    panel.y[3] =  8.0;

    panel.z[0] =  0.0;
    panel.z[1] =  0.0;
    panel.z[2] =  0.0;
    panel.z[3] =  0.0;

    panel.calculate_properties( );

    // Define integration function
    auto fcn_def    = []( cusfloat X, cusfloat Y, cusfloat ) -> cuscomplex
                       {
                           cuscomplex sol = ( pow2s( X ) + pow2s( Y ) ) + 0.0i;
                           
                           return sol;
                       };
    
    // Peform integration using the quadrature method
    cuscomplex  sol = adaptive_quadrature_panel(
                                                    &panel,
                                                    fcn_def,
                                                    1e-6,
                                                    9
                                                );
    
    // Check solution
    cuscomplex ref_sol = 32768.0/3.0 + 0.0i;
    if ( !assert_complex_equality( sol, ref_sol, EPS ) )
    {
        std::cerr << std::endl;
        std::cerr << "Test test_adaptive_quadrature/sub_test_square_1 failed!" << std::endl;
        std::cerr << " - Reference Solution: " << ref_sol << " - Solution: " << sol << std::endl;
        std::runtime_error( "" );
    }
}


void sub_test_square_2( void )
{
    // Define square panel
    PanelGeom panel;

    panel.num_nodes = 4;

    panel.x[0] = -1.0;
    panel.x[1] =  1.0;
    panel.x[2] =  1.0;
    panel.x[3] = -1.0;

    panel.y[0] = -1.0;
    panel.y[1] = -1.0;
    panel.y[2] =  1.0;
    panel.y[3] =  1.0;

    panel.z[0] =  0.0;
    panel.z[1] =  0.0;
    panel.z[2] =  0.0;
    panel.z[3] =  0.0;

    panel.calculate_properties( );

    // Define integration function
    auto fcn_def    = []( cusfloat X, cusfloat Y, cusfloat ) -> cuscomplex
                       {
                           cuscomplex sol = pow2s( std::cos( 8*PI*X ) * std::sin( 8*PI*Y ) ) + 0.0i;
                           
                           return sol;
                       };
    
    // Peform integration using the quadrature method
    cuscomplex  sol     = adaptive_quadrature_panel(
                                                        &panel,
                                                        fcn_def,
                                                        1e-6,
                                                        9
                                                    );
    
    // Check solution
    cuscomplex ref_sol = 1.0 + 0.0i;
    if ( !assert_complex_equality( sol, ref_sol, EPS ) )
    {
        std::cerr << std::endl;
        std::cerr << "Test test_adaptive_quadrature/sub_test_square_2 failed!" << std::endl;
        std::cerr << " - Reference Solution: " << ref_sol << " - Solution: " << sol << std::endl;
        std::runtime_error( "" );
    }

}


void sub_test_triangle_1( void )
{
    // Define square panel
    PanelGeom panel;

    panel.num_nodes = 3;

    panel.x[0] =  0.0;
    panel.x[1] =  8.0;
    panel.x[2] = -8.0;

    panel.y[0] =  0.0;
    panel.y[1] =  0.0;
    panel.y[2] =  8.0;

    panel.z[0] =  0.0;
    panel.z[1] =  0.0;
    panel.z[2] =  0.0;
    panel.z[3] =  0.0;

    panel.calculate_properties( );

    // Define integration function
    auto fcn_def    = []( cusfloat X, cusfloat Y, cusfloat ) -> cuscomplex
                       {
                           cuscomplex sol = ( pow2s( X ) + pow2s( Y ) ) + 0.0i;
                           
                           return sol;
                       };
    
    // Peform integration using the quadrature method
    cuscomplex  sol = adaptive_quadrature_panel(
                                                    &panel,
                                                    fcn_def,
                                                    1e-6,
                                                    9
                                                );
    
    // Check solution
    cuscomplex ref_sol = 2048.0/3.0 + 0.0i;
    if ( !assert_complex_equality( sol, ref_sol, EPS ) )
    {
        std::cerr << std::endl;
        std::cerr << "Test test_adaptive_quadrature/sub_test_triangle_1 failed!" << std::endl;
        std::cerr << " - Reference Solution: " << ref_sol << " - Solution: " << sol << std::endl;
        std::runtime_error( "" );
    }

}


void sub_test_triangle_2( void )
{
    // Define square panel
    PanelGeom panel;

    panel.num_nodes = 3;

    panel.x[0] =  0.0;
    panel.x[1] =  1.0;
    panel.x[2] = -1.0;

    panel.y[0] =  0.0;
    panel.y[1] =  0.0;
    panel.y[2] =  1.0;

    panel.z[0] =  0.0;
    panel.z[1] =  0.0;
    panel.z[2] =  0.0;

    panel.calculate_properties( );

    // Define integration function
    auto fcn_def    = []( cusfloat X, cusfloat Y, cusfloat ) -> cuscomplex
                       {
                           cuscomplex sol = pow2s( std::cos( 8*PI*X ) * std::sin( 8*PI*Y ) ) + 0.0i;
                           
                           return sol;
                       };
    
    // Peform integration using the quadrature method
    cuscomplex  sol     = adaptive_quadrature_panel(
                                                        &panel,
                                                        fcn_def,
                                                        1e-6,
                                                        9
                                                    );
    
    // Check solution
    cuscomplex ref_sol = 1.0/8.0 + 0.0i;
    if ( !assert_complex_equality( sol, ref_sol, EPS ) )
    {
        std::cerr << std::endl;
        std::cerr << "Test test_adaptive_quadrature/sub_test_triangle_2 failed!" << std::endl;
        std::cerr << " - Reference Solution: " << ref_sol << " - Solution: " << sol << std::endl;
        std::runtime_error( "" );
    }

}


int main( void )
{
    // First test:
    //  - Element type: Quadrilateral
    //  - Target function: x^2 + y^2
    sub_test_square_1( );

    // Second test:
    //  - Element type: Quadrilateral
    //  - Target function: ( cos(8*pi*x) * sin(8*pi*y) )^2
    sub_test_square_2( );

    // Third test:
    //  - Element type: Triangle
    //  - Target function: x^2 + y^2
    sub_test_triangle_1( );

    // Fourth test:
    //  - Element type: Quadrilateral
    //  - Target function: ( cos(8*pi*x) * sin(8*pi*y) )^2
    sub_test_triangle_2( );

    return 0;
}