
// Include general usage libraries
#include <iostream>

// Include local modules
#include "../../src/containers/panel_geom.hpp"
#include "../../src/math/integration.hpp"
#include "../../src/math/math_tools.hpp"

// Include namespaces
using namespace std::literals::complex_literals;


void sub_test_square( void )
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
    cuscomplex sol = quadrature_panel(
                                        &panel,
                                        fcn_def,
                                        9
                                    );
    cuscomplex sol2 = adaptive_quadrature_panel(
                                                    &panel,
                                                    fcn_def,
                                                    1e-6,
                                                    9
                                                );
    
    std::cout << "sol: " << sol << std::endl;
    std::cout << "sol2: " << sol2 << std::endl;
}


int main( void )
{
    // Launch test for the quadrilateral element
    sub_test_square( );

    return 0;
}