
// Include general usage libraries
#include <cassert>
#include <functional>
#include <iomanip>
#include <iostream>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/containers/panel_geom.hpp"
#include "../../src/math/gauss.hpp"
#include "../../src/math/chebyshev.hpp"
#include "../../src/math/math_interface.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/math/topology.hpp"
#include "../../src/green/pulsating_fin_depth.hpp"
#include "../../src/waves.hpp"


void sub_test_0( void )
{
    // Define geometry
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

    // Calculate local values
    panel.calculate_properties( );


    // Get gauss points for integration
    const int   gp_np = 4;
    cusfloat    gp_roots[gp_np]; clear_vector( gp_np, gp_roots );
    cusfloat    gp_weights[gp_np]; clear_vector( gp_np, gp_weights );
    get_gauss_legendre( gp_np, gp_roots, gp_weights );

    cusfloat    fcn_val         = 0.0;
    cusfloat    int_value       = 0.0;
    cusfloat    gp_global[3]    = { 0.0, 0.0, 0.0 };
    for ( int i=0; i<gp_np; i++ )
    {
        for ( int j=0; j<gp_np; j++ )
        {
            // Get global coordinates for the gauss points
            panel.local_to_global( 
                                    gp_roots[i],
                                    gp_roots[j],
                                    gp_global
                                );
            
            // Calculate target function value
            fcn_val = pow2s( gp_global[0] ) + pow2s( gp_global[1] );

            // Calculate function integral function value
            int_value += gp_weights[i]*gp_weights[j]*fcn_val*jacobi_det_2d( 
                                                                                panel.num_nodes,
                                                                                panel.x,
                                                                                panel.y,
                                                                                gp_roots[i],
                                                                                gp_roots[j]
                                                                            );

        }
    }

    // Check integral value
    cusfloat ref_value = std::pow( 8, 5 ) / 3.0;
    if ( !assert_scalar_equality( int_value, ref_value, 1e-9 ) )
    {
        std::cerr << std::endl;
        std::cerr << "test test_wave_term_john_quad/sub_test_0 failed!" << std::endl;
        std::cerr << " - Integral value: " << int_value << " - Reference Value: " << ref_value << std::endl;
        throw std::runtime_error( "" );
    }
    
}


void sub_test_1( void )
{
    // Define geometry
    PanelGeom panel;
    
    panel.num_nodes = 3;

    panel.x[0] =  0.0;
    panel.x[1] =  1.0;
    panel.x[2] = -1.0;

    panel.y[0] =  0.0;
    panel.y[1] =  0.0;
    panel.y[2] =  1.0;

    // Calculate local values
    panel.calculate_properties( );


    // Get gauss points for integration
    const int   gp_np = 4;
    cusfloat    gp_roots[gp_np]; clear_vector( gp_np, gp_roots );
    cusfloat    gp_weights[gp_np]; clear_vector( gp_np, gp_weights );
    get_gauss_legendre( gp_np, gp_roots, gp_weights );

    cusfloat    fcn_val         = 0.0;
    cusfloat    int_value       = 0.0;
    cusfloat    gp_global[3]    = { 0.0, 0.0, 0.0 };
    for ( int i=0; i<gp_np; i++ )
    {
        for ( int j=0; j<gp_np; j++ )
        {
            // Get global coordinates for the gauss points
            panel.local_to_global( 
                                    gp_roots[i],
                                    gp_roots[j],
                                    gp_global
                                );
            
            // Calculate target function value
            fcn_val = pow2s( gp_global[0] ) + pow2s( gp_global[1] );

            // Calculate function integral function value
            int_value += gp_weights[i]*gp_weights[j]*fcn_val*jacobi_det_2d( 
                                                                                panel.num_nodes,
                                                                                panel.x,
                                                                                panel.y,
                                                                                gp_roots[i],
                                                                                gp_roots[j]
                                                                            );

        }
    }

    // Check integral value
    cusfloat ref_value = 1.0 / 6.0;
    if ( !assert_scalar_equality( int_value, ref_value, 1e-9 ) )
    {
        std::cerr << std::endl;
        std::cerr << "test test_wave_term_john_quad/sub_test_1 failed!" << std::endl;
        std::cerr << " - Integral value: " << int_value << " - Reference Value: " << ref_value << std::endl;
        throw std::runtime_error( "" );
    }
    
}


cusfloat quadrature_unit_2d(
                            std::function <cusfloat(cusfloat,cusfloat)> f_def, 
                            int order
                            )
{
    // Check the order requested is in between the predefined 
    // limits
    const int max_order = 20;
    assert(order <= max_order && "Quadrature 2D maximum order reached!");

    // Get Gauss points
    cusfloat gp_roots[max_order];
    cusfloat gp_weights[max_order];
    get_gauss_legendre(order, gp_roots, gp_weights);
    // get_gauss_chebyshev(order, gp_weights, gp_roots);

    // Integrate function
    cusfloat sol = 0.0;
    for (int i=0; i<order; i++)
    {
        for (int j=0; j<order; j++)
        {
            sol += (
                    gp_weights[i]
                    *
                    gp_weights[j]
                    *
                    f_def(gp_roots[i], gp_roots[j])
                    );
        }
    }

    return sol;
}


int main(void)
{
    // Set test configuration
    const int N = 20;
    cusfloat int_values[N-1];

    // Define source input parameters
    cusfloat T = 0.01;
    cusfloat h = 100.0;
    cusfloat x = 0.0;
    cusfloat y = 0.0;
    cusfloat z = -h/100.0;

    cusfloat w0 = 2*PI/T;
    WaveDispersionData wd = WaveDispersionData(w0, 30, h, 9.81);
    wd.calculate_john_terms();
    cusfloat H = wd.nu*h;

    IntegralsDb idb = IntegralsDb();
    build_integrals_db(idb);
    idb.fold_h(H);

    // Launch test 0
    sub_test_0( );

    // Launch test 1
    sub_test_1( );

    // std::cout << "H: " << H << std::endl;
    // std::cout << "L1: " << idb.l1->get_value_abh(1.0, 0.0, H) << std::endl;
    // std::cout << "L2: " << idb.l2->get_value_abh(H) << std::endl;
    // std::cout << "L3: " << idb.l3->get_value_abh(1.0, 0.0, H) << std::endl;
    // cusfloat t0 = chebyshev_poly(11, 1.0);
    // cusfloat t1 = chebyshev_poly(0, 0.0);
    // cusfloat t2 = chebyshev_poly(3, H);
    // std::cout << "t0: " << t0 << std::endl;
    // std::cout << "t1: " << t1 << std::endl;
    // std::cout << "t2: " << t2 << std::endl;

    // Perform test
    auto f_def = [&wd]()->cusfloat {return john_series(1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, wd).real();};
    for (int i=1; i<N; i++)
    {
        int_values[i-1] = quadrature_unit_2d(
                                            [
                                                x,
                                                y,
                                                z,
                                                h,
                                                &wd,
                                                &idb
                                            ](
                                                cusfloat xi,
                                                cusfloat eta
                                                )-> cusfloat 
                                                {
                                                    // return john_series(
                                                    //                 x,
                                                    //                 y,
                                                    //                 z,
                                                    //                 xi,
                                                    //                 eta,
                                                    //                 -0.1,
                                                    //                 h,
                                                    //                 wd
                                                    //                 ).real();
                                                    return G_integral(
                                                                        x,
                                                                        y,
                                                                        z,
                                                                        xi,
                                                                        eta,
                                                                        -h/99,
                                                                        h,
                                                                        wd,
                                                                        idb
                                                                        ).real();
                                                    // return wave_term_fin_depth(
                                                    //                             x,
                                                    //                             y,
                                                    //                             z,
                                                    //                             xi,
                                                    //                             eta,
                                                    //                             -h/3,
                                                    //                             h,
                                                    //                             wd,
                                                    //                             idb
                                                    //                             ).real();
                                                    // return std::pow(xi, 2.0)+std::pow(eta, 2.0);
                                                }, 
                                                i
                                            );
    }

    // Print-out results
    for (int i=1; i<N; i++)
    {
        std::cout.precision(6);
        std::cout << "Gauss-Order: " << i << " - Integral Value: ";
        std::cout << std::scientific << int_values[i-1] << std::endl;
    }

    return 0;
}