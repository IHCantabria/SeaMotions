
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
#include "../../src/math/integration.hpp"
#include "../../src/math/math_interface.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/math/topology.hpp"
#include "../../src/green/pulsating_fin_depth.hpp"
#include "../../src/waves.hpp"


void sub_test_2( void )
{
    // Define source input parameters
    cusfloat T = 0.0001;
    cusfloat h = 20.0;
    cusfloat x = 0.0;
    cusfloat y = 0.0;
    cusfloat z = -h/1000.0;

    cusfloat w0 = 2*PI/T;
    WaveDispersionData wd = WaveDispersionData(w0, 30, h, 9.81);
    wd.calculate_john_terms();
    cusfloat H = wd.nu*h;

    IntegralsDb idb = IntegralsDb();
    build_integrals_db(idb);
    idb.fold_h(H);

    // Define panel
    PanelGeom panel;

    panel.num_nodes = 4;

    panel.num_nodes = 4;

    panel.x[0] = -1.0;
    panel.x[1] =  1.0;
    panel.x[2] =  1.0;
    panel.x[3] = -1.0;

    panel.y[0] = -1.0;
    panel.y[1] = -1.0;
    panel.y[2] =  1.0;
    panel.y[3] =  1.0;

    panel.z[0] = -h/999.99;
    panel.z[1] = -h/999.99;
    panel.z[2] = -h/999.99;
    panel.z[3] = -h/999.99;

    panel.calculate_properties( );

    // Define target function
    auto target_fcn = [
                            x,
                            y,
                            z,
                            h,
                            &wd,
                            &idb
                        ](
                            cusfloat xi,
                            cusfloat eta,
                            cusfloat zeta
                            )-> cuscomplex 
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
                                                    zeta,
                                                    h,
                                                    wd,
                                                    idb
                                                    );
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
                            };

    // Integrate using adaptative quadrature
    cuscomplex sol = adaptive_quadrature_panel(
                                                    &panel,
                                                    target_fcn,
                                                    1e-3,
                                                    5
                                                );

    std::cout << "sol: " << sol << std::endl;

}


int main(void)
{
    sub_test_2( );

    // // Set test configuration
    // const int N = 20;
    // cusfloat int_values[N-1];

    // // Define source input parameters
    // cusfloat T = 0.01;
    // cusfloat h = 100.0;
    // cusfloat x = 0.0;
    // cusfloat y = 0.0;
    // cusfloat z = -h/100.0;

    // cusfloat w0 = 2*PI/T;
    // WaveDispersionData wd = WaveDispersionData(w0, 30, h, 9.81);
    // wd.calculate_john_terms();
    // cusfloat H = wd.nu*h;

    // IntegralsDb idb = IntegralsDb();
    // build_integrals_db(idb);
    // idb.fold_h(H);

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

    // // Perform test
    // auto f_def = [&wd]()->cusfloat {return john_series(1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, wd).real();};
    // for (int i=1; i<N; i++)
    // {
    //     int_values[i-1] = quadrature_unit_2d(
    //                                         [
    //                                             x,
    //                                             y,
    //                                             z,
    //                                             h,
    //                                             &wd,
    //                                             &idb
    //                                         ](
    //                                             cusfloat xi,
    //                                             cusfloat eta
    //                                             )-> cusfloat 
    //                                             {
    //                                                 // return john_series(
    //                                                 //                 x,
    //                                                 //                 y,
    //                                                 //                 z,
    //                                                 //                 xi,
    //                                                 //                 eta,
    //                                                 //                 -0.1,
    //                                                 //                 h,
    //                                                 //                 wd
    //                                                 //                 ).real();
    //                                                 return G_integral(
    //                                                                     x,
    //                                                                     y,
    //                                                                     z,
    //                                                                     xi,
    //                                                                     eta,
    //                                                                     -h/99,
    //                                                                     h,
    //                                                                     wd,
    //                                                                     idb
    //                                                                     ).real();
    //                                                 // return wave_term_fin_depth(
    //                                                 //                             x,
    //                                                 //                             y,
    //                                                 //                             z,
    //                                                 //                             xi,
    //                                                 //                             eta,
    //                                                 //                             -h/3,
    //                                                 //                             h,
    //                                                 //                             wd,
    //                                                 //                             idb
    //                                                 //                             ).real();
    //                                                 // return std::pow(xi, 2.0)+std::pow(eta, 2.0);
    //                                             }, 
    //                                             i
    //                                         );
    // }

    // // Print-out results
    // for (int i=1; i<N; i++)
    // {
    //     std::cout.precision(6);
    //     std::cout << "Gauss-Order: " << i << " - Integral Value: ";
    //     std::cout << std::scientific << int_values[i-1] << std::endl;
    // }

    return 0;
}