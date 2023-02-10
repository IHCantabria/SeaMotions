
// Include general usage libraries
#include <cassert>
#include <functional>
#include <iomanip>
#include <iostream>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/math/gauss.hpp"
#include "../../src/math/chebyshev.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/green/pulsating_fin_depth.hpp"
#include "../../src/waves.hpp"


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
    cusfloat T = 10.0;
    cusfloat h = 10.0;
    cusfloat x = 10.0;
    cusfloat y = 0.0;
    cusfloat z = 0.0;

    cusfloat w0 = 2*PI/T;
    WaveDispersionData wd = WaveDispersionData(w0, 30, h, 9.81);
    wd.calculate_john_terms();

    IntegralsDb idb = IntegralsDb();
    build_integrals_db(idb);

    cusfloat H = wd.nu*h;
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
    // auto f_def = [&wd]()->cusfloat {return john_series(1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, wd).real();};
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
                                                    return john_series(
                                                                    x,
                                                                    y,
                                                                    z,
                                                                    xi,
                                                                    eta,
                                                                    0.0,
                                                                    h,
                                                                    wd
                                                                    ).real();
                                                    // return wave_term_fin_depth(
                                                    //                             x,
                                                    //                             y,
                                                    //                             z,
                                                    //                             xi,
                                                    //                             eta,
                                                    //                             0.0,
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