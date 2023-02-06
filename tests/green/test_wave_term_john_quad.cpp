
// Include general usage libraries
#include <cassert>
#include <functional>
#include <iostream>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/math/gauss.hpp"
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
    assert(order == max_order && "Quadrature 2D maximum order reached!");

    // Get Gauss points
    cusfloat gp_roots[max_order];
    cusfloat gp_weights[max_order];
    get_gauss_legendre(order, gp_weights, gp_roots);

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
    cusfloat T = 0.1;
    cusfloat h = 1.0;
    cusfloat x = 2.0;
    cusfloat y = 0.0;
    cusfloat z = 0.0;

    cusfloat w0 = 2*PI/T;
    WaveDispersionData wd = WaveDispersionData(w0, 30, h, 9.81);

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
                                                &wd
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
                                                }, 
                                            i
                                            );
    }

    return 0;
}