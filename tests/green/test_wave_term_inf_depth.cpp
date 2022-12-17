
// Include general usage libraries
#include <iomanip>
#include <iostream>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/green/pulsating.hpp"
#include "../../src/math_tools.hpp"
#include "../../src/special_math.hpp"


int main(void)
{
    cusfloat X = 10.0;
    cusfloat Y = 25.0;

    cusfloat f0 = wave_term_inf_depth_num(X, Y) + 1/sqrt(pow2s(X)+pow2s(Y));
    // cusfloat f0 = romberg_quadrature(
    //     [X, Y](cusfloat t)->cusfloat {return expint_inf_depth_num(X, Y, t);},
    //     0.0,
    //     Y,
    //     1e-12
    // );
    cusfloat f1 = wave_term_inf_depth(X, Y);

    std::cout << std::setprecision(15) << "f0: " << f0 << " - f1: " << f1;
    std::cout << " - Diff: " << std::abs(f1-f0) << std::endl;


    return 0;
}