
// Include general usage libraries
#include <iomanip>
#include <iostream>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/green/pulsating.hpp"
#include "../../src/math_tools.hpp"
#include "../../src/special_math.hpp"


void generate_domain(int N, cusfloat* X, cusfloat* Y, cusfloat x_min, cusfloat x_max,
    cusfloat y_min, cusfloat y_max)
{
    // Generate X coordinate points
    cusfloat dx = (x_max-x_min)/(N-1);
    for (int i=0; i<N; i++)
    {
        X[i] = x_min + i*dx;
    }

    // Generate Y coordinate points
    cusfloat dy = (y_max-y_min)/(N-1);
    for (int i=0; i<N; i++)
    {
        Y[i] = y_min + i*dy;
    }
}


bool launch_test(void)
{
    // Generate computational grid
    constexpr int N = 1000;
    cusfloat X[N], Y[N];
    generate_domain(N, X, Y, 0.001, 40.0, 0.001, 40.0);

    // Loop over grid points to check the accuracy of the model
    bool pass = true;
    cusfloat diff = 0.0;
    cusfloat f_ref = 0.0, f_mod = 0.0;
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            // Calculate reference value using Romberg integration method
            f_ref = wave_term_inf_depth_num(X[i], Y[j]) + 1/sqrt(pow2s(X[i])+pow2s(Y[j]));

            // Calculate model value
            f_mod = wave_term_inf_depth(X[i], Y[j]);

            // Evaluate difference
            diff = f_mod - f_ref;
            if (std::abs(diff) > 3e-6)
            {
                std::cerr << "X: " << X[i] << " - Y: " << Y[j] << " - f_ref: " << f_ref;
                std::cerr << " - f_mod: " << f_mod << " - diff: " << diff << std::endl;
                pass = false;
                break;
            }
        }
    }

    return pass;
}


int main(void)
{
    // Test modelling function over all XY plane
    int pass = launch_test();
    if (!pass)
    {
        std::cerr << "test_wave_term_inf_depth failed!" << std::endl;
        return 1;
    }

    return 0;
}