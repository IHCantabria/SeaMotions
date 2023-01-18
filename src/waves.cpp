
// Include general usage libraries
#include <iostream>

// Include general usage scientific libraries
#include <cmath>

// Inlude local modules
#include "config.hpp"
#include "./math/math_tools.hpp"


cusfloat dispersion_real_zero(cusfloat nu, cusfloat h, cusfloat k)
{
    return nu-k*std::tanh(k*h);
}


cusfloat dispersion_dk_real_zero(cusfloat nu, cusfloat h, cusfloat k)
{
    return -(std::tanh(k*h)+k*h*(1-pow2s(tanh(k*h))));
}


cusfloat w2k(cusfloat w, cusfloat h, cusfloat g)
{
    /**
     * @brief Calculate wave number from the dispersion relation.
     * 
     * This function calculates the real wave number which is the solution 
     * to the dispersion relation:
     *                 w^2/g=k·tanh(k·h)
     * 
     * \param w Angular frequency of the wave.
     * \param h Water depth
     * \return Wave number
     */

    // Calculate inteval bounds based on finite and 
    // infinite water depth aproximations of the 
    // dispersion relation
    cusfloat ki = pow2s(w)/g;
    cusfloat kf = std::sqrt(ki/h);
    cusfloat km = (ki<kf)?kf:ki;

    // Call bisection algorith to find zero of the dispersion
    // relation for this depth
    int info = -1;
    cusfloat k = -1;
    newton_raphson(
        [ki, h](cusfloat k){return dispersion_real_zero(ki, h, k);},
        [ki, h](cusfloat k){return dispersion_dk_real_zero(ki, h, k);},
        ki,
        1e-10,
        1e-6,
        20,
        false,
        k,
        info
    );

    // Check the stauts of the iterative search
    if (info != 0)
    {
        std::cerr << "Bisection method could not find the solution for the ";
        std::cerr << "dispersion equation with parameters:" << std::endl;
        std::cerr << "  -> T: " << 2*PI/w << std::endl;
        std::cerr << "  -> w: " << w << std::endl;
        std::cerr << "  -> nu: " << ki << std::endl;
        std::cerr << "  -> h: " << h << std::endl;
        std::cerr << std::endl;
    }

    return k;
}