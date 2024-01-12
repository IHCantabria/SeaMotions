
// Include general usage libraries
#include <iostream>

// Include general usage scientific libraries
#include <cmath>
#include <complex>

// Inlude local modules
#include "waves_common.hpp"

#include "../config.hpp"
#include "../math/math_tools.hpp"


cusfloat dispersion_imag_zero(cusfloat nu, cusfloat h, cusfloat k)
{
    return nu+k*std::tan(k*h);
}


cusfloat dispersion_real_zero(cusfloat nu, cusfloat h, cusfloat k)
{
    return nu-k*std::tanh(k*h);
}


cusfloat dispersion_real_zero_dk(cusfloat h, cusfloat k)
{
    return -(std::tanh(k*h)+k*h*(1-pow2s(tanh(k*h))));
}


cusfloat k2w( cusfloat k, cusfloat h, cusfloat g )
{
    return std::sqrt( g * k * std::tanh( k * h ) );
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

    // Call bisection algorith to find zero of the dispersion
    // relation for this depth
    int info = -1;
    cusfloat k = -1;
    newton_raphson(
        [ki, h](cusfloat k){return dispersion_real_zero(ki, h, k);},
        [h](cusfloat k){return dispersion_real_zero_dk(h, k);},
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
        std::string err_message = "w2k function could not find the real rool of the dispersion equation.";
        err_message += " See log file for more details.";
        throw std::runtime_error(err_message);
    }

    return k;
}


void w2ki(cusfloat w, cusfloat h, cusfloat g, int n, cusfloat* kn)
{
    /**
     * @brief Calculate imaginary wave numbers using the dispersion relation.
     * 
     * This function calculates the imaginary wave number which is the solution 
     * to the dispersion relation:
     *                 w^2/g=-k·tg(k·h)
     * 
     * \param w Angular frequency of the wave.
     * \param h Water depth
     * \param g Gravitational acceleration
     * \param n Number of imaginary solutions to calculate
     * \param kn Output channel to get the imaginary wave numbers
     */

    // Calculate inteval bounds based on finite and 
    // infinite water depth aproximations of the 
    // dispersion relation
    cusfloat nu = pow2s(w)/g;

    // Loop over the N first branches of the tangent function 
    // to find the N first imaginary solutions
    int info = -1;
    cusfloat ki = 0;
    for (int i=0; i<n; i++)
    {
        // Calculate solution at brank i+1
        bisection(
                [nu, h](cusfloat k){return dispersion_imag_zero(nu, h, k);},
                (((i+1)-0.5)*PI)/h+1e-12,
                (i+1)*PI/h,
                1e-10,
                1e-10,
                100,
                false,
                ki,
                info
                );
        
        // Check for correct convergence
        if (info != 0)
        {
            std::cerr << "Bisection method could not find the solution for the ";
            std::cerr << "dispersion equation with parameters:" << std::endl;
            std::cerr << "  -> T: " << 2*PI/w << std::endl;
            std::cerr << "  -> w: " << w << std::endl;
            std::cerr << "  -> nu: " << ki << std::endl;
            std::cerr << "  -> h: " << h << std::endl;
            std::cerr << "At branch : " << ((i+1)-0.5)*PI << " - " << (i+1)*PI << std::endl;
            std::cerr << std::endl;
            throw std::runtime_error("w2ki function had a convergence problem. Check log file for more details.");
        }

        // Storage the result
        kn[i] = ki;
    }
}


cuscomplex  wave_potential_fo_space( 
                                            cusfloat aw,
                                            cusfloat w,
                                            cusfloat k,
                                            cusfloat h,
                                            cusfloat g,
                                            cusfloat x,
                                            cusfloat y,
                                            cusfloat z,
                                            cusfloat mu
                                    )
{
    cuscomplex  exp_arg = cuscomplex( 0.0, k * ( x * std::cos( mu ) + y * std::sin( mu ) ) );
    cusfloat    phi_0   = aw * g * std::cosh( k * ( h + z ) ) / w / std::cosh( k * h );
    return cuscomplex( 0.0, -phi_0 ) * std::exp( exp_arg );
}


cuscomplex  wave_potential_fo_space_dx( 
                                            cusfloat aw,
                                            cusfloat w,
                                            cusfloat k,
                                            cusfloat h,
                                            cusfloat g,
                                            cusfloat x,
                                            cusfloat y,
                                            cusfloat z,
                                            cusfloat mu
                                        )
{
    cuscomplex  exp_arg = cuscomplex( 0.0, k * ( x * std::cos( mu ) + y * std::sin( mu ) ) );
    cusfloat    phi_0   = aw * g * std::cosh( k * ( h + z ) ) / w / std::cosh( k * h ) * k * std::cos( mu );
    return phi_0 * std::exp( exp_arg );
}


cuscomplex  wave_potential_fo_space_dy( 
                                            cusfloat aw,
                                            cusfloat w,
                                            cusfloat k,
                                            cusfloat h,
                                            cusfloat g,
                                            cusfloat x,
                                            cusfloat y,
                                            cusfloat z,
                                            cusfloat mu
                                        )
{
    cuscomplex  exp_arg = cuscomplex( 0.0, k * ( x * std::cos( mu ) + y * std::sin( mu ) ) );
    cusfloat    phi_0   = aw * g * std::cosh( k * ( h + z ) ) / w / std::cosh( k * h ) * k * std::sin( mu );
    return phi_0 * std::exp( exp_arg );
}


cuscomplex  wave_potential_fo_space_dz( 
                                            cusfloat aw,
                                            cusfloat w,
                                            cusfloat k,
                                            cusfloat h,
                                            cusfloat g,
                                            cusfloat x,
                                            cusfloat y,
                                            cusfloat z,
                                            cusfloat mu
                                        )
{
    cuscomplex  exp_arg = cuscomplex( 0.0, k * ( x * std::cos( mu ) + y * std::sin( mu ) ) );
    cusfloat    phi_0   = aw * g * k * std::sinh( k * ( h + z ) ) / w / std::cosh( k * h );
    return cuscomplex( 0.0, -phi_0 ) * std::exp( exp_arg );
}