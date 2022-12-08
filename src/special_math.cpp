
// Include general usage libraries
#include <iostream>
#include <string>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "math_tools.hpp"
#include "special_math.hpp"


//////////////////////////////////////////////
////// Bessel functions aproximations ////////
////// for simple precision           ////////
//////////////////////////////////////////////
#ifdef SIMPLE_PREC



//////////////////////////////////////////////
////// Bessel functions aproximations ////////
////// for double precision           ////////
//////////////////////////////////////////////
#else

#endif


cusfloat rational_fraction_f0(cusfloat x)
{
    // Define local constants
    cusfloat a0 = 0.79788454;
    cusfloat a1 = 5.46272781;
    cusfloat a2 = 3.02562477;

    cusfloat b1 = 6.90899779;
    cusfloat b2 = 4.12217805;

    // Calculate f0 function
    cusfloat x2 = x*x;
    cusfloat x4 = x2*x2;
    cusfloat f0 = (a0+a1/x2+a2/x4)/(1+b1/x2+b2/x4);

    return f0;
}


cusfloat rational_fraction_f1(cusfloat x)
{
    // Define local constants
    cusfloat a0 = 0.79788459;
    cusfloat a1 = 4.76650390;
    cusfloat a2 = 2.58896576;

    cusfloat b1 = 5.78645312;
    cusfloat b2 = 2.35033517;

    // Calculate f0 function
    cusfloat x2 = x*x;
    cusfloat x4 = x2*x2;
    cusfloat f1 = (a0+a1/x2+a2/x4)/(1+b1/x2+b2/x4);

    return f1;
}


cusfloat rational_fraction_th0(cusfloat x)
{
    // Define local constants
    cusfloat a0 = -0.12499967;
    cusfloat a1 = -1.07437411;
    cusfloat a2 = -0.75853664;

    cusfloat b1 = 9.11511321;
    cusfloat b2 = 9.19906287;

    // Calculate f0 function
    cusfloat x2 = x*x;
    cusfloat x3 = x*x2;
    cusfloat x4 = x2*x2;
    cusfloat th0 = (a0+a1/x2+a2/x4)/(x+b1/x+b2/x3) + x - PI/4.0;

    return th0;
}


cusfloat rational_fraction_th1(cusfloat x)
{
    // Define local constants
    cusfloat a0 = 0.37499947;
    cusfloat a1 = 2.77870488;
    cusfloat a2 = 1.39381402;

    cusfloat b1 = 7.84700458;
    cusfloat b2 = 6.19124657;

    // Calculate f0 function
    cusfloat x2 = x*x;
    cusfloat x3 = x*x2;
    cusfloat x4 = x2*x2;
    cusfloat th1 = (a0+a1/x2+a2/x4)/(x+b1/x+b2/x3) + x - 3.0*PI/4.0;

    return th1;
}


cusfloat besselj0(cusfloat x)
{
    // Check if the input is a real positive number
    if (x < 0.0)
    {
        std::string err_message("Bessel function of first kind and first order not defined for x < 0.0.");
        std::cerr << err_message << std::endl;
        throw std::runtime_error(err_message);
    }

    // Declare local variables
    cusfloat sol = 0.0;

    // Calculate polynomial approximation for 0 <= x <= 3
    if (x <= 3.0)
    {
        // Calculate polynomial coordinate to avoid
        // repeated calculation of it
        cusfloat xi = x/3.0;

        // Calculate power coefficients to economize
        // the computational effort
        cusfloat p2 = std::pow(xi, 2.0);
        cusfloat pc = p2;

        // Calculate polynomial approximation
        sol += 0.99999990;
        sol -= 2.24999239*pc;
        pc *= p2;
        sol += 1.26553572*pc;
        pc *= p2;
        sol -= 0.31602189*pc;
        pc *= p2;
        sol += 0.04374224*pc;
        pc *= p2;
        sol -= 0.00331563*pc;
    }
    else
    {
        sol = 1/std::sqrt(x)*rational_fraction_f0(x)*std::cos(rational_fraction_th0(x));
    }

    return sol;
}


cusfloat besselj1(cusfloat x)
{
    // Check if the input is a real positive number
    if (x < 0.0)
    {
        std::string err_message("Bessel function of first kind and second order not defined for x < 0.0.");
        std::cerr << err_message << std::endl;
        throw std::runtime_error(err_message);
    }

    // Declare local variables
    cusfloat sol = 0.0;

    // Calculate polynomial approximation for 0 <= x <= 3
    if (x <= 3.0)
    {
        // Calculate polynomial coordinate to avoid
        // repeated calculation of it
        cusfloat xi = x/3.0;

        // Calculate power coefficients to economize
        // the computational effort
        cusfloat p2 = std::pow(xi, 2.0);
        cusfloat pc = p2;

        // Calculate polynomial approximation
        sol += 0.49999999;
        sol -= 0.56249945*pc;
        pc *= p2;
        sol += 0.21093101*pc;
        pc *= p2;
        sol -= 0.03952287*pc;
        pc *= p2;
        sol += 0.00439494*pc;
        pc *= p2;
        sol -= 0.00028397*pc;
        sol *= x;
    }
    else
    {
        sol = 1/std::sqrt(x)*rational_fraction_f1(x)*std::cos(rational_fraction_th1(x));
    }

    return sol;
}


cusfloat bessely0(cusfloat x)
{
    // Check if the input is a real positive number
    if (x < 0.0)
    {
        std::string err_message("Bessel function of second kind and first order not defined for x < 0.0.");
        std::cerr << err_message << std::endl;
        throw std::runtime_error(err_message);
    }

    // Declare local variables
    cusfloat sol = 0.0;

    // Calculate polynomial approximation for 0 <= x <= 3
    if (x <= 3.0)
    {
        // Calculate polynomial coordinate to avoid
        // repeated calculation of it
        cusfloat xi = x/3.0;

        // Calculate power coefficients to economize
        // the computational effort
        cusfloat p2 = std::pow(xi, 2.0);
        cusfloat pc = p2;

        // Calculate polynomial approximation
        sol += (2.0/PI)*std::log(x/2.0)*besselj0(x);
        sol += 0.36746703;
        sol += 0.60558498*pc;
        pc *= p2;
        sol -= 0.74340225*pc;
        pc *= p2;
        sol += 0.25256673*pc;
        pc *= p2;
        sol -= 0.04177345*pc;
        pc *= p2;
        sol += 0.00353354*pc;
    }
    else
    {
        sol = 1/std::sqrt(x)*rational_fraction_f0(x)*std::sin(rational_fraction_th0(x));
    }

    return sol;
}


cusfloat bessely1(cusfloat x)
{
    // Check if the input is a real positive number
    if (x < 0.0)
    {
        std::string err_message("Bessel function of second kind and second order not defined for x < 0.0.");
        std::cerr << err_message << std::endl;
        throw std::runtime_error(err_message);
    }

    // Declare local variables
    cusfloat sol = 0.0;

    // Calculate polynomial approximation for 0 <= x <= 3
    if (x <= 3.0)
    {
        // Calculate polynomial coordinate to avoid
        // repeated calculation of it
        cusfloat xi = x/3.0;

        // Calculate power coefficients to economize
        // the computational effort
        cusfloat p2 = std::pow(xi, 2.0);
        cusfloat pc = xi;

        // Calculate polynomial approximation
        sol += (2.0/PI)*(std::log(x/2.0)*besselj1(x)-1/x);
        sol += 0.07373571*pc;
        pc *= p2;
        sol += 0.72276433*pc;
        pc *= p2;
        sol -= 0.43885620*pc;
        pc *= p2;
        sol += 0.10418264*pc;
        pc *= p2;
        sol -= 0.01340825*pc;
        pc *= p2;
        sol += 0.00094249*pc;
    }
    else
    {
        sol = 1/std::sqrt(x)*rational_fraction_f1(x)*std::sin(rational_fraction_th1(x));
    }

    return sol;
}