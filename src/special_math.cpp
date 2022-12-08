
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


cusfloat polynomial_f0(cusfloat x)
{
    // Define local variables
    cusfloat xi = 3.0/x;
    cusfloat p2 = xi*xi;
    cusfloat pc = p2;
    cusfloat f0 = 0.0;

    // Calculate polynomial expansion
    f0 += 0.79788454;
    f0 -= 0.00553897*pc;
    pc *= p2;
    f0 += 0.00099336*pc;
    pc *= p2;
    f0 -= 0.00044346*pc;
    pc *= p2;
    f0 += 0.00020445*pc;
    pc *= p2;
    f0 -= 0.00004959*pc;

    return f0;
}


cusfloat polynomial_f1(cusfloat x)
{
    // Define local variables
    cusfloat xi = 3.0/x;
    cusfloat p2 = xi*xi;
    cusfloat pc = p2;
    cusfloat f0 = 0.0;

    // Calculate polynomial expansion
    f0 += 0.79788459;
    f0 += 0.01662008*pc;
    pc *= p2;
    f0 -= 0.00187002*pc;
    pc *= p2;
    f0 += 0.00068519*pc;
    pc *= p2;
    f0 -= 0.00029440*pc;
    pc *= p2;
    f0 += 0.00006952*pc;

    return f0;
}


cusfloat polynomial_th0(cusfloat x)
{
    // Define local variables
    cusfloat xi = 3.0/x;
    cusfloat p2 = xi*xi;
    cusfloat pc = xi;
    cusfloat f0 = 0.0;

    // Calculate polynomial expansion
    f0 += x - PI/4.0;
    f0 -= -0.04166592*pc;
    pc *= p2;
    f0 += 0.00239399*pc;
    pc *= p2;
    f0 -= 0.00073984*pc;
    pc *= p2;
    f0 += 0.00031099*pc;
    pc *= p2;
    f0 -= 0.00007605*pc;

    return f0;
}


cusfloat polynomial_th1(cusfloat x)
{
    // Define local variables
    cusfloat xi = 3.0/x;
    cusfloat p2 = xi*xi;
    cusfloat pc = xi;
    cusfloat f0 = 0.0;

    // Calculate polynomial expansion
    f0 += x - 3.0*PI/4.0;
    f0 += 0.12499895*pc;
    pc *= p2;
    f0 -= 0.00605240*pc;
    pc *= p2;
    f0 += 0.00135825*pc;
    pc *= p2;
    f0 -= 0.00049616*pc;
    pc *= p2;
    f0 += 0.00011531*pc;

    return f0;
}


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


cusfloat rational_fraction_struve0(cusfloat x)
{
    // Define local constans
    cusfloat a0 = 0.99999906;
    cusfloat a1 = 4.77228920;
    cusfloat a2 = 3.85542044;
    cusfloat a3 = 0.32303607;

    cusfloat b1 = 4.88331068;
    cusfloat b2 = 4.28957333;
    cusfloat b3 = 0.52120508;

    // Compute struve factor
    cusfloat xi = 3.0/x;
    cusfloat x2 = xi*xi;
    cusfloat x4 = x2*x2;
    cusfloat x6 = x2*x4;
    cusfloat c0 = 2*(a0+a1*x2+a2*x4+a3*x6);
    cusfloat c1 = PI*x*(1+b1*x2+b2*x4+b3*x6);
    cusfloat sf = c0/c1;

    return sf;
}


cusfloat rational_fraction_struve1(cusfloat x)
{
    // Define local constans
    cusfloat a0 = 1.00000004;
    cusfloat a1 = 3.92205313;
    cusfloat a2 = 2.64893033;
    cusfloat a3 = 0.27450895;

    cusfloat b1 = 3.81095112;
    cusfloat b2 = 2.26216956;
    cusfloat b3 = 0.10885141;

    // Compute struve factor
    cusfloat xi = 3.0/x;
    cusfloat x2 = xi*xi;
    cusfloat x4 = x2*x2;
    cusfloat x6 = x2*x4;
    cusfloat c0 = 2*(a0+a1*x2+a2*x4+a3*x6);
    cusfloat c1 = PI*(1+b1*x2+b2*x4+b3*x6);
    cusfloat sf = c0/c1;

    return sf;
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


cusfloat struve0(cusfloat x)
{
    // Check if the input is a real positive number
    if (x < 0.0)
    {
        std::string err_message("Struve function of first order not defined for x < 0.0.");
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
        sol += 1.909859164*pc;
        pc *= p2;
        sol -= 1.909855001*pc;
        pc *= p2;
        sol += 0.687514637*pc;
        pc *= p2;
        sol -= 0.126164557*pc;
        pc *= p2;
        sol += 0.013828813*pc;
        pc *= p2;
        sol -= 0.000876918*pc;
    }
    else
    {
        sol = bessely0(x) + rational_fraction_struve0(x);
    }

    return sol;
}


cusfloat struve1(cusfloat x)
{
    // Check if the input is a real positive number
    if (x < 0.0)
    {
        std::string err_message("Struve function of second order not defined for x < 0.0.");
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
        sol += 1.909859286*pc;
        pc *= p2;
        sol -= 1.145914713*pc;
        pc *= p2;
        sol += 0.294656958*pc;
        pc *= p2;
        sol -= 0.042070508*pc;
        pc *= p2;
        sol += 0.003785727*pc;
        pc *= p2;
        sol -= 0.000207183*pc;
    }
    else
    {
        sol = bessely1(x) + rational_fraction_struve1(x);
    }

    return sol;
}