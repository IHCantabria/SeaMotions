
// Include general usage libraries
#include <iostream>
#include <string>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "math_tools.hpp"
#include "special_math.hpp"


cusfloat    besseli0( cusfloat x )
{
    // Check if the input is a real positive number
    if (x < -3.75)
    {
        std::string err_message("Modified Bessel function of first kind and first order not defined for x < -3.75.");
        std::cerr << err_message << std::endl;
        throw std::runtime_error(err_message);
    }

    // Declare local variables
    cusfloat sol = 0.0;

    // Calculate polynomial approximation for 0 <= x <= 3
    if (x <= 3.75)
    {
        // Calculate polynomial coordinate to avoid
        // repeated calculation of it
        cusfloat xi = x/3.75;

        // Calculate power coefficients to economize
        // the computational effort
        cusfloat p2 = xi*xi;
        cusfloat pc = p2;

        // Calculate polynomial approximation
        sol += 1.0;
        sol += 3.5156229*pc;
        pc *= p2;
        sol += 3.0899424*pc;
        pc *= p2;
        sol += 1.2067492*pc;
        pc *= p2;
        sol += 0.2659732*pc;
        pc *= p2;
        sol += 0.0360768*pc;
        pc *= p2;
        sol += 0.0045813*pc;
    }
    else
    {
        // Calculate polynomial coordinate to avoid
        // repeated calculation of it
        cusfloat xi = x/3.75;

        // Calculate power coefficients to economize
        // the computational effort
        cusfloat p1 = 1/xi;
        cusfloat pc = p1;

        // Calculate polynomial approximation
        sol += 0.39894228;
        sol += 0.01328592*pc;
        pc *= p1;
        sol += 0.00225319*pc;
        pc *= p1;
        sol -= 0.00157565*pc;
        pc *= p1;
        sol += 0.00916281*pc;
        pc *= p1;
        sol -= 0.02057706*pc;
        pc *= p1;
        sol += 0.02635537*pc;
        pc *= p1;
        sol -= 0.01647633*pc;
        pc *= p1;
        sol += 0.00392377*pc;
        sol *= std::exp(x)/std::sqrt(x);
    }

    return sol;
}


cusfloat    besseli1( cusfloat x )
{
    // Check if the input is a real positive number
    if (x < -3.75)
    {
        std::string err_message("Modified Bessel function of first kind and first order not defined for x < -3.75.");
        std::cerr << err_message << std::endl;
        throw std::runtime_error(err_message);
    }

    // Declare local variables
    cusfloat sol = 0.0;

    // Calculate polynomial approximation for 0 <= x <= 3
    if (x <= 3.75)
    {
        // Calculate polynomial coordinate to avoid
        // repeated calculation of it
        cusfloat xi = x/3.75;

        // Calculate power coefficients to economize
        // the computational effort
        cusfloat p2 = std::pow(xi, 2.0);
        cusfloat pc = p2;

        // Calculate polynomial approximation
        sol += 0.5;
        sol += 0.87890594*pc;
        pc *= p2;
        sol += 0.51498869*pc;
        pc *= p2;
        sol += 0.15084934*pc;
        pc *= p2;
        sol += 0.02658733*pc;
        pc *= p2;
        sol += 0.00301532*pc;
        pc *= p2;
        sol += 0.00032411*pc;
        sol *= x;
    }
    else
    {
        // Calculate polynomial coordinate to avoid
        // repeated calculation of it
        cusfloat xi = x/3.75;

        // Calculate power coefficients to economize
        // the computational effort
        cusfloat p1 = 1/xi;
        cusfloat pc = p1;

        // Calculate polynomial approximation
        sol += 0.39894228;
        sol -= 0.03988024*pc;
        pc *= p1;
        sol -= 0.00362018*pc;
        pc *= p1;
        sol += 0.00163801*pc;
        pc *= p1;
        sol -= 0.01031555*pc;
        pc *= p1;
        sol += 0.02282967*pc;
        pc *= p1;
        sol -= 0.02895312*pc;
        pc *= p1;
        sol += 0.01787654*pc;
        pc *= p1;
        sol -= 0.00420059*pc;
        sol *= std::exp(x)/std::sqrt(x);
    }

    return sol;
}


cusfloat    besselj0( cusfloat x )
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
        sol += 0.999999999;
        sol -= 2.249999879*pc;
        pc *= p2;
        sol += 1.265623060*pc;
        pc *= p2;
        sol -= 0.316394552*pc;
        pc *= p2;
        sol += 0.044460948*pc;
        pc *= p2;
        sol -= 0.003954479*pc;
        pc *= p2;
        sol += 0.000212950*pc;
    }
    else
    {
        sol = 1/std::sqrt(x)*polynomial_f0(x)*std::cos(polynomial_th0(x));
    }

    return sol;
}


cusfloat    besselj1( cusfloat x )
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
        sol += 0.500000000;
        sol -= 0.562499992*pc;
        pc *= p2;
        sol += 0.210937377*pc;
        pc *= p2;
        sol -= 0.039550040*pc;
        pc *= p2;
        sol += 0.004447331*pc;
        pc *= p2;
        sol -= 0.000330547*pc;
        pc *= p2;
        sol += 0.000015525*pc;
        sol *= x;
    }
    else
    {
        sol = 1/std::sqrt(x)*polynomial_f1(x)*std::cos(polynomial_th1(x));
    }

    return sol;
}


cusfloat    besseljn_cos_int( 
                               cusfloat alpha, 
                               cusfloat beta, 
                               cusfloat nu
                           )
{
    cusfloat int_value = 0.0;
    if ( beta < alpha )
    {
        int_value = (
                        std::cos( nu * std::asin( beta / alpha ) )
                        /
                        std::sqrt( pow2s( alpha ) - pow2s( beta ) )
                    );
    }
    else if ( beta > alpha )
    {
        cusfloat    f1 = std::sqrt( pow2s( beta ) - pow2s( alpha ) );
        int_value = (
                        - std::pow( alpha, nu ) * std::sin( nu * PI / 2.0 )
                        /
                        ( f1 * std::pow( beta + f1, nu ) )
                    );
    }
	
	return int_value;
}


cusfloat    besseljn_cos_kernel(
                                    cusfloat alpha,
                                    cusfloat beta,
                                    cusfloat nu,
                                    cusfloat x
                                )
{
    return cos( alpha * x ) * std::cyl_bessel_j( nu, beta * x );
}


cuscomplex  besseljn_expi_int( 
                                cusfloat alpha,
                                cusfloat beta,
                                cusfloat nu
                            )
{
    return cuscomplex( 
                            besseljn_cos_int( alpha, beta, nu ),
                            besseljn_sin_int( alpha, beta, nu )
                        );
}


cusfloat    besseljn_sin_int( 
                                cusfloat alpha, 
                                cusfloat beta, 
                                cusfloat nu
                            )
{
    cusfloat int_value = 0.0;
    if ( beta < alpha )
    {
        int_value = (
                        std::sin( nu * std::asin( beta / alpha ) )
                        /
                        std::sqrt( pow2s( alpha ) - pow2s( beta ) )
                    );
    }
    else if ( beta > alpha )
    {
        cusfloat    f1 = std::sqrt( pow2s( beta ) - pow2s( alpha ) );
        int_value = (
                        std::pow( alpha, nu ) * std::cos( nu * PI / 2.0 )
                        /
                        ( f1 * std::pow( beta + f1, nu ) )
                    );
    }
	
	return int_value;
}


cusfloat    besseljn_sin_kernel(
                                    cusfloat alpha,
                                    cusfloat beta,
                                    cusfloat nu,
                                    cusfloat x
                                )
{
    return sin( alpha * x ) * std::cyl_bessel_j( nu, beta * x );
}


cusfloat    besselk0( cusfloat x )
{
    // Check if the input is a real positive number
    if (!(x > 0.0))
    {
        std::string err_message("Modified Bessel function of second kind and first order not defined for x <= 0.0.");
        std::cerr << err_message << std::endl;
        throw std::runtime_error(err_message);
    }

    // Declare local variables
    cusfloat sol = 0.0;

    // Calculate polynomial approximation for 0 <= x <= 3
    if (x <= 2.0)
    {
        // Calculate polynomial coordinate to avoid
        // repeated calculation of it
        cusfloat xi = x/2.0;

        // Calculate power coefficients to economize
        // the computational effort
        cusfloat p2 = std::pow(xi, 2.0);
        cusfloat pc = p2;

        // Calculate polynomial approximation
        sol += -std::log(xi)*besseli0(x);
        sol -= 0.57721566;
        sol += 0.42278420*pc;
        pc *= p2;
        sol += 0.23069756*pc;
        pc *= p2;
        sol += 0.03488590*pc;
        pc *= p2;
        sol += 0.00262698*pc;
        pc *= p2;
        sol += 0.00010750*pc;
        pc *= p2;
        sol += 0.00000740*pc;
    }
    else
    {
        // Calculate polynomial coordinate to avoid
        // repeated calculation of it
        cusfloat xi = 2.0/x;

        // Calculate power coefficients to economize
        // the computational effort
        cusfloat pc = xi;

        // Calculate polynomial approximation
        sol += 1.25331414;
        sol -= 0.07832358*pc;
        pc *= xi;
        sol += 0.02189568*pc;
        pc *= xi;
        sol -= 0.01062446*pc;
        pc *= xi;
        sol += 0.00587872*pc;
        pc *= xi;
        sol -= 0.00251540*pc;
        pc *= xi;
        sol += 0.00053208*pc;
        sol *= 1.0/(std::exp(x)*std::sqrt(x));
    }

    return sol;
}


cusfloat    besselk1( cusfloat x )
{
    // Check if the input is a real positive number
    if (!(x > 0.0))
    {
        std::string err_message("Modified Bessel function of second kind and second order not defined for x > 0.0.");
        std::cerr << err_message << std::endl;
        throw std::runtime_error(err_message);
    }

    // Declare local variables
    cusfloat sol = 0.0;

    // Calculate polynomial approximation for 0 <= x <= 3
    if (x <= 2.0)
    {
        // Calculate polynomial coordinate to avoid
        // repeated calculation of it
        cusfloat xi = x/2.0;

        // Calculate power coefficients to economize
        // the computational effort
        cusfloat p2 = std::pow(xi, 2.0);
        cusfloat pc = p2;

        // Calculate polynomial approximation
        sol += x*std::log(xi)*besseli1(x)+1.0;
        sol += 0.15443144*pc;
        pc *= p2;
        sol -= 0.67278579*pc;
        pc *= p2;
        sol -= 0.18156897*pc;
        pc *= p2;
        sol -= 0.01919402*pc;
        pc *= p2;
        sol -= 0.00110404*pc;
        pc *= p2;
        sol -= 0.00004686*pc;
        sol /= x;
    }
    else
    {
        // Calculate polynomial coordinate to avoid
        // repeated calculation of it
        cusfloat xi = 2.0/x;

        // Calculate power coefficients to economize
        // the computational effort
        cusfloat pc = xi;

        // Calculate polynomial approximation
        sol += 1.25331414;
        sol += 0.23498619*pc;
        pc *= xi;
        sol -= 0.03655620*pc;
        pc *= xi;
        sol += 0.01504268*pc;
        pc *= xi;
        sol -= 0.00780353*pc;
        pc *= xi;
        sol += 0.00325614*pc;
        pc *= xi;
        sol -= 0.00068245*pc;
        sol /= std::exp(x)*std::sqrt(x);
    }

    return sol;
}


cusfloat    bessely0( cusfloat x )
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
        sol += 0.367466907;
        sol += 0.605593797*pc;
        pc *= p2;
        sol -= 0.743505078*pc;
        pc *= p2;
        sol += 0.253005481*pc;
        pc *= p2;
        sol -= 0.042619616*pc;
        pc *= p2;
        sol += 0.004285691*pc;
        pc *= p2;
        sol -= 0.000250716*pc;
    }
    else
    {
        sol = 1/std::sqrt(x)*polynomial_f0(x)*std::sin(polynomial_th0(x));
    }

    return sol;
}


cusfloat    bessely1( cusfloat x )
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
        sol += 0.073735531*pc;
        pc *= p2;
        sol += 0.722769344*pc;
        pc *= p2;
        sol -= 0.438896337*pc;
        pc *= p2;
        sol += 0.104320251*pc;
        pc *= p2;
        sol -= 0.013637596*pc;
        pc *= p2;
        sol += 0.001125970*pc;
        pc *= p2;
        sol -= 0.000056455*pc;
    }
    else
    {
        sol = 1/std::sqrt(x)*polynomial_f1(x)*std::sin(polynomial_th1(x));
    }

    return sol;
}


cusfloat    ep_n( int n )
{
    return ( n > 0 ) ? 2.0 : 1.0;
}


cusfloat    expint_i( cusfloat x )
{
    cusfloat ei = EULERGAMMA + std::log(x);
    cusfloat ei_old = ei;
    cusfloat dei = 0.0;
    cusfloat xn = x;
    cusfloat fn = 1.0;
    cusfloat n = 1.0;
    while (true)
    {
        // Calculate new term in the series
        dei = xn/(n*fn);
        ei += dei;
        if (std::abs(ei-ei_old) < EPS_PRECISION)
        {
            break;
        }
        
        // Advance iteration
        ei_old = ei;
        n += 1.0;
        fn *= n;
        xn *= x;
    }

    return ei;
}


cusfloat    legendre_poly_raw( 
                                int n, 
                                cusfloat x
                            )
{
    // Include name
    using namespace std;

    // Look for the Lengendre polynomial
    switch (n)
    {
        case 0:
        return 1.0;

        case 1:
        return x;

        case 2:
        return 0.5*(3.0*pow(x, 2.0)-1.0);

        case 3:
        return 0.5*(5.0*pow(x, 3.0)-3.0*x);

        case 4:
        return 0.125*(35.0*pow(x, 4.0)-30.0*pow(x, 2.0)+3.0);

        case 5:
        return 0.125*(63.0*pow(x, 5.0)-70.0*pow(x, 3.0)+15.0*x);

        case 6:
        return 0.0625*(231.0*pow(x, 6.0)-315.0*pow(x, 4.0)+105*pow(x, 2.0)-5.0);

        case 7:
        return 0.0625*(429.0*pow(x, 7.0)-693.0*pow(x, 5.0)+315.0*pow(x, 3.0)-35.0*x);

        case 8:
        return (1.0/128.0)*(6435.0*pow(x, 8.0)-12012*pow(x, 6.0)+6930.0*pow(x, 4.0)-1260.0*pow(x, 2.0)+35.0);

        case 9:
        return (1.0/128.0)*(12155.0*pow(x, 9.0)-25740.0*pow(x, 7.0)+18018.0*pow(x, 5.0)-4620.0*pow(x, 3.0)+315.0*x);

        case 10:
        return (1.0/256.0)*(46189.0*pow(x, 10.0)-109395.0*pow(x, 8.0)+90090*pow(x, 6.0)-30030*pow(x, 4.0)+3465.0*pow(x, 2.0)-63.0);

        default:
        std::string err_message("Raw Legendre polynomial defined for a degree en between 0 and 10.");
        std::cerr << err_message << std::endl;
        throw std::runtime_error(err_message);
    }
}


cusfloat    legendre_poly_der_raw(
                                    int n, 
                                    cusfloat x
                                )
{
    // Include name
    using namespace std;

    // Look for the Lengendre polynomial
    switch (n)
    {
        case 0:
        return 0.0;

        case 1:
        return 1.0;

        case 2:
        return 3*x;

        case 3:
        return 0.5*(15*pow(x, 2)-3);

        case 4:
        return 0.5*(35*pow(x, 3)-15*x);

        case 5:
        return (315*pow(x, 4)-210*pow(x, 2)+15)/8;

        case 6:
        return (693*pow(x, 5)-630*pow(x, 3)+105*x)/8;

        case 7:
        return (3003*pow(x, 6)-3465*pow(x, 4)+945*pow(x, 2)-35)/16;

        case 8:
        return (51480*pow(x, 7)-72072*pow(x, 5)+27720*pow(x, 3)-2520*x)/128;

        case 9:
        return (109395*pow(x, 8)-180180*pow(x, 6)+90090*pow(x, 4)-13860*pow(x, 2)+315)/128;

        case 10:
        return (461890*pow(x, 9)-875160*pow(x, 7)+540540*pow(x, 5)-120120*pow(x, 3)+6930*x)/256;

        default:
        std::string err_message("Raw Legendre polynomial defined for a degree en between 0 and 10.");
        std::cerr << err_message << std::endl;
        throw std::runtime_error(err_message);
    }
}


cusfloat    polynomial_f0(
                            cusfloat x
                        )
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


cusfloat    polynomial_f1( cusfloat x )
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


cusfloat    polynomial_th0( cusfloat x )
{
    // Define local variables
    cusfloat xi = 3.0/x;
    cusfloat p2 = xi*xi;
    cusfloat pc = xi;
    cusfloat f0 = 0.0;

    // Calculate polynomial expansion
    f0 += x - PI/4.0;
    f0 -= 0.04166592*pc;
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


cusfloat    polynomial_th1( cusfloat x )
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


cusfloat    psi_fun( int n )
{
    // Check for function domain bounds
    if (n <= 0)
    {
        std::string err_message("Psi function only defined for natural positive numbers.");
        std::cerr << err_message << std::endl;
        throw std::runtime_error(err_message);
    }

    // Calculate psi function values
    cusfloat sol = -EULERGAMMA;
    for (int i=1; i<n; i++)
    {
        sol += 1/static_cast<cusfloat>(i);
    }

    return sol;
}


cusfloat    rational_fraction_f0( cusfloat x )
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


cusfloat    rational_fraction_f1( cusfloat x )
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


cusfloat rational_fraction_struve0( cusfloat x )
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


cusfloat rational_fraction_struve1( cusfloat x )
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


cusfloat rational_fraction_th0( cusfloat x )
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


cusfloat rational_fraction_th1( cusfloat x )
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


cusfloat struve0( cusfloat x )
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


cusfloat struve1( cusfloat x )
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