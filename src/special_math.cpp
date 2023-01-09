
// Include general usage libraries
#include <iostream>
#include <string>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "math_tools.hpp"
#include "special_math.hpp"


cusfloat besseli0(cusfloat x)
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


cusfloat besseli1(cusfloat x)
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
        sol = 1/std::sqrt(x)*polynomial_f0(x)*std::cos(polynomial_th0(x));
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
        sol = 1/std::sqrt(x)*polynomial_f1(x)*std::cos(polynomial_th1(x));
    }

    return sol;
}


cusfloat besselk0(cusfloat x)
{
    // Check if the input is a real positive number
    if (!(x > 0.0))
    {
        std::string err_message("Modified Bessel function of second kind and first order not defined for x > 0.0.");
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


cusfloat besselk1(cusfloat x)
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
        sol = 1/std::sqrt(x)*polynomial_f0(x)*std::sin(polynomial_th0(x));
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
        sol = 1/std::sqrt(x)*polynomial_f1(x)*std::sin(polynomial_th1(x));
    }

    return sol;
}


cusfloat chebyshev_poly_raw(int n, cusfloat x)
{
    // Include name
    using namespace std;

    // Look for the Lengendre polynomial
    cusfloat tn = 0.0;
    cusfloat x2 = pow2s(x);
    cusfloat xi = 0.0;
    switch (n)
    {
        case 0:
        return 1.0;

        case 1:
        return x;

        case 2:
        tn = -1.0;
        xi = x2;
        tn += 2.0*xi;
        return tn;

        case 3:
        xi = x;
        tn = -3.0*xi;
        xi *= x2;
        tn += 4.0*xi;
        return tn;

        case 4:
        tn = 1.0;
        xi = x2;
        tn -= 8.0*xi;
        xi *= x2;
        tn += 8.0*xi;
        return tn;

        case 5:
        xi = x;
        tn = 5.0*xi;
        xi *= x2;
        tn -= 20.0*xi;
        xi *= x2;
        tn += 16.0*xi;
        return tn;

        case 6:
        tn = -1.0;
        xi = x2;
        tn += 18.0*xi;
        xi *= x2;
        tn -= 48.0*xi;
        xi *= x2;
        tn += 32.0*xi;
        return tn;

        case 7:
        xi = x;
        tn = -7.0*xi;
        xi *= x2;
        tn += 56.0*xi;
        xi *= x2;
        tn -= 112.0*xi;
        xi *= x2;
        tn += 64.0*xi;
        return tn;

        case 8:
        tn = 1.0;
        xi = x2;
        tn -= 32.0*xi;
        xi *= x2;
        tn += 160.0*xi;
        xi *= x2;
        tn -= 256.0*xi;
        xi *= x2;
        tn += 128.0*xi;
        return tn;

        case 9:
        xi = x;
        tn = 9.0*xi;
        xi *= x2;
        tn -= 120.0*xi;
        xi *= x2;
        tn += 432.0*xi;
        xi *= x2;
        tn -= 576.0*xi;
        xi *= x2;
        tn += 256.0*xi;
        return tn;

        case 10:
        tn = -1.0;
        xi = x2;
        tn += 50.0*xi;
        xi *= x2;
        tn -= 400.0*xi;
        xi *= x2;
        tn += 1120.0*xi;
        xi *= x2;
        tn -= 1280.0*xi;
        xi *= x2;
        tn += 512.0*xi;
        return tn;

        case 11:
        xi = x;
        tn = -11.0*xi;
        xi *= x2;
        tn += 220.0*xi;
        xi *= x2;
        tn -= 1232.0*xi;
        xi *= x2;
        tn += 2816.0*xi;
        xi *= x2;
        tn -= 2816.0*xi;
        xi *= x2;
        tn += 1024.0*xi;
        return tn;

        case 12:
        tn = 1.0;
        xi = x2;
        tn -= 72.0*xi;
        xi *= x2;
        tn += 840.0*xi;
        xi *= x2;
        tn -= 3584.0*xi;
        xi *= x2;
        tn += 6912.0*xi;
        xi *= x2;
        tn -= 6144.0*xi;
        xi *= x2;
        tn += 2048.0*xi;
        return tn;

        case 13:
        xi = x;
        tn = 13.0*xi;
        xi *= x2;
        tn -= 364.0*xi;
        xi *= x2;
        tn += 2912.0*xi;
        xi *= x2;
        tn -= 9984.0*xi;
        xi *= x2;
        tn += 16640.0*xi;
        xi *= x2;
        tn -= 13312.0*xi;
        xi *= x2;
        tn += 4096.0*xi;
        return tn;

        case 14:
        tn = -1.0;
        xi = x2;
        tn += 98.0*xi;
        xi *= x2;
        tn -= 1568.0*xi;
        xi *= x2;
        tn += 9408.0*xi;
        xi *= x2;
        tn -= 26880.0*xi;
        xi *= x2;
        tn += 39424.0*xi;
        xi *= x2;
        tn -= 28672.0*xi;
        xi *= x2;
        tn += 8192.0*xi;
        return tn;

        case 15:
        xi = x;
        tn = -15.0*xi;
        xi *= x2;
        tn += 560.0*xi;
        xi *= x2;
        tn -= 6048.0*xi;
        xi *= x2;
        tn += 28800.0*xi;
        xi *= x2;
        tn -= 70400.0*xi;
        xi *= x2;
        tn += 92160.0*xi;
        xi *= x2;
        tn -= 61440.0*xi;
        xi *= x2;
        tn += 16384.0*xi;
        return tn;

        case 16:
        tn = 1.0;
        xi = x2;
        tn -= 128.0*xi;
        xi *= x2;
        tn += 2688.0*xi;
        xi *= x2;
        tn -= 21504.0*xi;
        xi *= x2;
        tn += 84480.0*xi;
        xi *= x2;
        tn -= 180224.0*xi;
        xi *= x2;
        tn += 212992.0*xi;
        xi *= x2;
        tn -= 131072.0*xi;
        xi *= x2;
        tn += 32768.0*xi;
        return tn;

        case 17:
        xi = x;
        tn = 17.0*xi;
        xi *= x2;
        tn -= 816.0*xi;
        xi *= x2;
        tn += 11424.0*xi;
        xi *= x2;
        tn -= 71808.0*xi;
        xi *= x2;
        tn += 239360.0*xi;
        xi *= x2;
        tn -= 452608.0*xi;
        xi *= x2;
        tn += 487424.0*xi;
        xi *= x2;
        tn -= 278528.0*xi;
        xi *= x2;
        tn += 65536.0*xi;
        return tn;

        case 18:
        tn = -1.0;
        xi = x2;
        tn += 162.0*xi;
        xi *= x2;
        tn -= 4320.0*xi;
        xi *= x2;
        tn += 44352.0*xi;
        xi *= x2;
        tn -= 228096.0*xi;
        xi *= x2;
        tn += 658944.0*xi;
        xi *= x2;
        tn -= 1118208.0*xi;
        xi *= x2;
        tn += 1105920.0*xi;
        xi *= x2;
        tn -= 589824.0*xi;
        xi *= x2;
        tn += 131072.0*xi;
        return tn;

        case 19:
        xi = x;
        tn = -19.0*xi;
        xi *= x2;
        tn += 1140.0*xi;
        xi *= x2;
        tn -= 20064.0*xi;
        xi *= x2;
        tn += 160512.0*xi;
        xi *= x2;
        tn -= 695552.0*xi;
        xi *= x2;
        tn += 1770496.0*xi;
        xi *= x2;
        tn -= 2723840.0*xi;
        xi *= x2;
        tn += 2490368.0*xi;
        xi *= x2;
        tn -= 1245184.0*xi;
        xi *= x2;
        tn += 262144.0*xi;
        return tn;

        case 20:
        tn = 1.0;
        xi = x2;
        tn -= 200.0*xi;
        xi *= x2;
        tn += 6600.0*xi;
        xi *= x2;
        tn -= 84480.0*xi;
        xi *= x2;
        tn += 549120.0*xi;
        xi *= x2;
        tn -= 2050048.0*xi;
        xi *= x2;
        tn += 4659200.0*xi;
        xi *= x2;
        tn -= 6553600.0*xi;
        xi *= x2;
        tn += 5570560.0*xi;
        xi *= x2;
        tn -= 2621440.0*xi;
        xi *= x2;
        tn += 524288.0*xi;
        return tn;

        default:
        std::string err_message("Raw Chebyshev polynomial defined for a degree en between 0 and 10.");
        std::cerr << err_message << std::endl;
        throw std::runtime_error(err_message);
    }
}


cusfloat chebyshev_poly_der_raw(int n, cusfloat x)
{
    // Include name
    using namespace std;

    // Look for the Lengendre polynomial
    cusfloat tn = 0.0;
    cusfloat x2 = pow2s(x);
    cusfloat xi = 0.0;
    switch (n)
    {
        case 0:
        return 0.0;

        case 1:
        return 1.0;

        case 2:
        return 4*x;

        case 3:
        return 12*pow(x, 2)-3;

        case 4:
        xi = x;
        tn = -16.0*xi;
        xi *= x2;
        tn += 32.0*xi;
        return tn;

        case 5:
        tn = 5.0;
        xi = x2;
        tn -= 60.0*xi;
        xi *= x2;
        tn += 80.0*xi;
        return tn;

        case 6:
        xi = x;
        tn = 36.0*xi;
        xi *= x2;
        tn -= 192.0*xi;
        xi *= x2;
        tn += 192.0*xi;
        return tn;

        case 7:
        tn = -7.0;
        xi = x2;
        tn += 168.0*xi;
        xi *= x2;
        tn -= 560.0*xi;
        xi *= x2;
        tn += 448.0*xi;
        return tn;

        case 8:
        xi = x;
        tn = -64.0*xi;
        xi *= x2;
        tn += 640.0*xi;
        xi *= x2;
        tn -= 1536.0*xi;
        xi *= x2;
        tn += 1024.0*xi;
        return tn;

        case 9:
        tn = 9.0;
        xi = x2;
        tn -= 360.0*xi;
        xi *= x2;
        tn += 2160.0*xi;
        xi *= x2;
        tn -= 4032.0*xi;
        xi *= x2;
        tn += 2304.0*xi;
        return tn;

        case 10:
        xi = x;
        tn = 100.0*xi;
        xi *= x2;
        tn -= 1600.0*xi;
        xi *= x2;
        tn += 6720.0*xi;
        xi *= x2;
        tn -= 10240.0*xi;
        xi *= x2;
        tn += 5120.0*xi;
        return tn;

        case 11:
        tn = -11.0;
        xi = x2;
        tn += 660.0*xi;
        xi *= x2;
        tn -= 6160.0*xi;
        xi *= x2;
        tn += 19712.0*xi;
        xi *= x2;
        tn -= 25344.0*xi;
        xi *= x2;
        tn += 11264.0*xi;
        return tn;

        case 12:
        xi = x;
        tn = -144.0*xi;
        xi *= x2;
        tn += 3360.0*xi;
        xi *= x2;
        tn -= 21504.0*xi;
        xi *= x2;
        tn += 55296.0*xi;
        xi *= x2;
        tn -= 61440.0*xi;
        xi *= x2;
        tn += 24576.0*xi;
        return tn;

        case 13:
        tn = 13.0;
        xi = x2;
        tn -= 1092.0*xi;
        xi *= x2;
        tn += 14560.0*xi;
        xi *= x2;
        tn -= 69888.0*xi;
        xi *= x2;
        tn += 149760.0*xi;
        xi *= x2;
        tn -= 146432.0*xi;
        xi *= x2;
        tn += 53248.0*xi;
        return tn;

        case 14:
        xi = x;
        tn = 196.0*xi;
        xi *= x2;
        tn -= 6272.0*xi;
        xi *= x2;
        tn += 56448.0*xi;
        xi *= x2;
        tn -= 215040.0*xi;
        xi *= x2;
        tn += 394240.0*xi;
        xi *= x2;
        tn -= 344064.0*xi;
        xi *= x2;
        tn += 114688.0*xi;
        return tn;

        case 15:
        tn = -15.0;
        xi = x2;
        tn += 1680.0*xi;
        xi *= x2;
        tn -= 30240.0*xi;
        xi *= x2;
        tn += 201600.0*xi;
        xi *= x2;
        tn -= 633600.0*xi;
        xi *= x2;
        tn += 1013760.0*xi;
        xi *= x2;
        tn -= 798720.0*xi;
        xi *= x2;
        tn += 245760.0*xi;
        return tn;

        case 16:
        xi = x;
        tn = -256.0*xi;
        xi *= x2;
        tn += 10752.0*xi;
        xi *= x2;
        tn -= 129024.0*xi;
        xi *= x2;
        tn += 675840.0*xi;
        xi *= x2;
        tn -= 1802240.0*xi;
        xi *= x2;
        tn += 2555904.0*xi;
        xi *= x2;
        tn -= 1835008.0*xi;
        xi *= x2;
        tn += 524288.0*xi;
        return tn;

        case 17:
        tn = 17.0;
        xi = x2;
        tn -= 2448.0*xi;
        xi *= x2;
        tn += 57120.0*xi;
        xi *= x2;
        tn -= 502656.0*xi;
        xi *= x2;
        tn += 2154240.0*xi;
        xi *= x2;
        tn -= 4978688.0*xi;
        xi *= x2;
        tn += 6336512.0*xi;
        xi *= x2;
        tn -= 4177920.0*xi;
        xi *= x2;
        tn += 1114112.0*xi;
        return tn;

        case 18:
        xi = x;
        tn = 324.0*xi;
        xi *= x2;
        tn -= 17280.0*xi;
        xi *= x2;
        tn += 266112.0*xi;
        xi *= x2;
        tn -= 1824768.0*xi;
        xi *= x2;
        tn += 6589440.0*xi;
        xi *= x2;
        tn -= 13418496.0*xi;
        xi *= x2;
        tn += 15482880.0*xi;
        xi *= x2;
        tn -= 9437184.0*xi;
        xi *= x2;
        tn += 2359296.0*xi;
        return tn;

        case 19:
        tn = -19.0;
        xi = x2;
        tn += 3420.0*xi;
        xi *= x2;
        tn -= 100320.0*xi;
        xi *= x2;
        tn += 1123584.0*xi;
        xi *= x2;
        tn -= 6259968.0*xi;
        xi *= x2;
        tn += 19475456.0*xi;
        xi *= x2;
        tn -= 35409920.0*xi;
        xi *= x2;
        tn += 37355520.0*xi;
        xi *= x2;
        tn -= 21168128.0*xi;
        xi *= x2;
        tn += 4980736.0*xi;
        return tn;

        case 20:
        xi = x;
        tn = -400.0*xi;
        xi *= x2;
        tn += 26400.0*xi;
        xi *= x2;
        tn -= 506880.0*xi;
        xi *= x2;
        tn += 4392960.0*xi;
        xi *= x2;
        tn -= 20500480.0*xi;
        xi *= x2;
        tn += 55910400.0*xi;
        xi *= x2;
        tn -= 91750400.0*xi;
        xi *= x2;
        tn += 89128960.0*xi;
        xi *= x2;
        tn -= 47185920.0*xi;
        xi *= x2;
        tn += 10485760.0*xi;
        return tn;

        default:
        std::string err_message("Raw Chebyshev polynomial defined for a degree en between 0 and 10.");
        std::cerr << err_message << std::endl;
        throw std::runtime_error(err_message);
    }
}


cusfloat expint_i(cusfloat x)
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


cusfloat legendre_poly_raw(int n, cusfloat x)
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


cusfloat legendre_poly_der_raw(int n, cusfloat x)
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


cusfloat psi_fun(int n)
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