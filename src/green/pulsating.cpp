
// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/special_math.hpp"
#include "../../src/math_tools.hpp"


#ifdef SIMPLE_PREC
cusfloat PRECISION_ROMBERG = 1e-7;
#else
cusfloat PRECISION_ROMBERG = 1e-12;
#endif


cusfloat expint_inf_depth_num(cusfloat X, cusfloat Y, cusfloat t)
{
    return std::exp(t-Y)/std::sqrt(pow2s(X)+pow2s(t));
}


cusfloat expint_inf_depth_num_dx(cusfloat X, cusfloat Y, cusfloat t)
{
    return -X*std::exp(t-Y)/pow3s(std::sqrt(pow2s(X)+pow2s(t)));
}


cusfloat expint_inf_depth_num_dxt(cusfloat X, cusfloat Y, cusfloat t)
{
    return t*std::exp(t-Y)/std::sqrt(pow2s(X)+pow2s(t));
}


cusfloat wave_term_inf_depth(cusfloat X, cusfloat Y)
{
    cusfloat wave_term = 0.0;
    if ((X>8.0) && (Y>=20.0))
    {
        // Calculate formulation parameters
        cusfloat R = std::sqrt(pow2s(X)+pow2s(Y));
        cusfloat sf = Y/R;
        cusfloat rfs = 1/R;
        cusfloat rf = rfs;
        cusfloat pn = 1.0;
        for (int i=0; i<=4; i++)
        {
            // Add new term to the series
            wave_term += pn*legendre_poly_raw(i, sf)*rf;

            // Update formulation parameters
            rf *= rfs;
            pn *= static_cast<cusfloat>(i+1);
        }

        wave_term = rfs - PI*std::exp(-Y)*(struve0(X)+bessely0(X))-2.0*wave_term;
    }
    else if (Y > 2*X)
    {
        
    }

    return wave_term;
}


cusfloat wave_term_inf_depth_num(cusfloat X, cusfloat Y)
{
    cusfloat expint = romberg_quadrature(
        [X, Y](cusfloat t)->cusfloat {return expint_inf_depth_num(X, Y, t);},
        0,
        Y,
        PRECISION_ROMBERG
    );
    return -2*expint - PI*std::exp(-Y)*(bessely0(X)+struve0(X));
}


cusfloat wave_term_inf_depth_num_dx(cusfloat X, cusfloat Y)
{
    cusfloat expint = romberg_quadrature(
        [X, Y](cusfloat t)->cusfloat {return expint_inf_depth_num_dxt(X, Y, t);},
        0,
        Y,
        PRECISION_ROMBERG
    );
    cusfloat sol = -2.0/X*expint
                + 2.0*Y/(X*std::sqrt(pow2s(X)+pow2s(Y)))
                + PI*std::exp(-Y)*(bessely1(X)+struve1(X)-2.0/PI);
    return sol;
}