
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
        // Add radius term
        wave_term = 1/std::sqrt(pow2s(X)+pow2s(Y));

        // Define local variables
        cusfloat cumsum_n = 0.0;
        cusfloat cumsum_m = 0.0;
        cusfloat fx = -pow2s(X)/4.0;
        cusfloat fy = 1/Y;
        cusfloat tx = fx;
        cusfloat ty = fy*fy;
        cusfloat fn = 2.0;
        cusfloat fm = 2.0;
        cusfloat c_expi = std::exp(-Y)*expint_i(Y);

        // Add n=0 and n=1 terms
        cumsum_m = (fy+ty);
        cumsum_n = -c_expi + tx*(cumsum_m-c_expi);
        tx *= fx;
        ty *= fy;

        // Add terms up to n=9
        for (int n=2; n<=9; n++)
        {
            // Calculate m loop value
            for (int m=2*(n-1)+1; m<=2*n; m++)
            {
                // Upate m loop
                cumsum_m += fm*ty;

                // Calculate m loop variables
                fm *= m;
                ty *= fy;
            }

            // Calculate n loop contribution
            cumsum_n += tx/pow2s(fn)*(cumsum_m-c_expi);

            // Update n loop variables
            tx *= fx;
            fn *= n+1;
        }
        wave_term += 2.0*cumsum_n;
    }
    else if ((X >= 3.7) && (Y <= 0.25*X))
    {
        cusfloat fs = 1.0;
        cusfloat fn = 1.0;
        cusfloat fx = 1/pow2s(X);
        cusfloat fy = pow2s(Y);
        cusfloat tx = 1.0;
        cusfloat ty = 1.0;
        cusfloat fmult = 1.0;
        cusfloat i2n = 1-std::exp(-Y);
        for (int i=0; i<=3; i++)
        {
            // Add expansion series term
            wave_term += fs*tx/fn*fmult*i2n;

            // Update loop variables
            fs *= -1.0;
            fn *= (i+1);
            tx *= fx;
            ty *= fy;
            fmult *= (2*(i+1)-1.0)/2.0;
            i2n = ty - 2.0*(i+1.0)*ty/Y + 2.0*(i+1.0)*(2.0*(i+1.0)-1.0)*i2n;
        }
        wave_term /= X;
        wave_term = 1/std::sqrt(pow2s(X)+pow2s(Y)) - PI*std::exp(-Y)*(struve0(X)+bessely0(X)) - 2.0*wave_term;
    }
    else if ((X>0)&&(X<=3.7)&&(Y>0.0)&&(Y<=2.0))
    {
        // Define radial distance to be used along the module
        cusfloat R = std::sqrt(pow2s(X)+pow2s(Y));

        // Calculate double series first row of coefficients
        constexpr int N = 4;
        cusfloat cmn_coeffs[2*N+1];
        cusfloat fn = 2.0;
        for (int n=1; n<(2*N+1); n++)
        {
            cmn_coeffs[n-1] = 1/fn/(n+1);
            fn *= (n+2);
        }

        // Calculate double series of coefficients
        cusfloat fx = pow2s(X);
        cusfloat tx = fx;
        cusfloat ty = 0.0;
        cusfloat nf = 0.0;
        cusfloat cumsum_cmn = 0.0;
        for (int m=0; m<=N; m++)
        {
            ty = Y;
            for (int n=1; n<=(2*N+1-2*m); n++)
            {
                // Update n as float in order to perform the math operations
                // correctly
                nf = static_cast<cusfloat>(n);

                // Add new term to the series
                cumsum_cmn += cmn_coeffs[n-1]*tx*ty;

                // Update series coefficients
                cmn_coeffs[n-1] = -((nf+2.0)/(nf+1.0))*cmn_coeffs[n+1];
                ty *= Y;
            }
            // Update series coefficients
            tx *= fx;
        }
        cumsum_cmn *= R;

        // Addd rest of the function
        wave_term = 1/R - 2.0*std::exp(-Y)*(besselj0(X)*std::log(Y/X+std::sqrt(1.0+pow2s(Y/X)))
                                            + PI*bessely0(X)/2.0
                                            + PI*struve0(X)*R/2.0/X
                                            + cumsum_cmn
                                            );
    }
    else if ((X>=3.7)&&(Y<=2.0))
    {
        std::cout << "New term" << std::endl;
        cusfloat fs = 1.0;
        cusfloat fn = 1.0;
        cusfloat fx = 1/pow2s(X);
        cusfloat fy = pow2s(Y);
        cusfloat tx = 1.0;
        cusfloat ty = 1.0;
        cusfloat fmult = 1.0;
        cusfloat i2n = 1-std::exp(-Y);
        for (int i=0; i<=7; i++)
        {
            // Add expansion series term
            wave_term += fs*tx/fn*fmult*i2n;

            // Update loop variables
            fs *= -1.0;
            fn *= (i+1);
            tx *= fx;
            ty *= fy;
            fmult *= (2*(i+1)-1.0)/2.0;
            i2n = ty - 2.0*(i+1.0)*ty/Y + 2.0*(i+1.0)*(2.0*(i+1.0)-1.0)*i2n;
        }
        wave_term /= X;
        wave_term = 1/std::sqrt(pow2s(X)+pow2s(Y)) - PI*std::exp(-Y)*(struve0(X)+bessely0(X)) - 2.0*wave_term;
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