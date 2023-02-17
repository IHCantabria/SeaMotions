
// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "chebyshev_inf_depth.hpp"
#include "../../src/math/chebyshev.hpp"
#include "../../src/config.hpp"
#include "../../src/math/special_math.hpp"
#include "../../src/math/math_tools.hpp"


#ifdef SIMPLE_PREC
cusfloat PRECISION_ROMBERG = 1e-7;
#else
cusfloat PRECISION_ROMBERG = 1e-12;
#endif

///////////////////////////////////////////
///// Declare module local functions //////
///////////////////////////////////////////
void domain_inf_fit(cusfloat x, cusfloat y, cusfloat &xl, cusfloat &yl, cusfloat &jac);
cusfloat eval_chebyshev_fit(const int num_cheby, const cusfloat cheby_coeffs[num_cheby], const int cheby_order_0[num_cheby],
    const int cheby_order_1[num_cheby], cusfloat x, cusfloat y);
cusfloat eval_chebyshev_fit_dx(const int num_cheby, const cusfloat cheby_coeffs[num_cheby], const int cheby_order_0[num_cheby],
    const int cheby_order_1[num_cheby], cusfloat x, cusfloat y);
void get_inf_domain_bounds(cusfloat x, cusfloat y, cusfloat &x0, cusfloat &x1, cusfloat &y0, cusfloat &y1);


///////////////////////////////////////////
//////// Define module functions //////////
///////////////////////////////////////////
void domain_inf_fit(cusfloat x, cusfloat y, cusfloat &xl, cusfloat &yl, cusfloat &jac)
{
    // Calculate x local range
    cusfloat x0 = 0.0, x1 = 0.0;
    cusfloat y0 = 0.0, y1 = 0.0;
    get_inf_domain_bounds(x, y, x0, x1, y0, y1);
    cusfloat dx = x1-x0;
    xl = 2.0*(x-x0)/dx-1.0;
    jac = 2/dx;

    // Calculate y local range
    cusfloat dy = y1 - y0;
    yl = 2.0*(y-y0)/dy-1.0;

}


cusfloat eval_chebyshev_fit(const int num_cheby, const cusfloat cheby_coeffs[num_cheby], const int cheby_order_0[num_cheby],
    const int cheby_order_1[num_cheby], cusfloat x, cusfloat y)
{
    cusfloat sol = 0.0;
    cusfloat xl = 0.0, yl = 0.0, jac = 0.0;
    for (int i=0; i<num_cheby; i++)
    {
        domain_inf_fit(x, y, xl, yl, jac);
        sol += cheby_coeffs[i]*chebyshev_poly_raw(cheby_order_0[i], xl)*chebyshev_poly_raw(cheby_order_1[i], yl);
    }

    return sol;
}


cusfloat eval_chebyshev_fit_dxndim(const int num_cheby, const cusfloat cheby_coeffs[num_cheby], const int cheby_order_0[num_cheby],
    const int cheby_order_1[num_cheby], cusfloat x, cusfloat y)
{
    cusfloat sol = 0.0;
    cusfloat xl = 0.0, yl = 0.0, jac = 0.0;
    cusfloat ci = 0.0;
    for (int i=0; i<num_cheby; i++)
    {
        domain_inf_fit(x, y, xl, yl, jac);
        // std::cout << "x: " << x << " - y: " << y << " - xl: " << xl << " - yl: " << yl << std::endl;
        ci = cheby_coeffs[i]*chebyshev_poly_der_raw(cheby_order_0[i], xl)*chebyshev_poly_raw(cheby_order_1[i], yl)*jac;
        // std::cout << "i" << i << " - ci: " << ci << std::endl;
        sol += ci;
    }

    return sol;
}


cusfloat expint_inf_depth_num(cusfloat X, cusfloat Y)
{
    cusfloat expint = romberg_quadrature(
        [X, Y](cusfloat t)->cusfloat {return std::exp(t-Y)/std::sqrt(pow2s(X)+pow2s(t));},
        0,
        Y,
        PRECISION_ROMBERG
    );
    return expint;
}


cusfloat expint_inf_depth_num_dxndim(cusfloat X, cusfloat Y)
{
    cusfloat expint = romberg_quadrature(
        [X, Y](cusfloat t)->cusfloat {return -X*std::exp(t-Y)/pow3s(std::sqrt(pow2s(X)+pow2s(t)));},
        0,
        Y,
        PRECISION_ROMBERG
    );
    return expint;
}


cusfloat expint_inf_depth_num_dxtndim(cusfloat X, cusfloat Y)
{
    cusfloat expint = romberg_quadrature(
        [X, Y](cusfloat t)->cusfloat {return t*std::exp(t-Y)/std::sqrt(pow2s(X)+pow2s(t));},
        0,
        Y,
        PRECISION_ROMBERG
    );
    return expint;
}


void get_inf_domain_bounds(cusfloat x, cusfloat y, cusfloat &x0, cusfloat &x1, cusfloat &y0, cusfloat &y1)
{
    if (y<=4.0)
    {
        // Define X boundaries
        if (x <= 3.0)
        {
            x0 = y/2.0;
            x1 = 3.0;
        }
        else
        {
            x0 = 3.0;
            x1 = 4.0*y;
        }

        // Define Y boundaries
        y0 = 2.0;
        y1 = 4.0;
    }
    else if (y<=8.0)
    {
        // Define X boundaries
        if (x <= 8.0)
        {
            x0 = y/2.0;
            x1 = 8.0;
        }
        else
        {
            x0 = 8.0;
            x1 = 4.0*y;
        }

        // Define Y boundaries
        y0 = 4.0001;
        y1 = 8.0;
    }
    else
    {
        x0 = y/2.0;
        x1 = 4.0*y;
        y0 = 8.00001;
        y1 = 20.0;
    }
}


cusfloat wave_term_inf_depth_series(cusfloat X, cusfloat Y)
{
    // std::cout << "X: " << X << " - Y: " << Y << std::endl;
    cusfloat wave_term = 0.0;
    if ((X>=8.0) && (Y>=20.0))
    {
        // std::cout << "Region 0" << std::endl;

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

        wave_term = - PI*std::exp(-Y)*(struve0(X)+bessely0(X))-2.0*wave_term;
    }
    else if (Y > 2*X)
    {
        // std::cout << "Region 1" << std::endl;
        // Add radius term
        wave_term = 0.0;

        if (Y>100)
        {
            Y = 70;
        }

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
        // std::cout << "Region 2" << std::endl;

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
        wave_term = - PI*std::exp(-Y)*(struve0(X)+bessely0(X)) - 2.0*wave_term;
    }
    else if ((X>0)&&(X<=3.7)&&(Y>0.0)&&(Y<=2.0))
    {
        // std::cout << "Region 3" << std::endl;

        // Define radial distance to be used along the module
        cusfloat R = std::sqrt(pow2s(X)+pow2s(Y));

        // Define factor variables
        cusfloat fx = pow2s(X);
        cusfloat tx = fx;
        cusfloat ty = 0.0;

        // Calculate double series first row of coefficients
        constexpr int N = 9;
        cusfloat cmn_coeffs[2*N+1];
        cusfloat fn = 2.0;
        cusfloat cumsum_cmn = 0.0;
        ty = Y;
        for (int n=1; n<=(2*N+1); n++)
        {
            cmn_coeffs[n-1] = 1.0/fn/(n+1);
            cumsum_cmn += cmn_coeffs[n-1]*ty;
            fn *= (n+2);
            ty *= Y;
        }

        // Calculate double series of coefficients
        cusfloat nf = 0.0;
        for (int m=1; m<N; m++)
        {
            ty = Y;
            for (int n=1; n<=(2*N+1-2*m); n++)
            {
                // Update n as float in order to perform the math operations
                // correctly
                nf = static_cast<cusfloat>(n);

                // Update series coefficients
                cmn_coeffs[n-1] = -((nf+2.0)/(nf+1.0))*cmn_coeffs[n+1];
                
                // Add new term to the series
                cumsum_cmn += cmn_coeffs[n-1]*tx*ty;

                ty *= Y;
            }
            // Update series coefficients
            tx *= fx;
        }
        cumsum_cmn += cmn_coeffs[0]*tx*Y;
        cumsum_cmn *= R;

        // Addd rest of the function
        wave_term = - 2.0*std::exp(-Y)*(besselj0(X)*std::log(Y/X+std::sqrt(1.0+pow2s(Y/X)))
                                        + PI*bessely0(X)/2.0
                                        + PI*struve0(X)*R/2.0/X
                                        + cumsum_cmn
                                        );
    }
    else if ((X>=3.7)&&(Y<=2.0))
    {
        // std::cout << "Region 4" << std::endl;

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
        wave_term = - PI*std::exp(-Y)*(struve0(X)+bessely0(X)) - 2.0*wave_term;
    }
    else
    {
        // std::cout << "Region 5" << std::endl;

        // Calculate expint residual values - R(x,y)
        cusfloat rxy = 0.0;
        if (Y<=4.0)
        {
            if (X <= 3.0)
            {
                rxy = eval_chebyshev_fit(
                    chebyinf::num_cheby_a0,
                    chebyinf::cheby_coeff_a0,
                    chebyinf::cheby_order_0_a0,
                    chebyinf::cheby_order_1_a0,
                    X,
                    Y
                );
            }
            else
            {
                rxy = eval_chebyshev_fit(
                    chebyinf::num_cheby_a1,
                    chebyinf::cheby_coeff_a1,
                    chebyinf::cheby_order_0_a1,
                    chebyinf::cheby_order_1_a1,
                    X,
                    Y
                );
            }
        }
        else if (Y<=8.0)
        {
            if (X <= 8.0)
            {
                rxy = eval_chebyshev_fit(
                    chebyinf::num_cheby_b0,
                    chebyinf::cheby_coeff_b0,
                    chebyinf::cheby_order_0_b0,
                    chebyinf::cheby_order_1_b0,
                    X,
                    Y
                );
            }
            else
            {
                rxy = eval_chebyshev_fit(
                    chebyinf::num_cheby_b1,
                    chebyinf::cheby_coeff_b1,
                    chebyinf::cheby_order_0_b1,
                    chebyinf::cheby_order_1_b1,
                    X,
                    Y
                );
            }
        }
        else
        {
            rxy = eval_chebyshev_fit(
                chebyinf::num_cheby_c,
                chebyinf::cheby_coeff_c,
                chebyinf::cheby_order_0_c,
                chebyinf::cheby_order_1_c,
                X,
                Y
            );
        }

        // Add oscilatory bessel and struve terms to
        // calculate full expint integral
        cusfloat R = std::sqrt(pow2s(X)+pow2s(Y));
        cusfloat f1 = 1/R - std::exp(-Y)/X + Y/pow3s(R)*rxy;
        wave_term = - PI*std::exp(-Y)*(struve0(X)+bessely0(X)) - 2.0*f1;
    }

    return wave_term;
}


cusfloat wave_term_inf_depth_dxndim_series(cusfloat X, cusfloat Y)
{
    using namespace std;

    cusfloat wave_term = 0.0;
    if ((X>=8.0) && (Y>=20.0))
    {
        // Calculate formulation parameters
        cusfloat R = std::sqrt(pow2s(X)+pow2s(Y));
        cusfloat sf = Y/R;
        cusfloat rfs = 1/R;
        cusfloat rf0 = pow(rfs, 4.0);
        cusfloat rf1 = pow(rfs, 3.0);
        cusfloat xy = X*Y;
        cusfloat pn = 1.0;
        for (int i=0; i<=4; i++)
        {
            // Add new term to the series
            wave_term += pn*(
                -legendre_poly_der_raw(i, sf)*xy*rf0
                -legendre_poly_raw(i, sf)*(i+1.0)*X*rf1
            );


            // Update formulation parameters
            rf0 *= rfs;
            rf1 *= rfs;
            pn *= static_cast<cusfloat>(i+1);
        }

        wave_term = - PI*std::exp(-Y)*(2.0/PI-struve1(X)-bessely1(X))-2.0*wave_term;
    }
    else if (Y > 1.9*X)
    {
        // Add radius term
        wave_term = 0.0;

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
        cumsum_n = (cumsum_m-c_expi);
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
            cumsum_n += n*tx/pow2s(fn)*(cumsum_m-c_expi);

            // Update n loop variables
            tx *= fx;
            fn *= n+1;
        }
        wave_term += (-X)*cumsum_n;
    }
    else if ((X >= 3.7) && (Y <= 0.25*X))
    {
        cusfloat f0 = 0.0;
        cusfloat f1 = 0.0;
        cusfloat fs = 1.0;
        cusfloat fn = 1.0;
        cusfloat fx = 1/pow2s(X);
        cusfloat fx1 = 1/X;
        cusfloat fy = pow2s(Y);
        cusfloat tx = 1.0;
        cusfloat ty = 1.0;
        cusfloat fmult = 1.0;
        cusfloat i2n = 1-std::exp(-Y);
        for (int i=0; i<=3; i++)
        {
            // Add expansion series term
            f0 += fs*tx/fn*fmult*i2n;
            f1 += -1*fs*2*i*tx*fx1/fn*fmult*i2n;

            // Update loop variables
            fs *= -1.0;
            fn *= (i+1);
            tx *= fx;
            ty *= fy;
            fmult *= (2*(i+1)-1.0)/2.0;
            i2n = ty - 2.0*(i+1.0)*ty/Y + 2.0*(i+1.0)*(2.0*(i+1.0)-1.0)*i2n;
        }
        wave_term = -fx*f0 + fx1*f1;
        wave_term = - PI*std::exp(-Y)*(2.0/PI-struve1(X)-bessely1(X)) - 2.0*wave_term;
    }
    else if ((X>0)&&(X<=3.7)&&(Y>0.0)&&(Y<=2.05))
    {
        // Define radial distance to be used along the module
        cusfloat R = std::sqrt(pow2s(X)+pow2s(Y));

        // Define factor variables
        cusfloat fx = pow2s(X);
        cusfloat tx = fx;
        cusfloat ty = 0.0;

        // Calculate double series first row of coefficients
        constexpr int N = 10;
        cusfloat cmn_coeffs[2*N+1];
        cusfloat fn = 2.0;
        cusfloat cumsum_cmn = 0.0;
        cusfloat cumsum_cmn_der = 0.0;
        ty = Y;
        for (int n=1; n<=(2*N+1); n++)
        {
            cmn_coeffs[n-1] = 1.0/fn/(n+1);
            cumsum_cmn += cmn_coeffs[n-1]*ty;
            fn *= (n+2);
            ty *= Y;
        }

        // Calculate double series of coefficients
        cusfloat nf = 0.0;
        for (int m=1; m<N; m++)
        {
            ty = Y;
            for (int n=1; n<=(2*N+1-2*m); n++)
            {
                // Update n as float in order to perform the math operations
                // correctly
                nf = static_cast<cusfloat>(n);

                // Update series coefficients
                cmn_coeffs[n-1] = -((nf+2.0)/(nf+1.0))*cmn_coeffs[n+1];
                
                // Add new term to the series
                cumsum_cmn += cmn_coeffs[n-1]*tx*ty;
                cumsum_cmn_der += cmn_coeffs[n-1]*2*m*tx*ty/X;

                ty *= Y;
            }
            // Update series coefficients
            tx *= fx;
        }
        cumsum_cmn += cmn_coeffs[0]*tx*Y;
        cumsum_cmn_der += cmn_coeffs[0]*2*N*tx*Y/X;
        cumsum_cmn *= X/R;
        cumsum_cmn_der *= R;

        // Addd rest of the function
        cusfloat ix = 1/X;
        cusfloat ydx = Y*ix;
        cusfloat s0 = sqrt(1+pow2s(ydx));
        cusfloat s1 = ydx + s0;
        cusfloat f1 = -besselj1(X)*log(s1)
                    + besselj0(X)*((-ydx*ix-pow2s(ydx)*ix/s0)/s1);
        cusfloat f2 = -PI*bessely1(X)/2.0;
        cusfloat f3 = 0.5*PI*ix*(-ix*struve0(X)+2.0/PI-struve1(X))*R
                    + 0.5*PI*struve0(X)/R;
        wave_term = - 2.0*std::exp(-Y)*(f1
                                        + f2
                                        + f3
                                        + cumsum_cmn + cumsum_cmn_der
                                        );
    }
    else if ((X>=3.7)&&(Y<=2.05))
    {
        cusfloat f0 = 0.0;
        cusfloat f1 = 0.0;
        cusfloat fs = 1.0;
        cusfloat fn = 1.0;
        cusfloat fx = 1/pow2s(X);
        cusfloat fx1 = 1/X;
        cusfloat fy = pow2s(Y);
        cusfloat tx = 1.0;
        cusfloat ty = 1.0;
        cusfloat fmult = 1.0;
        cusfloat i2n = 1-std::exp(-Y);
        for (int i=0; i<=7; i++)
        {
            // Add expansion series term
            f0 += fs*tx/fn*fmult*i2n;
            f1 += -1*fs*2*i*tx*fx1/fn*fmult*i2n;

            // Update loop variables
            fs *= -1.0;
            fn *= (i+1);
            tx *= fx;
            ty *= fy;
            fmult *= (2*(i+1)-1.0)/2.0;
            i2n = ty - 2.0*(i+1.0)*ty/Y + 2.0*(i+1.0)*(2.0*(i+1.0)-1.0)*i2n;
        }
        wave_term = -fx*f0 + fx1*f1;
        wave_term = - PI*std::exp(-Y)*(2.0/PI-struve1(X)-bessely1(X)) - 2.0*wave_term;
    }
    else
    {
        // Calculate expint residual values - R(x,y)
        cusfloat rxy = 0.0;
        cusfloat rxy_dx = 0.0;
        if (Y<=4.0)
        {
            if (X <= 3.0)
            {
                rxy = eval_chebyshev_fit(
                    chebyinf::num_cheby_a0,
                    chebyinf::cheby_coeff_a0,
                    chebyinf::cheby_order_0_a0,
                    chebyinf::cheby_order_1_a0,
                    X,
                    Y
                );
                rxy_dx = eval_chebyshev_fit_dxndim(
                    chebyinf::num_cheby_a0,
                    chebyinf::cheby_coeff_a0,
                    chebyinf::cheby_order_0_a0,
                    chebyinf::cheby_order_1_a0,
                    X,
                    Y
                );
            }
            else
            {
                rxy = eval_chebyshev_fit(
                    chebyinf::num_cheby_a1,
                    chebyinf::cheby_coeff_a1,
                    chebyinf::cheby_order_0_a1,
                    chebyinf::cheby_order_1_a1,
                    X,
                    Y
                );
                rxy_dx = eval_chebyshev_fit_dxndim(
                    chebyinf::num_cheby_a1,
                    chebyinf::cheby_coeff_a1,
                    chebyinf::cheby_order_0_a1,
                    chebyinf::cheby_order_1_a1,
                    X,
                    Y
                );
            }
        }
        else if (Y<=8.0)
        {
            if (X <= 8.0)
            {
                rxy = eval_chebyshev_fit(
                    chebyinf::num_cheby_b0,
                    chebyinf::cheby_coeff_b0,
                    chebyinf::cheby_order_0_b0,
                    chebyinf::cheby_order_1_b0,
                    X,
                    Y
                );
                rxy_dx = eval_chebyshev_fit_dxndim(
                    chebyinf::num_cheby_b0,
                    chebyinf::cheby_coeff_b0,
                    chebyinf::cheby_order_0_b0,
                    chebyinf::cheby_order_1_b0,
                    X,
                    Y
                );
            }
            else
            {
                rxy = eval_chebyshev_fit(
                    chebyinf::num_cheby_b1,
                    chebyinf::cheby_coeff_b1,
                    chebyinf::cheby_order_0_b1,
                    chebyinf::cheby_order_1_b1,
                    X,
                    Y
                );
                rxy_dx = eval_chebyshev_fit_dxndim(
                    chebyinf::num_cheby_b1,
                    chebyinf::cheby_coeff_b1,
                    chebyinf::cheby_order_0_b1,
                    chebyinf::cheby_order_1_b1,
                    X,
                    Y
                );
            }
        }
        else
        {
            rxy = eval_chebyshev_fit(
                chebyinf::num_cheby_c,
                chebyinf::cheby_coeff_c,
                chebyinf::cheby_order_0_c,
                chebyinf::cheby_order_1_c,
                X,
                Y
            );
            rxy_dx = eval_chebyshev_fit_dxndim(
                chebyinf::num_cheby_c,
                chebyinf::cheby_coeff_c,
                chebyinf::cheby_order_0_c,
                chebyinf::cheby_order_1_c,
                X,
                Y
            );
        }

        // Add oscilatory bessel and struve terms to
        // calculate full expint integral
        cusfloat R = std::sqrt(pow2s(X)+pow2s(Y));
        cusfloat r3inv = X/pow3s(R);
        cusfloat f1 = -r3inv
                    + exp(-Y)/pow2s(X)
                    -3*r3inv*Y/pow2s(R)*rxy
                    +Y*r3inv*rxy_dx/X;
        wave_term = - PI*std::exp(-Y)*(2.0/PI-struve1(X)-bessely1(X)) - 2.0*f1;
    }

    return wave_term;
}


cusfloat wave_term_inf_depth_dyndim_series(cusfloat X, cusfloat Y)
{
    // Calculate derivative terms
    cusfloat rinv = 1/std::sqrt(pow2s(X)+pow2s(Y));
    cusfloat wave_term = -(2*rinv + wave_term_inf_depth_series(X, Y));

    return wave_term;
}


cusfloat wave_term_inf_depth_num(cusfloat X, cusfloat Y)
{
    return -2*expint_inf_depth_num(X, Y) - PI*std::exp(-Y)*(bessely0(X)+struve0(X));
}


cusfloat wave_term_inf_depth_num_dxndim(cusfloat X, cusfloat Y)
{
    cusfloat sol = -2.0/X*expint_inf_depth_num_dxtndim(X, Y)
                + 2.0*Y/(X*std::sqrt(pow2s(X)+pow2s(Y)))
                + PI*std::exp(-Y)*(bessely1(X)+struve1(X)-2.0/PI);
    return sol;
}


cusfloat wave_term_inf_depth_num_dyndim(cusfloat X, cusfloat Y)
{
    return -2/std::sqrt(pow2s(X)+pow2s(Y))-wave_term_inf_depth_num(X, Y);
}