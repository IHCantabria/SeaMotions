
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */


// Numerical integration for triple Hankel products over [1, inf).
#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <tuple>
#include <vector>

// Include local modules
#include "triple_hankel_far_field.hpp"
#include "../math/math_tools.hpp"


struct RotatedContourParameters {
                                    cusfloat sigma;
                                    cusfloat sigma_max;
                                    cusfloat segment_len;
                                    cusfloat x_max;
                                };

static void _validate_hankel_kinds( int hkind0, int hkind1 );


static RotatedContourParameters compute_rotated_contour_parameters(
                                                                        int         n1,
                                                                        int         n2,
                                                                        cusfloat    a,
                                                                        cusfloat    b,
                                                                        cusfloat    c,
                                                                        int         hkind0,
                                                                        int         hkind1,
                                                                        cusfloat    tail_cycles
                                                                    );


static cuscomplex hankel( int n, cusfloat x, int kind );


static cuscomplex romberg_complex(
                                    const std::function<cuscomplex(cusfloat)>& func,
                                    cusfloat    a,
                                    cusfloat    b,
                                    int         max_level,
                                    cusfloat    rtol,
                                    cusfloat    atol
                                );


static cusfloat romberg_real(
                                const std::function<cusfloat(cusfloat)>& func,
                                cusfloat    a,
                                cusfloat    b,
                                int         max_level,
                                cusfloat    rtol,
                                cusfloat    atol
                            );


static cuscomplex triple_hankel_integrand(
                                            cusfloat    r,
                                            int         n1,
                                            int         n2,
                                            int         n3,
                                            cusfloat    a,
                                            cusfloat    b,
                                            cusfloat    c,
                                            int         hkind0,
                                            int         hkind1
                                        );


static cuscomplex triple_hankel_integrand_asymp(
                                                    cusfloat y,
                                                    cusfloat x_max,
                                                    cusfloat sigma
                                                );


static void _validate_hankel_kinds( int hkind0, int hkind1 )
{
    if ( ( hkind0 != 1 && hkind0 != 2 ) 
         || 
         ( hkind1 != 1 && hkind1 != 2 ) 
        )
    {
        throw std::invalid_argument("hkind0 and hkind1 must be either 1 or 2.");
    }
}


static RotatedContourParameters compute_rotated_contour_parameters(
                                                                        int         n1,
                                                                        int         n2,
                                                                        cusfloat    a,
                                                                        cusfloat    b,
                                                                        cusfloat    c,
                                                                        int         hkind0,
                                                                        int         hkind1,
                                                                        cusfloat    tail_cycles
                                                                    )
{
    _validate_hankel_kinds( hkind0, hkind1 );

    cusfloat sign0          = hkind0 == 1 ? 1.0 : -1.0;
    cusfloat sign1          = hkind1 == 1 ? 1.0 : -1.0;
    cusfloat sigma          = sign0 * a + sign1 * b - c;
    cusfloat sigma_max      = std::max(std::abs(a) + std::abs(b), 1e-15);
    cusfloat segment_len    = 2.0 * PI / sigma_max / 10.0;

    cusfloat abs_sigma      = std::abs(sigma);
    cusfloat tail_scale     = abs_sigma > 1e-14 ? tail_cycles * 2.0 * PI / abs_sigma : tail_cycles * 2.0 * PI / sigma_max;
    cusfloat a_scale        = std::max(std::abs(a), 1e-15);
    cusfloat b_scale        = std::max(std::abs(b), 1e-15);
    cusfloat x_max          = std::max({static_cast<cusfloat>(n1) / a_scale, static_cast<cusfloat>(n2) / b_scale, tail_scale});
    x_max                   = x_max < 1.1 ? 1.0 : x_max;

    return { sigma, sigma_max, segment_len, x_max };

}


static cuscomplex hankel( int n, cusfloat x, int kind ) 
{
    if ( x == 0.0 ) 
    {
        return cuscomplex(0.0, 0.0);
    }

    cusfloat jn = std::cyl_bessel_j( static_cast<cusfloat>(n), x );
    cusfloat yn = std::cyl_neumann( static_cast<cusfloat>(n), x );

    if ( kind == 1 ) 
    {
        return cuscomplex( jn, yn );
    }

    if ( kind == 2 ) 
    {
        return cuscomplex( jn, -yn );
    }

    throw std::invalid_argument( "kind must be 1 or 2." );

}


static cuscomplex romberg_complex(
                                    const std::function<cuscomplex(cusfloat)>& func,
                                    cusfloat    a,
                                    cusfloat    b,
                                    int         max_level,
                                    cusfloat    rtol,
                                    cusfloat    atol
                                )
{
    cusfloat real_val = romberg_real(
                                        [&](cusfloat x) { return std::real(func(x)); },
                                        a,
                                        b,
                                        max_level,
                                        rtol,
                                        atol
                                    );

    cusfloat imag_val = romberg_real(
                                        [&](cusfloat x) { return std::imag(func(x)); },
                                        a,
                                        b,
                                        max_level,
                                        rtol,
                                        atol
                                    );

    return cuscomplex( real_val, imag_val );

}

static cusfloat romberg_real(
                                const std::function<cusfloat(cusfloat)>& func,
                                cusfloat    a,
                                cusfloat    b,
                                int         max_level,
                                cusfloat    rtol,
                                cusfloat    atol
                            )
{
    if ( a == b )
    {
        return 0.0;
    }

    std::vector<std::vector<cusfloat>> r(max_level + 1);
    r[0].resize(1);
    r[0][0] = 0.5 * (b - a) * (func(a) + func(b));

    for (int i = 1; i <= max_level; ++i)
    {
        int num_new = 1 << (i - 1);
        cusfloat h = (b - a) / static_cast<cusfloat>(1 << i);
        cusfloat sum = 0.0;
        for (int k = 1; k <= num_new; ++k) {
            cusfloat x = a + (2 * k - 1) * h;
            sum += func(x);
        }
        r[i].resize(i + 1);
        r[i][0] = 0.5 * r[i - 1][0] + h * sum;

        for (int j = 1; j <= i; ++j) {
            cusfloat factor = std::pow(4.0, j);
            r[i][j] = (factor * r[i][j - 1] - r[i - 1][j - 1]) / (factor - 1.0);
        }

        cusfloat err = std::abs(r[i][i] - r[i - 1][i - 1]);
        cusfloat tol = atol + rtol * std::abs(r[i][i]);
        if (err <= tol) {
            return r[i][i];
        }
    }

    return r[max_level][max_level];
}


static cuscomplex triple_hankel_integrand(
                                            cusfloat    r,
                                            int         n1,
                                            int         n2,
                                            int         n3,
                                            cusfloat    a,
                                            cusfloat    b,
                                            cusfloat    c,
                                            int         hkind0,
                                            int         hkind1
                                        )
{
    if ( r <= 0.0 )
    {
        return cuscomplex( 0.0, 0.0 );
    }

    cuscomplex val = r * hankel(n1, a * r, hkind0) * hankel(n2, b * r, hkind1) * hankel(n3, c * r, 2);

    return val;

}


static cuscomplex triple_hankel_integrand_asymp(
                                                    cusfloat y,
                                                    cusfloat x_max,
                                                    cusfloat sigma
                                                )
{
    cuscomplex z( x_max, y );
    if ( sigma >= 0.0 )
    {
        return cuscomplex( 0.0, 1.0 ) * std::exp( sigma * cuscomplex( -y, x_max ) ) / std::sqrt( z );
    }

    cuscomplex val = std::exp(sigma * cuscomplex(y, -x_max)) * std::conj(1.0 / std::sqrt(z));
    return cuscomplex( 0.0, -1.0 ) * std::conj( val );
}


cuscomplex integrate_triple_hankel(
                                        int                     n1,
                                        int                     n2,
                                        int                     n3,
                                        cusfloat                a,
                                        cusfloat                b,
                                        cusfloat                c,
                                        const TripleHankelIO&   options
                                    )
{
    _validate_hankel_kinds( options.hkind0, options.hkind1 );

    cusfloat k_max = std::max(std::abs(a), std::max(std::abs(b), std::abs(c)));
    if ( k_max == 0.0 )
    {
        return cuscomplex( 0.0, 0.0 );
    }

    cusfloat period         = 2.0 * PI / k_max;
    cusfloat segment_len    = std::max(period * options.segment_cycles, options.r_min);

    cuscomplex total(0.0, 0.0);
    cusfloat r_start        = std::max(options.r_min, 0.0);

    for ( int seg_idx = 0; seg_idx < options.max_segments; ++seg_idx )
    {
        cusfloat r_end = r_start + segment_len;
        cuscomplex seg_val = romberg_complex(
                                                [&](cusfloat r) {
                                                                    return triple_hankel_integrand(
                                                                                                        r, 
                                                                                                        n1, 
                                                                                                        n2, 
                                                                                                        n3, 
                                                                                                        a, 
                                                                                                        b, 
                                                                                                        c, 
                                                                                                        options.hkind0, 
                                                                                                        options.hkind1
                                                                                                    );
                                                                },
                                                r_start,
                                                r_end,
                                                options.romberg_max_level,
                                                options.rtol,
                                                options.atol
                                            );

        total += seg_val;

        if ( seg_idx + 1 >= options.min_segments )
        {
            if ( std::abs( seg_val ) <= options.atol + options.rtol * std::abs( total ) ) 
            {
                break;
            }
        }

        r_start = r_end;

    }

    return total;

}

cuscomplex integrate_triple_hankel_mod(
                                            int                     n1,
                                            int                     n2,
                                            int                     n3,
                                            cusfloat                a,
                                            cusfloat                b,
                                            cusfloat                c,
                                            const TripleHankelIO&   options
                                        )
{
    _validate_hankel_kinds( options.hkind0, options.hkind1 );

    cusfloat sign0 = options.hkind0 == 1 ? 1.0 : -1.0;
    cusfloat sign1 = options.hkind1 == 1 ? 1.0 : -1.0;
    RotatedContourParameters params = compute_rotated_contour_parameters(
                                                                            n1,
                                                                            n2,
                                                                            a,
                                                                            b,
                                                                            c,
                                                                            options.hkind0,
                                                                            options.hkind1,
                                                                            options.rotated_tail_cycles
                                                                        );

    cuscomplex finite_integral( 0.0, 0.0 );
    cusfloat head_start = std::max(options.r_min, 1.0);
    if ( params.x_max > head_start )
    {
        cusfloat    dx                  = params.x_max - head_start;
        int         segments            = static_cast<int>( std::ceil( dx / params.segment_len ) );
        cusfloat    head_segment_len    = dx / static_cast<cusfloat>( segments );

        for (int i = 0; i < segments; ++i)
        {
            cusfloat    start   = head_start + i * head_segment_len;
            cusfloat    end     = head_start + (i + 1) * head_segment_len;
            cuscomplex  segment = romberg_complex(
                                                    [&](cusfloat r) {
                                                        return triple_hankel_integrand(
                                                                                            r, 
                                                                                            n1, 
                                                                                            n2, 
                                                                                            n3, 
                                                                                            a, 
                                                                                            b, 
                                                                                            c, 
                                                                                            options.hkind0, 
                                                                                            options.hkind1
                                                                                        );
                                                    },
                                                    start,
                                                    end,
                                                    options.romberg_max_level,
                                                    options.rtol,
                                                    options.atol
                                                );
            finite_integral += segment;

            if ( options.verbose ) 
            {
                std::cout << "Head segment " << i << ": " << segment << " total: " << finite_integral
                          << " start: " << start << " end: " << end << "\n";
            }
        }
    }

    cuscomplex tail_integral( 0.0, 0.0 );
    cusfloat y_start = 0.0;
    for (int seg_idx = 0; seg_idx < options.max_segments; ++seg_idx) {
        cusfloat y_end = y_start + params.segment_len;
        cuscomplex segment = romberg_complex(
                                                [&](cusfloat y) { 
                                                                    return triple_hankel_integrand_asymp( y, params.x_max, params.sigma ); 
                                                                },
                                                y_start,
                                                y_end,
                                                options.romberg_max_level,
                                                options.rtol,
                                                options.atol
                                            );
        tail_integral += segment;

        if ( options.verbose ) 
        {
            std::cout << "Tail segment " << seg_idx << ": " << segment << "\n";
        }

        if ( seg_idx + 1 >= options.min_segments )
        {
            if ( std::abs( segment ) <= options.atol + options.rtol * std::abs( tail_integral ) ) 
            {
                break;
            }
        }

        y_start = y_end;
    }

    cusfloat    order_phase = -sign0 * n1 - sign1 * n2 + n3;
    cusfloat    one_phase   = -sign0 - sign1 + 1.0;
    cuscomplex  semi_scale  = (
                                    ( 2.0 / PI ) * std::sqrt( 2.0 / PI ) * std::sqrt( 1.0 / ( a * b * c ) ) 
                                    *
                                    std::exp( cuscomplex( 0.0, order_phase * PI / 2.0 + one_phase * PI / 4.0 ) )
                                );

    cuscomplex  total       = finite_integral + semi_scale * tail_integral;
    if ( options.verbose )
    {
        std::cout << finite_integral << " " << semi_scale * tail_integral << " " << total << "\n";
    }

    return total;
}


std::map<std::tuple<int, int, int>, cuscomplex> integrate_triple_hankel_orders(
                                                                                            const std::vector<int>& orders,
                                                                                            cusfloat                a,
                                                                                            cusfloat                b,
                                                                                            cusfloat                c,
                                                                                            const TripleHankelIO&   options
                                                                                        )
{
    std::map<std::tuple<int, int, int>, cuscomplex> results;
    for (int n1 : orders) {
        for (int n2 : orders) {
            for (int n3 : orders) {
                results[std::make_tuple(n1, n2, n3)] =
                    integrate_triple_hankel(n1, n2, n3, a, b, c, options);
            }
        }
    }
    return results;
}


cuscomplex triple_hankel_f(
                                int                     n1,
                                int                     n2,
                                int                     n3,
                                cusfloat                a,
                                cusfloat                b,
                                cusfloat                c,
                                TripleHankelIO&   options
                            )
{
    options.verbose = false;
    return integrate_triple_hankel_mod( n1, n2, n3, a, b, c, options );
}

cuscomplex triple_hankel_g(
                                int                 n1,
                                int                 n2,
                                int                 n3,
                                cusfloat            a,
                                cusfloat            b,
                                cusfloat            c,
                                TripleHankelIO&     options
                            )
{
    options.verbose = false;
    options.hkind0  = 1;
    cuscomplex val1 = integrate_triple_hankel_mod( n1, n2, n3, a, b, c, options );
    options.hkind0  = 2;
    cuscomplex val2 = integrate_triple_hankel_mod( n1, n2, n3, a, b, c, options );
    return ( val1 + val2 ) / 2.0;
}


cuscomplex triple_hankel_h(
                                int                 n1,
                                int                 n2,
                                int                 n3,
                                cusfloat            a,
                                cusfloat            b,
                                cusfloat            c,
                                TripleHankelIO&     options
                            )
{
    options.verbose = false;
    options.hkind0  = 2;
    options.hkind1  = 1;
    cuscomplex val1 = integrate_triple_hankel_mod( n1, n2, n3, a, b, c, options );
    options.hkind1  = 2;
    cuscomplex val2 = integrate_triple_hankel_mod( n1, n2, n3, a, b, c, options );
    return ( val1 + val2 ) / 2.0;
}


#ifdef SEAMOTIONS_INTEGRATE_RADIUS_FAR_FIELD_DEMO
int main() {
    int nu0 = 20;
    int nu1 = 20;
    int nu2 = 20;
    cuscomplex val = seamotions::integrate_triple_hankel(
        nu0, nu1, nu2, 1.0, 1.2, 0.8, {});
    std::cout << "Integral value: " << val << "\n";
    return 0;
}
#endif
