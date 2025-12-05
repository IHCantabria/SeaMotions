
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

#pragma once

// Include local modules
#include "chebyshev.hpp"


template<const std::size_t max_cheby_order, const std::size_t N>
void evaluate_chebyshev_polynomials_2d_vector( 
                                                cusfloat*   coeffs, 
                                                int*        ncx, 
                                                int*        ncy,
                                                int         np,
                                                cusfloat*   x,
                                                cusfloat*   y,
                                                cusfloat*   result
                                            )
{
    // Reset to zero in order to spurious values
    for ( int i=0; i<N; i++ )
    {
        result[i] = 0.0;
    }

    // Get chebyshev values up to the maximum order
    static cusfloat poly_x[ N * ( max_cheby_order + 1 ) ];
    static cusfloat poly_y[ N * ( max_cheby_order + 1 ) ];

    chebyshev_poly_upto_order<max_cheby_order, N>( x, poly_x );
    chebyshev_poly_upto_order<max_cheby_order, N>( y, poly_y );

    // Loop over chebyshev coefficients and their corresponding orders
    for ( std::size_t i=0; i<static_cast<std::size_t>( np ); i++ )
    {
        for ( int j=0; j<N; j++ )
        {
            result[j] += coeffs[i] * poly_x[ ncx[i]*N + j ] * poly_y[ ncy[i]*N + j ];
        }
    }
}


template<typename Derived, const std::size_t N, int mode_loop>
void evaluate_chebyshev_polynomials_2d_vector_t( 
                                                    const std::size_t   sp,
                                                    const std::size_t   np,
                                                    const std::size_t   nt,
                                                    const std::size_t   n,
                                                    cusfloat*           x,
                                                    cusfloat*           y,
                                                    cusfloat*           result
                                                )
{
    // Reset to zero in order to spurious values
    STATIC_LOOP( n, N, result[i] = 0.0; )

    // Get chebyshev values up to the maximum order
    static cusfloat poly_x[ N * ( Derived::max_cheby_order + 1 ) ];
    static cusfloat poly_y[ N * ( Derived::max_cheby_order + 1 ) ];

    chebyshev_poly_upto_order<N, mode_loop>( n, Derived::blocks_max_cheby_order[nt], x, poly_x );
    chebyshev_poly_upto_order<N, mode_loop>( n, Derived::blocks_max_cheby_order[nt], y, poly_y );

    // Loop over chebyshev coefficients and their corresponding orders
    for ( std::size_t j=0; j<static_cast<std::size_t>( np ); j++ )
    {
        STATIC_LOOP( n, N, result[i] += Derived::coeffs[sp+j] * poly_x[ Derived::ncx[sp+j]*n + i ] * poly_y[ Derived::ncy[sp+j]*n + i ]; )
    }
}


template<const std::size_t max_cheby_order>
void evaluate_chebyshev_polynomials_2d( 
                                            const cusfloat*     coeffs, 
                                            const std::size_t*  ncx, 
                                            const std::size_t*  ncy,
                                            std::size_t         np,
                                            cusfloat            x,
                                            cusfloat            y,
                                            cusfloat&           result,
                                            std::size_t&        count
                                        )
{
    // Reset to zero in order to spurious values
    result = 0.0;

    // Get chebyshev values up to the maximum order
    cusfloat poly_x[ ( max_cheby_order + 1 ) ];
    cusfloat poly_y[ ( max_cheby_order + 1 ) ];

    // std::cout << "chebyshev_poly_upto_order" << std::endl;
    chebyshev_poly_upto_order<max_cheby_order>( x, poly_x, count );
    chebyshev_poly_upto_order<max_cheby_order>( y, poly_y, count );

    // Loop over chebyshev coefficients and their corresponding orders
    count++;
    for ( std::size_t i=0; i<np; i++ )
    {
        result += coeffs[i] * poly_x[ ncx[i] ] * poly_y[ ncy[i] ];
        count++;
    }

}


template<typename Derived>
void evaluate_chebyshev_polynomials_2d_t( 
                                            const std::size_t   sp,
                                            const std::size_t   np,
                                            const std::size_t   nt,
                                            cusfloat            x,
                                            cusfloat            y,
                                            cusfloat&           result
                                        )
{
    // Reset to zero in order to spurious values
    result = 0.0;

    // Get chebyshev values up to the maximum order
    cusfloat poly_x[ ( Derived::max_cheby_order + 1 ) ];
    cusfloat poly_y[ ( Derived::max_cheby_order + 1 ) ];

    // std::cout << "chebyshev_poly_upto_order" << std::endl;
    chebyshev_poly_upto_order( Derived::blocks_max_cheby_order[nt], x, poly_x );
    chebyshev_poly_upto_order( Derived::blocks_max_cheby_order[nt], y, poly_y );
    
    // Loop over chebyshev coefficients and their corresponding orders
    for ( std::size_t i=0; i<np; i++ )
    {
        result += Derived::coeffs[sp+i] * poly_x[ Derived::ncx[sp+i] ] * poly_y[ Derived::ncy[sp+i] ];
    }
    
}
