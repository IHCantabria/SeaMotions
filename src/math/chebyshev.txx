
// Include general usage libraries
#include <cstddef>

// Include local modules
#include "../../src/config.hpp"


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

    chebyshev_poly_upto_order<MAX_ORDER, N>( x, poly_x );
    chebyshev_poly_upto_order<MAX_ORDER, N>( y, poly_y );

    // Loop over chebyshev coefficients and their corresponding orders
    for ( std::size_t i=0; i<static_cast<std::size_t>( np ); i++ )
    {
        for ( int j=0; j<N; j++ )
        {
            result[j] += coeffs[i] * poly_x[ ncx[i]*N + j ] * poly_y[ ncy[i]*N + j ];
        }
    }
}


template<typename Derived, const std::size_t N>
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
    for ( std::size_t i=0; i<n; i++ )
    {
        result[i] = 0.0;
    }

    // Get chebyshev values up to the maximum order
    static cusfloat poly_x[ N * ( Derived::max_cheby_order + 1 ) ];
    static cusfloat poly_y[ N * ( Derived::max_cheby_order + 1 ) ];

    chebyshev_poly_upto_order( n, Derived::blocks_max_cheby_order[nt], x, poly_x );
    chebyshev_poly_upto_order( n, Derived::blocks_max_cheby_order[nt], y, poly_y );

    // Loop over chebyshev coefficients and their corresponding orders
    for ( std::size_t i=0; i<static_cast<std::size_t>( np ); i++ )
    {
        for ( std::size_t j=0; j<n; j++ )
        {
            result[j] += Derived::coeffs[sp+i] * poly_x[ Derived::ncx[sp+i]*n + j ] * poly_y[ Derived::ncy[sp+i]*n + j ];
        }
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


template<const std::size_t max_order, const int N>
void chebyshev_poly_upto_order( cusfloat* x, cusfloat* results )
{
    if constexpr( max_order == 0 )
    {
        for ( int i=0; i<N; i++ )
        {
            results[i] = 1.0;
        }
        return;
    }
    if constexpr( max_order == 1 )
    {
        for ( int i=0; i<N; i++ )
        {
            results[i]              = 1.0;
            results[max_order*N+i]  = x[i];
        }
        return;
    }

    for ( int i=0; i<N; i++ )
    {
        results[i]      = 1.0;
        results[N+i]    = x[i];
    }

    for ( std::size_t i = 2; i <= max_order; i++ )
    {
        for ( int j=0; j<N; j++ )
        {
            results[i*N+j]  = 2.0 * x[j] * results[(i-1)*N+j] - results[(i-2)*N+j]; // Recurrence relation
        }
    }

}


template<const std::size_t max_order>
void chebyshev_poly_upto_order( cusfloat x, cusfloat* results )
{
    if constexpr( max_order == 0 )
    {
        results[0] = 1.0;
        return;
    }
    if constexpr( max_order == 1 )
    {
        results[0]  = 1.0;
        results[1]  = x;
        return;
    }

    results[0]  = 1.0;
    results[1]  = x;

    for ( std::size_t i = 2; i <= max_order; i++ )
    {
        results[i]  = 2.0 * x * results[(i-1)] - results[(i-2)]; // Recurrence relation
    }

}


inline void chebyshev_poly_upto_order( const std::size_t n, const std::size_t max_order, cusfloat* x, cusfloat* results )
{
    if ( max_order == 0 )
    {
        for ( std::size_t i=0; i<n; i++ )
        {
            results[i] = 1.0;
        }
        return;
    }
    if ( max_order == 1 )
    {
        for ( std::size_t i=0; i<n; i++ )
        {
            results[i]              = 1.0;
            results[max_order*n+i]  = x[i];
        }
        return;
    }

    for ( std::size_t i=0; i<n; i++ )
    {
        results[i]      = 1.0;
        results[n+i]    = x[i];
    }

    for ( std::size_t i = 2; i <= max_order; i++ )
    {
        for ( std::size_t j=0; j<n; j++ )
        {
            results[i*n+j]  = 2.0 * x[j] * results[(i-1)*n+j] - results[(i-2)*n+j]; // Recurrence relation
        }
    }

}


inline void chebyshev_poly_upto_order( std::size_t max_order, cusfloat x, cusfloat* results )
{
    if ( max_order == 0 )
    {
        results[0] = 1.0;
        return;
    }
    if ( max_order == 1 )
    {
        results[0]  = 1.0;
        results[1]  = x;
        return;
    }

    results[0]  = 1.0;
    results[1]  = x;

    for ( std::size_t i = 2; i <= max_order; i++ )
    {
        results[i]  = 2.0 * x * results[(i-1)] - results[(i-2)]; // Recurrence relation
    }

}