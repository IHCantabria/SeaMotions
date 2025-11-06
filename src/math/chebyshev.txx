
// Include general usage libraries
#include <cstddef>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/static_tools.hpp"




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


template< std::size_t N, int mode_loop>
inline void chebyshev_poly_upto_order( const std::size_t n, const std::size_t max_order, cusfloat* x, cusfloat* results )
{
    if ( max_order == 0 )
    {
        STATIC_LOOP( n, N, results[i] = 1.0; )
        return;
    }
    if ( max_order == 1 )
    {
        STATIC_LOOP( n, N, results[i]              = 1.0; )
        STATIC_LOOP( n, N, results[max_order*n+i]  = x[i]; )
        return;
    }

    // Calculate zero and first orders
    STATIC_LOOP( n, N, results[i]      = 1.0; )
    STATIC_LOOP( n, N, results[n+i]    = x[i]; )

    // Calculate rest of orders
    for ( std::size_t j = 2; j <= max_order; j++ )
    {
        STATIC_LOOP( n, N, results[j*n+i]  = 2.0 * x[i] * results[(j-1)*n+i] - results[(j-2)*n+i]; )
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