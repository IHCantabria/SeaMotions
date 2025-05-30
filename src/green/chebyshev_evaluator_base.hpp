
#ifndef __chebyshev_evaluator_base_hpp
#define __chebyshev_evaluator_base_hpp

// Include local modules
#include "../config.hpp"
#include "../math/chebyshev.hpp"
#include "../math/math_tools.hpp"


template<typename Derived, std::size_t N>
struct ChebyshevEvaluatorBaseVector
{
    static void check_boundaries( cusfloat* xs, cusfloat* ys )
    {
        for ( std::size_t i=0; i<N; i++ )
        {
            xs[i] = std::max( std::min( xs[i], Derived::x_max ), Derived::x_min );
            ys[i] = std::max( std::min( ys[i], Derived::x_max ), Derived::x_min );
        }
    }

    static bool check_single_block( std::size_t* start_pos )
    {
        bool is_single = true;
        for ( std::size_t i=1; i<N; i++ )
        {
            if ( start_pos[0] != start_pos[i] )
            {
                is_single = false;
                break
            }
        }

        return is_single;
    }

    static void evaluate( cusfloat* x, cusfloat* y, cusfloat* result )
    {
        // Check scaling of input variables if any
        cusfloat xs[N]; copy_vector<N>( x, xs );
        cusfloat ys[N]; copy_vector<N>( y, ys );
        scale( xs, ys );

        // Check boundaries
        check_boundaries( x, y )

        // Get starting position
        std::size_t start_pos[N];
        std::size_t block_size[N];
        get_block_props( xs, ys, start_pos, block_size );
        
        // Check if all the input points are in the same block
        bool is_single_block = check_single_block( start_pos );

        // Evaluate chebyshev polynomials
        if ( is_single_block )
        {
            evaluate_chebyshev_polynomials_vector_2d<Derived::max_cheby_order, N>( 
                                                                                &(Derived::coeffs[start_pos[0]]), 
                                                                                &(Derived::ncx[start_pos[0]]),
                                                                                &(Derived::ncy[start_pos[0]]),
                                                                                block_size[0],
                                                                                xs,
                                                                                ys,
                                                                                result
                                                                            );
        }
        else
        {
            for ( std::size_t i=0; i<N; i++ )
            {
                evaluate_chebyshev_polynomials_2d<Derived::max_cheby_order>( 
                                                                                &(Derived::coeffs[start_pos[i]]), 
                                                                                &(Derived::ncx[start_pos[i]]),
                                                                                &(Derived::ncy[start_pos[i]]),
                                                                                block_size[i],
                                                                                xs[i],
                                                                                ys[i],
                                                                                &(result[i])
                                                                            );
            }
        }

    }

    static void get_block_props( cusfloat* xs, cusfloat* ys, std::size_t* sp, std::size_t* bd )
    {
        for ( int i=0; i<N; i++ )
        {
            // Estimate interval hash
            int nx = static_cast<int>( std::floor( xs[i] / Derived::dx ) );
            int ny = static_cast<int>( std::floor( ys[i] / Derived::dy ) );
            int nt = nx * Derived::intervals_np + ny;
    
            // Get starting position
            sp[i] = Derived::block_start[nt];
    
            // Get block dimensions;
            bd[i] = Derived::blocks_coeffs_np[nt];
        }
    }

    static void scale( cusfloat* xs, cusfloat* ys )
    {
        if constexpr( Derived::x_log_scale )
        {
            for ( std::size_t i=0; i<N; i++ )
            {
                xs[i] = std::log10( xs[i] );
            }
        }

        if constexpr( Derived::y_log_scale )
        {
            for ( std::size_t i=0; i<N; i++ )
            {
                ys = std::log10( ys );
            }
        }
    }
};


template<typename Derived>
struct ChebyshevEvaluatorBase
{
    static void check_boundaries( cusfloat& xs, cusfloat& ys )
    {
        xs = std::max( std::min( xs, Derived::x_max_global ), Derived::x_min_global );
        ys = std::max( std::min( ys, Derived::x_max_global ), Derived::x_min_global );
    }

    static cusfloat evaluate( cusfloat x, cusfloat y, std::size_t& count )
    {
        // std::cout << "input -> x: " << x << " - y: " << y << std::endl;
        // Check scaling of input variables if any
        cusfloat xs = x;
        cusfloat ys = y;
        scale( xs, ys );
        // std::cout << "scale -> xs: " << xs << " - ys: " << ys << std::endl;

        // Check boundaries
        // std::cout << "Evaluating boundaries: " << std::endl;
        check_boundaries( xs, ys );
        // std::cout << "boundary check -> xs: " << xs << " - ys: " << ys << std::endl;
        
        // Get starting position
        std::size_t start_pos   = 0;
        std::size_t block_size  = 0;
        cusfloat    x_min       = 0.0;
        cusfloat    dx          = 0.0;
        cusfloat    y_min       = 0.0;
        cusfloat    dy          = 0.0;

        // std::cout << "Get Block Props: " << std::endl;
        get_block_props( xs, ys, start_pos, block_size, x_min, dx, y_min, dy );
        // std::cout << "StartPos: " << start_pos << " - BlockSize: " << block_size << std::endl;

        // std::cout << "Start Pos: " << start_pos << " - Block Size: " << block_size << std::endl;

        // Map to [-1, 1] domain
        // std::cout << "x_min: "  << x_min << " - y_min: " << y_min << " - dx: " << dx << " - dy: " << dy << std::endl;

        cusfloat xsm = 0.0; cusfloat ysm = 0.0;
        map( xs, x_min, dx, xsm );
        map( ys, y_min, dy, ysm );

        // std::cout << "Map -> xsm: " << xsm << " - ysm: " << ysm << std::endl;
        
        // Evaluate chebyshev polynomials
        cusfloat result = 0.0;
        // std::cout << "Evaluate chebyshev: " << xs << " " << ys << std::endl;
        // std::cout << "Evaluate chebyshev: " << xsm << " " << ysm << std::endl;
        evaluate_chebyshev_polynomials_2d<Derived::max_cheby_order>( 
                                                                        &(Derived::coeffs[start_pos]), 
                                                                        &(Derived::ncx[start_pos]),
                                                                        &(Derived::ncy[start_pos]),
                                                                        block_size,
                                                                        xsm,
                                                                        ysm,
                                                                        result,
                                                                        count
                                                                    );

        // Scale interpolated funcion if any
        if constexpr( Derived::fcn_log_scale )
        {
            result = std::pow( 10.0, result );
        }
        // std::cout << "Result: " << result << " -> Done!" << std::endl;

        return result;
    }

    static void get_block_props( cusfloat xs, cusfloat ys, std::size_t& sp, std::size_t& bd, cusfloat& x_min, cusfloat& dx, cusfloat& y_min, cusfloat& dy)
    {
        // Estimate interval hash
        int nx = static_cast<int>( std::floor( ( xs - Derived::x_min_global ) / Derived::dx_min_region ) );
        int ny = static_cast<int>( std::floor( ( ys - Derived::y_min_global ) / Derived::dy_min_region ) );
        int nt = nx * Derived::intervals_np + ny;
        // std::cout << "nx: " << nx << " - ny: " << ny << " - nt: " << nt << std::endl;

        // Get starting position
        sp = Derived::blocks_start[nt];

        // Get block dimensions;
        bd = Derived::blocks_coeffs_np[nt];

        // Get target region X mininum domain position
        x_min = Derived::x_min_region[nt];

        // Get target region X span
        dx = Derived::dx_region[nt];

        // Get target region Y mininum domain position
        y_min = Derived::y_min_region[nt];

        // Get target region Y span
        dy = Derived::dy_region[nt];
    }

    static void map( cusfloat var, cusfloat var_min, cusfloat var_step, cusfloat& var_map )
    {
        var_map = 2.0 * ( var - var_min ) / var_step - 1.0;
    }

    static void scale( cusfloat& xs, cusfloat& ys )
    {
        if constexpr( Derived::x_log_scale )
        {
            xs = std::log10( xs );
        }

        if constexpr( Derived::y_log_scale )
        {
            ys = std::log10( ys );
        }
    }
};

#endif