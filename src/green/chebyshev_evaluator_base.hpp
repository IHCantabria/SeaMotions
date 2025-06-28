
#ifndef __chebyshev_evaluator_base_hpp
#define __chebyshev_evaluator_base_hpp

// Include local modules
#include "../config.hpp"
#include "../math/chebyshev.hpp"
#include "../math/math_tools.hpp"


/**************************************/
/******* Define Module Macros *********/
/**************************************/
#define MAP_LOOP(x)                                                                                         \
for (std::size_t i = 0; i < N; i++) {                                                                       \
    x##m[i] = 2.0 * (x##s[i] - Derived::x##_min_region[nt[i]]) / Derived::d##x##_region[nt[i]] - 1.0;       \
}                                                                                                           \


template<typename Derived, std::size_t N>
struct ChebyshevEvaluatorBaseVector
{
    static void check_boundaries( cusfloat* xs, cusfloat* ys )
    {
        for ( std::size_t i=0; i<N; i++ )
        {
            xs[i] = std::max( std::min( xs[i], Derived::x_max_global ), Derived::x_min_global );
            ys[i] = std::max( std::min( ys[i], Derived::y_max_global ), Derived::x_min_global );
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
                break;
            }
        }

        return is_single;
    }

    static void evaluate( cusfloat* x, cusfloat* y, cusfloat* result )
    {
        // Check scaling of input variables if any
        cusfloat xs[N]; copy_vector<cusfloat, N>( x, xs );
        cusfloat ys[N]; copy_vector<cusfloat, N>( y, ys );
        scale( xs, ys );

        std::cout << "SCALE:" << std::endl;
        for ( std::size_t i=0; i<N; i++ )
        {
            std::cout << "Index: " << i;
            std::cout << " - xs: " << xs[i];
            std::cout << " - ys: " << ys[i] << std::endl;
        }

        // Check boundaries
        check_boundaries( xs, ys );

        // Get starting position
        std::size_t start_pos[N];
        std::size_t block_size[N];
        std::size_t nt[N];
        get_block_props( xs, ys, start_pos, block_size, nt );

        std::cout << "BLOCK PROPERTIES" << std::endl;
        for ( std::size_t i=0; i<N; i++ )
        {
            std::cout << "Index: " << i;
            std::cout << " - StartPos: " << start_pos[i];
            std::cout << " - BlockSize: " << block_size[i];
            std::cout << " - Nt: " << nt[i] << std::endl;
        }

        std::cout << "Derived Coeffs:" << std::endl;
        for ( int i=0; i<5; i++ )
        {
            std::cout << "Index: " << i << " - Coeffs: " << Derived::coeffs[i] << std::endl;
        }

        // Map coordinates
        cusfloat xm[N];
        cusfloat ym[N];
        
        MAP_LOOP( x )
        MAP_LOOP( y )
        
        // Check if all the input points are in the same block
        bool is_single_block = check_single_block( start_pos );

        // Evaluate chebyshev polynomials
        if ( is_single_block )
        {
            // evaluate_chebyshev_polynomials_2d_vector<Derived::max_cheby_order, N>( 
            //                                                                     &(Derived::coeffs[start_pos[0]]), 
            //                                                                     &(Derived::ncx[start_pos[0]]),
            //                                                                     &(Derived::ncy[start_pos[0]]),
            //                                                                     block_size[0],
            //                                                                     xm,
            //                                                                     ym,
            //                                                                     result
            //                                                                 );
            evaluate_chebyshev_polynomials_2d_vector_t<Derived, N>( 
                                                                        start_pos[0],
                                                                        block_size[0],
                                                                        xm,
                                                                        ym,
                                                                        result
                                                                    );
        }
        else
        {
            for ( std::size_t i=0; i<N; i++ )
            {
                evaluate_chebyshev_polynomials_2d_t<Derived>( 
                                                                start_pos[i],
                                                                block_size[i],
                                                                xm[i],
                                                                ym[i],
                                                                result[i]
                                                            );
            }
        }

    }

    static void get_block_props( cusfloat* xs, cusfloat* ys, std::size_t* sp, std::size_t* bd, std::size_t* ntv )
    {
        constexpr cusfloat dx = ( Derived::x_max_global - Derived::x_min_global ) / Derived::intervals_np;
        constexpr cusfloat dy = ( Derived::y_max_global - Derived::y_min_global ) / Derived::intervals_np;
        for ( std::size_t i=0; i<N; i++ )
        {
            // Estimate interval hash
            std::size_t nx = static_cast<std::size_t>( std::floor( ( xs[i] - Derived::x_min_global ) / dx ) );
            std::size_t ny = static_cast<std::size_t>( std::floor( ( ys[i] - Derived::y_min_global ) / dy ) );
            std::size_t nt = nx * Derived::intervals_np + ny;
    
            // Get starting position
            sp[i] = Derived::blocks_start[nt];
    
            // Get block dimensions;
            bd[i] = Derived::blocks_coeffs_np[nt];

            // Storage position in virtual grid
            ntv[i] = nt;
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
        evaluate_chebyshev_polynoials_2d<Derived::max_cheby_order>( 
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