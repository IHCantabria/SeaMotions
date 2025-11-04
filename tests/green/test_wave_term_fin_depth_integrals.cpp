
// Include general usage libraries
#include <fstream>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/green/chebyshev_traits.hpp"
#include "../../src/green/chebyshev_evaluator_interface.hpp"
#include "../../src/green/integrals_database.hpp"
#include "../../src/green/fin_depth_coeffs/L1.hpp"


/************ Define module macros ******************/
#define CHECK_GLOBAL_ERROR( err_thr, N, err, err_count, db_name, fcn_type )                                         \
{                                                                                                                   \
    if( ( err_count / N ) > 0.01 )                                                                                  \
    {                                                                                                               \
        std::cerr << "TEST FAILED for " << db_name << " and " << fcn_type << " function" << std::endl;               \
        std::cerr << "--> Global Analysis: test exceeds the number of faults over the threshold." << std::endl;     \
        throw std::exception( );                                                                                    \
    }                                                                                                               \
                                                                                                                    \
    if ( err > 10 * err_thr )                                                                                       \
    {                                                                                                               \
        std::cerr << "TEST FAILED for " << db_name << " and" << fcn_type << " function" << std::endl;               \
        std::cerr << "--> Global Analysis: test exceeds the error threshold value in a factor of 10." << std::endl; \
        throw std::exception( );                                                                                    \
    }                                                                                                               \
}                                                                                                                   \

#define COND_LINE( expr, line )     \
{                                   \
    if constexpr( expr )            \
    {                               \
        line;                       \
    }                               \
}                                   \

#define UPDATE_GLOBAL_ERROR( r0, r1, max_err, count_thr )   \
{                                                           \
    cusfloat err = std::abs( r0 - r1 );                     \
    if ( err > max_err )                                    \
    {                                                       \
        max_err = err;                                      \
    }                                                       \
    count_thr++;                                            \
}                                                           \


enum class ErrorControl{ Local, Global };


cusfloat EPS_PREC = 5E-4;


template<int N, int M>
struct RefData
{
public:
    cusfloat*   X0          = nullptr;
    cusfloat*   X1          = nullptr;
    cusfloat*   X2          = nullptr;
    cusfloat*   G           = nullptr;
    cusfloat*   G_dx0       = nullptr;
    cusfloat*   G_dx1       = nullptr;
    int         num_x0      = 0;
    int         num_x1      = 0;
    int         num_x2      = 0;
    int         total_size  = 0;

    static_assert( N>0 || N<4 );
    static_assert( M>0 || M<4 );

    RefData( std::string fipath )
    {
        std::string aux_str;

        // Open file unit
        std::ifstream infile( fipath );

        // Read number of X0 values
        infile >> aux_str >> this->num_x0;
        this->X0 = generate_empty_vector<cusfloat>( this->num_x0 );
        for ( int i=0; i<this->num_x0; i++ )
        {
            infile >> this->X0[i];
        }

        if constexpr( N > 1 )
        {
            // Read number of X1 values
            infile >> aux_str >> this->num_x1;
            this->X1 = generate_empty_vector<cusfloat>( this->num_x1 );
            for ( int i=0; i<this->num_x1; i++ )
            {
                infile >> this->X1[i];
            }

            // Read number of X2 values
            if constexpr( N > 2 )
            {
                infile >> aux_str >> this->num_x2;
                this->X2 = generate_empty_vector<cusfloat>( this->num_x2 );
                for ( int i=0; i<this->num_x2; i++ )
                {
                    infile >> this->X2[i];
                }
            }
        }

        // Allocate space to storage space for functions
        COND_LINE( N > 0, this->total_size  = this->num_x0 );
        COND_LINE( N > 1, this->total_size *= this->num_x1 );
        COND_LINE( N > 2, this->total_size *= this->num_x2 );

        COND_LINE( M > 0, this->G       = generate_empty_vector<cusfloat>( this->total_size ) );
        COND_LINE( M > 1, this->G_dx0   = generate_empty_vector<cusfloat>( this->total_size ) );
        COND_LINE( M > 2, this->G_dx1   = generate_empty_vector<cusfloat>( this->total_size ) );

        // Read functions header
        COND_LINE( M > 0, infile >> aux_str );
        COND_LINE( M > 1, infile >> aux_str );
        COND_LINE( M > 2, infile >> aux_str );

        // Read functions value
        for ( int i=0; i<this->total_size; i++ )
        {
            COND_LINE( M > 0, infile >> this->G[i] );
            COND_LINE( M > 1, infile >> this->G_dx0[i] );
            COND_LINE( M > 2, infile >> this->G_dx1[i] );
        }

        // Close file unit
        infile.close( );
    }

    void print( void )
    {
        // Print out integral domain extents
        COND_LINE( N > 0, std::cout << "Num.X0: " << this->num_x0 << std::endl );
        COND_LINE( N > 0, print_vector( this->num_x0, this->X0, 0, 6 ) );
        COND_LINE( N > 1, std::cout << "Num.X1: " << this->num_x1 << std::endl );
        COND_LINE( N > 1, print_vector( this->num_x1, this->X1, 0, 6 ) );
        COND_LINE( N > 2, std::cout << "Num.X2: " << this->num_x2 << std::endl );
        COND_LINE( N > 2, print_vector( this->num_x2, this->X2, 0, 6 ) );

        // Print out function values
        COND_LINE( M > 0, std::cout << "G" );
        COND_LINE( M > 1, std::cout << "    G_dx0" );
        COND_LINE( M > 2, std::cout << "    G_dx1" );
        COND_LINE( M > 0, std::cout << std::endl );

        for ( int i=0; i<this->total_size; i++ )
        {
            COND_LINE( M > 0, std::cout << this->G[i] );
            COND_LINE( M > 1, std::cout << "    " << this->G_dx0[i] );
            COND_LINE( M > 2, std::cout << "    " << this->G_dx1[i] );
            COND_LINE( M > 0, std::cout << std::endl );
        }   
    }

    ~RefData( )
    {
        COND_LINE( M > 0, mkl_free( this->X0 ) );
        COND_LINE( M > 1, mkl_free( this->X1 ) );
        COND_LINE( M > 2, mkl_free( this->X2 ) );

        COND_LINE( M > 0, mkl_free( this->G ) );
        COND_LINE( M > 1, mkl_free( this->G_dx0 ) );
        COND_LINE( M > 2, mkl_free( this->G_dx1 ) );
    }

};


template<typename T, ErrorControl EC>
void compare_1d_database( std::string fipath, std::string db_name )
{
    std::cout << "COMPARE " << db_name << " DATABASE" << std::endl;
    // Define chebyshev evaluator functions
    using TEV       = ChebyshevTraits<T>;

    // Load data
    RefData<1, 1> ref_data( fipath );

    // Prepare data variables to global error control if any
    cusfloat    g_max_err       = 0.0;
    int         g_count_thr     = 0;

    // Loop over database to check its value
    for ( int i=0; i<ref_data.num_x0; i++ )
    {
        // Fold database
        fold_database_1d<T>( ref_data.X0[i] );

        if ( !assert_scalar_equality( TEV::coeffs, ref_data.G[i], EPS_PREC ) )
        {
            if constexpr( EC == ErrorControl::Local )
            {
                std::cerr << "TEST FAILED for " << db_name << " and G function" << std::endl;
                std::cerr << "Expected: " << ref_data.G[i] << " - Calculated: " << TEV::coeffs << std::endl;
                std::cerr << "H[" << i << "]: " << ref_data.X0[i] << std::endl;
                throw std::exception( );
            }
            else if constexpr( EC == ErrorControl::Global )
            {
                UPDATE_GLOBAL_ERROR( TEV::coeffs, ref_data.G[i], g_max_err, g_count_thr )
            }
        }
    }

    // Check global error if selected
    if constexpr( EC == ErrorControl::Global )
    {
        int N2 = ref_data.num_x0;
        CHECK_GLOBAL_ERROR( EPS_PREC, N2, g_max_err,     g_count_thr,     db_name, "G"     )
    }
    std::cout << " -> DONE" << std::endl;
}


template<typename T, typename T_dX0, ErrorControl EC>
void compare_2d_database( std::string fipath, std::string db_name )
{
    std::cout << "COMPARE " << db_name << " DATABASE" << std::endl;
    // Define local variables
    constexpr int   N = 1;
    int             index = 0;
    cusfloat        result[N];
    cusfloat        result_dx0[N];

    // Define chebyshev evaluator functions
    using TEV       = ChebyshevEvaluatorBaseVector<ChebyshevTraits<T>, N, STATIC_LOOP_ON>;
    using T_dX0EV   = ChebyshevEvaluatorBaseVector<ChebyshevTraits<T_dX0>, N, STATIC_LOOP_ON>;

    // Load data
    RefData<2, 2> ref_data( fipath );

    // Prepare data variables to global error control if any
    cusfloat    g_max_err       = 0.0;
    int         g_count_thr     = 0;
    cusfloat    g_dx0_max_err   = 0.0;
    int         g_dx0_count_thr = 0;

    // Loop over database to check its value
    for ( int i=0; i<ref_data.num_x0; i++ )
    {
        for ( int j=0; j<ref_data.num_x1; j++ )
        {
            index = ( 
                            i * ( ref_data.num_x0 )
                            +
                            j
                        );
            
            // Interpolate in the database
            TEV::evaluate( N, &(ref_data.X0[i]), &(ref_data.X1[j]), result );
            T_dX0EV::evaluate( N, &(ref_data.X0[i]), &(ref_data.X1[j]), result_dx0 );

            if ( !assert_scalar_equality( result[0], ref_data.G[index], EPS_PREC ) )
            {
                if constexpr( EC == ErrorControl::Local )
                {
                    std::cerr << "TEST FAILED for " << db_name << " and G function" << std::endl;
                    std::cerr << "Expected: " << ref_data.G[index] << " - Calculated: " << result[0] << std::endl;
                    std::cerr << "A[" << i  << "]: " << ref_data.X0[i] << " - B[" << j << "]: " << ref_data.X1[j] << std::endl;
                    throw std::exception( );
                }
                else if constexpr( EC == ErrorControl::Global )
                {
                    UPDATE_GLOBAL_ERROR( result[0], ref_data.G[index], g_max_err, g_count_thr )
                }
            }
            
            if ( !assert_scalar_equality( result_dx0[0], ref_data.G_dx0[index], EPS_PREC ) )
            {
                if constexpr( EC == ErrorControl::Local )
                {
                    std::cerr << "TEST FAILED for " << db_name << " and G_dA function" << std::endl;
                    std::cerr << "Expected: " << ref_data.G_dx0[index] << " - Calculated: " << result_dx0[0] << std::endl;
                    std::cerr << "A[" << i  << "]: " << ref_data.X0[i] << " - B[" << j << "]: " << ref_data.X1[j] << " - H[" << j << "]: " << std::endl;
                    throw std::exception( );
                }
                else if constexpr( EC == ErrorControl::Global )
                {
                    UPDATE_GLOBAL_ERROR( result_dx0[0], ref_data.G_dx0[index], g_dx0_max_err, g_dx0_count_thr )
                }
            }
        }

    }

    // Check global error if selected
    if constexpr( EC == ErrorControl::Global )
    {
        int N2 = ref_data.num_x0 * ref_data.num_x1;
        CHECK_GLOBAL_ERROR( EPS_PREC, N2, g_max_err,     g_count_thr,     db_name, "G"     )
        CHECK_GLOBAL_ERROR( EPS_PREC, N2, g_dx0_max_err, g_dx0_count_thr, db_name, "G_dx0" )
    }
    std::cout << " -> DONE" << std::endl;
}


template<typename T, typename T_dX0, typename T_dX1, ErrorControl EC>
void compare_3d_database( std::string fipath, std::string db_name )
{
    std::cout << "COMPARE " << db_name << " DATABASE" << std::endl;
    // Define local variables
    constexpr int   N = 1;
    int             index = 0;
    cusfloat        result[N];
    cusfloat        result_dx0[N];
    cusfloat        result_dx1[N];

    // Define chebyshev evaluator functions
    using TEV       = ChebyshevEvaluatorBaseVector<ChebyshevTraits<T>, N, STATIC_LOOP_ON>;
    using T_dX0EV   = ChebyshevEvaluatorBaseVector<ChebyshevTraits<T_dX0>, N, STATIC_LOOP_ON>;
    using T_dX1EV   = ChebyshevEvaluatorBaseVector<ChebyshevTraits<T_dX1>, N, STATIC_LOOP_ON>;

    // Load data
    RefData<3, 3> ref_data( fipath );

    // Prepare data variables to global error control if any
    cusfloat    g_max_err       = 0.0;
    int         g_count_thr     = 0;
    cusfloat    g_dx0_max_err   = 0.0;
    int         g_dx0_count_thr = 0;
    cusfloat    g_dx1_max_err   = 0.0;
    int         g_dx1_count_thr = 0;

    // Loop over database to check its value
    for ( int i=0; i<ref_data.num_x2; i++ )
    {
        // std::cout << "X2: " << ref_data.X2[i] << std::endl;
        // Fold database
        fold_database_3d<T>( ref_data.X2[i] );
        fold_database_3d<T_dX0>( ref_data.X2[i] );
        fold_database_3d<T_dX1>( ref_data.X2[i] );

        for ( int j=0; j<ref_data.num_x0; j++ )
        {
            for ( int k=0; k<ref_data.num_x1; k++ )
            {
                // Interpolate in the database
                TEV::evaluate( N, &(ref_data.X0[j]), &(ref_data.X1[k]), result );
                T_dX0EV::evaluate( N, &(ref_data.X0[j]), &(ref_data.X1[k]), result_dx0 );
                T_dX1EV::evaluate( N, &(ref_data.X0[j]), &(ref_data.X1[k]), result_dx1 );

                // Compare result with the target
                index = ( 
                            i * ( ref_data.num_x0 * ref_data.num_x1 )
                            +
                            j * ( ref_data.num_x1 )
                            +
                            k
                        );
                // std::cout << "Index: " << index << std::endl;
                
                if ( !assert_scalar_equality( result[0], ref_data.G[index], EPS_PREC ) )
                {
                    if constexpr( EC == ErrorControl::Local )
                    {
                        std::cerr << "TEST FAILED for " << db_name << " and G function" << std::endl;
                        std::cerr << "Expected: " << ref_data.G[index] << " - Calculated: " << result[0] << std::endl;
                        std::cerr << "A[" << j  << "]: " << ref_data.X0[j] << " - B[" << k << "]: " << ref_data.X1[k] << " - H[" << i << "]: " << ref_data.X2[i] << std::endl;
                        throw std::exception( );
                    }
                    else if constexpr( EC == ErrorControl::Global )
                    {
                        UPDATE_GLOBAL_ERROR( result[0], ref_data.G[index], g_max_err, g_count_thr )
                    }
                }
                
                if ( !assert_scalar_equality( result_dx0[0], ref_data.G_dx0[index], EPS_PREC ) )
                {
                    if constexpr( EC == ErrorControl::Local )
                    {
                        std::cerr << "TEST FAILED for " << db_name << " and G_dA function" << std::endl;
                        std::cerr << "Expected: " << ref_data.G_dx0[index] << " - Calculated: " << result_dx0[0] << std::endl;
                        std::cerr << "A[" << j  << "]: " << ref_data.X0[j] << " - B[" << k << "]: " << ref_data.X1[k] << " - H[" << i << "]: " << ref_data.X2[i] << std::endl;
                        throw std::exception( );
                    }
                    else if constexpr( EC == ErrorControl::Global )
                    {
                        UPDATE_GLOBAL_ERROR( result_dx0[0], ref_data.G_dx0[index], g_dx0_max_err, g_dx0_count_thr )
                    }
                }
                
                if ( !assert_scalar_equality( result_dx1[0], ref_data.G_dx1[index], EPS_PREC ) )
                {
                    if constexpr( EC == ErrorControl::Local )
                    {
                        std::cerr << "TEST FAILED for " << db_name << " and G_dB function" << std::endl;
                        std::cerr << "Expected: " << ref_data.G_dx1[index] << " - Calculated: " << result_dx1[0] << std::endl;
                        std::cerr << "A[" << j  << "]: " << ref_data.X0[j] << " - B[" << k << "]: " << ref_data.X1[k] << " - H[" << i << "]: " << ref_data.X2[i] << std::endl;
                        throw std::exception( );
                    }
                    else if constexpr( EC == ErrorControl::Global )
                    {
                        UPDATE_GLOBAL_ERROR( result_dx1[0], ref_data.G_dx1[index], g_dx1_max_err, g_dx1_count_thr )
                    }
                }
            }
        }
    }

    // Check global error if selected
    if constexpr( EC == ErrorControl::Global )
    {
        int N2 = ref_data.num_x0 * ref_data.num_x1 * ref_data.num_x2;
        CHECK_GLOBAL_ERROR( EPS_PREC, N2, g_max_err,     g_count_thr,     db_name, "G"     )
        CHECK_GLOBAL_ERROR( EPS_PREC, N2, g_dx0_max_err, g_dx0_count_thr, db_name, "G_dx0" )
        CHECK_GLOBAL_ERROR( EPS_PREC, N2, g_dx1_max_err, g_dx1_count_thr, db_name, "G_dx1" )
    }

    std::cout << " -> DONE" << std::endl;
}


int main( int argc, char* argv[] )
{
    // Read command line arguments
    if ( !check_num_cmd_args( argc, 7 ) )
    {
        return 1;
    }

    std::string l1_fipath(argv[1]);
    std::string l2_fipath(argv[2]);
    std::string l3_fipath(argv[3]);
    std::string m1_fipath(argv[4]);
    std::string m2_fipath(argv[5]);
    std::string m3_fipath(argv[6]);
    std::string r11_fipath(argv[7]);
    
    // Launch database comparison
    compare_3d_database<L1C, L1_dAC, L1_dBC, ErrorControl::Local>( l1_fipath, "L1" );
    compare_1d_database<L2C, ErrorControl::Local>( l2_fipath, "L2" );
    compare_3d_database<L3C, L3_dAC, L3_dBC, ErrorControl::Local>( l3_fipath, "L3" );
    compare_3d_database<M1C, M1_dAC, M1_dBC, ErrorControl::Local>( m1_fipath, "M1" );
    compare_1d_database<M2C, ErrorControl::Local>( m2_fipath, "M2" );
    compare_3d_database<M3C, M3_dAC, M3_dBC, ErrorControl::Local>( m3_fipath, "M3" );
    compare_2d_database<R11C, R11_dXC, ErrorControl::Local>( r11_fipath, "R11" );

    return 0;
}
