
// Include general usage libraries
#include <fstream>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/green/chebyshev_traits.hpp"
#include "../../src/green/chebyshev_evaluator_interface.hpp"
#include "../../src/green/integrals_database.hpp"
#include "../../src/green/fin_depth_coeffs/L1.hpp"


#define COND_LINE( expr, line )     \
{                                   \
    if constexpr( expr )            \
    {                               \
        line;                       \
    }                               \
}                                   \


cusfloat EPS_PREC = 1E-4;


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


template<typename T>
void compare_1d_database( std::string fipath, std::string db_name )
{
    std::cout << "COMPARE " << db_name << " DATABASE" << std::endl;
    // Define chebyshev evaluator functions
    using TEV       = ChebyshevTraits<T>;

    // Load data
    RefData<1, 1> ref_data( fipath );

    // Loop over database to check its value
    for ( int i=0; i<ref_data.num_x0; i++ )
    {
        // Fold database
        fold_database_1d<T>( ref_data.X0[i] );

        if ( !assert_scalar_equality( TEV::coeffs, ref_data.G[i], EPS_PREC ) )
        {
            std::cerr << "TEST FAILED for " << db_name << " and G function" << std::endl;
            std::cerr << "Expected: " << ref_data.G[i] << " - Calculated: " << TEV::coeffs << std::endl;
            std::cerr << "H[" << i << "]: " << ref_data.X0[i] << std::endl;
            throw std::exception( );
        }
    }
    std::cout << " -> DONE" << std::endl;
}


template<typename T, typename T_dX0>
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

            // result[0] *= 2.0;
            // result_dx0[0] *= 2.0;
            // std::cout << "X: " << ref_data.X0[i] << " - Y: " << ref_data.X1[j] << std::endl;
            // std::cout << "\t->result: " << result[0] << " - " << ref_data.G[index] << " - " << result[0]-ref_data.G[index] << std::endl;
            // std::cout << "\t->result_dx0: " << result_dx0[0] << " - " << ref_data.G_dx0[index] << " - " << result_dx0[0]-ref_data.G_dx0[index] << std::endl;

            if ( !assert_scalar_equality( result[0], ref_data.G[index], EPS_PREC ) )
            {
                std::cerr << "TEST FAILED for " << db_name << " and G function" << std::endl;
                std::cerr << "Expected: " << ref_data.G[index] << " - Calculated: " << result[0] << std::endl;
                std::cerr << "A[" << i  << "]: " << ref_data.X0[i] << " - B[" << j << "]: " << ref_data.X1[j] << std::endl;
                throw std::exception( );
            }
            
            if ( !assert_scalar_equality( result_dx0[0], ref_data.G_dx0[index], EPS_PREC ) )
            {
                std::cerr << "TEST FAILED for " << db_name << " and G_dA function" << std::endl;
                std::cerr << "Expected: " << ref_data.G_dx0[index] << " - Calculated: " << result_dx0[0] << std::endl;
                std::cerr << "A[" << i  << "]: " << ref_data.X0[i] << " - B[" << j << "]: " << ref_data.X1[j] << " - H[" << j << "]: " << std::endl;
                throw std::exception( );
            }
        }

    }
    std::cout << " -> DONE" << std::endl;
}


template<typename T, typename T_dX0, typename T_dX1>
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
                    std::cerr << "TEST FAILED for " << db_name << " and G function" << std::endl;
                    std::cerr << "Expected: " << ref_data.G[index] << " - Calculated: " << result[0] << std::endl;
                    std::cerr << "A[" << j  << "]: " << ref_data.X0[j] << " - B[" << k << "]: " << ref_data.X1[k] << " - H[" << i << "]: " << ref_data.X2[i] << std::endl;
                    throw std::exception( );
                }
                
                if ( !assert_scalar_equality( result_dx0[0], ref_data.G_dx0[index], EPS_PREC ) )
                {
                    std::cerr << "TEST FAILED for " << db_name << " and G_dA function" << std::endl;
                    std::cerr << "Expected: " << ref_data.G_dx0[index] << " - Calculated: " << result_dx0[0] << std::endl;
                    std::cerr << "A[" << j  << "]: " << ref_data.X0[j] << " - B[" << k << "]: " << ref_data.X1[k] << " - H[" << i << "]: " << ref_data.X2[i] << std::endl;
                    throw std::exception( );
                }
                
                if ( !assert_scalar_equality( result_dx1[0], ref_data.G_dx1[index], EPS_PREC ) )
                {
                    std::cerr << "TEST FAILED for " << db_name << " and G_dB function" << std::endl;
                    std::cerr << "Expected: " << ref_data.G_dx1[index] << " - Calculated: " << result_dx1[0] << std::endl;
                    std::cerr << "A[" << j  << "]: " << ref_data.X0[j] << " - B[" << k << "]: " << ref_data.X1[k] << " - H[" << i << "]: " << ref_data.X2[i] << std::endl;
                    throw std::exception( );
                }
            }
        }
    }
    std::cout << " -> DONE" << std::endl;
}


int main( void )
{

    std::string l1_fipath   = "D:/sergio/developments/SeaMotions/aux_tools/0_databases/L1.dat";
    std::string l2_fipath   = "D:/sergio/developments/SeaMotions/aux_tools/0_databases/L2.dat";
    std::string l3_fipath   = "D:/sergio/developments/SeaMotions/aux_tools/0_databases/L3.dat";
    std::string m1_fipath   = "D:/sergio/developments/SeaMotions/aux_tools/0_databases/M1.dat";
    std::string m2_fipath   = "D:/sergio/developments/SeaMotions/aux_tools/0_databases/M2.dat";
    std::string m3_fipath   = "D:/sergio/developments/SeaMotions/aux_tools/0_databases/M3.dat";
    std::string r11_fipath  = "D:/sergio/developments/SeaMotions/aux_tools/0_databases/R11.dat";
    // // std::string l1_fipath = "C:/Users/feruanos/Downloads/GreenFunction-master (1)/GreenFunction-master/L1.dat";
    // std::string l2_fipath = "C:/Users/feruanos/Downloads/GreenFunction-master (1)/GreenFunction-master/L2.dat";
    // std::string l3_fipath = "C:/Users/feruanos/Downloads/GreenFunction-master (1)/GreenFunction-master/L3.dat";
    // std::string m1_fipath = "C:/Users/feruanos/Downloads/GreenFunction-master (1)/GreenFunction-master/M1.dat";
    // std::string m2_fipath = "C:/Users/feruanos/Downloads/GreenFunction-master (1)/GreenFunction-master/M2.dat";
    // std::string m3_fipath = "C:/Users/feruanos/Downloads/GreenFunction-master (1)/GreenFunction-master/M3.dat";
    
    compare_3d_database<L1C, L1_dAC, L1_dBC>( l1_fipath, "L1" );
    compare_1d_database<L2C>( l2_fipath, "L2" );
    compare_3d_database<L3C, L3_dAC, L3_dBC>( l3_fipath, "L3" );
    compare_3d_database<M1C, M1_dAC, M1_dBC>( m1_fipath, "M1" );
    compare_1d_database<M2C>( m2_fipath, "M2" );
    compare_3d_database<M3C, M3_dAC, M3_dBC>( m3_fipath, "M3" );
    compare_2d_database<R11C, R11_dXC>( r11_fipath, "R11" );

    return 0;
}


// #include <iostream>
// // #include "base_struct.hpp"
// struct MyStruct
// {
//     static double my_vec[10];
// };

// double MyStruct::my_vec[10];

// template<typename T>
// struct MyStruct2
// {
//     inline static double* my_vec = T::my_vec;
// };

// using MS2 = MyStruct2<MyStruct>;


// int main( void )
// {
//     MS2::my_vec[0] = 0.0;

//     std::cout << "Program finished!" << std::endl;

//     return 0;
// }