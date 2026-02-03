
// Include general usage modules
#include <cstdlib>
#include <iostream>
#include <vector>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/math/cusvector.hpp"
#include "../../src/math/math_tools.hpp"


// Macro definition
#define IS_PASS_TEST( test_name, pass )                                                         \
    if ( !pass )                                                                                \
    {                                                                                           \
        std::cerr << "ERROR - Test for operator: " << test_name << " failed!" << std::endl;     \
        return 1;                                                                               \
    }                                                                                           \



class C
{
public:
    double a = 0.0;

    C( double a_in ): a( a_in )
    {
        std::cout << "Constructor a: " << this->a << std::endl; 
    };

    C( C& c_in )
    {
        std::cout << "Activating copy constructor" << std::endl;
        this->a = c_in.a;
        std::cout << "Done!" << std::endl; 
    }

    friend C& operator + ( const double& b, C& other )
    {
        std::cout << "Performing LHS operation..." << std::endl;
        other.a += b;
        std::cout << "Other: " << other.a << std::endl;
        std::cout << "Done!" << std::endl;

        return other;
    }
};


// C& operator + ( const double& b, C& other )
// {
//     std::cout << "Performing LHS operation..." << std::endl;
//     other.a += b;
//     std::cout << "Other: " << other.a << std::endl;
//     std::cout << "Done!" << std::endl;

//     return other;
// }


// Define auxiliar funcions
template<typename T, typename U>
void compare_results( T& vec, U& cusvec )
{
    bool pass = true;
    for ( int i=0; i<vec.size( ); i++ )
    {
        if (  )
    }
}


template<typename T, typename U>
void copy_vec_to_cusvec( T& vec, U& cusvec )
{
    for ( size_t i=0; i<vec.size( ); i++ )
    {
        cusvec[i] = vec[i];
    }
}


// Function test definition
bool test_div_lhs( void )
{
    // Create reference vectors
    int         N       = 100;
    cusfloat    s       = 2.365;
    std::vector<cusfloat> v0( N );
    std::vector<cusfloat> v2( N );
    
    for ( int i=0; i<N; i++ )
    {
        v0[i] = std::rand( );
        v2[i] = s / v0[i];
    }

    // Create Custom Vectors
    CusVector<cusfloat> cv0( N );

    copy_vec_to_cusvec( v0, cv0 );

    // Sum up Custom vectors
    CusVector<cusfloat> cv2 = s / cv0 ;

    // Compare Vector results with CusVector ones
    bool pass = assert_vector_equality( 
                                            v2.size( ) ,
                                            v2.data( ),
                                            cv2.data( ),
                                            EPS_PRECISION
                                        );
    std::cout << "pass: " << pass << std::endl;
    return pass;
}

bool test_div_rhs( void )
{
    // Create reference vectors
    int N = 100;
    std::vector<cusfloat> v0( N );
    std::vector<cusfloat> v1( N );
    std::vector<cusfloat> v2( N );
    
    for ( int i=0; i<N; i++ )
    {
        v0[i] = std::rand( );
        v1[i] = std::rand( );
        v2[i] = v0[i] / v1[i];
    }

    // Create Custom Vectors
    CusVector<cusfloat> cv0( N );
    CusVector<cusfloat> cv1( N );

    copy_vec_to_cusvec( v0, cv0 );
    copy_vec_to_cusvec( v1, cv1 );

    // Sum up Custom vectors
    cv0 /= cv1;

    // Compare Vector results with CusVector ones
    bool pass = assert_vector_equality( 
                                            v2.size( ) ,
                                            v2.data( ),
                                            cv0.data( ),
                                            EPS_PRECISION
                                        );

    return pass;
}


bool test_sub_lhs( void )
{
    // Create reference vectors
    int         N       = 100;
    cusfloat    s       = 2.365;
    std::vector<cusfloat> v0( N );
    std::vector<cusfloat> v2( N );
    
    for ( int i=0; i<N; i++ )
    {
        v0[i] = std::rand( );
        v2[i] = s - v0[i];
    }

    // Create Custom Vectors
    CusVector<cusfloat> cv0( N );

    copy_vec_to_cusvec( v0, cv0 );

    // Sum up Custom vectors
    CusVector<cusfloat> cv2 = s - cv0 ;

    // Compare Vector results with CusVector ones
    bool pass = assert_vector_equality( 
                                            v2.size( ) ,
                                            v2.data( ),
                                            cv2.data( ),
                                            EPS_PRECISION
                                        );

    return pass;
}


bool test_sub_rhs( void )
{
    // Create reference vectors
    int N = 100;
    std::vector<cusfloat> v0( N );
    std::vector<cusfloat> v1( N );
    std::vector<cusfloat> v2( N );
    
    for ( int i=0; i<N; i++ )
    {
        v0[i] = std::rand( );
        v1[i] = std::rand( );
        v2[i] = v0[i] - v1[i];
    }

    // Create Custom Vectors
    CusVector<cusfloat> cv0( N );
    CusVector<cusfloat> cv1( N );

    copy_vec_to_cusvec( v0, cv0 );
    copy_vec_to_cusvec( v1, cv1 );

    // Sum up Custom vectors
    cv0 -= cv1;

    // Compare Vector results with CusVector ones
    bool pass = assert_vector_equality( 
                                            v2.size( ) ,
                                            v2.data( ),
                                            cv0.data( ),
                                            EPS_PRECISION
                                        );

    return pass;
}


bool test_sum_lhs( void )
{
    // Create reference vectors
    int         N       = 100;
    cusfloat    s       = 2.365;
    std::vector<cusfloat> v0( N );
    std::vector<cusfloat> v2( N );
    
    for ( int i=0; i<N; i++ )
    {
        v0[i] = std::rand( );
        v2[i] = v0[i] + s;
    }

    // Create Custom Vectors
    CusVector<cusfloat> cv0( N );

    copy_vec_to_cusvec( v0, cv0 );

    // Sum up Custom vectors
    CusVector<cusfloat> cv2 = s + cv0 ;

    // Compare Vector results with CusVector ones
    bool pass = assert_vector_equality( 
                                            v2.size( ) ,
                                            v2.data( ),
                                            cv2.data( ),
                                            EPS_PRECISION
                                        );

    return pass;
}


bool test_sum_rhs( void )
{
    // Create reference vectors
    int N = 100;
    std::vector<cusfloat> v0( N );
    std::vector<cusfloat> v1( N );
    std::vector<cusfloat> v2( N );
    
    for ( int i=0; i<N; i++ )
    {
        v0[i] = std::rand( );
        v1[i] = std::rand( );
        v2[i] = v0[i] + v1[i];
    }

    // Create Custom Vectors
    CusVector<cusfloat> cv0( N );
    CusVector<cusfloat> cv1( N );

    copy_vec_to_cusvec( v0, cv0 );
    copy_vec_to_cusvec( v1, cv1 );

    // Sum up Custom vectors
    cv0 += cv1;

    // Compare Vector results with CusVector ones
    bool pass = assert_vector_equality( 
                                            v2.size( ) ,
                                            v2.data( ),
                                            cv0.data( ),
                                            EPS_PRECISION
                                        );

    return pass;
}


bool test_mult_lhs( void )
{
    // Create reference vectors
    int         N       = 100;
    cusfloat    s       = 2.365;
    std::vector<cusfloat> v0( N );
    std::vector<cusfloat> v2( N );
    
    for ( int i=0; i<N; i++ )
    {
        v0[i] = std::rand( );
        v2[i] = s * v0[i];
    }

    // Create Custom Vectors
    CusVector<cusfloat> cv0( N );

    copy_vec_to_cusvec( v0, cv0 );

    // Sum up Custom vectors
    CusVector<cusfloat> cv2 = s * cv0 ;

    // Compare Vector results with CusVector ones
    bool pass = assert_vector_equality( 
                                            v2.size( ) ,
                                            v2.data( ),
                                            cv2.data( ),
                                            EPS_PRECISION
                                        );

    return pass;
}


bool test_mult_rhs( void )
{
    // Create reference vectors
    int N = 100;
    std::vector<cusfloat> v0( N );
    std::vector<cusfloat> v1( N );
    std::vector<cusfloat> v2( N );
    
    for ( int i=0; i<N; i++ )
    {
        v0[i] = std::rand( );
        v1[i] = std::rand( );
        v2[i] = v0[i] * v1[i];
    }

    // Create Custom Vectors
    CusVector<cusfloat> cv0( N );
    CusVector<cusfloat> cv1( N );

    copy_vec_to_cusvec( v0, cv0 );
    copy_vec_to_cusvec( v1, cv1 );

    // Sum up Custom vectors
    cv0 *= cv1;

    // Compare Vector results with CusVector ones
    bool pass = assert_vector_equality( 
                                            v2.size( ) ,
                                            v2.data( ),
                                            cv0.data( ),
                                            EPS_PRECISION
                                        );

    return pass;
}


bool test_pow_rhs( void )
{
    // Create reference vectors
    int N = 100;
    std::vector<cusfloat> v0( N );
    std::vector<cusfloat> v1( N );
    std::vector<cusfloat> v2( N );
    
    for ( int i=0; i<N; i++ )
    {
        v0[i] = std::rand( );
        v1[i] = std::rand( );
        v2[i] = std::pow( v0[i], v1[i] );
    }

    // Create Custom Vectors
    CusVector<cusfloat> cv0( N );
    CusVector<cusfloat> cv1( N );

    copy_vec_to_cusvec( v0, cv0 );
    copy_vec_to_cusvec( v1, cv1 );

    // Sum up Custom vectors
    cv0 ^ cv1;

    // Compare Vector results with CusVector ones
    bool pass = assert_vector_equality( 
                                            v2.size( ) ,
                                            v2.data( ),
                                            cv0.data( ),
                                            EPS_PRECISION
                                        );

    return pass;
}


bool test_square_bracket( void )
{
    // Define test properties
    int N = 100;

    // Create Vector containers
    CusVector<cusfloat>     my_vector( N );
    std::vector<cusfloat>   vec( N );

    // Fill in vector containers
    for ( int i=0; i<N; i++ )
    {
        my_vector[i]    = i;
        vec[i]          = i;
    }

    for ( int i=0; i<N; i++ )
    {
        my_vector[i]    += i;
        vec[i]          += i;
    }

    // Check if the vectors are equal
    bool pass = true;
    for ( int i=0; i<N; i++ )
    {
        if ( std::abs( my_vector[i] - vec[i] ) > ZEROTH_EPS )
        {
            pass = false;
            break;
        }

    }

    return pass;
}


int main( void )
{
    IS_PASS_TEST( "DIV_LHS", test_div_lhs( ) );
    std::cout << "P1" << std::endl;
    IS_PASS_TEST( "DIV_RHS", test_div_rhs( ) );
    std::cout << "P2" << std::endl;
    IS_PASS_TEST( "SUB_LHS", test_sub_lhs( ) );
    std::cout << "P3" << std::endl;
    IS_PASS_TEST( "SUB_RHS", test_sub_rhs( ) );
    std::cout << "P4" << std::endl;
    IS_PASS_TEST( "SUM_LHS", test_sum_lhs( ) );
    std::cout << "P5" << std::endl;
    IS_PASS_TEST( "SUM_RHS", test_sum_rhs( ) );
    std::cout << "P6" << std::endl;
    IS_PASS_TEST( "MULT_LHS", test_mult_lhs( ) );
    std::cout << "P7" << std::endl;
    IS_PASS_TEST( "MULT_RHS", test_mult_rhs( ) );
    std::cout << "P8" << std::endl;
    IS_PASS_TEST( "POW_RHS", test_pow_rhs( ) );
    std::cout << "P9" << std::endl;
    IS_PASS_TEST( "SQUARE_BRACKET", test_square_bracket( ) );

    // cusfloat a = 2.0;
    // cusfloat b = 3.0;
    // CusVector<cusfloat> c0( 10 );
    // CusVector<cusfloat> c1 = ( a + c0 ) + b;

    // std::cout << c0[0] << " - " << c0[1] << " - " << c0[2] << std::endl;


    std::cout << "Test Completed!" << std::endl;

    return 0;
}