
// Include local modules
#include "../../src/green/integrals_db.hpp"
#include "../../src/math/math_tools.hpp"


void    calculate_l1( 
                            IntegralsDb*    int_db,
                            cusfloat        A, 
                            cusfloat        B 
                    )
{
    int_db->l1->get_value_ab( A, B );
}


void    calculate_l1_da( 
                            IntegralsDb*    int_db,
                            cusfloat        A, 
                            cusfloat        B 
                        )
{
    int_db->l1_da->get_value_ab( A, B );
}


void    calculate_l1_db( 
                            IntegralsDb*    int_db,
                            cusfloat        A, 
                            cusfloat        B 
                        )
{
    int_db->l1_db->get_value_ab( A, B );
}


void    calculate_l2(
                            IntegralsDb*    int_db
                    )
{
    int_db->l2->int_1d;
}


void    calculate_l3(
                            IntegralsDb*    int_db,
                            cusfloat        A,
                            cusfloat        B
                    )
{
    int_db->l3->get_value_ab( A, B );
}


void    calculate_l3_da(
                            IntegralsDb*    int_db,
                            cusfloat        A,
                            cusfloat        B
                        )
{
    int_db->l3_da->get_value_ab( A, B );
}


void    calculate_l3_db(
                            IntegralsDb*    int_db,
                            cusfloat        A,
                            cusfloat        B
                        )
{
    int_db->l3_db->get_value_ab( A, B );
}


void    calculate_m1(
                            IntegralsDb*    int_db,
                            cusfloat        A,
                            cusfloat        B
                    )
{
    int_db->m1->get_value_ab( A, B );
}


void    calculate_m1_da(
                            IntegralsDb*    int_db,
                            cusfloat        A,
                            cusfloat        B
                        )
{
    int_db->m1_da->get_value_ab( A, B );
}


void    calculate_m1_db(
                            IntegralsDb*    int_db,
                            cusfloat        A,
                            cusfloat        B
                        )
{
    int_db->m1_db->get_value_ab( A, B );
}


void    calculate_m2(
                            IntegralsDb*    int_db
                    )
{
    int_db->m2->int_1d;
}


void    calculate_m3(
                            IntegralsDb*    int_db,
                            cusfloat        A,
                            cusfloat        B
                    )
{
    int_db->m3->get_value_ab( A, B );
}


void    calculate_m3_da(
                            IntegralsDb*    int_db,
                            cusfloat        A,
                            cusfloat        B
                        )
{
    int_db->m3_da->get_value_ab( A, B );
}


void    calculate_m3_db(
                            IntegralsDb*    int_db,
                            cusfloat        A,
                            cusfloat        B
                        )
{
    int_db->m3_db->get_value_ab( A, B );
}


int main( void )
{
    // Define integrals parameters
    cusfloat A              = 1.0;
    cusfloat ang_freq       = 2 * PI / 62.0;
    cusfloat B              = 1.0;
    cusfloat grav_acc       = 9.81;
    cusfloat water_depth    = 1000.0;

    // Create integrals database
    IntegralsDb* int_db = new IntegralsDb( );

    // Fold for current frequency and water depth
    cusfloat H = pow2s( ang_freq ) * water_depth / grav_acc;
    int_db->fold_h( H );

    // Iterative loop to have repetitions and 
    // to have an accurate measurement of the time
    for ( int i=0; i<1000000; i++ )
    {
        // std::cout << "Iteration: " << i << std::endl;
        // Calculate L1
        calculate_l1( int_db, A, B );

        // Calculate L1_dA
        calculate_l1_da( int_db, A, B );

        // Calculate L1_dB
        calculate_l1_db( int_db, A, B );

        // Calculate L2
        calculate_l2( int_db );

        // Calculate L3
        calculate_l3( int_db, A, B );

        // Calculate L3_dA
        calculate_l3_da( int_db, A, B );

        // Calculate L3_dB
        calculate_l3_db( int_db, A, B );

        // Calculate M1
        calculate_m1( int_db, A, B );

        // Calcualte M1_dA
        calculate_m1_da( int_db, A, B );

        // Calculate M1_dB
        calculate_m1_db( int_db, A, B );

        // Calculate M2
        calculate_m2( int_db );

        // Calculate M3
        calculate_m3( int_db, A, B );

        // Calculate M3_dA
        calculate_m3_da( int_db, A, B );

        // Calculate M3_dB
        calculate_m3_db( int_db, A, B );
    }

    return 0;
}
