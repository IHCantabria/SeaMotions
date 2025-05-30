
#pragma once

// Include local modules
#include "../config.hpp"
#include "math_tools.hpp"


// Define powers size
constexpr std::size_t   NPOWER = 14;


class BesselFactory
{
private:
    // Define private attributes
    cusfloat _f0                = 0.0;
    cusfloat _f1                = 0.0;
    cusfloat _sf0               = 0.0;
    cusfloat _sf1               = 0.0;
    cusfloat _th0               = 0.0;
    cusfloat _th1               = 0.0;
    cusfloat _powers[NPOWER];

    // Define private methods
    void _calculate_j0_ascending( void )
    {
        this->j0 = 0.999999999;
        this->j0 -= 2.249999879 * this->_powers[2];
        this->j0 += 1.265623060 * this->_powers[4];
        this->j0 -= 0.316394552 * this->_powers[6];
        this->j0 += 0.044460948 * this->_powers[8];
        this->j0 -= 0.003954479 * this->_powers[10];
        this->j0 += 0.000212950 * this->_powers[12];
    }

    void _calculate_j0_inverse( const cusfloat x )
    {
        this->j0 = 1 / std::sqrt( x ) * this->_f0 * std::cos( this->_th0 );
    }

    void _calculate_j1_ascending( const cusfloat x )
    {
        this->j1 = 0.500000000;
        this->j1 -= 0.562499992 * this->_powers[2];
        this->j1 += 0.210937377 * this->_powers[4];
        this->j1 -= 0.039550040 * this->_powers[6];
        this->j1 += 0.004447331 * this->_powers[8];
        this->j1 -= 0.000330547 * this->_powers[10];
        this->j1 += 0.000015525 * this->_powers[12];
        this->j1 *= x;
    }

    void _calculate_j1_inverse( const cusfloat x )
    {
        this->j1 = 1 / std::sqrt( x ) * this->_f1 * std::cos( this->_th1 );
    }

    void _calculate_y0_ascending( const cusfloat x )
    {
        this->y0 = ( 2.0 / PI ) * std::log( x / 2.0 ) * this->j0;
        this->y0 += 0.367466907;
        this->y0 += 0.605593797 * this->_powers[2];
        this->y0 -= 0.743505078 * this->_powers[4];
        this->y0 += 0.253005481 * this->_powers[6];
        this->y0 -= 0.042619616 * this->_powers[8];
        this->y0 += 0.004285691 * this->_powers[10];
        this->y0 -= 0.000250716 * this->_powers[12];
    }

    void _calculate_y0_inverse( const cusfloat x )
    {
        this->y0 = 1 / std::sqrt( x ) * this->_f0 * std::sin( this->_th0 );
    }

    void _calculate_y1_ascending( const cusfloat x )
    {
        this->y1 = ( 2.0 / PI ) * ( std::log( x / 2.0 ) * this->j1 - 1 / x );
        this->y1 += 0.073735531 * this->_powers[1];
        this->y1 += 0.722769344 * this->_powers[3];
        this->y1 -= 0.438896337 * this->_powers[5];
        this->y1 += 0.104320251 * this->_powers[7];
        this->y1 -= 0.013637596 * this->_powers[9];
        this->y1 += 0.001125970 * this->_powers[11];
        this->y1 -= 0.000056455 * this->_powers[13];
    }

    void _calculate_y1_inverse( const cusfloat x )
    {
        this->y1 = 1 / std::sqrt( x ) * this->_f1 * std::sin( this->_th1 );
    }

    void _calculate_struve0_ascending( void )
    {
        this->struve0 = 1.909859164  * this->_powers[1];
        this->struve0 -= 1.909855001 * this->_powers[3];
        this->struve0 += 0.687514637 * this->_powers[5];
        this->struve0 -= 0.126164557 * this->_powers[7];
        this->struve0 += 0.013828813 * this->_powers[9];
        this->struve0 -= 0.000876918 * this->_powers[11];
    }

    void _calculate_struve0_inverse( const cusfloat x )
    {
        this->struve0 = this->y0 + this->_sf0;
    }

    void _calculate_struve1_ascending( void )
    {
        this->struve1 = 1.909859286  * this->_powers[2];
        this->struve1 -= 1.145914713 * this->_powers[4];
        this->struve1 += 0.294656958 * this->_powers[6];
        this->struve1 -= 0.042070508 * this->_powers[8];
        this->struve1 += 0.003785727 * this->_powers[10];
        this->struve1 -= 0.000207183 * this->_powers[12];
    }

    void _calculate_struve1_inverse( const cusfloat x )
    {
        this->struve1 = this->y1 + this->_sf1;
    }

    void _calculate_polynomial_f0_th0( const cusfloat x )
    {
        this->_f0 = 0.79788454;
        this->_f0 -= 0.00553897 * this->_powers[2];
        this->_f0 += 0.00099336 * this->_powers[4];
        this->_f0 -= 0.00044346 * this->_powers[6];
        this->_f0 += 0.00020445 * this->_powers[8];
        this->_f0 -= 0.00004959 * this->_powers[10];

        this->_th0 = x - PI/4.0;
        this->_th0 -= 0.04166592 * this->_powers[1];
        this->_th0 += 0.00239399 * this->_powers[3];
        this->_th0 -= 0.00073984 * this->_powers[5];
        this->_th0 += 0.00031099 * this->_powers[7];
        this->_th0 -= 0.00007605 * this->_powers[9];
    }

    void _calculate_polynomial_f1_th1( const cusfloat x )
    {
        this->_f1 += 0.79788459;
        this->_f1 += 0.01662008 * this->_powers[2];
        this->_f1 -= 0.00187002 * this->_powers[4];
        this->_f1 += 0.00068519 * this->_powers[6];
        this->_f1 -= 0.00029440 * this->_powers[8];
        this->_f1 += 0.00006952 * this->_powers[10];

        this->_th1 += x - 3.0*PI/4.0;
        this->_th1 += 0.12499895 * this->_powers[1];
        this->_th1 -= 0.00605240 * this->_powers[3];
        this->_th1 += 0.00135825 * this->_powers[5];
        this->_th1 -= 0.00049616 * this->_powers[7];
        this->_th1 += 0.00011531 * this->_powers[9];
    }

    void _calculate_ascending_powers( const cusfloat x )
    {
        cusfloat pc = x / 3.0;

        this->_powers[0] = 1.0;
        this->_powers[1] = pc;

        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            this->_powers[i] = this->_powers[i-1] * pc;
        }
    }

    void _calculate_inverse_powers( const cusfloat x )
    {
        cusfloat pc = 3.0 / x;

        this->_powers[0] = 1.0;
        this->_powers[1] = pc;

        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            this->_powers[i] = this->_powers[i-1] * pc;
        }
    }

    void _calculate_rational_fraction_struve0( const cusfloat x )
    {
        // Define local constans
        cusfloat a0 = 0.99999906;
        cusfloat a1 = 4.77228920;
        cusfloat a2 = 3.85542044;
        cusfloat a3 = 0.32303607;

        cusfloat b1 = 4.88331068;
        cusfloat b2 = 4.28957333;
        cusfloat b3 = 0.52120508;

        // Compute struve factor
        cusfloat c0 = 2 * ( a0 + a1 * this->_powers[2] +a2 * this->_powers[4] + a3 * this->_powers[6] );
        cusfloat c1 = PI * x * ( 1 + b1 * this->_powers[2] + b2 * this->_powers[4] + b3 * this->_powers[6] );
        this->_sf0  = c0 / c1;

    }

    void _calculate_rational_fraction_struve1( const cusfloat x )
    {
        // Define local constans
        cusfloat a0 = 1.00000004;
        cusfloat a1 = 3.92205313;
        cusfloat a2 = 2.64893033;
        cusfloat a3 = 0.27450895;

        cusfloat b1 = 3.81095112;
        cusfloat b2 = 2.26216956;
        cusfloat b3 = 0.10885141;

        // Compute struve factor
        cusfloat c0 = 2 * ( a0 + a1 * this->_powers[2] + a2 * this->_powers[4] + a3 * this->_powers[6] );
        cusfloat c1 = PI * ( 1 + b1 * this->_powers[2] + b2 * this->_powers[4] + b3 * this->_powers[6] );
        this->_sf1  = c0 / c1;

    }

public:
    // Define class public attributes
    cusfloat j0         = 0.0;
    cusfloat j1         = 0.0;
    cusfloat y0         = 0.0;
    cusfloat y1         = 0.0;
    cusfloat struve0    = 0.0;
    cusfloat struve1    = 0.0;

    // Define class constructors
    BesselFactory( void) = default;

    // Define class methods
    void calculate_cheby( const cusfloat x )
    {
        
        // Check branch selection 
        if ( x < 3 )
        {
            // Calculate powers
            this->_calculate_ascending_powers( x );

            // Calculate first kind bessel functions
            this->_calculate_j0_ascending( );
            this->_calculate_j1_ascending( x );

            // Calculate second king bessel functions
            this->_calculate_y0_ascending( x );
            this->_calculate_y1_ascending( x );

            // Calculate struve functions
            this->_calculate_struve0_ascending( );
            this->_calculate_struve1_ascending( );
        }
        else
        {
            // Calculate powers
            this->_calculate_inverse_powers( x );

            // Calculate polynomial expansions
            this->_calculate_polynomial_f0_th0( x );
            this->_calculate_polynomial_f1_th1( x );
            this->_calculate_rational_fraction_struve0( x );
            this->_calculate_rational_fraction_struve1( x );
            
            // Calculate first kind bessel functions
            this->_calculate_j0_inverse( x );
            this->_calculate_j1_inverse( x );

            // Calculate second king bessel functions
            this->_calculate_y0_inverse( x );
            this->_calculate_y1_inverse( x );

            // Calculate struve functions
            this->_calculate_struve0_inverse( x );
            this->_calculate_struve0_inverse( x );

        }

    }

};