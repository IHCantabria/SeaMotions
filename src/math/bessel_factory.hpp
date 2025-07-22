
#pragma once


// Include local modules
#include "../config.hpp"
#include "math_tools.hpp"
#include "math_interface.hpp"


// Define powers size
constexpr std::size_t   NPOWER  = 14;
constexpr cusfloat      BRST    = 3.0;
constexpr cusfloat      BRMI    = 3.75;
constexpr cusfloat      BRMK    = 2.0;


// Define macros to evaluate ascending and descending masks
#define _ASDMASK( x ) ( 1 + x ) / 2.0
#define _INVMASK( x ) ( 1 - x ) / 2.0


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
    cusfloat _powers_modi[NPOWER];
    cusfloat _powers_modk[NPOWER];

    // Define private methods
    void _calculate_bessel_standard( const cusfloat x )
    {
        // Check branch selection 
        if ( x < BRST )
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
            this->_calculate_rational_fraction_struve1( );
            
            // Calculate first kind bessel functions
            this->_calculate_j0_inverse( x );
            this->_calculate_j1_inverse( x );

            // Calculate second king bessel functions
            this->_calculate_y0_inverse( x );
            this->_calculate_y1_inverse( x );

            // Calculate struve functions
            this->_calculate_struve0_inverse( );
            this->_calculate_struve1_inverse( );

        }

    }

    void _calculate_bessel_modified( const cusfloat x )
    {
        // Calculate first order modified
        if ( x < BRMI )
        {
            // Calculate ascending powers
            this->_calculate_ascending_powers_modi( x );

            // Calculate function values
            this->_calculate_i0_ascending( );
            this->_calculate_i1_ascending( x );

        }
        else
        {
            // Calculate inverse powers
            this->_calculate_inverse_powers_modi( x );

            // Calculate function values
            this->_calculate_i0_inverse( x );
            this->_calculate_i1_inverse( x );

        }

        // Calculate second order bessel modified
        if ( x < BRMK )
        {
            // Calculate ascending powers
            this->_calculate_ascending_powers_modk( x );

            // Calculate function values
            this->_calculate_k0_ascending( );
            this->_calculate_k1_ascending( x );

        }
        else
        {
            // Calculate inverse powers
            this->_calculate_inverse_powers_modk( x );

            // Calculate function values
            this->_calculate_k0_inverse( x );
            this->_calculate_k1_inverse( x );

        }

    }

    void _calculate_i0_ascending( void )
    {
        this->i0 = 1.0;
        this->i0 += 3.5156229 * this->_powers_modi[2];
        this->i0 += 3.0899424 * this->_powers_modi[4];
        this->i0 += 1.2067492 * this->_powers_modi[6];
        this->i0 += 0.2659732 * this->_powers_modi[8];
        this->i0 += 0.0360768 * this->_powers_modi[10];
        this->i0 += 0.0045813 * this->_powers_modi[12];
    }

    void _calculate_i0_inverse( const cusfloat x )
    {
        this->i0 = 0.39894228;
        this->i0 += 0.01328592 * this->_powers_modi[1];
        this->i0 += 0.00225319 * this->_powers_modi[2];
        this->i0 -= 0.00157565 * this->_powers_modi[3];
        this->i0 += 0.00916281 * this->_powers_modi[4];
        this->i0 -= 0.02057706 * this->_powers_modi[5];
        this->i0 += 0.02635537 * this->_powers_modi[6];
        this->i0 -= 0.01647633 * this->_powers_modi[7];
        this->i0 += 0.00392377 * this->_powers_modi[8];
        this->i0 *= std::exp( x ) / std::sqrt( x );
    }

    void _calculate_i1_ascending( const cusfloat x )
    {
        this->i1 = 0.5;
        this->i1 += 0.87890594 * this->_powers_modi[2];
        this->i1 += 0.51498869 * this->_powers_modi[4];
        this->i1 += 0.15084934 * this->_powers_modi[6];
        this->i1 += 0.02658733 * this->_powers_modi[8];
        this->i1 += 0.00301532 * this->_powers_modi[10];
        this->i1 += 0.00032411 * this->_powers_modi[12];
        this->i1 *= x;
    }

    void _calculate_i1_inverse( const cusfloat x )
    {
        this->i1 = 0.39894228;
        this->i1 -= 0.03988024 * this->_powers_modi[1];
        this->i1 -= 0.00362018 * this->_powers_modi[2];
        this->i1 += 0.00163801 * this->_powers_modi[3];
        this->i1 -= 0.01031555 * this->_powers_modi[4];
        this->i1 += 0.02282967 * this->_powers_modi[5];
        this->i1 -= 0.02895312 * this->_powers_modi[6];
        this->i1 += 0.01787654 * this->_powers_modi[7];
        this->i1 -= 0.00420059 * this->_powers_modi[8];
        this->i1 *= std::exp( x ) / std::sqrt( x );
    }

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

    void _calculate_k0_ascending( void )
    {
        this->k0 = -std::log( this->_powers_modk[1] ) * this->i0;
        this->k0 -= 0.57721566;
        this->k0 += 0.42278420 * this->_powers_modk[2];
        this->k0 += 0.23069756 * this->_powers_modk[4];
        this->k0 += 0.03488590 * this->_powers_modk[6];
        this->k0 += 0.00262698 * this->_powers_modk[8];
        this->k0 += 0.00010750 * this->_powers_modk[10];
        this->k0 += 0.00000740 * this->_powers_modk[12];
    }

    void _calculate_k0_inverse( const cusfloat x )
    {
        this->k0 = 1.25331414;
        this->k0 -= 0.07832358 * this->_powers_modk[1];
        this->k0 += 0.02189568 * this->_powers_modk[2];
        this->k0 -= 0.01062446 * this->_powers_modk[3];
        this->k0 += 0.00587872 * this->_powers_modk[4];
        this->k0 -= 0.00251540 * this->_powers_modk[5];
        this->k0 += 0.00053208 * this->_powers_modk[6];
        this->k0 *= 1.0 / ( std::exp( x ) * std::sqrt( x ) );
    }

    void _calculate_k1_ascending( const cusfloat x )
    {
        this->k1 = x * std::log( this->_powers_modk[1] ) * this->i1 + 1.0;
        this->k1 += 0.15443144 * this->_powers_modk[2];
        this->k1 -= 0.67278579 * this->_powers_modk[4];
        this->k1 -= 0.18156897 * this->_powers_modk[6];
        this->k1 -= 0.01919402 * this->_powers_modk[8];
        this->k1 -= 0.00110404 * this->_powers_modk[10];
        this->k1 -= 0.00004686 * this->_powers_modk[12];
        this->k1 /= x;
    }

    void _calculate_k1_inverse( const cusfloat x )
    {
        this->k1 = 1.25331414;
        this->k1 += 0.23498619 * this->_powers_modk[1];
        this->k1 -= 0.03655620 * this->_powers_modk[2];
        this->k1 += 0.01504268 * this->_powers_modk[3];
        this->k1 -= 0.00780353 * this->_powers_modk[4];
        this->k1 += 0.00325614 * this->_powers_modk[5];
        this->k1 -= 0.00068245 * this->_powers_modk[6];
        this->k1 /= std::exp( x ) * std::sqrt( x );
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

    void _calculate_struve0_inverse( void )
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

    void _calculate_struve1_inverse( void )
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
        this->_f1 = 0.79788459;
        this->_f1 += 0.01662008 * this->_powers[2];
        this->_f1 -= 0.00187002 * this->_powers[4];
        this->_f1 += 0.00068519 * this->_powers[6];
        this->_f1 -= 0.00029440 * this->_powers[8];
        this->_f1 += 0.00006952 * this->_powers[10];

        this->_th1 = x - 3.0*PI/4.0;
        this->_th1 += 0.12499895 * this->_powers[1];
        this->_th1 -= 0.00605240 * this->_powers[3];
        this->_th1 += 0.00135825 * this->_powers[5];
        this->_th1 -= 0.00049616 * this->_powers[7];
        this->_th1 += 0.00011531 * this->_powers[9];
    }

    void _calculate_ascending_powers( const cusfloat x )
    {
        cusfloat pc = x / BRST;

        this->_powers[0] = 1.0;
        this->_powers[1] = pc;

        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            this->_powers[i] = this->_powers[i-1] * pc;
        }
    }

    void _calculate_ascending_powers_modi( const cusfloat x )
    {
        cusfloat pc = x / BRMI;

        this->_powers_modi[0] = 1.0;
        this->_powers_modi[1] = pc;

        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            this->_powers_modi[i] = this->_powers_modi[i-1] * pc;
        }
    }

    void _calculate_ascending_powers_modk( const cusfloat x )
    {
        cusfloat pc = x / BRMK;

        this->_powers_modk[0] = 1.0;
        this->_powers_modk[1] = pc;

        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            this->_powers_modk[i] = this->_powers_modk[i-1] * pc;
        }
    }

    void _calculate_inverse_powers( const cusfloat x )
    {
        cusfloat pc = BRST / x;

        this->_powers[0] = 1.0;
        this->_powers[1] = pc;

        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            this->_powers[i] = this->_powers[i-1] * pc;
        }
    }

    void _calculate_inverse_powers_modi( const cusfloat x )
    {
        cusfloat pc = BRMI / x;

        this->_powers_modi[0] = 1.0;
        this->_powers_modi[1] = pc;

        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            this->_powers_modi[i] = this->_powers_modi[i-1] * pc;
        }
    }
    
    void _calculate_inverse_powers_modk( const cusfloat x )
    {
        cusfloat pc = BRMK / x;

        this->_powers_modk[0] = 1.0;
        this->_powers_modk[1] = pc;

        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            this->_powers_modk[i] = this->_powers_modk[i-1] * pc;
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

    void _calculate_rational_fraction_struve1( void )
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
    cusfloat i0         = 0.0;
    cusfloat i1         = 0.0;
    cusfloat j0         = 0.0;
    cusfloat j1         = 0.0;
    cusfloat k0         = 0.0;
    cusfloat k1         = 0.0;
    cusfloat y0         = 0.0;
    cusfloat y1         = 0.0;
    cusfloat struve0    = 0.0;
    cusfloat struve1    = 0.0;

    // Define class constructors
    BesselFactory( void ) = default;

    // Define class methods
    void calculate_cheby( const cusfloat x )
    {
        this->_calculate_bessel_standard( x );
    }

    void calculate_series( const cusfloat x )
    {
        this->_calculate_bessel_standard( x );
        this->_calculate_bessel_modified( x );
    }

};


template<int N>
class BesselFactoryVec
{
private:
    // Define private attributes
    MEMALINGR   cusfloat    _f0[N];
    MEMALINGR   cusfloat    _f1[N];
    MEMALINGR   cusfloat    _fcn_costh0[N];
    MEMALINGR   cusfloat    _fcn_costh1[N];
    MEMALINGR   cusfloat    _fcn_expx[N];
    MEMALINGR   cusfloat    _fcn_logpowers[N];
    MEMALINGR   cusfloat    _fcn_logx2[N];
    MEMALINGR   cusfloat    _fcn_sinth0[N];
    MEMALINGR   cusfloat    _fcn_sinth1[N];
    MEMALINGR   cusfloat    _fcn_sqrtx_inv[N];
    MEMALINGR   cusfloat    _sf0[N];
    MEMALINGR   cusfloat    _sf1[N];
    MEMALINGR   cusfloat    _th0[N];
    MEMALINGR   cusfloat    _th1[N];
    MEMALINGR   cusfloat    _mask_std[N];
    MEMALINGR   cusfloat    _mask_modi[N];
    MEMALINGR   cusfloat    _mask_modk[N];
                int         _na = 0;
                int         _ni = 0;
    MEMALINGR   cusfloat    _powers[N*NPOWER];
    MEMALINGR   cusfloat    _powers_inv[N*NPOWER];
    MEMALINGR   cusfloat    _powers_modi[N*NPOWER];
    MEMALINGR   cusfloat    _powers_modi_inv[N*NPOWER];
    MEMALINGR   cusfloat    _powers_modk[N*NPOWER];
    MEMALINGR   cusfloat    _powers_modk_inv[N*NPOWER];
    MEMALINGR   cusfloat    _xd2[N];
    MEMALINGR   cusfloat    _x_inv[N];

    // Define private methods
    void _calculate_bessel_standard( cusfloat* x )
    {
        // Calculate base powers
        this->_calculate_ascending_powers( x );
        this->_calculate_inverse_powers( x );

        // Calculate polynomial expansions
        this->_calculate_polynomial_f0_th0( x );
        this->_calculate_polynomial_f1_th1( x );
        this->_calculate_rational_fraction_struve0( x );
        this->_calculate_rational_fraction_struve1( );

        // Calculate standard auxiliar functions
        this->_calculate_standard_aux_fcns( x );

        // Calculate first kind bessel functions
        this->_calculate_j0( x );
        this->_calculate_j1( x );

        // Calculate second king bessel functions
        this->_calculate_y0( x );
        this->_calculate_y1( x );

        // Calculate struve functions
        this->_calculate_struve0( x );
        this->_calculate_struve1( x );

    }

    void _calculate_bessel_modified( cusfloat* x )
    {
        // Calculate powers
        this->_calculate_ascending_powers_modi( x );
        this->_calculate_ascending_powers_modk( x );
        this->_calculate_inverse_powers_modi( x );
        this->_calculate_inverse_powers_modk( x );
        
        // Calculate standard auxiliar functions
        this->_calculate_modified_aux_fcns( x );

        // Calculate first order modified
        this->_calculate_i0( x );
        this->_calculate_i1( x );

        this->_calculate_k0( x );
        this->_calculate_k1( x );

    }

    void _calculate_i0( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->i0[i] =   _ASDMASK( this->_mask_modi[i] )
                            *
                            (
                                + 1.0       * this->_powers_modi[0*N+i]
                                + 3.5156229 * this->_powers_modi[2*N+i]
                                + 3.0899424 * this->_powers_modi[4*N+i]
                                + 1.2067492 * this->_powers_modi[6*N+i]
                                + 0.2659732 * this->_powers_modi[8*N+i]
                                + 0.0360768 * this->_powers_modi[10*N+i]
                                + 0.0045813 * this->_powers_modi[12*N+i]
                            )
                            +
                            _INVMASK( this->_mask_modi[i] )
                            *
                            (
                                + 0.39894228 * this->_powers_modi_inv[0*N+i]
                                + 0.01328592 * this->_powers_modi_inv[1*N+i]
                                + 0.00225319 * this->_powers_modi_inv[2*N+i]
                                - 0.00157565 * this->_powers_modi_inv[3*N+i]
                                + 0.00916281 * this->_powers_modi_inv[4*N+i]
                                - 0.02057706 * this->_powers_modi_inv[5*N+i]
                                + 0.02635537 * this->_powers_modi_inv[6*N+i]
                                - 0.01647633 * this->_powers_modi_inv[7*N+i]
                                + 0.00392377 * this->_powers_modi_inv[8*N+i]
                            ) * this->_fcn_expx[i] * this->_fcn_sqrtx_inv[i];

        }
    }

    void _calculate_i1( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->i1[i] =   _ASDMASK( this->_mask_modi[i] )
                            *
                            (
                                + 0.5        * this->_powers_modi[0*N+i]
                                + 0.87890594 * this->_powers_modi[2*N+i]
                                + 0.51498869 * this->_powers_modi[4*N+i]
                                + 0.15084934 * this->_powers_modi[6*N+i]
                                + 0.02658733 * this->_powers_modi[8*N+i]
                                + 0.00301532 * this->_powers_modi[10*N+i]
                                + 0.00032411 * this->_powers_modi[12*N+i]
                            ) * x[i]
                            +
                            _INVMASK( this->_mask_modi[i] )
                            *
                            (
                                + 0.39894228 * this->_powers_modi_inv[0*N+i]
                                - 0.03988024 * this->_powers_modi_inv[1*N+i]
                                - 0.00362018 * this->_powers_modi_inv[2*N+i]
                                + 0.00163801 * this->_powers_modi_inv[3*N+i]
                                - 0.01031555 * this->_powers_modi_inv[4*N+i]
                                + 0.02282967 * this->_powers_modi_inv[5*N+i]
                                - 0.02895312 * this->_powers_modi_inv[6*N+i]
                                + 0.01787654 * this->_powers_modi_inv[7*N+i]
                                - 0.00420059 * this->_powers_modi_inv[8*N+i]
                                
                            ) * this->_fcn_expx[i] * this->_fcn_sqrtx_inv[i];

        }
    }

    void _calculate_j0( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->j0[i] =  _ASDMASK( this->_mask_std[i] )
                            * 
                            (
                                + 0.999999999 * this->_powers[0*N+i]
                                - 2.249999879 * this->_powers[2*N+i]
                                + 1.265623060 * this->_powers[4*N+i]
                                - 0.316394552 * this->_powers[6*N+i]
                                + 0.044460948 * this->_powers[8*N+i]
                                - 0.003954479 * this->_powers[10*N+i]
                                + 0.000212950 * this->_powers[12*N+i]
                            )
                            +
                            _INVMASK( this->_mask_std[i] )
                            *
                            (
                                this->_f0[i] * this->_fcn_costh0[i] * this->_fcn_sqrtx_inv[i]
                            );

        }
    }

    void _calculate_j1( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->j1[i]     =   _ASDMASK( this->_mask_std[i] )
                                *
                                (
                                    + 0.500000000 * this->_powers[0*N+i]
                                    - 0.562499992 * this->_powers[2*N+i]
                                    + 0.210937377 * this->_powers[4*N+i]
                                    - 0.039550040 * this->_powers[6*N+i]
                                    + 0.004447331 * this->_powers[8*N+i]
                                    - 0.000330547 * this->_powers[10*N+i]
                                    + 0.000015525 * this->_powers[12*N+i]
                                ) * x[i]
                                +
                                _INVMASK( this->_mask_std[i] )
                                *
                                (
                                    this->_f1[i] * this->_fcn_costh1[i] * this->_fcn_sqrtx_inv[i]
                                );

        }
    }

    void _calculate_k0( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->k0[i]     =   _ASDMASK( this->_mask_modk[i] )
                                *
                                (
                                    - this->_fcn_logpowers[i] * this->i0[i]
                                    - 0.57721566 * this->_powers_modk[0*N+i]
                                    + 0.42278420 * this->_powers_modk[2*N+i]
                                    + 0.23069756 * this->_powers_modk[4*N+i]
                                    + 0.03488590 * this->_powers_modk[6*N+i]
                                    + 0.00262698 * this->_powers_modk[8*N+i]
                                    + 0.00010750 * this->_powers_modk[10*N+i]
                                    + 0.00000740 * this->_powers_modk[12*N+i]
                                )
                                +
                                _INVMASK( this->_mask_modk[i] )
                                *
                                (
                                    + 1.25331414 * this->_powers_modk_inv[0*N+i]
                                    - 0.07832358 * this->_powers_modk_inv[1*N+i]
                                    + 0.02189568 * this->_powers_modk_inv[2*N+i]
                                    - 0.01062446 * this->_powers_modk_inv[3*N+i]
                                    + 0.00587872 * this->_powers_modk_inv[4*N+i]
                                    - 0.00251540 * this->_powers_modk_inv[5*N+i]
                                    + 0.00053208 * this->_powers_modk_inv[6*N+i]
                                ) / ( this->_fcn_expx[i] * this->_fcn_sqrtx_inv[i] );

        }
    }

    void _calculate_k1( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->k1[i]     = 
                                _ASDMASK( this->_mask_modk[i] )
                                *
                                (
                                    x[i] * this->_fcn_logpowers[i] * this->i1[i] + 1.0
                                    + 0.15443144 * this->_powers_modk[2*N+i]
                                    - 0.67278579 * this->_powers_modk[4*N+i]
                                    - 0.18156897 * this->_powers_modk[6*N+i]
                                    - 0.01919402 * this->_powers_modk[8*N+i]
                                    - 0.00110404 * this->_powers_modk[10*N+i]
                                    - 0.00004686 * this->_powers_modk[12*N+i]
                                ) * this->_x_inv[i]
                                +
                                _INVMASK( this->_mask_modk[i] )
                                *
                                (
                                    1.25331414
                                    + 0.23498619 * this->_powers_modk_inv[1*N+i]
                                    - 0.03655620 * this->_powers_modk_inv[2*N+i]
                                    + 0.01504268 * this->_powers_modk_inv[3*N+i]
                                    - 0.00780353 * this->_powers_modk_inv[4*N+i]
                                    + 0.00325614 * this->_powers_modk_inv[5*N+i]
                                    - 0.00068245 * this->_powers_modk_inv[6*N+i]
                                ) / ( this->_fcn_expx[i] * this->_fcn_sqrtx_inv[i] );

        }
    }

    void _calculate_y0( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->y0[i] =   _ASDMASK( this->_mask_std[i] )
                            *
                            (
                                ( 2.0 / PI ) * this->_fcn_logx2[i] * this->j0[i]
                                + 0.367466907 * this->_powers[0*N+i]
                                + 0.605593797 * this->_powers[2*N+i]
                                - 0.743505078 * this->_powers[4*N+i]
                                + 0.253005481 * this->_powers[6*N+i]
                                - 0.042619616 * this->_powers[8*N+i]
                                + 0.004285691 * this->_powers[10*N+i]
                                - 0.000250716 * this->_powers[12*N+i]
                            )
                            +
                            _INVMASK( this->_mask_std[i] )
                            *
                            (
                                this->_fcn_sqrtx_inv[i] * this->_f0[i] * this->_fcn_sinth0[i]
                            );
        }
    }

    void _calculate_y1( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->y1[i] =   _ASDMASK( this->_mask_std[i] )
                            *
                            (
                                ( 2.0 / PI ) * ( this->_fcn_logx2[i] * this->j1[i] - this->_x_inv[i] )
                                + 0.073735531 * this->_powers[1*N+i]
                                + 0.722769344 * this->_powers[3*N+i]
                                - 0.438896337 * this->_powers[5*N+i]
                                + 0.104320251 * this->_powers[7*N+i]
                                - 0.013637596 * this->_powers[9*N+i]
                                + 0.001125970 * this->_powers[11*N+i]
                                - 0.000056455 * this->_powers[13*N+i]
                            )
                            +
                            _INVMASK( this->_mask_std[i] )
                            *
                            (
                                this->_fcn_sqrtx_inv[i] * this->_f1[i] * this->_fcn_sinth1[i]
                            );

        }
    }

    void _calculate_struve0( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->struve0[i]    =   _ASDMASK( this->_mask_std[i] )
                                    *
                                    (
                                        + 1.909859164 * this->_powers[1*N+i]
                                        - 1.909855001 * this->_powers[3*N+i]
                                        + 0.687514637 * this->_powers[5*N+i]
                                        - 0.126164557 * this->_powers[7*N+i]
                                        + 0.013828813 * this->_powers[9*N+i]
                                        - 0.000876918 * this->_powers[11*N+i]
                                    )
                                    +
                                    _INVMASK( this->_mask_std[i] )
                                    *
                                    (
                                        this->y0[i] + this->_sf0[i]
                                    );

        }
    }

    void _calculate_struve1( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->struve1[i] =  _ASDMASK( this->_mask_std[i] )
                                *
                                (
                                    + 1.909859286 * this->_powers[2*N+i]
                                    - 1.145914713 * this->_powers[4*N+i]
                                    + 0.294656958 * this->_powers[6*N+i]
                                    - 0.042070508 * this->_powers[8*N+i]
                                    + 0.003785727 * this->_powers[10*N+i]
                                    - 0.000207183 * this->_powers[12*N+i]
                                )
                                +
                                _INVMASK( this->_mask_std[i] )
                                *
                                (
                                    this->y1[i] + this->_sf1[i]
                                );

        }
    }

    void _calculate_polynomial_f0_th0( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->_f0[i]    =   (
                                    + 0.79788454 * this->_powers_inv[0*N+i]
                                    - 0.00553897 * this->_powers_inv[2*N+i]
                                    + 0.00099336 * this->_powers_inv[4*N+i]
                                    - 0.00044346 * this->_powers_inv[6*N+i]
                                    + 0.00020445 * this->_powers_inv[8*N+i]
                                    - 0.00004959 * this->_powers_inv[10*N+i]
                                );
    
        }

        for ( int i=0; i<N; i++ )
        {
            this->_th0[i] =     (
                                    x[i]         * this->_powers_inv[0*N+i]
                                    - PI/4.0     * this->_powers_inv[0*N+i]
                                    - 0.04166592 * this->_powers_inv[1*N+i]
                                    + 0.00239399 * this->_powers_inv[3*N+i]
                                    - 0.00073984 * this->_powers_inv[5*N+i]
                                    + 0.00031099 * this->_powers_inv[7*N+i]
                                    - 0.00007605 * this->_powers_inv[9*N+i]
                                );
        }
    }

    void _calculate_polynomial_f1_th1( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->_f1[i]    =  (
                                    + 0.79788459 * this->_powers_inv[0*N+i]
                                    + 0.01662008 * this->_powers_inv[2*N+i]
                                    - 0.00187002 * this->_powers_inv[4*N+i]
                                    + 0.00068519 * this->_powers_inv[6*N+i]
                                    - 0.00029440 * this->_powers_inv[8*N+i]
                                    + 0.00006952 * this->_powers_inv[10*N+i]
                                );
        }

        for ( int i=0; i<N; i++ )
        {
            this->_th1[i]   =   (
                                    x[i]         * this->_powers_inv[0*N+i]
                                    - 3.0*PI/4.0 * this->_powers_inv[0*N+i]
                                    + 0.12499895 * this->_powers_inv[1*N+i]
                                    - 0.00605240 * this->_powers_inv[3*N+i]
                                    + 0.00135825 * this->_powers_inv[5*N+i]
                                    - 0.00049616 * this->_powers_inv[7*N+i]
                                    + 0.00011531 * this->_powers_inv[9*N+i]
                                );
        }
    }

    void _calculate_ascending_powers( cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<N; i++ )
        {
            pc[i] = x[i] / BRST;
        }

        for ( int i=0; i<N; i++ )
        {
            this->_powers[i]        = 1.0;
            this->_powers[1*N+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<N; j++ )
            {
                this->_powers[i*N+j] = this->_powers[(i-1)*N+j] * pc[j];
            }
        }
    }

    void _calculate_ascending_powers_modi( cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<N; i++ )
        {
            pc[i] = x[i] / BRMI;
        }

        for ( int i=0; i<N; i++ )
        {
            this->_powers_modi[i]        = 1.0;
            this->_powers_modi[1*N+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<N; j++ )
            {
                this->_powers_modi[i*N+j] = this->_powers_modi[(i-1)*N+j] * pc[j];
            }
        }
    }

    void _calculate_ascending_powers_modk( cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<N; i++ )
        {
            pc[i] = x[i] / BRMK;
        }

        for ( int i=0; i<N; i++ )
        {
            this->_powers_modk[i]        = 1.0;
            this->_powers_modk[1*N+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<N; j++ )
            {
                this->_powers_modk[i*N+j] = this->_powers_modk[(i-1)*N+j] * pc[j];
            }
        }
    }

    void _calculate_inverse_powers( cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<N; i++ )
        {
            pc[i] = BRST / x[i];
        }

        for ( int i=0; i<N; i++ )
        {
            this->_powers_inv[i]        = 1.0;
            this->_powers_inv[1*N+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<N; j++ )
            {
                this->_powers_inv[i*N+j] = this->_powers_inv[(i-1)*N+j] * pc[j];
            }
        }

    }

    void _calculate_inverse_powers_modi( cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<N; i++ )
        {
            pc[i] = BRMI / x[i];
        }

        for ( int i=0; i<N; i++ )
        {
            this->_powers_modi_inv[i]        = 1.0;
            this->_powers_modi_inv[1*N+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<N; j++ )
            {
                this->_powers_modi_inv[i*N+j] = this->_powers_modi_inv[(i-1)*N+j] * pc[j];
            }
        }
    }
    
    void _calculate_inverse_powers_modk( cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<N; i++ )
        {
            pc[i] = BRMK / x[i];
        }

        for ( int i=0; i<N; i++ )
        {
            this->_powers_modk_inv[i]        = 1.0;
            this->_powers_modk_inv[1*N+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<N; j++ )
            {
                this->_powers_modk_inv[i*N+j] = this->_powers_modk_inv[(i-1)*N+j] * pc[j];
            }
        }
    }

    void _calculate_rational_fraction_struve0( cusfloat* x )
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
        cusfloat c0 = 0.0;
        cusfloat c1 = 0.0;

        for ( int i=0; i<N; i++ )
        {
            c0              = 2 * ( a0 * this->_powers_inv[0*N+i] + a1 * this->_powers_inv[2*N+i] +a2 * this->_powers_inv[4*N+i] + a3 * this->_powers_inv[6*N+i] );
            c1              = PI * x[i] * ( 1 * this->_powers_inv[0*N+i] + b1 * this->_powers_inv[2*N+i] + b2 * this->_powers_inv[4*N+i] + b3 * this->_powers_inv[6*N+i] );
            this->_sf0[i]   = c0 / c1;

        }

    }

    void _calculate_rational_fraction_struve1( void )
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
        cusfloat c0 = 0.0;
        cusfloat c1 = 0.0;
        for ( int i=0; i<N; i++ )
        {
            c0              = 2 * ( a0 * this->_powers_inv[0*N+i] + a1 * this->_powers_inv[2*N+i] + a2 * this->_powers_inv[4*N+i] + a3 * this->_powers_inv[6*N+i] );
            c1              = PI * ( 1 * this->_powers_inv[0*N+i] + b1 * this->_powers_inv[2*N+i] + b2 * this->_powers_inv[4*N+i] + b3 * this->_powers_inv[6*N+i] );
            this->_sf1[i]   = c0 / c1;
        }

    }

    void _calculate_standard_aux_fcns( cusfloat* x )
    {
        // Pre-process x for functions
        for ( int i=0; i<N; i++ )
        {
            this->_x_inv[i] = 1 / x[i];
        }

        for ( int i=0; i<N; i++ )
        {
            this->_xd2[i] = x[i] / 2.0;
        }

        // Calculate function values
        lv_cos<cusfloat>( N, this->_th0, this->_fcn_costh0 );
        lv_cos<cusfloat>( N, this->_th1, this->_fcn_costh1 );
        lv_exp<cusfloat>( N, x, this->_fcn_expx );
        lv_log<cusfloat>( N, this->_xd2, this->_fcn_logx2 );
        lv_sin<cusfloat>( N, this->_th0, this->_fcn_sinth0 );
        lv_sin<cusfloat>( N, this->_th1, this->_fcn_sinth1 );
        lv_sqrt<cusfloat>( N, this->_x_inv, this->_fcn_sqrtx_inv );
    }

    void _calculate_modified_aux_fcns( cusfloat* x )
    {
        // Calculate function values
        lv_log<cusfloat>( N, &(this->_powers_modk[1*N]), this->_fcn_logpowers );
    }

    void _distribute_data_standard( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->_mask_std[i] = ( x[i] < BRST ) ? 1: -1;
        }
    }

    void _distribute_data_modified_i( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->_mask_modi[i] = ( x[i] < BRMI ) ? 1: -1;
        }
    }

    void _distribute_data_modified_k( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->_mask_modk[i] = ( x[i] < BRMK ) ? 1: -1;
        }
    }

public:
    // Define class public attributes
    cusfloat i0[N];
    cusfloat i1[N];
    cusfloat j0[N];
    cusfloat j1[N];
    cusfloat k0[N];
    cusfloat k1[N];
    cusfloat y0[N];
    cusfloat y1[N];
    cusfloat struve0[N];
    cusfloat struve1[N];

    // Define class constructors
    BesselFactoryVec( void ) = default;

    // Define class methods
    void calculate_cheby( cusfloat* x )
    {
        this->_distribute_data_standard( x );
        this->_calculate_bessel_standard( x );
    }

    void calculate_series( cusfloat* x )
    {
        this->_distribute_data_standard( x );
        this->_distribute_data_modified_i( x );
        this->_distribute_data_modified_k( x );

        this->_calculate_bessel_standard( x );
        this->_calculate_bessel_modified( x );
    }

};


template<int N>
class BesselFactoryVecUpTo
{
private:
    // Define private attributes
    MEMALINGR   cusfloat    _f0[N];
    MEMALINGR   cusfloat    _f1[N];
    MEMALINGR   cusfloat    _fcn_costh0[N];
    MEMALINGR   cusfloat    _fcn_costh1[N];
    MEMALINGR   cusfloat    _fcn_expx[N];
    MEMALINGR   cusfloat    _fcn_logpowers[N];
    MEMALINGR   cusfloat    _fcn_logx2[N];
    MEMALINGR   cusfloat    _fcn_sinth0[N];
    MEMALINGR   cusfloat    _fcn_sinth1[N];
    MEMALINGR   cusfloat    _fcn_sqrtx_inv[N];
    MEMALINGR   cusfloat    _sf0[N];
    MEMALINGR   cusfloat    _sf1[N];
    MEMALINGR   cusfloat    _th0[N];
    MEMALINGR   cusfloat    _th1[N];
    MEMALINGR   cusfloat    _mask_std[N];
    MEMALINGR   cusfloat    _mask_modi[N];
    MEMALINGR   cusfloat    _mask_modk[N];
                int         _na = 0;
                int         _ni = 0;
    MEMALINGR   cusfloat    _powers[N*NPOWER];
    MEMALINGR   cusfloat    _powers_inv[N*NPOWER];
    MEMALINGR   cusfloat    _powers_modi[N*NPOWER];
    MEMALINGR   cusfloat    _powers_modi_inv[N*NPOWER];
    MEMALINGR   cusfloat    _powers_modk[N*NPOWER];
    MEMALINGR   cusfloat    _powers_modk_inv[N*NPOWER];
    MEMALINGR   cusfloat    _xd2[N];
    MEMALINGR   cusfloat    _x_inv[N];

    // Define private methods
    void _calculate_bessel_standard( int n, cusfloat* x )
    {
        // Calculate base powers
        this->_calculate_ascending_powers( n, x );
        this->_calculate_inverse_powers( n, x );

        // Calculate polynomial expansions
        this->_calculate_polynomial_f0_th0( n, x );
        this->_calculate_polynomial_f1_th1( n, x );
        this->_calculate_rational_fraction_struve0( n, x );
        this->_calculate_rational_fraction_struve1( n );

        // Calculate standard auxiliar functions
        this->_calculate_standard_aux_fcns( n, x );

        // Calculate first kind bessel functions
        this->_calculate_j0( n, x );
        this->_calculate_j1( n, x );

        // Calculate second king bessel functions
        this->_calculate_y0( n, x );
        this->_calculate_y1( n, x );

        // Calculate struve functions
        this->_calculate_struve0( n, x );
        this->_calculate_struve1( n, x );

    }

    void _calculate_bessel_modified( int n, cusfloat* x )
    {
        // Calculate powers
        this->_calculate_ascending_powers_modi( n, x );
        this->_calculate_ascending_powers_modk( n, x );
        this->_calculate_inverse_powers_modi( n, x );
        this->_calculate_inverse_powers_modk( n, x );

        // Calculate standard auxiliar functions
        this->_calculate_modified_aux_fcns( n, x );

        // Calculate first order modified
        this->_calculate_i0( n, x );
        this->_calculate_i1( n, x );

        this->_calculate_k0( n, x );
        this->_calculate_k1( n, x );

    }

    void _calculate_i0( int n, cusfloat* x )
    {
        for ( int i=0; i<n; i++ )
        {
            this->i0[i] =   _ASDMASK( this->_mask_modi[i] )
                            *
                            (
                                + 1.0       * this->_powers_modi[0*n+i]
                                + 3.5156229 * this->_powers_modi[2*n+i]
                                + 3.0899424 * this->_powers_modi[4*n+i]
                                + 1.2067492 * this->_powers_modi[6*n+i]
                                + 0.2659732 * this->_powers_modi[8*n+i]
                                + 0.0360768 * this->_powers_modi[10*n+i]
                                + 0.0045813 * this->_powers_modi[12*n+i]
                            )
                            +
                            _INVMASK( this->_mask_modi[i] )
                            *
                            (
                                + 0.39894228 * this->_powers_modi_inv[0*n+i]
                                + 0.01328592 * this->_powers_modi_inv[1*n+i]
                                + 0.00225319 * this->_powers_modi_inv[2*n+i]
                                - 0.00157565 * this->_powers_modi_inv[3*n+i]
                                + 0.00916281 * this->_powers_modi_inv[4*n+i]
                                - 0.02057706 * this->_powers_modi_inv[5*n+i]
                                + 0.02635537 * this->_powers_modi_inv[6*n+i]
                                - 0.01647633 * this->_powers_modi_inv[7*n+i]
                                + 0.00392377 * this->_powers_modi_inv[8*n+i]
                            ) * this->_fcn_expx[i] * this->_fcn_sqrtx_inv[i];

        }
    }

    void _calculate_i1( int n, cusfloat* x )
    {
        for ( int i=0; i<n; i++ )
        {
            this->i1[i] =   _ASDMASK( this->_mask_modi[i] )
                            *
                            (
                                + 0.5        * this->_powers_modi[0*n+i]
                                + 0.87890594 * this->_powers_modi[2*n+i]
                                + 0.51498869 * this->_powers_modi[4*n+i]
                                + 0.15084934 * this->_powers_modi[6*n+i]
                                + 0.02658733 * this->_powers_modi[8*n+i]
                                + 0.00301532 * this->_powers_modi[10*n+i]
                                + 0.00032411 * this->_powers_modi[12*n+i]
                            ) * x[i]
                            +
                            _INVMASK( this->_mask_modi[i] )
                            *
                            (
                                + 0.39894228 * this->_powers_modi_inv[0*n+i]
                                - 0.03988024 * this->_powers_modi_inv[1*n+i]
                                - 0.00362018 * this->_powers_modi_inv[2*n+i]
                                + 0.00163801 * this->_powers_modi_inv[3*n+i]
                                - 0.01031555 * this->_powers_modi_inv[4*n+i]
                                + 0.02282967 * this->_powers_modi_inv[5*n+i]
                                - 0.02895312 * this->_powers_modi_inv[6*n+i]
                                + 0.01787654 * this->_powers_modi_inv[7*n+i]
                                - 0.00420059 * this->_powers_modi_inv[8*n+i]
                                
                            ) * this->_fcn_expx[i] * this->_fcn_sqrtx_inv[i];

        }
    }

    void _calculate_j0( int n, cusfloat* x )
    {
        for ( int i=0; i<n; i++ )
        {
            this->j0[i] =  _ASDMASK( this->_mask_std[i] )
                            * 
                            (
                                + 0.999999999 * this->_powers[0*n+i]
                                - 2.249999879 * this->_powers[2*n+i]
                                + 1.265623060 * this->_powers[4*n+i]
                                - 0.316394552 * this->_powers[6*n+i]
                                + 0.044460948 * this->_powers[8*n+i]
                                - 0.003954479 * this->_powers[10*n+i]
                                + 0.000212950 * this->_powers[12*n+i]
                            )
                            +
                            _INVMASK( this->_mask_std[i] )
                            *
                            (
                                this->_f0[i] * this->_fcn_costh0[i] * this->_fcn_sqrtx_inv[i]
                            );

        }
    }

    void _calculate_j1( int n, cusfloat* x )
    {
        for ( int i=0; i<n; i++ )
        {
            this->j1[i]     =   _ASDMASK( this->_mask_std[i] )
                                *
                                (
                                    + 0.500000000 * this->_powers[0*n+i]
                                    - 0.562499992 * this->_powers[2*n+i]
                                    + 0.210937377 * this->_powers[4*n+i]
                                    - 0.039550040 * this->_powers[6*n+i]
                                    + 0.004447331 * this->_powers[8*n+i]
                                    - 0.000330547 * this->_powers[10*n+i]
                                    + 0.000015525 * this->_powers[12*n+i]
                                ) * x[i]
                                +
                                _INVMASK( this->_mask_std[i] )
                                *
                                (
                                    this->_f1[i] * this->_fcn_costh1[i] * this->_fcn_sqrtx_inv[i]
                                );

        }
    }

    void _calculate_k0( int n, cusfloat* x )
    {
        for ( int i=0; i<n; i++ )
        {
            this->k0[i]     =   _ASDMASK( this->_mask_modk[i] )
                                *
                                (
                                    - this->_fcn_logpowers[i] * this->i0[i]
                                    - 0.57721566 * this->_powers_modk[0*n+i]
                                    + 0.42278420 * this->_powers_modk[2*n+i]
                                    + 0.23069756 * this->_powers_modk[4*n+i]
                                    + 0.03488590 * this->_powers_modk[6*n+i]
                                    + 0.00262698 * this->_powers_modk[8*n+i]
                                    + 0.00010750 * this->_powers_modk[10*n+i]
                                    + 0.00000740 * this->_powers_modk[12*n+i]
                                )
                                +
                                _INVMASK( this->_mask_modk[i] )
                                *
                                (
                                    + 1.25331414 * this->_powers_modk_inv[0*n+i]
                                    - 0.07832358 * this->_powers_modk_inv[1*n+i]
                                    + 0.02189568 * this->_powers_modk_inv[2*n+i]
                                    - 0.01062446 * this->_powers_modk_inv[3*n+i]
                                    + 0.00587872 * this->_powers_modk_inv[4*n+i]
                                    - 0.00251540 * this->_powers_modk_inv[5*n+i]
                                    + 0.00053208 * this->_powers_modk_inv[6*n+i]
                                ) * this->_fcn_sqrtx_inv[i] / this->_fcn_expx[i];
        }
    }

    void _calculate_k1( int n, cusfloat* x )
    {
        for ( int i=0; i<n; i++ )
        {
            this->k1[i]     = 
                                _ASDMASK( this->_mask_modk[i] )
                                *
                                (
                                    x[i] * this->_fcn_logpowers[i] * this->i1[i] + 1.0
                                    + 0.15443144 * this->_powers_modk[2*n+i]
                                    - 0.67278579 * this->_powers_modk[4*n+i]
                                    - 0.18156897 * this->_powers_modk[6*n+i]
                                    - 0.01919402 * this->_powers_modk[8*n+i]
                                    - 0.00110404 * this->_powers_modk[10*n+i]
                                    - 0.00004686 * this->_powers_modk[12*n+i]
                                ) * this->_x_inv[i]
                                +
                                _INVMASK( this->_mask_modk[i] )
                                *
                                (
                                    1.25331414
                                    + 0.23498619 * this->_powers_modk_inv[1*n+i]
                                    - 0.03655620 * this->_powers_modk_inv[2*n+i]
                                    + 0.01504268 * this->_powers_modk_inv[3*n+i]
                                    - 0.00780353 * this->_powers_modk_inv[4*n+i]
                                    + 0.00325614 * this->_powers_modk_inv[5*n+i]
                                    - 0.00068245 * this->_powers_modk_inv[6*n+i]
                                ) * this->_fcn_sqrtx_inv[i] / this->_fcn_expx[i];

        }
    }

    void _calculate_y0( int n, cusfloat* x )
    {
        for ( int i=0; i<n; i++ )
        {
            this->y0[i] =   _ASDMASK( this->_mask_std[i] )
                            *
                            (
                                ( 2.0 / PI ) * this->_fcn_logx2[i] * this->j0[i]
                                + 0.367466907 * this->_powers[0*n+i]
                                + 0.605593797 * this->_powers[2*n+i]
                                - 0.743505078 * this->_powers[4*n+i]
                                + 0.253005481 * this->_powers[6*n+i]
                                - 0.042619616 * this->_powers[8*n+i]
                                + 0.004285691 * this->_powers[10*n+i]
                                - 0.000250716 * this->_powers[12*n+i]
                            )
                            +
                            _INVMASK( this->_mask_std[i] )
                            *
                            (
                                this->_fcn_sqrtx_inv[i] * this->_f0[i] * this->_fcn_sinth0[i]
                            );
        }
    }

    void _calculate_y1( int n, cusfloat* x )
    {
        for ( int i=0; i<n; i++ )
        {
            this->y1[i] =   _ASDMASK( this->_mask_std[i] )
                            *
                            (
                                ( 2.0 / PI ) * ( this->_fcn_logx2[i] * this->j1[i] - this->_x_inv[i] )
                                + 0.073735531 * this->_powers[1*n+i]
                                + 0.722769344 * this->_powers[3*n+i]
                                - 0.438896337 * this->_powers[5*n+i]
                                + 0.104320251 * this->_powers[7*n+i]
                                - 0.013637596 * this->_powers[9*n+i]
                                + 0.001125970 * this->_powers[11*n+i]
                                - 0.000056455 * this->_powers[13*n+i]
                            )
                            +
                            _INVMASK( this->_mask_std[i] )
                            *
                            (
                                this->_fcn_sqrtx_inv[i] * this->_f1[i] * this->_fcn_sinth1[i]
                            );

        }
    }

    void _calculate_struve0( int n, cusfloat* x )
    {
        for ( int i=0; i<n; i++ )
        {
            this->struve0[i]    =   _ASDMASK( this->_mask_std[i] )
                                    *
                                    (
                                        + 1.909859164 * this->_powers[1*n+i]
                                        - 1.909855001 * this->_powers[3*n+i]
                                        + 0.687514637 * this->_powers[5*n+i]
                                        - 0.126164557 * this->_powers[7*n+i]
                                        + 0.013828813 * this->_powers[9*n+i]
                                        - 0.000876918 * this->_powers[11*n+i]
                                    )
                                    +
                                    _INVMASK( this->_mask_std[i] )
                                    *
                                    (
                                        this->y0[i] + this->_sf0[i]
                                    );

        }
    }

    void _calculate_struve1( int n, cusfloat* x )
    {
        for ( int i=0; i<n; i++ )
        {
            this->struve1[i] =  _ASDMASK( this->_mask_std[i] )
                                *
                                (
                                    + 1.909859286 * this->_powers[2*n+i]
                                    - 1.145914713 * this->_powers[4*n+i]
                                    + 0.294656958 * this->_powers[6*n+i]
                                    - 0.042070508 * this->_powers[8*n+i]
                                    + 0.003785727 * this->_powers[10*n+i]
                                    - 0.000207183 * this->_powers[12*n+i]
                                )
                                +
                                _INVMASK( this->_mask_std[i] )
                                *
                                (
                                    this->y1[i] + this->_sf1[i]
                                );

        }
    }

    void _calculate_polynomial_f0_th0( int n, cusfloat* x )
    {
        for ( int i=0; i<n; i++ )
        {
            this->_f0[i]    =   (
                                    + 0.79788454 * this->_powers_inv[0*n+i]
                                    - 0.00553897 * this->_powers_inv[2*n+i]
                                    + 0.00099336 * this->_powers_inv[4*n+i]
                                    - 0.00044346 * this->_powers_inv[6*n+i]
                                    + 0.00020445 * this->_powers_inv[8*n+i]
                                    - 0.00004959 * this->_powers_inv[10*n+i]
                                );
    
        }

        for ( int i=0; i<n; i++ )
        {
            this->_th0[i] =     (
                                    x[i]         * this->_powers_inv[0*n+i]
                                    - PI/4.0     * this->_powers_inv[0*n+i]
                                    - 0.04166592 * this->_powers_inv[1*n+i]
                                    + 0.00239399 * this->_powers_inv[3*n+i]
                                    - 0.00073984 * this->_powers_inv[5*n+i]
                                    + 0.00031099 * this->_powers_inv[7*n+i]
                                    - 0.00007605 * this->_powers_inv[9*n+i]
                                );
        }
    }

    void _calculate_polynomial_f1_th1( int n, cusfloat* x )
    {
        for ( int i=0; i<n; i++ )
        {
            this->_f1[i]    =  (
                                    + 0.79788459 * this->_powers_inv[0*n+i]
                                    + 0.01662008 * this->_powers_inv[2*n+i]
                                    - 0.00187002 * this->_powers_inv[4*n+i]
                                    + 0.00068519 * this->_powers_inv[6*n+i]
                                    - 0.00029440 * this->_powers_inv[8*n+i]
                                    + 0.00006952 * this->_powers_inv[10*n+i]
                                );
        }

        for ( int i=0; i<n; i++ )
        {
            this->_th1[i]   =   (
                                    x[i]         * this->_powers_inv[0*n+i]
                                    - 3.0*PI/4.0 * this->_powers_inv[0*n+i]
                                    + 0.12499895 * this->_powers_inv[1*n+i]
                                    - 0.00605240 * this->_powers_inv[3*n+i]
                                    + 0.00135825 * this->_powers_inv[5*n+i]
                                    - 0.00049616 * this->_powers_inv[7*n+i]
                                    + 0.00011531 * this->_powers_inv[9*n+i]
                                );
        }
    }

    void _calculate_ascending_powers( int n, cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<n; i++ )
        {
            pc[i] = x[i] / BRST;
        }

        for ( int i=0; i<n; i++ )
        {
            this->_powers[i]        = 1.0;
            this->_powers[1*n+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<n; j++ )
            {
                this->_powers[i*n+j] = this->_powers[(i-1)*n+j] * pc[j];
            }
        }
    }

    void _calculate_ascending_powers_modi( int n, cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<n; i++ )
        {
            pc[i] = x[i] / BRMI;
        }

        for ( int i=0; i<n; i++ )
        {
            this->_powers_modi[i]        = 1.0;
            this->_powers_modi[1*n+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<n; j++ )
            {
                this->_powers_modi[i*n+j] = this->_powers_modi[(i-1)*n+j] * pc[j];
            }
        }
    }

    void _calculate_ascending_powers_modk( int n, cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<n; i++ )
        {
            pc[i] = x[i] / BRMK;
        }

        for ( int i=0; i<n; i++ )
        {
            this->_powers_modk[i]        = 1.0;
            this->_powers_modk[1*n+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<n; j++ )
            {
                this->_powers_modk[i*n+j] = this->_powers_modk[(i-1)*n+j] * pc[j];
            }
        }
    }

    void _calculate_inverse_powers( int n, cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<n; i++ )
        {
            pc[i] = BRST / x[i];
        }

        for ( int i=0; i<n; i++ )
        {
            this->_powers_inv[i]        = 1.0;
            this->_powers_inv[1*n+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<n; j++ )
            {
                this->_powers_inv[i*n+j] = this->_powers_inv[(i-1)*n+j] * pc[j];
            }
        }

    }

    void _calculate_inverse_powers_modi( int n, cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<n; i++ )
        {
            pc[i] = BRMI / x[i];
        }

        for ( int i=0; i<n; i++ )
        {
            this->_powers_modi_inv[i]        = 1.0;
            this->_powers_modi_inv[1*n+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<n; j++ )
            {
                this->_powers_modi_inv[i*n+j] = this->_powers_modi_inv[(i-1)*n+j] * pc[j];
            }
        }
    }
    
    void _calculate_inverse_powers_modk( int n, cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<n; i++ )
        {
            pc[i] = BRMK / x[i];
        }

        for ( int i=0; i<n; i++ )
        {
            this->_powers_modk_inv[i]        = 1.0;
            this->_powers_modk_inv[1*n+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<n; j++ )
            {
                this->_powers_modk_inv[i*n+j] = this->_powers_modk_inv[(i-1)*n+j] * pc[j];
            }
        }
    }

    void _calculate_rational_fraction_struve0( int n, cusfloat* x )
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
        cusfloat c0 = 0.0;
        cusfloat c1 = 0.0;

        for ( int i=0; i<n; i++ )
        {
            c0              = 2 * ( a0 * this->_powers_inv[0*n+i] + a1 * this->_powers_inv[2*n+i] +a2 * this->_powers_inv[4*n+i] + a3 * this->_powers_inv[6*n+i] );
            c1              = PI * x[i] * ( 1 * this->_powers_inv[0*n+i] + b1 * this->_powers_inv[2*n+i] + b2 * this->_powers_inv[4*n+i] + b3 * this->_powers_inv[6*n+i] );
            this->_sf0[i]   = c0 / c1;

        }

    }

    void _calculate_rational_fraction_struve1( int n )
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
        cusfloat c0 = 0.0;
        cusfloat c1 = 0.0;
        for ( int i=0; i<n; i++ )
        {
            c0              = 2 * ( a0 * this->_powers_inv[0*n+i] + a1 * this->_powers_inv[2*n+i] + a2 * this->_powers_inv[4*n+i] + a3 * this->_powers_inv[6*n+i] );
            c1              = PI * ( 1 * this->_powers_inv[0*n+i] + b1 * this->_powers_inv[2*n+i] + b2 * this->_powers_inv[4*n+i] + b3 * this->_powers_inv[6*n+i] );
            this->_sf1[i]   = c0 / c1;
        }

    }

    void _calculate_standard_aux_fcns( int n, cusfloat* x )
    {
        // Pre-process x for functions
        for ( int i=0; i<n; i++ )
        {
            this->_x_inv[i] = 1 / x[i];
        }

        for ( int i=0; i<n; i++ )
        {
            this->_xd2[i] = x[i] / 2.0;
        }

        // Calculate function values
        lv_cos<cusfloat>( n, this->_th0, this->_fcn_costh0 );
        lv_cos<cusfloat>( n, this->_th1, this->_fcn_costh1 );
        lv_exp<cusfloat>( n, x, this->_fcn_expx );
        lv_log<cusfloat>( n, this->_xd2, this->_fcn_logx2 );
        lv_sin<cusfloat>( n, this->_th0, this->_fcn_sinth0 );
        lv_sin<cusfloat>( n, this->_th1, this->_fcn_sinth1 );
        lv_sqrt<cusfloat>( n, this->_x_inv, this->_fcn_sqrtx_inv );
    }

    void _calculate_modified_aux_fcns( int n, cusfloat* x )
    {
        // Calculate function values
        lv_log<cusfloat>( n, &(this->_powers_modk[1*n]), this->_fcn_logpowers );
    }

    void _distribute_data_standard( int n, cusfloat* x )
    {
        for ( int i=0; i<n; i++ )
        {
            this->_mask_std[i] = ( x[i] < BRST ) ? 1: -1;
        }
    }

    void _distribute_data_modified_i( int n, cusfloat* x )
    {
        for ( int i=0; i<n; i++ )
        {
            this->_mask_modi[i] = ( x[i] < BRMI ) ? 1: -1;
        }
    }

    void _distribute_data_modified_k( int n, cusfloat* x )
    {
        for ( int i=0; i<n; i++ )
        {
            this->_mask_modk[i] = ( x[i] < BRMK ) ? 1: -1;
        }
    }

public:
    // Define class public attributes
    cusfloat i0[N];
    cusfloat i1[N];
    cusfloat j0[N];
    cusfloat j1[N];
    cusfloat k0[N];
    cusfloat k1[N];
    cusfloat y0[N];
    cusfloat y1[N];
    cusfloat struve0[N];
    cusfloat struve1[N];

    // Define class constructors
    BesselFactoryVecUpTo( ) = default;

    // Define class methods
    void calculate_cheby( int n, cusfloat* x )
    {
        this->_distribute_data_standard( n, x );
        this->_calculate_bessel_standard( n, x );
    }

    void calculate_series( int n, cusfloat* x )
    {
        this->_distribute_data_standard( n, x );
        this->_distribute_data_modified_i( n, x );
        this->_distribute_data_modified_k( n, x );

        this->_calculate_bessel_standard( n, x );
        this->_calculate_bessel_modified( n, x );
    }

    void calculate_modified( int n, cusfloat* x )
    {
        // Check limits
        cusfloat xl[N];
        for ( int i=0; i<n; i++ )
        {
            xl[i] = std::min( x[i], 50.0 );
        }

        // Pre-calcualte modified bessel function paramets
        this->_distribute_data_modified_i( n, xl );
        this->_distribute_data_modified_k( n, xl );

        // Calculate auxiliar function. In this case is included the "standard ones" also
        // beacause some of them are appearing in modified bessel formulations
        this->_calculate_standard_aux_fcns( n, xl );

        // Calculate modified bessel function values
        this->_calculate_bessel_modified( n, xl );
    }

};


template<int N>
class BesselFactoryBranch
{
private:
    // Define private attributes
    MEMALINGR   cusfloat    _f0[N];
    MEMALINGR   cusfloat    _f1[N];
    MEMALINGR   cusfloat    _fcn_costh0[N];
    MEMALINGR   cusfloat    _fcn_costh1[N];
    MEMALINGR   cusfloat    _fcn_expx[N];
    MEMALINGR   cusfloat    _fcn_logpowers[N];
    MEMALINGR   cusfloat    _fcn_logx2[N];
    MEMALINGR   cusfloat    _fcn_sinth0[N];
    MEMALINGR   cusfloat    _fcn_sinth1[N];
    MEMALINGR   cusfloat    _fcn_sqrtx_inv[N];
    MEMALINGR   cusfloat    _sf0[N];
    MEMALINGR   cusfloat    _sf1[N];
    MEMALINGR   cusfloat    _th0[N];
    MEMALINGR   cusfloat    _th1[N];
    MEMALINGR   cusfloat     _mask_std[N];
    MEMALINGR   cusfloat     _mask_modi[N];
    MEMALINGR   cusfloat     _mask_modk[N];
                int         _na = 0;
                int         _ni = 0;
    MEMALINGR   cusfloat    _powers[N*NPOWER];
    MEMALINGR   cusfloat    _powers_inv[N*NPOWER];
    MEMALINGR   cusfloat    _powers_modi[N*NPOWER];
    MEMALINGR   cusfloat    _powers_modi_inv[N*NPOWER];
    MEMALINGR   cusfloat    _powers_modk[N*NPOWER];
    MEMALINGR   cusfloat    _powers_modk_inv[N*NPOWER];
    MEMALINGR   cusfloat    _xd2[N];
    MEMALINGR   cusfloat    _x_inv[N];

    // Define private methods
    void _calculate_bessel_standard( cusfloat* x )
    {
        // Calculate base powers
        this->_calculate_ascending_powers( x );
        this->_calculate_inverse_powers( x );

        // Calculate polynomial expansions
        this->_calculate_polynomial_f0_th0( x );
        this->_calculate_polynomial_f1_th1( x );
        this->_calculate_rational_fraction_struve0( x );
        this->_calculate_rational_fraction_struve1( );

        // Calculate standard auxiliar functions
        this->_calculate_standard_aux_fcns( x );

        // Calculate first kind bessel functions
        this->_calculate_j0( x );
        this->_calculate_j1( x );

        // Calculate second king bessel functions
        this->_calculate_y0( x );
        this->_calculate_y1( x );

        // Calculate struve functions
        this->_calculate_struve0( x );
        this->_calculate_struve1( x );

    }

    void _calculate_bessel_modified( cusfloat* x )
    {
        // Calculate powers
        this->_calculate_ascending_powers_modi( x );
        this->_calculate_ascending_powers_modk( x );
        this->_calculate_inverse_powers_modi( x );
        this->_calculate_inverse_powers_modk( x );

        // Calculate standard auxiliar functions
        this->_calculate_modified_aux_fcns( x );

        // Calculate first order modified
        this->_calculate_i0( x );
        this->_calculate_i1( x );

        this->_calculate_k0( x );
        this->_calculate_k1( x );

    }

    void _calculate_i0( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            if ( this->_mask_modi[i] > 0 )
            {
                this->i0[i] =   (
                                    + 1.0       * this->_powers_modi[0*N+i]
                                    + 3.5156229 * this->_powers_modi[2*N+i]
                                    + 3.0899424 * this->_powers_modi[4*N+i]
                                    + 1.2067492 * this->_powers_modi[6*N+i]
                                    + 0.2659732 * this->_powers_modi[8*N+i]
                                    + 0.0360768 * this->_powers_modi[10*N+i]
                                    + 0.0045813 * this->_powers_modi[12*N+i]
                                );

            }
            else
            {
                this->i0[i] =   (
                                    + 0.39894228 * this->_powers_modi_inv[0*N+i]
                                    + 0.01328592 * this->_powers_modi_inv[1*N+i]
                                    + 0.00225319 * this->_powers_modi_inv[2*N+i]
                                    - 0.00157565 * this->_powers_modi_inv[3*N+i]
                                    + 0.00916281 * this->_powers_modi_inv[4*N+i]
                                    - 0.02057706 * this->_powers_modi_inv[5*N+i]
                                    + 0.02635537 * this->_powers_modi_inv[6*N+i]
                                    - 0.01647633 * this->_powers_modi_inv[7*N+i]
                                    + 0.00392377 * this->_powers_modi_inv[8*N+i]
                                ) * this->_fcn_expx[i] * this->_fcn_sqrtx_inv[i];
            }

        }
    }

    void _calculate_i1( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            if ( this->_mask_modi[i] > 0 )
            {
                this->i1[i] =   (
                                    + 0.5        * this->_powers_modi[0*N+i]
                                    + 0.87890594 * this->_powers_modi[2*N+i]
                                    + 0.51498869 * this->_powers_modi[4*N+i]
                                    + 0.15084934 * this->_powers_modi[6*N+i]
                                    + 0.02658733 * this->_powers_modi[8*N+i]
                                    + 0.00301532 * this->_powers_modi[10*N+i]
                                    + 0.00032411 * this->_powers_modi[12*N+i]
                                ) * this->_x_inv[i];

            }
            else
            {
                this->i1[i] = (
                                    + 0.39894228 * this->_powers_modi_inv[0*N+i]
                                    - 0.03988024 * this->_powers_modi_inv[1*N+i]
                                    - 0.00362018 * this->_powers_modi_inv[2*N+i]
                                    + 0.00163801 * this->_powers_modi_inv[3*N+i]
                                    - 0.01031555 * this->_powers_modi_inv[4*N+i]
                                    + 0.02282967 * this->_powers_modi_inv[5*N+i]
                                    - 0.02895312 * this->_powers_modi_inv[6*N+i]
                                    + 0.01787654 * this->_powers_modi_inv[7*N+i]
                                    - 0.00420059 * this->_powers_modi_inv[8*N+i]
                                    
                                ) * this->_fcn_expx[i] * this->_fcn_sqrtx_inv[i];

            }

        }
    }

    void _calculate_j0( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            if ( this->_mask_std[i] > 0 )
            {
                this->j0[i] =  (
                                    + 0.999999999 * this->_powers[0*N+i]
                                    - 2.249999879 * this->_powers[2*N+i]
                                    + 1.265623060 * this->_powers[4*N+i]
                                    - 0.316394552 * this->_powers[6*N+i]
                                    + 0.044460948 * this->_powers[8*N+i]
                                    - 0.003954479 * this->_powers[10*N+i]
                                    + 0.000212950 * this->_powers[12*N+i]
                                );

            }
            else
            {
                this->j0[i] =   (
                                    this->_f0[i] * this->_fcn_costh0[i] * this->_fcn_sqrtx_inv[i]
                                );

            }

        }
    }

    void _calculate_j1( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            if ( this->_mask_std[i] > 0 )
            {
                this->j1[i]     =   (
                                        + 0.500000000 * this->_powers[0*N+i]
                                        - 0.562499992 * this->_powers[2*N+i]
                                        + 0.210937377 * this->_powers[4*N+i]
                                        - 0.039550040 * this->_powers[6*N+i]
                                        + 0.004447331 * this->_powers[8*N+i]
                                        - 0.000330547 * this->_powers[10*N+i]
                                        + 0.000015525 * this->_powers[12*N+i]
                                    ) * this->_x_inv[i];

            }
            else
            {
                this->j1[i]     =   (
                                        this->_f1[i] * this->_fcn_costh1[i] * this->_fcn_sqrtx_inv[i]
                                    );

            }

        }
    }

    void _calculate_k0( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            if ( this->_mask_modk[i] > 0 )
            {
                this->k0[i]     =   (
                                        -this->_fcn_logpowers[i] * this->i0[i]
                                        - 0.57721566 * this->_powers_modk[0*N+i]
                                        + 0.42278420 * this->_powers_modk[2*N+i]
                                        + 0.23069756 * this->_powers_modk[4*N+i]
                                        + 0.03488590 * this->_powers_modk[6*N+i]
                                        + 0.00262698 * this->_powers_modk[8*N+i]
                                        + 0.00010750 * this->_powers_modk[10*N+i]
                                        + 0.00000740 * this->_powers_modk[12*N+i]
                                    );

            }
            else
            {
                this->k0[i]     =   (
                                        + 1.25331414 * this->_powers_modk_inv[0*N+i]
                                        - 0.07832358 * this->_powers_modk_inv[1*N+i]
                                        + 0.02189568 * this->_powers_modk_inv[2*N+i]
                                        - 0.01062446 * this->_powers_modk_inv[3*N+i]
                                        + 0.00587872 * this->_powers_modk_inv[4*N+i]
                                        - 0.00251540 * this->_powers_modk_inv[5*N+i]
                                        + 0.00053208 * this->_powers_modk_inv[6*N+i]
                                    ) / ( this->_fcn_expx[i] * this->_fcn_sqrtx_inv[i] );

            }

        }
    }

    void _calculate_k1( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            if ( this->_mask_modk[i] > 0 )
            {
                this->k1[i]     =   (
                                        x[i] * this->_fcn_logpowers[i] * this->i1[i] + 1.0
                                        + 0.15443144 * this->_powers_modk[2*N+i]
                                        - 0.67278579 * this->_powers_modk[4*N+i]
                                        - 0.18156897 * this->_powers_modk[6*N+i]
                                        - 0.01919402 * this->_powers_modk[8*N+i]
                                        - 0.00110404 * this->_powers_modk[10*N+i]
                                        - 0.00004686 * this->_powers_modk[12*N+i]
                                    ) * this->_x_inv[i];

            }
            else
            {
                this->k1[i]     =   (
                                        + 1.25331414 * this->_powers_modk_inv[0*N+i]
                                        + 0.23498619 * this->_powers_modk_inv[1*N+i]
                                        - 0.03655620 * this->_powers_modk_inv[2*N+i]
                                        + 0.01504268 * this->_powers_modk_inv[3*N+i]
                                        - 0.00780353 * this->_powers_modk_inv[4*N+i]
                                        + 0.00325614 * this->_powers_modk_inv[5*N+i]
                                        - 0.00068245 * this->_powers_modk_inv[6*N+i]
                                    ) / ( this->_fcn_expx[i] * this->_fcn_sqrtx_inv[i] );

            }
            
        }
    }

    void _calculate_y0( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            if ( this->_mask_std[i] > 0 )
            {
                this->y0[i] =   (
                                    ( 2.0 / PI ) * this->_fcn_logx2[i] * this->j0[i]
                                    + 0.367466907 * this->_powers[0*N+i]
                                    + 0.605593797 * this->_powers[2*N+i]
                                    - 0.743505078 * this->_powers[4*N+i]
                                    + 0.253005481 * this->_powers[6*N+i]
                                    - 0.042619616 * this->_powers[8*N+i]
                                    + 0.004285691 * this->_powers[10*N+i]
                                    - 0.000250716 * this->_powers[12*N+i]
                                );
            }
            else
            {
                this->y0[i] =   (
                                    this->_f0[i] * this->_fcn_sinth0[i] * this->_fcn_sqrtx_inv[i]
                                );
            }
        }
    }

    void _calculate_y1( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            if ( this->_mask_std[i] > 0 )
            {
                this->y1[i] =   (
                                    ( 2.0 / PI ) * ( this->_fcn_logx2[i] * this->j1[i] - this->_x_inv[i] )
                                    + 0.073735531 * this->_powers[1*N+i]
                                    + 0.722769344 * this->_powers[3*N+i]
                                    - 0.438896337 * this->_powers[5*N+i]
                                    + 0.104320251 * this->_powers[7*N+i]
                                    - 0.013637596 * this->_powers[9*N+i]
                                    + 0.001125970 * this->_powers[11*N+i]
                                    - 0.000056455 * this->_powers[13*N+i]
                                );

            }
            else
            {
                this->y1[i] =   (
                                    this->_f1[i] * this->_fcn_sinth1[i] * this->_fcn_sqrtx_inv[i]
                                );
            }

        }
    }

    void _calculate_struve0( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            if ( this->_mask_std[i] > 0 )
            {
                this->struve0[i]    =   (
                                            + 1.909859164 * this->_powers[1*N+i]
                                            - 1.909855001 * this->_powers[3*N+i]
                                            + 0.687514637 * this->_powers[5*N+i]
                                            - 0.126164557 * this->_powers[7*N+i]
                                            + 0.013828813 * this->_powers[9*N+i]
                                            - 0.000876918 * this->_powers[11*N+i]
                                        );
            }
            else
            {
                this->struve0[i]    =   (
                                            this->y0[i] + this->_sf0[i]
                                        );

            }

        }
    }

    void _calculate_struve1( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            if ( this->_mask_std[i] > 0 )
            {
                this->struve1[i] =  (
                                        + 1.909859286 * this->_powers[2*N+i]
                                        - 1.145914713 * this->_powers[4*N+i]
                                        + 0.294656958 * this->_powers[6*N+i]
                                        - 0.042070508 * this->_powers[8*N+i]
                                        + 0.003785727 * this->_powers[10*N+i]
                                        - 0.000207183 * this->_powers[12*N+i]
                                    );
            }
            else
            {
                this->struve1[i] = (
                                        this->y1[i] + this->_sf1[i]
                                    );
            }
        }
    }

    void _calculate_polynomial_f0_th0( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            if ( this->_mask_std[i] < 0 )
            {
                this->_f0[i]    =   (
                                        + 0.79788454 * this->_powers_inv[0*N+i]
                                        - 0.00553897 * this->_powers_inv[2*N+i]
                                        + 0.00099336 * this->_powers_inv[4*N+i]
                                        - 0.00044346 * this->_powers_inv[6*N+i]
                                        + 0.00020445 * this->_powers_inv[8*N+i]
                                        - 0.00004959 * this->_powers_inv[10*N+i]
                                    );
            }

            
        }

        for ( int i=0; i<N; i++ )
        {
            if ( this->_mask_std[i] < 0 )
            {
                this->_th0[i] =     (
                                        x[i] - PI/4.0
                                        - 0.04166592 * this->_powers_inv[1*N+i]
                                        + 0.00239399 * this->_powers_inv[3*N+i]
                                        - 0.00073984 * this->_powers_inv[5*N+i]
                                        + 0.00031099 * this->_powers_inv[7*N+i]
                                        - 0.00007605 * this->_powers_inv[9*N+i]
                                    );
            }
        }
    }

    void _calculate_polynomial_f1_th1( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            if ( this->_mask_std[i] < 0 )
            {
                this->_f1[i]    =  (
                                        + 0.79788459 * this->_powers_inv[0*N+i]
                                        + 0.01662008 * this->_powers_inv[2*N+i]
                                        - 0.00187002 * this->_powers_inv[4*N+i]
                                        + 0.00068519 * this->_powers_inv[6*N+i]
                                        - 0.00029440 * this->_powers_inv[8*N+i]
                                        + 0.00006952 * this->_powers_inv[10*N+i]
                                    );
            }
        }

        for ( int i=0; i<N; i++ )
        {
            if ( this->_mask_std[i] < 0 )
            {
                this->_th1[i]   =   (
                                        x[i] - 3.0*PI/4.0
                                        + 0.12499895 * this->_powers_inv[1*N+i]
                                        - 0.00605240 * this->_powers_inv[3*N+i]
                                        + 0.00135825 * this->_powers_inv[5*N+i]
                                        - 0.00049616 * this->_powers_inv[7*N+i]
                                        + 0.00011531 * this->_powers_inv[9*N+i]
                                    );
            }
        }
    }

    void _calculate_ascending_powers( cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<N; i++ )
        {
            pc[i] = x[i] / BRST;
        }

        for ( int i=0; i<N; i++ )
        {
            this->_powers[i]        = 1.0;
            this->_powers[1*N+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<N; j++ )
            {
                this->_powers[i*N+j] = this->_powers[(i-1)*N+j] * pc[j];
            }
        }
    }

    void _calculate_ascending_powers_modi( cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<N; i++ )
        {
            pc[i] = x[i] / BRMI;
        }

        for ( int i=0; i<N; i++ )
        {
            this->_powers_modi[i]        = 1.0;
            this->_powers_modi[1*N+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<N; j++ )
            {
                this->_powers_modi[i*N+j] = this->_powers_modi[(i-1)*N+j] * pc[j];
            }
        }
    }

    void _calculate_ascending_powers_modk( cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<N; i++ )
        {
            pc[i] = x[i] / BRMK;
        }

        for ( int i=0; i<N; i++ )
        {
            this->_powers_modk[i]        = 1.0;
            this->_powers_modk[1*N+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<N; j++ )
            {
                this->_powers_modk[i*N+j] = this->_powers_modk[(i-1)*N+j] * pc[j];
            }
        }
    }

    void _calculate_inverse_powers( cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<N; i++ )
        {
            pc[i] = BRST / x[i];
        }

        for ( int i=0; i<N; i++ )
        {
            this->_powers_inv[i]        = 1.0;
            this->_powers_inv[1*N+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<N; j++ )
            {
                this->_powers_inv[i*N+j] = this->_powers_inv[(i-1)*N+j] * pc[j];
            }
        }
    }

    void _calculate_inverse_powers_modi( cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<N; i++ )
        {
            pc[i] = BRMI / x[i];
        }

        for ( int i=0; i<N; i++ )
        {
            this->_powers_modi_inv[i]        = 1.0;
            this->_powers_modi_inv[1*N+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<N; j++ )
            {
                this->_powers_modi_inv[i*N+j] = this->_powers_modi_inv[(i-1)*N+j] * pc[j];
            }
        }
    }
    
    void _calculate_inverse_powers_modk( cusfloat* x )
    {
        cusfloat pc[N];
        for ( int i=0; i<N; i++ )
        {
            pc[i] = BRMK / x[i];
        }

        for ( int i=0; i<N; i++ )
        {
            this->_powers_modk_inv[i]        = 1.0;
            this->_powers_modk_inv[1*N+i]    = pc[i];
        }
        
        for ( std::size_t i=2; i<NPOWER; i++ )
        {
            for ( int j=0; j<N; j++ )
            {
                this->_powers_modk_inv[i*N+j] = this->_powers_modk_inv[(i-1)*N+j] * pc[j];
            }
        }
    }

    void _calculate_rational_fraction_struve0( cusfloat* x )
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
        cusfloat c0 = 0.0;
        cusfloat c1 = 0.0;

        for ( int i=0; i<N; i++ )
        {
            c0              = 2 * ( a0 * this->_powers_inv[0*N+i] + a1 * this->_powers_inv[2*N+i] +a2 * this->_powers_inv[4*N+i] + a3 * this->_powers_inv[6*N+i] );
            c1              = PI * x[i] * ( 1 * this->_powers_inv[0*N+i] + b1 * this->_powers_inv[2*N+i] + b2 * this->_powers_inv[4*N+i] + b3 * this->_powers_inv[6*N+i] );
            this->_sf0[i]   = c0 / c1;

        }

    }

    void _calculate_rational_fraction_struve1( void )
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
        cusfloat c0 = 0.0;
        cusfloat c1 = 0.0;
        for ( int i=0; i<N; i++ )
        {
            c0              = 2 * ( a0 * this->_powers_inv[0*N+i] + a1 * this->_powers_inv[2*N+i] + a2 * this->_powers_inv[4*N+i] + a3 * this->_powers_inv[6*N+i] );
            c1              = PI * ( 1 * this->_powers_inv[0*N+i] + b1 * this->_powers_inv[2*N+i] + b2 * this->_powers_inv[4*N+i] + b3 * this->_powers_inv[6*N+i] );
            this->_sf1[i]   = c0 / c1;
        }

    }

    void _calculate_standard_aux_fcns( cusfloat* x )
    {
        // Pre-process x for functions
        for ( int i=0; i<N; i++ )
        {
            this->_x_inv[i] = 1 / x[i];
        }

        for ( int i=0; i<N; i++ )
        {
            this->_xd2[i] = x[i] / 2.0;
        }

        // Calculate function values
        lv_cos<cusfloat>( N, this->_th0, this->_fcn_costh0 );
        lv_cos<cusfloat>( N, this->_th1, this->_fcn_costh1 );
        lv_exp<cusfloat>( N, x, this->_fcn_expx );
        lv_log<cusfloat>( N, this->_xd2, this->_fcn_logx2 );
        lv_sin<cusfloat>( N, this->_th0, this->_fcn_sinth0 );
        lv_sin<cusfloat>( N, this->_th1, this->_fcn_sinth1 );
        lv_sqrt<cusfloat>( N, this->_x_inv, this->_fcn_sqrtx_inv );
    }

    void _calculate_modified_aux_fcns( cusfloat* x )
    {
        // Calculate function values
        lv_log<cusfloat>( N, &(this->_powers_modk[1*N]), this->_fcn_logpowers );
    }

    void _distribute_data_standard( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->_mask_std[i] = ( x[i] < BRST ) ? 1: -1;
        }
    }

    void _distribute_data_modified_i( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->_mask_modi[i] = ( x[i] < BRMI ) ? 1: -1;
        }
    }

    void _distribute_data_modified_k( cusfloat* x )
    {
        for ( int i=0; i<N; i++ )
        {
            this->_mask_modk[i] = ( x[i] < BRMK ) ? 1: -1;
        }
    }

public:
    // Define class public attributes
    cusfloat i0[N];
    cusfloat i1[N];
    cusfloat j0[N];
    cusfloat j1[N];
    cusfloat k0[N];
    cusfloat k1[N];
    cusfloat y0[N];
    cusfloat y1[N];
    cusfloat struve0[N];
    cusfloat struve1[N];

    // Define class constructors
    BesselFactoryBranch( void ) = default;

    // Define class methods
    void calculate_cheby( cusfloat* x )
    {
        this->_distribute_data_standard( x );
        this->_calculate_bessel_standard( x );
    }

    void calculate_series( cusfloat* x )
    {
        this->_distribute_data_standard( x );
        this->_distribute_data_modified_i( x );
        this->_distribute_data_modified_k( x );

        this->_calculate_bessel_standard( x );
        this->_calculate_bessel_modified( x );
    }

};