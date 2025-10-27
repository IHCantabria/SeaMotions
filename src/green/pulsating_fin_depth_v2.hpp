
#pragma once

// Include general usage libraries
#include <tuple>

// Include local modules
#include "../config.hpp"
#include "../math/bessel_factory.hpp"
#include "../math/special_math.hpp"
#include "../waves/wave_dispersion_fo.hpp"


// Define module macros


// Define module constants
#define     GREEN_ZEROTH_DR     1e-4
#define     GREEN_DR_EPS        1e-6

// Declare module functions
template<
            std::size_t N, 
            int         mode_loop, 
            int         mode_f, 
            int         mode_dfdr, 
            int         mode_dfdz
        >                           void        Fxy(
                                                                        const std::size_t           n,
                                                                        cusfloat*                   X,
                                                                        cusfloat*                   Y,
                                                                        BesselFactoryVecUpTo<N>&    bessel_factory,
                                                                        cusfloat*                   results,
                                                                        cusfloat*                   results_dx,
                                                                        cusfloat*                   results_dy
                                                    );


template<
            std::size_t N, 
            int         mode_loop, 
            int         mode_f, 
            int         mode_dfdr, 
            int         mode_dfdz
        >                           void        G1_Hlt1(
                                                                        const std::size_t   n,
                                                                        cusfloat*           A,
                                                                        cusfloat*           B,
                                                                        cusfloat*           results,
                                                                        cusfloat*           results_da,
                                                                        cusfloat*           results_db
                                                        );

                                            
template<
            std::size_t N, 
            int         mode_loop, 
            int         mode_f, 
            int         mode_dfdr, 
            int         mode_dfdz
        >                           void        G1_Hgt1(
                                                                        const std::size_t   n,
                                                                        cusfloat*           A,
                                                                        cusfloat*           B,
                                                                        cusfloat*           results,
                                                                        cusfloat*           results_da,
                                                                        cusfloat*           results_db
                                                        );


template<
            std::size_t N, 
            int         mode_loop, 
            int         mode_f, 
            int         mode_dfdr, 
            int         mode_dfdz
        >                           void        G2_Hlt1(
                                                                        const std::size_t   n,
                                                                        cusfloat*           A,
                                                                        cusfloat*           B,
                                                                        cusfloat*           results,
                                                                        cusfloat*           results_da,
                                                                        cusfloat*           results_db
                                                        );


template<
            std::size_t N, 
            int         mode_loop, 
            int         mode_f, 
            int         mode_dfdr, 
            int         mode_dfdz
        >                           void        G2_Hgt1(
                                                                        const std::size_t   n,
                                                                        cusfloat*           A,
                                                                        cusfloat*           B,
                                                                        cusfloat*           results,
                                                                        cusfloat*           results_da,
                                                                        cusfloat*           results_db
                                                );


template<
            std::size_t N, 
            int         mode_f, 
            int         mode_dfdr, 
            int         mode_dfdz
        >                           void        john_series(
                                                                        cusfloat*               R,
                                                                        cusfloat*               z,
                                                                        cusfloat*               zeta,
                                                                        cusfloat                h,
                                                                        BesselFactoryVecUpTo<N> &bessel_factory,
                                                                        WaveDispersionFONK      &wave_data,
                                                                        cuscomplex*             G,
                                                                        cuscomplex*             G_dr,
                                                                        cuscomplex*             G_dz,
                                                                        cuscomplex*             G_dzeta
                                                                );


template<
            std::size_t N, 
            int         mode_f, 
            int         mode_dfdr, 
            int         mode_dfdz,
            int         mode_fslid
        >                           void        wave_term_integral(
                                                                        cusfloat*                   R,
                                                                        cusfloat*                   z,
                                                                        cusfloat*                   zeta,
                                                                        cusfloat                    h,
                                                                        BesselFactoryVecUpTo<N>     &bessel_factory,
                                                                        WaveDispersionFONK          &wave_data,
                                                                        cuscomplex*                 G,
                                                                        cuscomplex*                 G_dr,
                                                                        cuscomplex*                 G_dz,
                                                                        cuscomplex*                 G_dzeta
                                                                    );


inline cuscomplex  john_series(
                                                cusfloat R,
                                                cusfloat z, 
                                                cusfloat zeta, 
                                                cusfloat h,
                                                WaveDispersionFONK &wave_data
                        )
{
    /**
     * @brief John series representation of the finite water depth Green function
     * 
     * This is an eigenfunction expasion generated by Fritz John in the article
     * "On the Motion of Floating Bodies II". In the current implementation the 
     * series has been modified to work withou hyperbolic cosines in order to 
     * reduce numerical problems with big number. It is explained at:
     * "Consistent expression for the free-surfae Green function in finite water
     * depth" - Ed Mackay.
     * 
     * \param R Horizontal distance between source and the field point
     * \param z Vertical coordinate of the field point
     * \param zeta Vertical coordinate of the source point
     * \param h Water depth
     * \param wave_data Wave dispersion data object initialized with John's constants
     * \return value of the Jonh series
     */

    // Check if John's coeffcients were precalculated
    // assert(wave_data.is_john);

    // Get some data from wave data to have shorter names
    cusfloat k0 = wave_data.k0;

    // Calculate vertical distances between the source and
    // the field point
    cusfloat v3 = z+zeta;
    cusfloat v4 = z-zeta+2*h;
    cusfloat v5 = zeta-z+2*h;
    cusfloat v6 = z+zeta+4*h;

    // Calcuate real root series part
    cusfloat k0r = k0*R;
    cusfloat expsum = (
                        + exp(k0*v3)
                        + exp(-k0*v4)
                        + exp(-k0*v5)
                        + exp(-k0*v6)
                    );
    cusfloat    c_real = -wave_data.k0nu*expsum;
    cuscomplex  sol( bessely0(k0r), -besselj0(k0r) );
                sol = c_real*sol;

    // Calculate imag root series part
    cusfloat c_imag = 0.0;
    cusfloat ci = 0.0;
    int count_k = 0;
    cusfloat kni = 0.0;
    cusfloat zetah = zeta+h;
    cusfloat zh = z+h;
    while (true)
    {
        // Calculate i term of the series
        kni = wave_data.kn[count_k];
        ci = 4*wave_data.knnu[count_k]*cos(kni*zh)*cos(kni*zetah)*besselk0(kni*R);
        c_imag += ci;
        // std::cout << "knnu: " << wave_data.knnu[count_k] << " - cos(kni*zh): " << cos(kni*zh) << " - cos(kni*zetah): " << cos(kni*zetah) << " - knir: " << kni*R << " - besselk0(kni*R): " << besselk0(kni*R) << " - Cimag: " << c_imag << std::endl;

        // Check for convergence
        if (abs(ci)<EPS_PRECISION)
        {
            break;
        }

        // Check for the limit in imaginary roots
        if (count_k > (wave_data.num_kn-2))
        {
            std::cerr << "Jonh series could not converge up to the precision required with" << std::endl;
            std::cerr << "Input parameters:"  << std::endl;
            std::cerr << "  - R/h: " << R/h << std::endl;
            std::cerr << "  - R: " << R << std::endl;
            std::cerr << "  - z: " << z << std::endl;
            std::cerr << "  - zeta: " << zeta << std::endl;
            std::cerr << "  - h: " << h << std::endl;
            std::cerr << "  - nu: " << wave_data.nu << std::endl;
            std::cerr << "  - k0: " << k0 << std::endl;
            std::cerr << "  - num_kn: " << wave_data.num_kn << std::endl;
            std::cerr << "* Current term value is: " << ci << std::endl;
            throw std::runtime_error("Jonh series value could not converge. See log file for details.");
        }

        // Update counter
        count_k++;
    }
    sol += c_imag;

    return sol;
}



// Include definition module
#include "pulsating_fin_depth_v2.txx"