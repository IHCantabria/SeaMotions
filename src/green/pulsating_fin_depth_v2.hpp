
#pragma once

// Include general usage libraries
#include <tuple>

// Include local modules
#include "../config.hpp"
#include "../math/bessel_factory.hpp"
#include "../waves/wave_dispersion_fo.hpp"


// Define module constants
#define     GREEN_ZEROTH_DR     1e-4
#define     GREEN_DR_EPS        1e-6

// Declare module functions
template<std::size_t N>     void         Fxy(
                                                                        const std::size_t           n,
                                                                        cusfloat*                   X,
                                                                        cusfloat*                   Y,
                                                                        BesselFactoryVecUpTo<N>&    bessel_factory,
                                                                        cusfloat*                   results,
                                                                        cusfloat*                   results_dx,
                                                                        cusfloat*                   results_dy
                                            );


template<std::size_t N>     void         G1_Hlt1(
                                                                        const std::size_t   n,
                                                                        cusfloat*           A,
                                                                        cusfloat*           B,
                                                                        cusfloat*           results,
                                                                        cusfloat*           results_da,
                                                                        cusfloat*           results_db
                                                );

                                            
template<std::size_t N>     void         G1_Hgt1(
                                                                        const std::size_t   n,
                                                                        cusfloat*           A,
                                                                        cusfloat*           B,
                                                                        cusfloat*           results,
                                                                        cusfloat*           results_da,
                                                                        cusfloat*           results_db
                                                );


template<std::size_t N>     void         G2_Hlt1(
                                                                        const std::size_t   n,
                                                                        cusfloat*           A,
                                                                        cusfloat*           B,
                                                                        cusfloat*           results,
                                                                        cusfloat*           results_da,
                                                                        cusfloat*           results_db
                                                );


template<std::size_t N>     void         G2_Hgt1(
                                                                        const std::size_t   n,
                                                                        cusfloat*           A,
                                                                        cusfloat*           B,
                                                                        cusfloat*           results,
                                                                        cusfloat*           results_da,
                                                                        cusfloat*           results_db
                                                );


template<std::size_t N>     void        john_series(
                                                                        const std::size_t       n,
                                                                        cusfloat*               R,
                                                                        cusfloat*               z,
                                                                        cusfloat*               zeta,
                                                                        cusfloat                h,
                                                                        BesselFactoryVecUpTo<N> &bessel_factory,
                                                                        WaveDispersionFO        &wave_data,
                                                                        cuscomplex*             G,
                                                                        cuscomplex*             G_dr,
                                                                        cuscomplex*             G_dz
                                                    );


template<std::size_t N>     void        wave_term_fin_depth(
                                                                        cusfloat*               R,
                                                                        cusfloat*               z,
                                                                        cusfloat*               zeta,
                                                                        cusfloat                h,
                                                                        BesselFactoryVecUpTo<N> &bessel_factory,
                                                                        WaveDispersionFO        &wave_data,
                                                                        cuscomplex*             G,
                                                                        cuscomplex*             G_dr,
                                                                        cuscomplex*             G_dz
                                                        );


template<std::size_t N>     void        custom_template(
                                                                        cusfloat*               R
                                                        );


template<std::size_t N>     void        wave_term_fin_depth_integral(
                                                                        const std::size_t           n,
                                                                        cusfloat*                   R,
                                                                        cusfloat*                   z,
                                                                        cusfloat*                   zeta,
                                                                        cusfloat                    h,
                                                                        BesselFactoryVecUpTo<N>     &bessel_factory,
                                                                        WaveDispersionFO            &wave_data,
                                                                        cuscomplex*                 G,
                                                                        cuscomplex*                 G_dr,
                                                                        cuscomplex*                 G_dz
                                                                    );


// Include definition module
#include "pulsating_fin_depth_v2.txx"