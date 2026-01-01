
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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
            int         mode_f, 
            int         mode_dfdr, 
            int         mode_dfdz, 
            int         mode_fslid
        >                           void        fd_green_rad_asymptotic(
                                                                                cusfloat*                   Ri,
                                                                                cusfloat*                   zi,
                                                                                cusfloat*                   zeta,
                                                                                cusfloat                    h,
                                                                                BesselFactoryVecUpTo<N>     &bessel_factory,
                                                                                WaveDispersionFONK          &wave_data,
                                                                                cuscomplex*                 G,
                                                                                cuscomplex*                 G_dr,
                                                                                cuscomplex*                 G_dz,
                                                                                cuscomplex*                 G_dzeta
                                                                        );


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


inline                              cuscomplex  john_series(
                                                                        cusfloat R,
                                                                        cusfloat z, 
                                                                        cusfloat zeta, 
                                                                        cusfloat h,
                                                                        WaveDispersionFONK &wave_data
                                                            );


template<
            std::size_t N
        >                   void        wave_term_fin_depth(
                                                                        cusfloat*               R,
                                                                        cusfloat*               z,
                                                                        cusfloat*               zeta,
                                                                        cusfloat                h,
                                                                        BesselFactoryVecUpTo<N> &bessel_factory,
                                                                        WaveDispersionFONK      &wave_data,
                                                                        cuscomplex*             G,
                                                                        cuscomplex*             G_dr,
                                                                        cuscomplex*             G_dz
                                                        );


template<
            std::size_t N
        >                   void        wave_term_fin_depth_integral(
                                                                        const std::size_t           n,
                                                                        cusfloat*                   R,
                                                                        cusfloat*                   z,
                                                                        cusfloat*                   zeta,
                                                                        cusfloat                    h,
                                                                        BesselFactoryVecUpTo<N>     &bessel_factory,
                                                                        WaveDispersionFONK          &wave_data,
                                                                        cuscomplex*                 G,
                                                                        cuscomplex*                 G_dr,
                                                                        cuscomplex*                 G_dz
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


template<
            std::size_t N, 
            int         mode_f, 
            int         mode_dfdr, 
            int         mode_dfdz, 
            int         mode_fslid
        >                           void        wave_term_integral_inf_freq(
                                                                                    cusfloat*                   Ri,
                                                                                    cusfloat*                   zi,
                                                                                    cusfloat*                   zeta,
                                                                                    cusfloat                    h,
                                                                                    BesselFactoryVecUpTo<N>     ,
                                                                                    WaveDispersionFONK          ,
                                                                                    cuscomplex*                 G,
                                                                                    cuscomplex*                 G_dr,
                                                                                    cuscomplex*                 G_dz,
                                                                                    cuscomplex*                 G_dzeta
                                                                            );

                                        
template<
            std::size_t N, 
            int         mode_f, 
            int         mode_dfdr, 
            int         mode_dfdz, 
            int         mode_fslid
        >
                                    void        wave_term_integral_zero_freq(
                                                                                            cusfloat*                   Ri,
                                                                                            cusfloat*                   zi,
                                                                                            cusfloat*                   zeta,
                                                                                            cusfloat                    h,
                                                                                            BesselFactoryVecUpTo<N>     ,
                                                                                            WaveDispersionFONK          ,
                                                                                            cuscomplex*                 G,
                                                                                            cuscomplex*                 G_dr,
                                                                                            cuscomplex*                 G_dz,
                                                                                            cuscomplex*                 G_dzeta
                                                                            );


// Include definition module
#include "pulsating_fin_depth.txx"