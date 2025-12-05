
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

// Include general usage libraries
#include <cassert>
#include <iostream>

// Include general usage scientific libraries
#include <cmath>
#include <complex>

// Inlude local modules
#include "waves_common.hpp"

#include "../config.hpp"
#include "../math/math_tools.hpp"
#include "wave_dispersion_so.hpp"


cusfloat    cosh_cosh_factor( 
                                            cusfloat k,
                                            cusfloat h,
                                            cusfloat z
                            )
{
    cusfloat    kz  = k * z;
    cusfloat    kh  = k * h;
    return ( std::exp( kz ) + std::exp( -2.0 * kh ) * std::exp( -kz ) ) / ( 1.0 + std::exp( -2.0 * kh ) );
}


cusfloat    sinh_cosh_factor( 
                                            cusfloat k,
                                            cusfloat h,
                                            cusfloat z
                            )
{
    cusfloat    kz  = k * z;
    cusfloat    kh  = k * h;
    return ( std::exp( kz ) - std::exp( -2.0 * kh ) * std::exp( -kz ) ) / ( 1.0 + std::exp( -2.0 * kh ) );
}


cuscomplex  wave_potential_fo_space( 
                                            cusfloat aw,
                                            cusfloat w,
                                            cusfloat k,
                                            cusfloat h,
                                            cusfloat g,
                                            cusfloat x,
                                            cusfloat y,
                                            cusfloat z,
                                            cusfloat mu
                                    )
{
    cuscomplex  exp_arg = cuscomplex( 0.0, k * ( x * std::cos( mu ) + y * std::sin( mu ) ) );
    cusfloat    phi_0   = aw * g / w * cosh_cosh_factor( k, h, z );
    return cuscomplex( 0.0, -phi_0 ) * std::exp( exp_arg );
}


cuscomplex  wave_potential_fo_space_dx( 
                                            cusfloat aw,
                                            cusfloat w,
                                            cusfloat k,
                                            cusfloat h,
                                            cusfloat g,
                                            cusfloat x,
                                            cusfloat y,
                                            cusfloat z,
                                            cusfloat mu
                                        )
{
    cuscomplex  exp_arg = cuscomplex( 0.0, k * ( x * std::cos( mu ) + y * std::sin( mu ) ) );
    cusfloat    phi_0   = aw * g / w * cosh_cosh_factor( k, h, z ) * k * std::cos( mu );
    return phi_0 * std::exp( exp_arg );
}


cuscomplex  wave_potential_fo_space_dy( 
                                            cusfloat aw,
                                            cusfloat w,
                                            cusfloat k,
                                            cusfloat h,
                                            cusfloat g,
                                            cusfloat x,
                                            cusfloat y,
                                            cusfloat z,
                                            cusfloat mu
                                        )
{
    cuscomplex  exp_arg = cuscomplex( 0.0, k * ( x * std::cos( mu ) + y * std::sin( mu ) ) );
    cusfloat    phi_0   = aw * g / w * cosh_cosh_factor( k, h, z ) * k * std::sin( mu );
    return phi_0 * std::exp( exp_arg );
}


cuscomplex  wave_potential_fo_space_dz( 
                                            cusfloat aw,
                                            cusfloat w,
                                            cusfloat k,
                                            cusfloat h,
                                            cusfloat g,
                                            cusfloat x,
                                            cusfloat y,
                                            cusfloat z,
                                            cusfloat mu
                                        )
{
    cuscomplex  exp_arg = cuscomplex( 0.0, k * ( x * std::cos( mu ) + y * std::sin( mu ) ) );
    cusfloat    phi_0   = aw * g * k / w * sinh_cosh_factor( k, h, z );
    return cuscomplex( 0.0, -phi_0 ) * std::exp( exp_arg );
}


cuscomplex  wave_potential_so_space( 
                                        cusfloat            x,
                                        cusfloat            y,
                                        cusfloat            z,
                                        WaveDispersionSO*   wd,
                                        int                 qtf_type
                                    )
{
    // Check if there is difference or summation term
    cusfloat    k_ds_mod    = 0.0;
    cusfloat    k_ds_vec[2] = { 0.0, 0.0 };
    cusfloat    s0          = 0.0;
    cusfloat    s1          = 0.0;
    cusfloat    w_ds        = 0.0;
    cusfloat    w_k_ds      = 0.0;

    assert( ( ( qtf_type == 0 ) || ( qtf_type == 1 ) ) && "qtf_type should be 0 or 1" );
    if ( qtf_type == 0 )
    {
        k_ds_mod        = wd->k_diff_mod;
        k_ds_vec[0]     = wd->k_diff_vec[0];
        k_ds_vec[1]     = wd->k_diff_vec[1];
        s0              = 1.0;
        s1              = -1.0;
        w_ds            = wd->w_diff;
        w_k_ds          = wd->w_k_diff;
    }
    else if ( qtf_type == 1 )
    {
        k_ds_mod        = wd->k_sum_mod;
        k_ds_vec[0]     = wd->k_sum_vec[0];
        k_ds_vec[1]     = wd->k_sum_vec[1];
        s0              = -1.0;
        s1              = 1.0;
        w_ds            = wd->w_sum;
        w_k_ds          = wd->w_k_sum;
    }

    // Calculate first term
    cuscomplex  a0  =   (
                            cuscomplex( 0.0, 1.0 )
                            *
                            wd->a0
                            *
                            wd->a1
                            *
                            pow2s( wd->grav_acc )
                            *
                            std::exp( cuscomplex( 0.0, 1.0 ) * ( k_ds_vec[0] * x + k_ds_vec[1] * y ) )
                            *
                            wave_vertical_profile_fo( k_ds_mod, wd->water_depth, z )
                        );
    
    cuscomplex  a1  =   (
                            -pow2s( w_ds )
                            +
                            pow2s( w_k_ds )
                        );

    cuscomplex  b   =   (
                            w_ds / wd->w0 / wd->w1
                            *
                            ( 
                                wd->k0_vec[0] * wd->k1_vec[0]
                                +
                                wd->k0_vec[1] * wd->k1_vec[1]
                                +
                                s0 * wd->k0_deep_water * wd->k1_deep_water
                            )
                        );

    cuscomplex  c   =   0.5 * (
                            ( pow2s( wd->k0 ) - pow2s( wd->k0_deep_water ) ) / wd->w0
                            +
                            s1 * ( pow2s( wd->k1 ) - pow2s( wd->k1_deep_water ) ) / wd->w1
                        );

    // Calculate potential
    cuscomplex  pot2    = a0/a1 * ( b + c );

    return pot2;
}


cuscomplex  wave_potential_so_space_dx( 
                                        cusfloat            x,
                                        cusfloat            y,
                                        cusfloat            z,
                                        WaveDispersionSO*   wd,
                                        int                 qtf_type
                                    )
{
    // Check if there is difference or summation term
    cusfloat    k_ds_vec[2] = { 0.0, 0.0 };

    assert( ( ( qtf_type == 0 ) || ( qtf_type == 1 ) ) && "qtf_type should be 0 or 1" );
    if ( qtf_type == 0 )
    {
        k_ds_vec[0]     = wd->k_diff_vec[0];
        k_ds_vec[1]     = wd->k_diff_vec[1];
    }
    else if ( qtf_type == 1 )
    {
        k_ds_vec[0]     = wd->k_sum_vec[0];
        k_ds_vec[1]     = wd->k_sum_vec[1];
    }

    // Calculate potential
    cuscomplex  pot2    =   k_ds_vec[0] * wave_potential_so_space( 
                                                                        x,
                                                                        y,
                                                                        z,
                                                                        wd,
                                                                        qtf_type
                                                                );

    return pot2;
}


cuscomplex  wave_potential_so_space_dy( 
                                        cusfloat            x,
                                        cusfloat            y,
                                        cusfloat            z,
                                        WaveDispersionSO*   wd,
                                        int                 qtf_type
                                    )
{
    // Check if there is difference or summation term
    cusfloat    k_ds_vec[2] = { 0.0, 0.0 };

    assert( ( ( qtf_type == 0 ) || ( qtf_type == 1 ) ) && "qtf_type should be 0 or 1" );
    if ( qtf_type == 0 )
    {
        k_ds_vec[0]     = wd->k_diff_vec[0];
        k_ds_vec[1]     = wd->k_diff_vec[1];
    }
    else if ( qtf_type == 1 )
    {
        k_ds_vec[0]     = wd->k_sum_vec[0];
        k_ds_vec[1]     = wd->k_sum_vec[1];
    }

    // Calculate potential
    cuscomplex  pot2    =   k_ds_vec[1] * wave_potential_so_space( 
                                                                        x,
                                                                        y,
                                                                        z,
                                                                        wd,
                                                                        qtf_type
                                                                    );

    return pot2;
}


cuscomplex  wave_potential_so_space_dz( 
                                        cusfloat            x,
                                        cusfloat            y,
                                        cusfloat            z,
                                        WaveDispersionSO*   wd,
                                        int                 qtf_type
                                    )
{
    // Check if there is difference or summation term
    cusfloat    k_ds_mod    = 0.0;

    assert( ( ( qtf_type == 0 ) || ( qtf_type == 1 ) ) && "qtf_type should be 0 or 1" );
    if ( qtf_type == 0 )
    {
        k_ds_mod        = wd->k_diff_mod;
    }
    else if ( qtf_type == 1 )
    {
        k_ds_mod        = wd->k_sum_mod;
    }

    // Calculate scale factor
    cuscomplex  scale   = (
                                wave_vertical_profile_fo_dz( k_ds_mod, wd->water_depth, z )
                                /
                                wave_vertical_profile_fo( k_ds_mod, wd->water_depth, z )
                            );

    // Calculate potential

    cuscomplex  pot2    =   scale * wave_potential_so_space( 
                                                                x,
                                                                y,
                                                                z,
                                                                wd,
                                                                qtf_type
                                                            );

    return pot2;
}


cusfloat    wave_vertical_profile_fo(
                                            cusfloat    k,
                                            cusfloat    h,
                                            cusfloat    z
                                    )
{
    return std::cosh( k * ( h + z ) ) / std::cosh( k * h );
}


cusfloat    wave_vertical_profile_fo_dz(
                                            cusfloat    k,
                                            cusfloat    h,
                                            cusfloat    z
                                        )
{
    return k * std::sinh( k * ( h + z ) ) / std::cosh( k * h );
}


cusfloat    wave_vertical_profile_mod_fo(
                                            cusfloat    k,
                                            cusfloat    h,
                                            cusfloat    z
                                        )
{
    return wave_vertical_profile_fo( k, h, z ) / ( k * h * ( 1 - pow2s( std::tanh( k * h ) ) ) + std::tanh( k * h ) );
}