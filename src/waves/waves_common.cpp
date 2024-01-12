
// Include general usage libraries
#include <iostream>

// Include general usage scientific libraries
#include <cmath>
#include <complex>

// Inlude local modules
#include "waves_common.hpp"

#include "../config.hpp"
#include "../math/math_tools.hpp"
#include "wave_dispersion_so.hpp"


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
    cusfloat    phi_0   = aw * g * std::cosh( k * ( h + z ) ) / w / std::cosh( k * h );
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
    cusfloat    phi_0   = aw * g * std::cosh( k * ( h + z ) ) / w / std::cosh( k * h ) * k * std::cos( mu );
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
    cusfloat    phi_0   = aw * g * std::cosh( k * ( h + z ) ) / w / std::cosh( k * h ) * k * std::sin( mu );
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
    cusfloat    phi_0   = aw * g * k * std::sinh( k * ( h + z ) ) / w / std::cosh( k * h );
    return cuscomplex( 0.0, -phi_0 ) * std::exp( exp_arg );
}


cuscomplex  wave_potential_so_space( 
                                        cusfloat            x,
                                        cusfloat            y,
                                        cusfloat            z,
                                        WaveDispersionSO*   wd,
                                        bool                is_diff
                                    )
{
    // Check if there is difference or summation term
    cusfloat    k_ds_mod    = 0.0;
    cusfloat    k_ds_vec[2] = { 0.0, 0.0 };
    cusfloat    s0          = 0.0;
    cusfloat    s1          = 0.0;
    cusfloat    w_ds        = 0.0;
    cusfloat    w_k_ds      = 0.0;

    if ( is_diff )
    {
        k_ds_mod        = wd->k_diff_mod;
        k_ds_vec[0]     = wd->k_diff_vec[0];
        k_ds_vec[1]     = wd->k_diff_vec[1];
        s0              = 1.0;
        s1              = -1.0;
        w_ds            = wd->w_diff;
        w_k_ds          = wd->w_k_diff;
    }
    else
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
                                        bool                is_diff
                                    )
{
    // Check if there is difference or summation term
    cusfloat    k_ds_vec[2] = { 0.0, 0.0 };

    if ( is_diff )
    {
        k_ds_vec[0]     = wd->k_diff_vec[0];
        k_ds_vec[1]     = wd->k_diff_vec[1];
    }
    else
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
                                                                        is_diff
                                                                );

    return pot2;
}


cuscomplex  wave_potential_so_space_dy( 
                                        cusfloat            x,
                                        cusfloat            y,
                                        cusfloat            z,
                                        WaveDispersionSO*   wd,
                                        bool                is_diff
                                    )
{
    // Check if there is difference or summation term
    cusfloat    k_ds_vec[2] = { 0.0, 0.0 };

    if ( is_diff )
    {
        k_ds_vec[0]     = wd->k_diff_vec[0];
        k_ds_vec[1]     = wd->k_diff_vec[1];
    }
    else
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
                                                                        is_diff
                                                                    );

    return pot2;
}


cuscomplex  wave_potential_so_space_dz( 
                                        cusfloat            x,
                                        cusfloat            y,
                                        cusfloat            z,
                                        WaveDispersionSO*   wd,
                                        bool                is_diff
                                    )
{
    // Check if there is difference or summation term
    cusfloat    k_ds_mod    = 0.0;

    if ( is_diff )
    {
        k_ds_mod        = wd->k_diff_mod;
    }
    else
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
                                                                is_diff
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