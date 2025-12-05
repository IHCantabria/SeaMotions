
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

// Include local modules
#include "../config.hpp"


struct WaveDispersionSO
{
private:
    // Define private class methods
    void    _update_status(
                                cusfloat    w0_in,
                                cusfloat    w1_in,
                                cusfloat    head_0_in,
                                cusfloat    head_1_in
                            );

public:
    // Define public class attributes
    cusfloat    a0              = 0.0;
    cusfloat    a1              = 0.0;
    cusfloat    grav_acc        = 0.0;
    cusfloat    head_0          = 0.0;
    cusfloat    head_1          = 0.0;
    cusfloat    k0              = 0.0;
    cusfloat    k0_deep_water   = 0.0;
    cusfloat    k0_vec[2]       = { 0.0, 0.0 };
    cusfloat    k1              = 0.0;
    cusfloat    k1_deep_water   = 0.0;
    cusfloat    k1_vec[2]       = { 0.0, 0.0 };
    cusfloat    k_diff_mod      = 0.0;
    cusfloat    k_diff_vec[2]   = { 0.0, 0.0 };
    cusfloat    k_sum_mod       = 0.0;
    cusfloat    k_sum_vec[2]    = { 0.0, 0.0 };
    cusfloat    w0              = 0.0;
    cusfloat    w1              = 0.0;
    cusfloat    water_depth     = 0.0;
    cusfloat    w_diff          = 0.0;
    cusfloat    w_k_diff        = 0.0;
    cusfloat    w_k_sum         = 0.0;
    cusfloat    w_sum           = 0.0;

    // Define class constructors and destructor
    WaveDispersionSO( ) = default;

    WaveDispersionSO(
                        cusfloat    a0_in,
                        cusfloat    a1_in,
                        cusfloat    w0_in,
                        cusfloat    w1_in,
                        cusfloat    head_0_in,
                        cusfloat    head_1_in,
                        cusfloat    water_depth,
                        cusfloat    grav_acc
                    );

    // Define class methods
    cusfloat    get_w_ds(
                                int     qtf_type
                        );

    void        set_new_data(
                                cusfloat w0_in,
                                cusfloat w1_in,
                                cusfloat head_0_in,
                                cusfloat head_1_in
                        );

};
