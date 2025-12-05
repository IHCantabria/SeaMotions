
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
#include "../math/math_tools.hpp"
#include "../waves/wave_dispersion_so.hpp"


// Declare module functions
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
                                    );

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
                                        );

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
                                        );

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
                                        );

cuscomplex  wave_potential_so_space( 
                                            cusfloat            x,
                                            cusfloat            y,
                                            cusfloat            z,
                                            WaveDispersionSO*   wd,
                                            int                 qtf_type
                                    );

cuscomplex  wave_potential_so_space_dx( 
                                            cusfloat            x,
                                            cusfloat            y,
                                            cusfloat            z,
                                            WaveDispersionSO*   wd,
                                            int                 qtf_type
                                    );

cuscomplex  wave_potential_so_space_dy( 
                                            cusfloat            x,
                                            cusfloat            y,
                                            cusfloat            z,
                                            WaveDispersionSO*   wd,
                                            int                 qtf_type
                                    );

cuscomplex  wave_potential_so_space_dz( 
                                            cusfloat            x,
                                            cusfloat            y,
                                            cusfloat            z,
                                            WaveDispersionSO*   wd,
                                            int                 qtf_type
                                    );

cusfloat    wave_vertical_profile_fo(
                                            cusfloat    k,
                                            cusfloat    h,
                                            cusfloat    z
                                    );

cusfloat    wave_vertical_profile_fo_dz(
                                            cusfloat    k,
                                            cusfloat    h,
                                            cusfloat    z
                                        );

cusfloat    wave_vertical_profile_mod_fo(
                                            cusfloat    k,
                                            cusfloat    h,
                                            cusfloat    z
                                        );
