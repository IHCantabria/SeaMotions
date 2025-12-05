
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
#include "integrals_db.hpp"

cusfloat calculate_dr_dx(cusfloat R, cusfloat dx);
cusfloat calculate_r(
                    cusfloat x,
                    cusfloat y,
                    cusfloat xi,
                    cusfloat eta
                    );

cusfloat wave_term_inf_depth(
                            cusfloat x_ndim,
                            cusfloat y_ndim,
                            IntegralsDb &idb
                            );
std::tuple<cusfloat, cusfloat> wave_term_inf_depth_dhoriz(
                                                            cusfloat x,
                                                            cusfloat y,
                                                            cusfloat z,
                                                            cusfloat xi,
                                                            cusfloat eta,
                                                            cusfloat zeta,
                                                            cusfloat nu,
                                                            IntegralsDb &idb
                                                            );
cusfloat wave_term_inf_depth_dx(
                                cusfloat x,
                                cusfloat y,
                                cusfloat z,
                                cusfloat xi,
                                cusfloat eta,
                                cusfloat zeta,
                                cusfloat nu,
                                IntegralsDb &idb
                                );
cusfloat wave_term_inf_depth_dxndim(
                                cusfloat x_ndim,
                                cusfloat y_ndim,
                                IntegralsDb &idb
                                );
cusfloat wave_term_inf_depth_dy(
                                cusfloat x,
                                cusfloat y,
                                cusfloat z,
                                cusfloat xi,
                                cusfloat eta,
                                cusfloat zeta,
                                cusfloat nu,
                                IntegralsDb &idb
                                );
cusfloat wave_term_inf_depth_dyndim(
                                    cusfloat x_ndim,
                                    cusfloat y_ndim,
                                    IntegralsDb &idb
                                    );
cusfloat wave_term_inf_depth_dz(
                                cusfloat x_ndim,
                                cusfloat z,
                                cusfloat zeta,
                                cusfloat nu,
                                IntegralsDb &idb
                                );
cusfloat wave_term_inf_depth_dz(
                                cusfloat x,
                                cusfloat y,
                                cusfloat z,
                                cusfloat xi,
                                cusfloat eta,
                                cusfloat zeta,
                                cusfloat nu,
                                IntegralsDb &idb
                                );
