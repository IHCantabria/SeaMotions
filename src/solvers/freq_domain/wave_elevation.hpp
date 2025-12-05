
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
#include "../../config.hpp"
#include "../../containers/matlin_group.hpp"
#include "../../inout/input.hpp"
#include "../../containers/mpi_config.hpp"


void    calculate_relative_wave_elevation_lin(
                                                Input*          input,
                                                MLGCmpx*        pot_gp,
                                                cuscomplex*     potpanel_total,
                                                cusfloat        ang_freq,
                                                cuscomplex*     raos,
                                                cuscomplex*     rel_wave_elevation
                                            );


void    calculate_wave_elevation_lin(
                                                cuscomplex*     pot_total,
                                                int             pot_total_np,
                                                cusfloat        ang_freq,
                                                cusfloat        grav_acc,
                                                cuscomplex*     wave_elevation
                                    );
