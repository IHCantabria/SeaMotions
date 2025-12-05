
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
#include "../../containers/mpi_config.hpp"
#include "../../containers/simulation_data.hpp"
#include "../../inout/input.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_raddif_velocity_mat_steady(
                                                    Input*          input,
                                                    MeshGroup*      mesh_gp,
                                                    MLGCmpx*        vel_x_gp,
                                                    MLGCmpx*        vel_y_gp,
                                                    MLGCmpx*        vel_z_gp
                                            );


void    calculate_raddif_velocity_mat_steady_nlin(
                                                    Input*          input,
                                                    MeshGroup*      mesh_gp,
                                                    MLGCmpx*        vel_x_gp,
                                                    MLGCmpx*        vel_y_gp,
                                                    MLGCmpx*        vel_z_gp
                                                );


void    calculate_raddif_velocity_mat_wave(
                                                    Input*          input,
                                                    MeshGroup*      mesh_gp,
                                                    cusfloat        ang_freq,
                                                    MLGCmpx*        vel_x_gp,
                                                    MLGCmpx*        vel_y_gp,
                                                    MLGCmpx*        vel_z_gp
                                            );


void    calculate_velocities_total(
                                                    Input*          input,
                                                    MpiConfig*      mpi_config,
                                                    MeshGroup*      mesh_gp,
                                                    cusfloat        ang_freq,
                                                    cuscomplex*     intensities,
                                                    cuscomplex*     raos,
                                                    MLGCmpx*        vel_x_gp,
                                                    MLGCmpx*        vel_y_gp,
                                                    MLGCmpx*        vel_z_gp,
                                                    cuscomplex*     vel_x_total,
                                                    cuscomplex*     vel_y_total,
                                                    cuscomplex*     vel_z_total
                                    );
